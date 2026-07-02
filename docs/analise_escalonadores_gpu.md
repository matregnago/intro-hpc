# Análise comparativa dos escalonadores e da camada de GPU: PaRSEC vs StarPU

Documento de referência consolidando todas as análises discutidas sobre:
(i) os escalonadores `gd` e `lfq` do PaRSEC;
(ii) como a decisão CPU-vs-GPU é tomada no PaRSEC (e por que algumas gemms
acabam indo para CPU mesmo com GPU disponível);
(iii) comparação da nossa versão (PaRSEC v3) com o upstream (v4.0.2411);
(iv) o que o StarPU `dmda`/`dmdas` faz de diferente (com performance model
aprendido em runtime + prefetch).

---

## Parte 1 — Os escalonadores `gd` e `lfq` do PaRSEC

### O ponto-chave que muda tudo

`gd` e `lfq` são **escalonadores de CPU**. Eles decidem *qual thread de
CPU (execution stream) pega qual tarefa*. Eles **NÃO decidem se uma tarefa
vai para GPU** — isso é responsabilidade de uma camada totalmente separada,
o **device layer CUDA** (`parsec/parsec/devices/cuda/dev_cuda.c`), invocada
por `parsec_gpu_get_best_device()`.

Para uma tarefa do Chameleon, o ciclo de vida é:

```
           [escalona quem pega a tarefa]      [decide GPU vs CPU]       [executa]
DAG pronto ────────────────────────────► delegated task ──────────► device layer ─► kernel
                  lfq ou gd                                    parsec_gpu_get_best_device
```

O `lfq`/`gd` só atua quando a tarefa vai rodar em **CPU**. Para tarefas
marcadas como delegáveis a uma GPU, o runtime pula esses escalonadores e
chama o `parsec_gpu_kernel_scheduler` direto.

Também documentado nos próprios scripts do repositório:

- `scripts/doe/doe_gpu_tile.r:30` — *"a escolha CPU-vs-GPU do PaRSEC está
  no device layer (`parsec_gpu_get_best_device`), igual p/ lfq/gd."*
- `scripts/gpu_sweep.sh:40-43` — *"PaRSEC has no perf-model scheduler;
  lfq is its standard ready-task queue and GPU placement is decided by
  PaRSEC's heterogeneous device heuristic regardless."*

### 1.1 `gd` — Global Dequeue (o "simples")

Arquivo: `parsec/parsec/mca/sched/gd/sched_gd_module.c` (191 linhas).
Cabeçalho (`sched_gd.h:14-16`): *"Global Dequeue Scheduler"*.

#### Estrutura de dados
Existe **uma única fila global**, compartilhada por **todos** os threads
de todos os VPs (virtual processes). É um `parsec_dequeue_t` — uma deque
*lock-free* (sem mutex) operável pelas duas pontas.

A inicialização em `flow_gd_init` (linhas 77-123) é esclarecedora:

- O **primeiro** execution stream de cada VP **cria** a deque
  (`OBJ_NEW(parsec_dequeue_t)`, linhas 90/92).
- Todos os **outros** streams, após uma barreira, simplesmente recebem a
  **mesma referência** (linhas 98-106):
  `es->scheduler_object = LOCAL_SCHED_OBJECT(vp->execution_streams[0])`.
- Ou seja: **1 fila, N consumidores**.

#### Inserção: `sched_gd_schedule` (linhas 143-166)

```c
if( (new_context->task_class->flags & PARSEC_HIGH_PRIORITY_TASK) &&
    (0 == distance) ) {
    parsec_dequeue_chain_front( dq, (parsec_list_item_t*)new_context);
} else {
    parsec_dequeue_chain_back( dq, (parsec_list_item_t*)new_context);
}
```

- Pode inserir uma **cadeia inteira** de tarefas (`chain` = encadeia lista
  ligada).
- Alta prioridade + distância 0 → insere pela frente (consumida primeiro).
- Caso contrário → fila FIFO normal (inserção no fim).

#### Seleção: `sched_gd_select` (linhas 125-141)

```c
parsec_task_t *context =
    (parsec_task_t*)parsec_dequeue_try_pop_front( sd );
*distance = 0;
return context;
```

Simplesmente tenta **popar da frente** da deque global. `try_pop_front`
não bloqueia; retorna NULL se vazia. Se voltar NULL, aquele thread fica
ocioso.

#### Características práticas do `gd`

- **Sem noção de localidade**: nada sabe se a tarefa roda melhor no mesmo
  NUMA/socket que produziu seus dados. É só uma fila FIFO global.
- **Sem prioridade entre tarefas comuns**: só o flag
  `PARSEC_HIGH_PRIORITY_TASK` vai para a frente; o resto é estritamente
  FIFO.
- **Sem hierarquia de filas**: menor overhead por tarefa (não há deques
  locais, nem hbbuffer, nem hwloc), **mas** a deque é ponto quente de
  contenção quando há muitos threads (24 no caso do poti) competindo.
- **Balanceamento de carga perfeito**: todos puxam da mesma fila, nenhum
  thread fica ocioso enquanto há trabalho.

> Em uma frase: o `gd` é uma **fila FIFO global lock-free**. Cada worker
> pega a próxima tarefa disponível, ponto final. É o escalonador mais
> barato e mais justo, mas ignora localidade de memória.

### 1.2 `lfq` — Locality First Queue (o "esperto")

Arquivo: `parsec/parsec/mca/sched/lfq/sched_lfq_module.c` (248 linhas).
Cabeçalho: *"Local Flat Queue Scheduler"*.

Aqui aparece a **hierarquia de filas** baseada em `hwloc` (a topologia
real do hardware — núcleos, sockets, NUMA). A estrutura do estado
por-thread está em `sched_local_queues_utils.h:24-34`:

```c
typedef struct {
    parsec_dequeue_t   *system_queue;      // fila de overflow, compartilhada (1 por VP)
    parsec_hbbuffer_t  *task_queue;        // fila LOCAL, privada deste thread
    int                 nb_hierarch_queues;
    parsec_hbbuffer_t **hierarch_queues;   // lista das filas dos vizinhos, da mais próxima à mais distante
} parsec_mca_sched_local_queues_scheduler_object_t;
```

Elementos-chave:

- `parsec_hbbuffer_t` = "**hierarchical bounded buffer**". É um buffer de
  cada worker com tamanho fixo (`queue_size = vp->nb_cores * 4`, linha
  82). Quando enche, transborda (*overflow*) para a fila de sistema.
- A fila de sistema é **uma deque global por VP** (criada só pelo thread
  0, linhas 76-77, e compartilhada, linha 89) — análoga à deque do `gd`,
  mas usada só como "ralo" quando o buffer local enche.

#### Inicialização: `flow_lfq_init` (linhas 64-169)

Cada thread:

1. Cria sua **própria** fila local privada (`parsec_hbbuffer_new`,
   linha 92), com `push_in_system_queue_wrapper` como callback de
   overflow.
2. Constrói uma **lista ordenada das filas dos outros threads**
   (`hierarch_queues`), da **mais próxima** à **mais distante** em
   termos de distância de hardware (linhas 100-138):
   - Se `hwloc` está disponível, usa
     `parsec_hwloc_distance(es->th_id, id)` (linhas 118-135) e agrupa
     por "nível" (linha 124: `d == 2*level || d == 2*level + 1`).
     Resultado: o thread vê primeiro as filas dos colegas do mesmo
     núcleo físico / mesmo L1/L2 / mesmo L3 / mesmo NUMA, e por último
     os do outro socket.
   - Se `hwloc` não está disponível, faz Round-Robin simples (linhas
     110-113).

#### Seleção: `sched_lfq_select` (linhas 171-201) — o coração do "Locality First"

```c
// 1. Tenta minha própria fila local
task = parsec_hbbuffer_pop_best(task_queue, priority_comparator);
if(task) { *distance = 0; return task; }

// 2. Tenta as filas dos vizinhos, da mais próxima à mais distante
for(i = 0; i < nb_hierarch_queues; i++) {
    task = parsec_hbbuffer_pop_best(hierarch_queues[i], priority_comparator);
    if(task) { *distance = i + 1; return task; }
}

// 3. Último recurso: fila de sistema global
task = parsec_mca_sched_pop_from_system_queue_wrapper(...);
if(task) *distance = 1 + nb_hierarch_queues;  // distância máxima
return task;
```

A ordem é o ponto:

1. **Próprio buffer local** → máxima localidade (o dado provavelmente foi
   produzido neste mesmo núcleo).
2. **Filas dos vizinhos, em ordem crescente de distância hwloc** → rouba
   trabalho (*work stealing*) do vizinho mais barato primeiro.
3. **Fila de sistema** → quando tudo mais está vazio, busca do ralo
   global compartilhado.

O `*distance` retornado indica "a que distância" a tarefa foi encontrada
(0 = local, crescente = cada vez mais longe). Essa distância é passada
de volta para `sched_lfq_schedule` quando a tarefa gera sucessoras, para
que elas **sejam inseridas na mesma vizinhança** (afinidade de
produtor-consumidor).

#### Inserção: `sched_lfq_schedule` (linhas 203-211)

```c
parsec_hbbuffer_push_all(task_queue, new_context, distance);
```

Insere no **próprio buffer local**, passando a `distance` herdada. Se o
buffer encher, o callback
`parsec_mca_sched_push_in_system_queue_wrapper` (em
`sched_local_queues_utils.h:72-82`) faz o overflow cair na **deque de
sistema**.

> Em uma frase: o `lfq` é um **work-stealing hierárquico com consciência
> de localidade** — cada worker tem seu buffer local; ao esvaziar, ele
> rouba do vizinho mais próximo (topologicamente) primeiro, e só vai ao
> ralo global em último caso.

### 1.3 `gd` vs `lfq` lado a lado

| Aspecto | `gd` (Global Dequeue) | `lfq` (Locality First Queue) |
|---|---|---|
| Nº de filas | **1** (deque global) | **N+1**: 1 buffer local por thread + 1 deque de sistema/overflow |
| Consciência de topologia | Nenhuma | **Sim**, via `hwloc` (sockets, NUMA, caches) |
| Estratégia de select | pop_front da fila global | local → vizinhos próximos → ... → ralo global |
| Work stealing | Não (todo mundo puxa do mesmo lugar) | **Sim**, hierárquico por distância |
| Custo por tarefa | muito baixo (uma operação lock-free) | maior (buffer local + hierarquia, mas sem lock na maioria dos casos) |
| Contenção sob 24 threads | **alta** (todos competem pela mesma deque) | **baixa** (cada thread passa a maior parte do tempo no seu buffer privado) |
| Bom para | workloads com muitas tarefas curtas e sem reuso de dado | workloads com **afinidade de dado** (Cholesky/QR, onde tiles são iterativamente reusados) |
| Equivalente StarPU | `lws` (local work stealing) — embora lfq seja ainda mais localizado | próximo a `ws`/`dmdar` em espírito, **sem** performance model |

Para a workload deste trabalho (Cholesky/QR Chameleon, 24 threads,
tiles de 480-8000), o `lfq` tende a vencer o `gd` porque:

- A fatoração tem **forte reuso de tiles** (cada tile é lido/escrito por
  várias tarefas consecutivas);
- Mantê-las no cache do mesmo núcleo (produtor → consumidor no mesmo
  buffer) reduz tráfego de memória;
- 24 threads competindo por uma única deque (`gd`) gera contenção
  mensurável.

---

## Parte 2 — Como fica o caminho quando há uma GPU (RTX 4070)

Quando `PARSEC_MCA_device_cuda_enabled=1` (configurado em
`gpu_doe_sweep.sh:161`), o registro de dispositivos em `parsec_gpu_init`
(`dev_cuda.c:504-609`) cria **uma `gpu_device_t` por GPU física**, e a
registra via `parsec_devices_add` (linha 609). Cada dispositivo é
identificado por um `device_index`.

### 2.1 A decisão CPU-vs-GPU: `parsec_gpu_get_best_device` (`dev_cuda.c:1323-1362`)

```c
int parsec_gpu_get_best_device(parsec_task_t* this_task, double ratio) {
    // Passo 1: achamos o flow de SAÍDA em modo WRITE da tarefa
    for(i = 0; i < nb_flows; i++)
        if(out[i] & FLOW_ACCESS_WRITE) {
            dev_index = data_in->original->owner_device;   // quem "possui" o dado hoje
            if(dev_index > 1) break;   // já está em uma GPU, mantém lá
        }

    // Passo 2: se é a 1ª vez (dono é CPU, index 0), escolhemos a melhor GPU
    if(dev_index <= 1) {
        best_weight = parsec_device_load[0] + ratio * parsec_device_sweight[0];   // custo CPU
        for(dev_index = 2; dev_index < n_devices; dev_index++) {
            if(!(tp->devices_index_mask & (1<<dev_index))) continue;   // GPU desabilitada p/ este taskpool
            weight = parsec_device_load[dev] + ratio * parsec_device_sweight[dev];   // custo GPU
            if(best_weight > weight) { best_index = dev; best_weight = weight; }    // menor ganha
        }
        dev_index = best_index;
    }
    return dev_index;
}
```

O **algoritmo de colocação** é:

1. **Afinidade de dado primeiro**: se a tarefa grava em um tile que **já
   foi alocado numa GPU**, ela vai para essa mesma GPU (mantém o dado
   local, evita D2H+H2D).
2. **Caso contrário**, equilíbrio de carga ponderado: cada dispositivo
   tem dois contadores:
   - `parsec_device_load[dev]` — a soma da carga (flops/ciclos) das
     tarefas atualmente em voo nele; **dinâmico**.
   - `parsec_device_sweight[dev]` — peso **estático** do dispositivo =
     `(SM count) × (throughput rate) × (clock)`, computado uma vez em
     `parsec_gpu_init` (linhas 575-578). Para a CPU vale
     `parsec_device_load[0] + ratio × sweight[0]`.
   - `ratio` é um fator de compensação que o Chameleon passa ao chamar;
     pondera entre "tarefa já mapeada" e "tarefa pronta para a GPU".
3. O dispositivo com **menor peso** (carga atual + ajuste estático)
   ganha.

**Crucial**: `lfq` e `gd` **não participam** dessa decisão. A função é
chamada pela engine quando a `task_class` tem `PARSEC_FOLLOW_GPU` ou
está registrada como CUDA-capable. Para essas tarefas, o caminho de
escalonamento é completamente diferente:

```
tarefa GPU-ready ─► engine delega ─► parsec_gpu_kernel_scheduler(es, gpu_task, which_gpu=best)
                                              │
                          (cf. dev_cuda.c:1967-2124)
```

### 2.2 As 6 streams por GPU e o que elas fazem

`PARSEC_MAX_STREAMS = 6` (`dev_cuda.h:26`). Para **cada GPU** são
criadas 6 *CUDA streams* (`dev_cuda.c:508-557`):

| `j` | nome | papel |
|---|---|---|
| **0** | `h2d(0)` | **Host → Device** (movein). `prof_event_key=movein_key`. |
| **1** | `d2h(1)` | **Device → Host** (moveout). `prof_event_key=moveout_key`. |
| **2 a 5** | `cuda(2..5)` | **execução de kernels** (4 streams concorrentes). |

Cada stream tem um buffer circular de até
`PARSEC_MAX_EVENTS_PER_STREAM = 4` eventos em voo (`dev_cuda.h:27`). A
sincronização é feita por **`cudaEvent_t`** gravados na stream após
cada operação — nunca blocking (`cudaEventQuery`).

O `exec_stream` usado para executar é escolhido em **round-robin**
(`dev_cuda.c:2038`):

```c
exec_stream = (exec_stream + 1) % (gpu_device->max_exec_streams - 2);  // max_streams=6 - 2 = 4 streams de execução
```

ou seja, dentre os índices 2..5.

### 2.3 O pipeline em `parsec_gpu_kernel_scheduler` (`dev_cuda.c:1967-2124`)

Pseudocódigo simplificado:

```
1. Toma o mutex da GPU (atomic). Se != 0, outro thread está mexendo:
   empurra a tarefa em gpu_device->pending (FIFO) e sai (RETURN_ASYNC).
2. cudaSetDevice(gpu)
3. check_in_deps:
     progress_stream(h2d,  parsec_gpu_kernel_push,   gpu_task)   ← movein
        · chama parsec_gpu_data_reserve_device_space (verifica memória; se faltar, cria tarefa W2R p/ liberar LRU)
        · para cada flow com READ/OUTPUT: parsec_gpu_data_stage_in:
            · se o dado está em versão diferente da cópia device → cudaMemcpyAsync H2D
            · marca data_transfer_status = UNDER_TRANSFER
            · grava cudaEvent no stream h2d
     · exec_stream = (exec_stream+1) % 4   ← escolhe cuda(2..5) em RR
     progress_stream(cuda, NULL, gpu_task)                          ← execução
        · chama submit do Chameleon (parsec_gpu_task_t.submit)
        · grava cudaEvent
4. get_data_out_of_device:
     progress_stream(d2h,  parsec_gpu_kernel_pop,    progress_task)  ← moveout
        · para cada flow WRITE com pushout=1: cudaMemcpyAsync D2H
        · dados só-LEITURA voltam para gpu_mem_lru (cache de device)
        · decrementa parsec_device_load[dev] em gpu_task->load
5. fetch_task_from_shared_queue:
     pega mais uma tarefa de gpu_device->pending (FIFO) e volta a check_in_deps
6. complete_task:
     __parsec_complete_execution(es, task)      ← avança dependências do DAG
```

**Insight fundamental**: a GPU é operada por **um único "dono"** por
vez (o mutex atomic em linha 1991). Quem chega primeiro marcha pelas
três fases (movein → exec → moveout), e enquanto isso coloca as
próximas tarefas na **`gpu_device->pending`** (uma FIFO, linhas 1993 e
2101). Ou seja: **a partir do momento em que a tarefa é delegada à GPU,
o escalonador lfq/gd perdeu a jurisdição** — a partir daí é a
`gpu_device->pending` que dita a ordem, e a ordem é **FIFO** (com lapso
de `PARSEC_PUSH_TASK` ordenado por prioridade só na
`exec_stream[0].fifo_pending`, linhas 1427/1366-1371).

### 2.4 A transferência host ↔ device em detalhe

**H2D — `parsec_gpu_data_stage_in` (`dev_cuda.c:1065-1159+`)**:

- `parsec_data_transfer_ownership_to_copy(original, device_index, type)`
  decide se precisa transferir comparando **versões** do dado no host
  (`in_elem->version`) e no device (`gpu_elem->version`). Se iguais →
  *coherent*, **nenhum byte transferido** (linha 1155+).
- Se precisa: `cudaMemcpyAsync(dst=gpu_elem->device_private,
  src=in_elem->device_private, size=nb_elts, cudaMemcpyHostToDevice,
  h2d_stream->cuda_stream)` (linhas 1139-1142).
- Marca o dado como `DATA_STATUS_UNDER_TRANSFER` e registra qual tarefa
  iniciou o push (`gpu_elem->push_task = gpu_task->ec`). Só essa tarefa
  pode baixar o status para COMPLETE (em
  `parsec_gpu_callback_complete_push`, linhas 1377-1426).
- Contabiliza em `gpu_device->super.transferred_data_in += nb_elts`
  (linha 1146) — número que aparece em
  `PARSEC_MCA_device_show_statistics=1`.

**D2H — `parsec_gpu_kernel_pop` (`dev_cuda.c:1681-1808`)**:

- Só faz D2H para flows que têm `FLOW_ACCESS_WRITE` **e** o bit
  `gpu_task->pushout & (1<<i)` está setado (linha 1765). Ou seja: o
  runtime tem **controle granular** sobre quais resultados voltam — só
  escreve de volta o que uma tarefa consumidora futura realmente
  precisará na CPU.
- Para dados só-LEITURA, não há cópia de volta: a cópia device vai para
  `gpu_mem_lru` (linha 1744) como **cache** (afinidade de device — é o
  cachê LRU `gpu_device->gpu_mem_lru`).
- Usa o stream `d2h(1)` (`gpu_stream->cuda_stream` passado é o stream
  1, linhagem 2076-2078).
- `transferred_data_out += nb_elts` (linha 1790).

**A tarefa "W2R" (Write-to-Read)**: quando não há memória livre no
device (`parsec_gpu_data_reserve_device_space` retorna negativo, linhas
1634-1636), o scheduler cria uma tarefa sintética
(`parsec_gpu_create_W2R_task`, linha 2031) do tipo
`GPU_TASK_TYPE_D2HTRANSFER` (constante 111, `dev_cuda.h:45`). Essa
tarefa **força** um D2H de algum bloco sujo da `gpu_mem_owned_lru` para
liberar espaço no device — é a **evicção controlada**. É o "thrashing"
que `docs/plano_experimentos_gpu.md:17` cita como métrica causal nos
experimentos (*"tráfego H2D/D2H na GPU - thrashing do pushout do DTD do
PaRSEC"*) — o `lfq`/`gd` não conseguem reduzir, só uma política melhor
de pushout poderia.

### 2.5 Síntese do caminho completo de UMA tarefa GPU

```
DAG: tarefa fica pronta (deps satisfeitas)
   │
   ├─ se task_class é CUDA-capable ─► engine chama parsec_gpu_get_best_device
   │     (decisão: CPU porque é o melhor? GPU? qual GPU?)   ← AQUI lfq/gd não opinam
   │
   ├─ se CPU: tarefa é inserida no escalonador (lfq: buffer local; gd: deque global)
   │     e pega um thread worker normal.
   │
   └─ se GPU: tarefa é encapsulada em parsec_gpu_task_t e:
         engine chama parsec_gpu_kernel_scheduler(es, gpu_task, which_gpu)
            (1) mutex atômico na GPU. Se livre, sigo. Se não, push em gpu_device->pending.
            (2) H2D: p/ cada flow READ/OUTPUT desatualizado ─► cudaMemcpyAsync H2D no stream h2d(0)
                    (se faltar memória: evict W2R do LRU antes)
            (3) EXEC: escolhe cuda(2..5) por round-robin, chama gpu_task->submit (=kernel Chameleon)
            (4) D2H: p/ cada flow WRITE com pushout ─► cudaMemcpyAsync D2H no stream d2h(1)
                    (dados só-leitura ficam na cache device LRU p/ a próxima)
            (5) complete: __parsec_complete_execution (destrava sucessoras no DAG)
            (6) loop: pega próxima tarefa em gpu_device->pending (FIFO)
```

---

## Parte 3 — O que define CPU-vs-GPU no PaRSEC: a peça que faltava

Os três mecanismos que, juntos, definem o caminho:

### 3.1 Registro das codelets — "*esta tarefa pode rodar em GPU*"

Cada `parsec_task_class_t` tem uma lista de **incarnations** (em DTD, é o
array `incarnations[]` modificado em
`parsec/parsec/interfaces/superscalar/insert_function.c`). Uma
incarnation é um par `(device_type, função)`. Para o `dgemm` do
Chameleon há duas entradas:

```c
{.type = PARSEC_DEV_CPU,     .hook = cl_dgemm_cpu_func },
{.type = PARSEC_DEV_CUDA,    .hook = cl_dgemm_cuda_func },
```

É o Chameleon (via `parsec_dtd_taskpool_insert_task_cuda` no
runtime/backend PaRSEC) que declara essas incarnations ao registrar cada
codelet. **Se a codelet só tem a incarnation CPU, a tarefa nunca vai
para GPU** — a engine nunca a considera como candidata a GPU.

Em `chameleon/runtime/parsec/codelets/codelet_zgemm.c:126`:

```c
parsec_dtd_taskpool_insert_task_cuda(
    PARSEC_dtd_taskpool,
    CORE_zgemm_parsec,        // CPU function
    CORE_zgemm_parsec_cuda,   // CUDA function
    ...);
```

A versão `_cuda` do insert registra **duas incarnations** para a mesma
task class: `(PARSEC_DEV_CPU, CORE_zgemm_parsec)` e
`(PARSEC_DEV_CUDA, CORE_zgemm_parsec_cuda)`. Então **a tarefa é
GPU-capable**.

### 3.2 A máscara de dispositivos por taskpool

Quando o Chameleon cria o taskpool, ele registra quais dispositivos
aceita. Em `parsec_gpu_get_best_device` (`dev_cuda.c:1350`):

```c
if(!(tp->devices_index_mask & (1 << dev_index))) continue;
```

Mesmo que a codelet tenha incarnation CUDA, **se o taskpool não listar
aquela GPU na sua máscara**, ela fica de fora.

### 3.3 O algoritmo `parsec_gpu_get_best_device` detalhado

Quando uma tarefa GPU-ready fica pronta e a engine precisa decidir para
onde mandá-la, chama-se `parsec_gpu_get_best_device(task, ratio)`
(`dev_cuda.c:1323-1362`). Há **apenas dois critérios**:

**Critério A — Afinidade de dado (dono atual do output write)** (linhas
1329-1338): se o dado que esta tarefa vai escrever **já está alocado numa
GPU** (porque uma tarefa anterior o colocou lá), a tarefa **vai para a
mesma GPU**. É a regra de *data affinity*: minimiza transfers.

**Critério B — Se é a primeira vez (dono é CPU, index 0): minimize carga
ponderada** (linhas 1342-1358):

```c
best_weight = parsec_device_load[0] + ratio * parsec_device_sweight[0];   // custo CPU
for(dev_index = 2; dev_index < n_devices; dev_index++) {
    weight = parsec_device_load[dev] + ratio * parsec_device_sweight[dev]; // custo GPU
    if(best_weight > weight) { best_index = dev; best_weight = weight; }   // menor ganha
}
```

Onde:
- `parsec_device_load[dev]` — soma das `gpu_task->load`s atualmente em
  voo naquele dispositivo. **Load dinâmico, mas `gpu_task->load` é
  declarado pela codelet (estático), não medido.**
- `parsec_device_sweight[dev]` — peso estático do hardware, computado
  uma única vez em `parsec_gpu_init` (linhas 575-578):
  `SM_count × throughput_rate × clock`.
- `ratio` — escalar passado pelo Chameleon; pondera a importância
  relativa do peso estático vs. load atual.

**Crucial**: isso é **estático**. Não há "tempo previsto de execução
aprendido em runtime" (como o
`starpu_task_worker_expected_length` do StarPU). A `load` da codelet é
um escalar fixo declarado pelo Chameleon ao registrar a task_class — não
é atualizado por medições.

### 3.4 Onde o `lfq`/`gd` entra? Só em paralelo, como "coadjuvantes"

Para tarefas **não-GPU** (todas as demais do Cholesky/QR: `dpotrf`,
`dtrsm`, `dsyrk` que você não marcou como CUDA-capable), o `lfq`/`gd`
decide qual thread de CPU executa. É a única jurisdição deles.

---

## Parte 4 — Como `device_sweight` é realmente computado (e por que gemms vão para CPU)

Em `parsec/parsec/devices/device.c:314-363`, a função
`parsec_devices_freeze` faz uma normalização **inesperada** dos pesos:

1. **Primeiro copia reto** os `device_sweight[i]` brutos para o array
   global (linhas 335-338) e soma em `total_sperf`.
2. **Depois inverte** (linhas 352-355) — substitui cada peso por
   **`total / peso`**:

```c
parsec_device_sweight[i] = (total_sperf / parsec_device_sweight[i]);
```

Ou seja, o sweight final de cada dispositivo é
**`total_de_todos / peso_bruto`** — uma **razão relativa inversa**, não
o throughput:

- O dispositivo **mais rápido** fica com **menor sweight** (porque
  `total/grande` = pequeno).
- O dispositivo **mais lento** fica com **sweight maior**.

### 4.1 O `device_sweight` da CPU não é zero

Em `parsec/parsec/devices/device.c:489`:

```c
device->device_sweight = nstreams * fp_ipc * freq;
```

Onde (linhas 432-448):
- `fp_ipc` (single-precision IPC): **64** se AVX-512, **32** se AVX2,
  **16** se AVX, **8** senão.
- `dp_ipc` (double-precision): 32 / 16 / 8 / 4 respectivamente.
- `nstreams` = número de streams PaRSEC (= número de threads de
  computação, 24 no poti).
- `freq` = clock do CPU em GHz (lido em `/proc/cpuinfo` —
  `clock : XMHz`, ou do model string `@ XGHz`).

Para o poti (Xeon com AVX2 típico), algo como `freq ≈ 2.5 GHz`,
`fp_ipc = 32`, `nstreams = 24`:

- `device_sweight[CPU] = 24 × 32 × 2.5 = **1920**` (Gflops/s teóricos
  single)

Para a RTX 4070 (`dev_cuda.c:575-578`):

- `device_sweight[GPU] = streaming_multiprocessor × srate × clockRate ×
  2e-3`
- SM=46, srate depende de CUDA (FP32 ~ 256 ops/clock/SM típico de Ada),
  clockRate ≈ 2.5 GHz:
- valor bruto ~ 46 × 256 × 2500 × 2e-3 ≈ **59000** Gflops/s

### 4.2 A inversão e o que isso significa

Para `total_sperf ≈ 60920`:

| Dispositivo | sweight bruto | sweight após inversão (`total/bruto`) |
|---|---|---|
| CPU | ~1920 | `60920 / 1920` ≈ **31.7** |
| GPU | ~59000 | `60920 / 59000` ≈ **1.03** |

Com `load = 0` (sem tarefas em voo), `ratio = 1`:

- `weight[CPU] = 0 + 1 × 31.7 = 31.7`
- `weight[GPU] = 0 + 1 × 1.03 = 1.03`
- `weight[CPU] > weight[GPU]` → GPU ganha → tarefa vai pra GPU. ✓

### 4.3 Por que gemms vão para CPU quando o sistema está em movimento

O `parsec_device_load[gpu]` cresce a cada tarefa submetida na GPU, e só
é decrementado ao final (`dev_cuda.c:2126`:
`parsec_device_load[gpu] -= gpu_task->load`). O `gpu_task->load` é
declarado pela codelet; para o Chameleon, provavelmente `1` (default
comum).

Cada gemm submetida à GPU faz `parsec_device_load[gpu] += 1`. Após **~30
gemms em voo**:

- `weight[GPU] = 30 + 1.03 = 31.03`
- `weight[CPU] = 0 + 31.7 = 31.7`
- `weight[GPU] ≈ weight[CPU]` → aproxima do empate.

Um pouquinho mais (31 gemms pending):

- `weight[GPU] = 31 + 1.03 = 32.03 > 31.7 = weight[CPU]`
- **CPU ganha** → próxima gemm vai pra CPU.

Em outras palavras: o PaRSEC **deliberadamente distribui** gemms entre
CPU e GPU quando tem ~30 tarefas na fila de pending da GPU. Não é bug —
é balanceamento de carga com base em contagem de tarefas. Mas o problema
é:

1. **A `load` da codelet não reflete a duração real** (`load=1` para uma
   GEMM N=960 é o mesmo que `load=1` para N=480).
2. **Sem previsão de transferência** — se a GPU está congestionada em
   H2D (só 1 stream `h2d(0)` para todas as tarefas!), o
   `device_load[gpu]` não captura o custo do bottleneck de
   transferência, cresce só pelo número de kernels em voo.
3. **Sem previsão de speedup real** — o `device_sweight` teórico
   despreza JIT, memory bandwidth e a real diferença assintótica:
   gemm_N=960 na CPU pode levar ~50 ms; na GPU ~0.5 ms (100×). Mas o
   modelo entende como se fosse só a razão de TFLOPS pico (~30×).

### 4.4 Os três motivos juntos, resumidos

```
Soma dos load's estáticos das tarefas pending
    (cada gemm conta +1 — não sabe quão grande é)

      weight[GPU] = device_load[gpu] + ratio × sweight_inverso[gpu]
                    ↑ cresce a cada gemm submetida      ↑ constante baixa (~1)
                    └ só desce quando a tarefa termina

      weight[CPU] = device_load[cpu] + ratio × sweight_inverso[cpu]
                    ↑ estável (CPU não fica "load submitted")   ↑ constante alta (~32)

      escolha = argmin(weight)
      CPU ganha quando: device_load[gpu] + 1 > 0 + 32
                       ≈ 31 gemms pending na GPU (ou menos, se load != 1)
```

Quando isso acontece:

- 30+ gemms em voo na GPU (porque a GPU não conseguiu drenar tão rápido
  quanto o DAG liberou).
- **merecedor de CPU**: o algoritmo decide mandar a próxima pra CPU,
  porque "a GPU está mais cheia".
- Mas a GPU, mesmo com 30 pending, terminaria cada gemm em ~0.5 ms —
  enquanto a CPU levaria 50 ms em uma só. Resultado: **a tarefa demora
  muito mais** na CPU do que se tivesse esperado a fila da GPU.

Contraste com StarPU `dmda`:

- `exp_end[CPU][nimpl]` é **previsto** pela perfmodel: ~50 ms (aprendido
  em runtime).
- `exp_end[GPU][nimpl]` é **previsto**: considerando a fila de pending,
  ~5 ms.
- `fitness[CPU] = alpha × (5 - 5) + ...` é bom se CPU estiver livre.
- `fitness[GPU] = alpha × (30 × 0.5 - 5) + ... = 10 + …` pior que CPU?
  Sim! Mas `local_data_penalty[CPU]` seria enorme (precisa trazer os
  tiles de device pra host — `beta × transfer_large`), e o fitness final
  preferiria esperar a GPU.

No PaRSEC o algoritmo não capta isso porque:

- Não há `local_data_penalty` (sem data-aware modeling, só
  `device_load`).
- Não há perfmodel aprendido.
- A "carga" é só contagem nominal.

---

## Parte 5 — Comparação com o upstream PaRSEC (`parsec-atual/`, v4.0.2411)

A versão upstream passou por uma **reforma arquitetural massiva**, não
apenas adições incrementais. O CHANGELOG `v4.0.2411` confirma: migra
para "PaRSEC API 4.0", remove a API 3.0, e reescreve a camada de
dispositivos como framework MCA.

### 5.1 Escalonadores (mca/sched/)

#### Novo escalonador: `llp`
Upstream adicionou `mca/sched/llp/` (657 linhas). É o **"Local LIFO with
support for Priorities"** — cada thread mantém uma pilha (LIFO, via
`parsec_lifo`) ordenada por prioridade; work-stealing entre threads.
Prioridade 2 (mais baixa, no mesmo nível do `ll`).

Nossa árvore tem 11 escalonadores (`ap gd ip lfq lhq ll ltq pbq rnd
spq` + default); upstream tem **12** (acrescenta `llp`).

#### Mudanças mecânicas em todos os escalonadores existentes
Nenhuma assinatura de API mudou. Mas o código foi reescrito:

| Mudança | Impacto |
|---|---|
| `return 0` → `return PARSEC_SUCCESS` | semântico (`PARSEC_SUCCESS` agora vale `0` no upstream, era `-1` na v3) |
| `OBJ_NEW`/`OBJ_RETAIN`/`OBJ_RELEASE`/`OBJ_DESTRUCT` → `PARSEC_OBJ_*` | renomeação de namespace |
| `papi_sde_register_fp_counter(handle, …)` → `parsec_papi_sde_register_fp_counter(…)` (sem handle global) | wrapper PaRSEC ao redor do PAPI; prefixo dos contadores `"PARSEC::SCHEDULER::…"` → `"SCHEDULER::…"` |
| Macro `TAKE_TIME(…)` removida de `lfq`/`ll` | era no-op em build sem prof trace |
| Copyright atualizado; `lfq` agora traz *"Copyright (c) 2024 NVIDIA Corporation"* | sinal de que NVIDIA está usando lfq como backend GPU-CPU |

#### Bug fix importante em `sched_local_queues_utils.h`
`parsec_mca_sched_list_local_counter_pop_back()` na nossa versão chama
**`parsec_list_pop_front`** (bug — diz "pop back" mas faz "pop front").
Upstream corrigiu para `parsec_list_pop_back`. **Isso afeta os
escalonadores `ltq`, `lhq`, `ap` que usam esse utilitário.**

#### Remoção do enum hard-coded
Nosso `parsec/dplasma/testing/common.c:31-41` tem
`char *PARSEC_SCHED_NAME[] = {"", "lfq", "ltq", ...}`. Upstream
**deletou `dplasma/` inteiro**; os testes agora iteram `${MCA_sched}`
descoberto pelo CMake e passam `--mca mca_sched <name>` em runtime.

#### Mudanças no despachante (`parsec/scheduling.c`)

- `current_scheduler` → `parsec_current_scheduler`.
- `__parsec_schedule` agora é `inline` e ganhou
  **`PARSEC_PINS(local_es, SCHEDULE_BEGIN/END, tasks_ring)`** — hooks
  de instrumentação PINS em volta de cada invocação do escalonador
  (novo na v4).
- Novo helper `__parsec_schedule_vp(…)` para submeter array de anéis com
  uma entrada por VP.
- `task->data[i].data_repo` → `task->data[i].source_repo_entry` (rename).
- `FLOW_ACCESS_*` → `PARSEC_FLOW_ACCESS_*`.

#### Prioridades (component_query)
**Não mudaram**: `lfq` continua prioridade 20 (default, mais alta); `gd`
continua 10; tudo igual.

### 5.2 Camada de dispositivos — reforma arquitetural completa

Esta é a maior mudança. Upstream **migrada para um framework MCA** com
backends plugin.

#### Reestruturação de paths

```
NOSSA VERSÃO (A)                      UPSTREAM (B)
parsec/parsec/devices/device.c (646)  → parsec/mca/device/device.c (1229)
parsec/parsec/devices/device.h (156)  → parsec/mca/device/device.h (381)
parsec/parsec/devices/cuda/dev_cuda.c (2149) → parsec/mca/device/device_gpu.c (3613) + cuda/device_cuda_module.c (717) + cuda/device_cuda_component.c (297)
parsec/parsec/devices/cuda/dev_cuda.h (242)  → parsec/mca/device/device_gpu.h (452) + cuda/device_cuda.h (82) + cuda/device_cuda_internal.h (35)
parsec/parsec/devices/cuda/transfer.c (315)  → parsec/mca/device/transfer_gpu.c (362)
```

#### Três novos backends de acelerador

- **`mca/device/level_zero/`** — Intel GPU via oneAPI Level Zero, **com
  bridge SYCL/DPC++** (`device_level_zero_dpcpp_interface.cpp`). 1284
  linhas em 7 arquivos. Implementa os mesmos callbacks do CUDA.
- **`mca/device/hip/`** — AMD ROCm (por enquanto placeholder, só
  `ValidateModule.CMake`).
- **`mca/device/template/`** — esqueleto documentado para novos
  backends.

Bits de tipo de dispositivo:

- `PARSEC_DEV_INTEL_PHI` (1<<3) → **`PARSEC_DEV_HIP`**
- `PARSEC_DEV_OPENCL` (1<<4) → **`PARSEC_DEV_LEVEL_ZERO`**
- **Novo** `PARSEC_DEV_TEMPLATE` (1<<7)
- **Novo** `PARSEC_DEV_CHORE_ALLOW_BATCH` (0x100) — habilita batched
  submission
- **Novo** `PARSEC_DEV_GPU_MASK = CUDA|HIP|LEVEL_ZERO` e helper
  `PARSEC_DEV_IS_GPU(t)`

#### Abstração da backend via vtable
A struct `parsec_device_gpu_module_s` (em `device_gpu.h:246-281`) agora
carrega **function pointers** que cada backend preenche:

```c
parsec_device_set_device_fn_t       set_device;
parsec_device_memcpy_async_fn_t     memcpy_async;
parsec_device_event_query_fn_t      event_query;
parsec_device_event_record_fn_t     event_record;
parsec_device_memory_info_fn_t      memory_info;
parsec_device_memory_allocate_fn_t  memory_allocate;
parsec_device_memory_free_fn_t      memory_free;
parsec_device_find_incarnation_fn_t find_incarnation;
```

Na nossa versão, o engine GPU chamava
`cudaMemcpyAsync`/`cudaEventRecord`/`cudaSetDevice`/etc. **direto**.
Agora cada backend implementa os equivalentes (CUDA em
`device_cuda_module.c:318-397`; Level Zero em
`device_level_zero_module.c`).

#### Streams sem fields CUDA
`parsec_gpu_exec_stream_t` (nossa `dev_cuda.h:89-113`) embute
`cudaEvent_t *events` e `cudaStream_t cuda_stream` **diretamente na
struct**. Upstream (`device_gpu.h:283-298`) **remove esses campos
CUDA** — eles vivem numa subclasse backend-específica. Stream comum
ganha `parsec_info_object_array_t infos` (per-stream info/MCA params).

`PARSEC_MAX_STREAMS` foi renomeado para **`PARSEC_GPU_MAX_STREAMS`**
(mesmo valor: 6). `PARSEC_MAX_EVENTS_PER_STREAM` e
`PARSEC_GPU_MAX_WORKSPACE` inalterados.

#### API de device publicamente renomeada
Toda a família `parsec_devices_*` → `parsec_mca_device_*`, e
`parsec_device_t` → `parsec_device_module_t` (agora é `parsec_object_t`,
ref-counted). Funções novas:

- `parsec_mca_device_is_gpu(devindex)` — novo
- `parsec_devices_release_memory(void)` — novo
- `parsec_devices_reset_load(context)` — novo
- `parsec_devices_save_statistics / _free_statistics / _print_statistics`
  — split de `dump_and_reset_statistics`
- `parsec_advise_data_on_device(data, device, advice)` — **novo
  mecanismo de advice** (PREFETCH/PREFERRED_DEVICE/WARMUP)
- `parsec_select_best_device(task)` — substitui
  `parsec_gpu_get_best_device(task, ratio)` (sem `ratio` arg, movida
  para `device.c`)

### 5.3 GPU engine (`parsec_gpu_*` → `parsec_device_*`)

Todas as funções centrais mudaram de nome e assinatura:

| Nossa versão | Upstream | Mudança |
|---|---|---|
| `parsec_gpu_get_best_device(task, ratio)` | `parsec_select_best_device(task)` | sem `ratio`; movida pra `device.c` |
| `parsec_gpu_kernel_scheduler(es, gpu_task, which_gpu)` | `parsec_device_kernel_scheduler(module, es, void*)` | `which_gpu` vira `module` (device abstraído) |
| `parsec_gpu_kernel_push(gpu_device*, task, stream)` | `parsec_device_kernel_push(gpu_module*, es, task, stream)` | adicionou `es` |
| `parsec_gpu_kernel_pop(...)` | `parsec_device_kernel_pop(...)` | mesma reestruturação |
| `parsec_gpu_data_stage_in` | `parsec_device_data_stage_in` | type reescrito |
| `parsec_gpu_create_W2R_task` | `parsec_gpu_create_w2r_task` | snake_case |
| `parsec_gpu_W2R_task_fini` | `parsec_gpu_complete_w2r_task` | rename |
| `parsec_gpu_sort_pending_list` | `parsec_device_sort_pending_list` | rename |
| `parsec_gpu_push_workspace`/`pop_workspace`/`free_workspace` | `parsec_device_*` | rename |

Compute capabilities suportadas — CUDA: nossa versão
`{10,11,12,13,20,21,30,32,35,37,50,52,53,60,61,62,70}`; upstream
`{35,37,50,52,53,60,61,62,70,72,75,80,90}` — **dropou as antigas,
adicionou 7.2/7.5/8.0/9.0** (Hopper).

#### Mutex do kernel_scheduler
Nossa versão: `parsec_atomic_fetch_inc_int32(&gpu_device->mutex)` — uma
tentativa só, se != 0 desiste (ASYNC). Upstream: **busy-wait CAS loop**
que distingue `rc < 0` de `rc >= 0` — tenta mais agressivamente,
eventualmente bloqueia.

### 5.4 Transferências host↔device — várias mudanças

#### Constantes de task type
Nossa versão tem só `GPU_TASK_TYPE_D2HTRANSFER = 111`. Upstream refez
como **bitmask 4-bit**:

```c
PARSEC_GPU_TASK_TYPE_KERNEL       0x0000
PARSEC_GPU_TASK_TYPE_D2HTRANSFER  0x1000
PARSEC_GPU_TASK_TYPE_PREFETCH     0x2000   // NOVO
PARSEC_GPU_TASK_TYPE_WARMUP       0x4000   // NOVO
PARSEC_GPU_TASK_TYPE_D2D_COMPLETE 0x8000   // NOVO — device-to-device
PARSEC_GPU_TASK_TYPE_INVALID      0xf000
```

#### P2P: Transferências D2D peer-to-peer (NOVO)
Upstream introduz um **enum de direção** que não existe na nossa versão:

```c
typedef enum {
    parsec_device_gpu_transfer_direction_h2d,
    parsec_device_gpu_transfer_direction_d2h,
    parsec_device_gpu_transfer_direction_d2d   // NOVO
} parsec_device_transfer_direction_t;
```

Em `parsec_default_gpu_stage_in/out`, o engine agora verifica
`src_dev->peer_access_mask & (1 << dst_dev->device_index)` e, se peer
access está on, emite `memcpy_async(... d2d)` direto — **sem bounce
buffer pelo host**. Isso é enorme para multi-GPU.

Há também um novo tipo de tarefa sintética
`PARSEC_GPU_TASK_TYPE_D2D_COMPLETE` que avisa o dispositivo destino que
a transferência D2D pousou.

#### Multi-GPU / NVLink seletivo
Nossa versão habilita peer-access para **todos** os pares possíveis,
incondicionalmente. Upstream adiciona um MCA param:

- **`device_cuda_nvlink_mask`** (default `0xffffffff`) — bitmask por
  source GPU, permite desligar peer access seletivamente.
- O loop de peer-access sai de `parsec_gpu_init` e vira um hook
  `parsec_cuda_all_devices_attached(device)` — só roda depois que
  **todos** devices (incluindo não-CUDA) estão registrados, permitindo
  matriz de interconexão heterogênea.

#### Prefetch explícito (NOVO — 50 referências onde tínhamos 0)
Mecanismo completo para empurrar dados ao device antes de qualquer
tarefa precisar:

```c
parsec_advise_data_on_device(data, device, PARSEC_DEV_DATA_ADVICE_PREFETCH);
```

Isso cria uma fake `parsec_gpu_task_t` com
`task_type = PARSEC_GPU_TASK_TYPE_PREFETCH` e joga na
`gpu_device->pending`. Ela atravessa o mesmo pipeline h2d/exec/d2h mas
sem corpo de kernel. Há task class dedicada
`parsec_device_data_prefetch_tc`, flow dedicado, profiling key dedicado
`"prefetch"` (`fill:#66ff66`), e função de release dedicada.

Também há `PARSEC_DEV_DATA_ADVICE_WARMUP` — warmup do kernel CUDA (carrega
modules, JIT). `PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE` — seta
`data->preferred_device` (novo campo em `parsec_data_t`) para influenciar
`parsec_select_best_device`.

#### Batching de tarefas GPU (NOVO)
Nova API:

```c
int parsec_gpu_task_collect_batch(stream, batch_head, callback, callback_data);
```

O callback decide se cada tarefa pendente pode ser anexada ao batch
atual na stream. Gated por:

- `PARSEC_DEV_CHORE_ALLOW_BATCH` flag no chore type
- `parsec_mca_device_type_supports_batch(device_type)`
- `parsec_device_enable_batching` global (default se construir com
  capability)

Permite ao Chameleon/CuBLAS fazer `cublasGemmBatched`/
`cublasGemmGroupedBatched` transparentemente.

#### Flags de profiling CUDA→GPU
`PARSEC_PROFILE_CUDA_TRACK_*` → `PARSEC_PROFILE_GPU_TRACK_*`. Adicionadas
duas:

- **`PARSEC_PROFILE_GPU_TRACK_MEM_USE = 0x0010`** — rastreia
  alocação/liberação/uso de memória device
- **`PARSEC_PROFILE_GPU_TRACK_PREFETCH = 0x0020`** — rastreia tarefas
  de prefetch

Novas chaves de profiling: `parsec_gpu_allocate_memory_key`,
`parsec_gpu_free_memory_key`, `parsec_gpu_use_memory_key_start/end`,
`parsec_gpu_prefetch_key_start/end`.

O `trackable_events` antes era global (`parsec_cuda_trackable_events`);
agora é **per-device** (campo em `parsec_device_gpu_module_t`).

#### Variável `parsec_device_skip_empty_events` (NOVO)
Default on — otimização que pula eventos de stream input/output
estritamente desnecessários, reduzindo overhead de profiling.

#### `parsec_gpu_d2h_max_flows` (NOVO)
Limita quantos flows uma tarefa W2R (eviction) agrega.

### 5.5 Camada de dados (`parsec/data.{h,c}`, `arena.*`)

#### `parsec_data_t`
Novos campos struct:

- **`parsec_atomic_lock_t lock;`** — cada data tem seu próprio atomic
  lock (granularidade fina, em vez do lock global)
- **`int8_t preferred_device;`** — hint via ADVICE_PREFERRED_DEVICE
- **`int32_t nb_copies;`**

Outros:

- `parsec_data_key_t` alargado de `uint32_t` → **`uint64_t`** (suporta
  taskpools enormes)
- `nb_elts` renomeado para **`span`** (semanticamente: span em bytes,
  não só elementos)
- `device_copies[1]` → `device_copies[]` flexible array member
- Removido `extern uint32_t parsec_supported_number_of_devices;`

#### `parsec_data_copy_t`
Novos campos:

- `parsec_datatype_t dtt;` — datatype aware
- `parsec_data_copy_alloc_cb *alloc_cb;` e
  `parsec_data_copy_release_cb *release_cb;` — hooks de alloc/release
  por cópia
- `PARSEC_DATA_FLAG_CPU_MIRROR_PROTECTED (1<<2)`,
  `PARSEC_DATA_FLAG_EVICTED (1<<5)`,
  `PARSEC_DATA_FLAG_PARSEC_MANAGED (1<<6)`,
  `PARSEC_DATA_FLAG_PARSEC_OWNED (1<<7)` — novo modelo de ownership
- Sentinel `PARSEC_DATA_CREATE_ON_DEMAND`

#### Novo modelo de ownership (doxygen novo em `data.h:35-55`)

- **PARSEC-owned e managed** — vêm de arenas, GPU copies; runtime sabe
  mover/alloc/free
- **Não-owned mas managed** — data collections; runtime duplica mas não
  libera `device_private`
- **Não-owned não-managed** — user-managed copies; usuário cuida do
  ciclo de vida

#### Novas APIs de data

- `parsec_data_start_transfer_ownership_to_copy(...)` /
  `_end_transfer_ownership_to_copy(...)` — lifecycle explícito
- `parsec_data_release_self_contained_data(data)` — libera mirror CPU
  self-contained para GPU-aware remote deps
- `parsec_data_protect_cpu_mirror(data)` — protege mirror CPU durante
  propagação GPU
- `parsec_data_copy_new(data, device, dtt, flags)` — adicionou args
  `dtt` e `flags`

#### Arena
`parsec_arena_s` agora herda de **`parsec_object_t`** (refcounted).
`parsec_arena_construct(arena, elem_size, alignment)` perdeu o arg
`opaque_dtt`. Nova API pública:

- `parsec_arena_datatype_new()`,
  `parsec_arena_datatype_set_type()`, `parsec_arena_datatype_release()`
  — anteriormente `parsec_dtd_arena_datatypes[]` era array global;
  agora objetos vivem dentro do context.

### 5.6 Comm engine / remote deps (novo grande subsistema)

Upstream tem **três** arquivos novos que não existem na nossa versão:

- **`parsec/parsec_comm_engine.{c,h}`** (58 + 207 linhas) — abstração
  sobre MPI
- **`parsec/parsec_mpi_funnelled.{c,h}`** (83 + 1617 linhas) —
  implementação funnelled-MPI active-message que roda no thread pool
  PaRSEC
- **`parsec/parsec_reshape.c`** (786 linhas) — redistribution via
  futures/promises (`parsec_base_future_t`, estados
  `PARSEC_UNFULFILLED_RESHAPE_PROMISE`/
  `PARSEC_FULFILLED_RESHAPE_PROMISE`)

Isso, combinado com o mirror-CPU-protected
`PARSEC_DATA_FLAG_CPU_MIRROR_PROTECTED`, permite **GPU-aware remote
deps** — uma tarefa remota pode receber dado que está numa GPU local
sem bounce eager para host.

### 5.7 Profiling / tracing

#### Novo: `mca/pins/alperf/`
Módulo PINS que inspeciona o profiling dictionary em cada evento e
avalia função associada à propriedade, escrevendo em área de shared
memory. Agregação algébrica online.

#### `parsec/dictionary.{c,h}` (361 + 947 linhas) — NOVO
Dictionary de profiling runtime; propriedades requisitadas via MCA
`profiling_properties=dgemm_NN_summa:GEMM:flops;...`.

#### `parsec/profiling_nvtx.{c,h}` (34 + 365 linhas) — NOVO
Adapter NVTX (NVIDIA Tools Extension) — exporta eventos PaRSEC para
`nsys`/`nsight-systems`. **Não existe na nossa versão.**

#### `parsec/termdet.c` (127 linhas) + `mca/termdet/` framework (NOVO)
Detecção de terminação plugável. Componentes: `fourcounter`, `local`,
`user_trigger`.

#### `parsec/compound.c` (135 linhas) — NOVO
Mecanismo de tarefas compostas (fusão de tasks).

#### `h5totrace` substitui `dbp2paje`
O upstream removeu `dbp2paje` (que era o que
`scripts/analysis/parsec_phase1.sh` usa para o caminho "oficial" de
phase-1). A troca é `h5totrace`. Isso afeta diretamente o pipeline de
análise StarVZ — o comentário no CLAUDE.md sobre "Official
(`scripts/analysis/trace_phase1.sh`): `dbp2paje` + `starvz -1 -t`"
deixa de funcionar sem mudança.

### 5.8 O que NÃO mudou (importante mencionar)

- **`PARSEC_GPU_MAX_STREAMS = 6`** e
  **`PARSEC_MAX_EVENTS_PER_STREAM = 4`** — unchanged.
- **Hierarquia h2d(0) / d2h(1) / cuda(2..5)** — mesma estrutura, ainda
  1 stream H2D + 1 stream D2H + 4 streams de exec.
- **Sem `cudaMemPool` / `cudaMallocAsync`** em nenhuma das duas árvores
  (0 referências nas duas). PaRSEC não usa CUDA stream-ordered memory
  pools.
- **Sem "AMM" / "lookahead"** em nenhuma versão.
- **Algoritmo base de `parsec_gpu_get_best_device`/
  `select_best_device`** — mesma lógica (afinidade de dado WRITE
  primeiro, depois balance por `device_load + ratio*sweight`), só perdeu
  o `ratio` externo e ganhou o `data->preferred_device` como input
  adicional.
- **lfq** e **gd** ainda funcionam igual como escalonadores CPU-only;
  não ganham vs. GPU novos papéis. A divisão "CPU scheduler (`lfq`/`gd`)
  não decide GPU" continua valendo.

### 5.9 Implicações práticas para o trabalho

| Você usa hoje | Upstream tem | Impacto |
|---|---|---|
| `mca_sched=lfq` (v3) | lfq ainda existe, mesma prioridade, mesmo algoritmo | code portable direto; ganho marginal (bugfix `pop_back`) |
| `mca_sched=gd` (v3) | gd idem | sem mudança funcional |
| `PARSEC_MCA_device_cuda_enabled=1` | ainda funciona | mas não há `parsec_gpu_init`/`fini`; lifecycle é MCA component |
| `dbp2paje` (seu `trace_phase1.sh`) | **removido** | teria que migrar para `h5totrace` |
| 1 GPU (RTX 4070) | D2D peer-to-peer é irrelevante (1 device só) | sem ganho de D2D |
| — | Prefetch mechanism | oportunidade: prefetchar tiles de coluna antes da POTRF chegar |
| — | Batching (`parsec_gpu_task_collect_batch`) | se o Chameleon adotar, `dgemm` GPU batchéavel |
| `CMAKE_CUDA_COMPILER=nvcc` | compute caps agora exigem ≥3.5 | RTX 4070 = sm 8.9 → OK |
| `trackable_events` global | vira per-device | código que lê `parsec_cuda_trackable_events` quebra |
| `papi_sde_register_fp_counter(handle, …)` | `parsec_papi_sde_*` (sem handle) | se você instrumenta PAPI direto, muda |
| `task->data[i].data_repo` | renomeado `source_repo_entry` | se seu código de análise tocar nesses campos |
| `parsec_data_key_t` 32-bit | 64-bit | geralmente compatível, mas pode expor bugs de cast em código antigo |

### 5.10 Resumo em uma frase

Upstream v4.0 refatorou **a camada de dispositivos inteira em framework
MCA plugin** (CUDA/HIP/Level-Zero/template) com vtable, adicionou
**prefetch explícito, batching, D2D peer-to-peer, NVLink mask**,
reformulou o **data layer** para ter locks e datatypes por data,
reescreveu o **comm engine**, e adicionou **NVTX tracing + dictionary de
profiling + termdet MCA + alperf PINS**. Os escalonadores (`lfq`, `gd`)
só tiveram mudanças cosméticas + 1 bugfix (`pop_back`) + 1 novo entry
(`llp`); **a decisão CPU-vs-GPU ainda é responsabilidade de
`parsec_select_best_device` (ex-`parsec_gpu_get_best_device`), e nem
`lfq` nem `gd` opinam nela, exatamente como na nossa versão**.

---

## Parte 6 — StarPU `dmda`/`dmdas`: o oposto, com performance model

StarPU faz o caminho inverso do PaRSEC: o **próprio escalonador** decide
qual worker (CPU ou GPU) pega a tarefa, usando previsões de duração e
custo de transferência. O arquivo
`starpu/src/sched_policies/deque_modeling_policy_data_aware.c` (1153
linhas) define uma família de 5 escalonadores derivados do mesmo código
base:

| Nome político | struct no arquivo | descrição | `push_task` | `pop_task` |
|---|---|---|---|---|
| `dm` | `_starpu_sched_dm_policy` (1047) | "performance model" | `dm_push_task` (sem data-aware) | `dmda_pop_task` |
| **`dmda`** | `_starpu_sched_dmda_policy` (1065) | **"data-aware performance model"** | `dmda_push_task` | `dmda_pop_task` |
| `dmdap` | `_starpu_sched_dmda_prio_policy` (1083) | data-aware + prioridade | `dmda_push_sorted_task` | `dmda_pop_task` |
| **`dmdas`** | `_starpu_sched_dmda_sorted_policy` (1101) | **"data-aware performance model (sorted)"** | `dmda_push_sorted_task` | `dmda_pop_ready_task` (só pega ready!) |
| `dmdasd` | `_starpu_sched_dmda_sorted_decision_policy` (1119) | sorted + decisão ordenada | `dmda_push_sorted_decision_task` | `dmda_pop_ready_task` |
| `dmdar` | `_starpu_sched_dmda_ready_policy` (1137) | data-aware ready | `dmda_push_task` | `dmda_pop_ready_task` |

### 6.1 Estado interno (`deque_modeling_policy_data_aware.h:25-39`)

Cada contexto de scheduling mantém:

```c
struct _starpu_dmda_data {
    double alpha;                                       // peso computação
    double beta;                                         // peso transferência
    double _gamma;                                       // peso energia
    double idle_power;                                   // idle consumption
    struct starpu_st_fifo_taskq queue_array[STARPU_NMAXWORKERS];   // 1 FIFO por worker!
    long int total_task_cnt, ready_task_cnt, eager_task_cnt;       // contadores
    int num_priorities;
};
```

Cada worker tem sua própria fila Fifo (`queue_array[w]`). Knobs
`alpha/beta/gamma` têm defaults `1.0 / 1.0 / 1000.0` (linha 156-158) e
são ajustáveis via env vars `STARPU_SCHED_ALPHA/BETA/GAMMA`. Há um
**`starpu_perf_knob_group`** que permite tunar esses knobs em runtime,
**por scheduler**.

### 6.2 O algoritmo (`_dmda_push_task`, linha 604-751) — onde o performance model entra

Quando uma tarefa fica pronta, o StarPU chama `dmda_push_task` →
`_dmda_push_task`:

**Passo 1: `compute_all_performance_predictions` (linha 411-602)** —
para **cada worker** e **cada implementação**
(`STARPU_MAXIMPLEMENTATIONS` do codelet), computar:

- `local_task_length[w][nimpl]` via
  `starpu_task_worker_expected_length(task, workerid, sched_ctx_id,
  nimpl)` (linha 507) — **consulta o performance model aprendido do
  codelet** para aquele (arquitetura, implementação). Se não tem modelo
  ou está calibrando, retorna `NaN`.
- `local_data_penalty[w][nimpl]` via
  `starpu_task_expected_data_transfer_time_for(task, workerid)` (linha
  509) — **estima o tempo de H2D/D2H** baseado em onde os handles estão
  agora e para onde irão.
- `local_energy[w][nimpl]` via
  `starpu_task_worker_expected_energy(...)` (511) — energy model.
- `exp_end[w][nimpl] = exp_start + prev_exp_len` (489) — tempo previsto
  de término neste worker (fila + transferência + duração).
- `ntasks_end = fifo_ntasks / relative_speedup(perf_arch)` (516) —
  heurística de greedy fallback para quando não há modelo.

Detecta `calibrating`/`unknown` (NaN) e resolve para o `ntasks_best`
(linha 531-551) — se a tarefa não tem modelo ou está em calibração,
**usa greedy** baseado em contagem de tarefas × speedup relativo, e
seleciona workers com NaN para calibrar lá.

**Passo 2: calcular `fitness` (linha 654-710)** — só se não houve greedy
forçado. A função fitness da Linha 680:

```c
fitness[w][nimpl] = dt->alpha * (exp_end[w][nimpl] - min_exp_end_of_task)     // Δ finish time
                  + dt->beta  * (local_data_penalty[w][nimpl])                   // custo transfer
                  + dt->_gamma * (local_energy[w][nimpl]);                       // energia
```

Mais (linha 691): se este worker ficaria mais ocupado que
`max_exp_end_of_workers`, adiciona penalidade de **idle power** dos
outros cores:

```c
fitness += dt->_gamma * dt->idle_power * (exp_end[w] - max_exp_end_of_workers) / 1e6
```

Escolhe o `(w, nimpl)` com **menor fitness** (linha 696-706).

**Passo 3: `push_task_on_best_worker` (linha 271-408)** — coloca a task
na Fifo do worker escolhido, atualiza o estado de overlapping de
pipeline da Fifo, e dispara
**`starpu_prefetch_task_input_for(task, best_workerid)`** (linha 355) —
ou seja, **o escalonador dispara explicitamente o prefetch dos dados
para o dispositivo escolhido**. Esse é o equivalente direto do "advice
prefetch" que o PaRSEC v4.0 acabou de adicionar — mas que StarPU tem há
muito tempo, e que **o dmda/dmdas usa ativamente**.

### 6.3 Diferenças específicas `dmda` vs. `dmdas`

O **único diff estrutural** entre os dois:

| Comportamento | `dmda` | `dmdas` |
|---|---|---|
| Ordem de Push | FIFO (insere no **fim** da fila do worker, linha 395) | **Prioridade** (insere ordenado na fila do worker, linha 377) |
| Ordem de Pop | `dmda_pop_task` = pop_first_ready_task ou pop_local_task | `dmda_pop_ready_task` = **só pop_first_ready_task** (linha 261) |
| `simule_push` | `_dmda_push_task(prio=0, da=1, simulate=0, sorted=0)` | `_dmda_push_task(prio=1, da=1, simulate=0, sorted=0)` |
| Indexação de prioridades | Não usa `num_priorities` efetivamente | **Usa** `exp_len_per_priority[]` e `ntasks_per_priority[]` |

Em outras palavras:

- **`dmda`** = data-aware modeling, **sem prioridade entre tarefas**.
  Tarefas entram numa FIFO pelo worker escolhido pelo model. Bom para um
  workload homogêneo (só GEMM, mesmo tamanho).
- **`dmdas`** = mesma modelagem + ordenação por prioridade (tanto no
  push quanto no pop — só retira tarefas marcadas como ready). Adequado
  para **workloads com dependências que têm caminho crítico** (Cholesky,
  QR): prioridades das tarefas críticas (POTRF, TRSM) são explicitamente
  honradas — o Chameleon já seta essas prioridades; `dmdas` respeita,
  `dmda` não.

### 6.4 Performance model: como StarPU consegue prever duração

`starpu_task_worker_expected_length` é a chave. StarPU mantém, por
codelet, por arquétipo de worker (`CPU/CUDA` + compute capability) e por
implementação, uma tabela de durações observadas (`struct
starpu_perfmodel` no `core/perfmodel/`). Cada codelet do Chameleon
(`cl_zgemm`, `cl_zpotrf`, etc.) declara um `model` que referencia um
símbolo (string). À medida que tarefas executam, o StarPU atualiza as
estatísticas (média, sigma, nsample). Quando se faz
`STARPU_CALIBRATE=1`, ele recalibra do zero; com `STARPU_CALIBRATE=0`, ele
lê de `$HOME/.starpu/sampling/`.

Essa é a característica fundamental que separa StarPU do PaRSEC:

- StarPU: codelet carrega um perfmodel declarado como string na codelet;
  scheduler consulta em runtime.
- PaRSEC: a task_class tem um escalar `load` (no caso GPU é
  `gpu_task->load`), declarado pelo desenvolvedor ou pelo Chameleon ao
  registrar a codelet — sem histórico em runtime.

### 6.5 Comparação lado-a-lado do fluxo de decisão

Para uma tarefa de Chameleon (`dgemm`) que pode rodar em CPU ou GPU:

```
=== PaRSEC (lfq/gd + dev_cuda) ===
1. DAG libera tarefa
2. se a task_class é CUDA-capable:
     engine chama parsec_gpu_get_best_device(task, ratio)          [device layer]
     dev_index = argmin(load[dev] + ratio*sweight[dev])
     engine chama parsec_gpu_kernel_scheduler(es, gpu_task, dev_index)
3. se não é CUDA-capable:
     engine chama sched_X_schedule(es, task, distance)             [scheduler]
     alguma CPU peça task via sched_X_select                       [scheduler]

=== StarPU (dmdas) ===
1. DAG libera tarefa
2. engine chama _starpu_push_task → dmda_push_sorted_task           [scheduler!]
3. scheduler:
     para cada worker w (CPU E CUDA), cada implementação nimpl:
         exp_end[w][nimpl] = exp_start[w] + fifo[w].exp_len +
                            + data_transfer_penalty(task, w)         [perf model]
     fitness = alpha*Δfinish + beta*transfer + gamma*energy (+idle_power)
     best = argmin fitness[w][nimpl]
4. push_task_on_best_worker(task, best, ...)
     → starpu_prefetch_task_input_for(task, best)                   [prefetch]
```

Note o quanto o **dmda/dmdas unifica CPU e GPU numa mesma filosofia** —
ambos são workers, ambos entram na mesma matriz de fitness. O
`lfq`/`gd` do PaRSEC **só opera CPU**; quando a tarefa é GPU-able, ele
é completamente bypassado.

---


## Parte 8 — Referências de arquivo

### Escalonadores (PaRSEC v3)

- **Enum dos escalonadores**: `parsec/dplasma/testing/common.c:31-41`
  (lfq=1, gd=5).
- **gd**: `parsec/parsec/mca/sched/gd/sched_gd_module.c` (init :77,
  schedule :143, select :125). Cabeçalho: `sched_gd.h:14-16`.
- **lfq**: `parsec/parsec/mca/sched/lfq/sched_lfq_module.c` (init :64
  com hierarquia hwloc, schedule :203, select :171). Estrutura de dados
  em `sched_local_queues_utils.h:24-34`.

### Device layer (PaRSEC v3)

- **Decisão GPU**: `parsec/parsec/devices/cuda/dev_cuda.c:1323`
  (`parsec_gpu_get_best_device`).
- **Streams/init**: `dev_cuda.c:504-557` (6 streams: h2d/d2h/cuda×4 com
  eventos), `dev_cuda.h:26-27`.
- **Pipeline GPU**: `dev_cuda.c:1967` (`parsec_gpu_kernel_scheduler`);
  fases H2D/exec/D2H em linhas 2015, 2044, 2076.
- **H2D**: `dev_cuda.c:1065` (`parsec_gpu_data_stage_in`);
  `cudaMemcpyAsync` em :1139.
- **D2H**: `dev_cuda.c:1681` (`parsec_gpu_kernel_pop`);
  `cudaMemcpyAsync` em :1700 e :1780.
- **Evicção (W2R)**: `dev_cuda.c:2031` (`parsec_gpu_create_W2R_task`);
  marca `GPU_TASK_TYPE_D2HTRANSFER = 111` em `dev_cuda.h:45`.
- **Peso estático do device (GPU)**: `dev_cuda.c:575-578`
  (`device_sweight` para CPU vs GPU).
- **Peso estático do device (CPU)**: `parsec/parsec/devices/device.c:489`
  (`nstreams × fp_ipc × freq`).
- **Inversão de pesos**: `parsec/parsec/devices/device.c:314-363`
  (`parsec_devices_freeze`, em particular linhas 352-355).

### StarPU dmda/dmdas

- **Implementação**:
  `starpu/src/sched_policies/deque_modeling_policy_data_aware.c` (1153
  linhas).
- **Estado interno**:
  `starpu/src/sched_policies/deque_modeling_policy_data_aware.h:25-39`.
- **Dispatch de prefetch**: linha 354-355
  (`starpu_prefetch_task_input_for(task, best_workerid)`).
- **Streams CUDA do StarPU**:
  `starpu/src/drivers/cuda/driver_cuda.c:104-105`
  (`in_transfer_streams[STARPU_MAXCUDADEVS]`,
  `out_transfer_streams[STARPU_MAXCUDADEVS]`).

### Análise de traces (no repositório)

- **Reconstrução dos traces PaRSEC**:
  `scripts/analysis/parsec_tasks_to_parquet.r:53-66` (categorias),
  :161-171 (CPU<N> vs CUDA<dev>_<stream>), :195-226 (de-duplicação do
  perfil duplo de kernels GPU).
- **Comentário do plano**: `docs/plano_experimentos_gpu.md:17` (métrica
  de thrashing H2D/D2H).

### Scripts do runner

- `scripts/run.sh:77-84` (parsec arm) e `scripts/gpu_doe_sweep.sh:159-185`
  (GPU parsec arm, com `device_cuda_enabled` e
  `mca_pins=task_profiler`).
