# Escalonadores `gd` e `lfq` do PaRSEC — como funcionam na prática

Documento de referência para o trabalho de Introdução ao HPC. Descreve
os dois escalonadores do PaRSEC usados nos experimentos (`lfq` e `gd`),
com foco em (1) a estrutura de dados e o алгоритmo de seleção/inserção
de cada um, e (2) o caminho completo de uma tarefa quando há uma GPU
envolvida — incluindo a decisão CPU-vs-GPU, o pipeline de execução no
device e as transferências H2D/D2H.

## O ponto-chave que muda tudo

Antes de mais nada: **`gd` e `lfq` são escalonadores de CPU**. Eles
decidem *qual thread de CPU (execution stream) pega qual tarefa*. **Eles
NÃO decidem se uma tarefa vai para GPU** — isso é responsabilidade de
uma camada totalmente separada, o **device layer CUDA**
(`parsec/parsec/devices/cuda/dev_cuda.c`), invocada por
`parsec_gpu_get_best_device()`. Está documentado nos comentários dos
próprios scripts do repositório:

- `scripts/doe/doe_gpu_tile.r:30` — *"a escolha CPU-vs-GPU do PaRSEC está
  no device layer (`parsec_gpu_get_best_device`), igual p/ lfq/gd."*
- `scripts/gpu_sweep.sh:40-43` — *"PaRSEC has no perf-model scheduler;
  lfq is its standard ready-task queue and GPU placement is decided by
  PaRSEC's heterogeneous device heuristic regardless."*

Para uma tarefa do Chameleon, o ciclo de vida é:

```
           [escalona quem pega a tarefa]      [decide GPU vs CPU]       [executa]
DAG pronto ────────────────────────────► delegated task ──────────► device layer ─► kernel
                  lfq ou gd                                    parsec_gpu_get_best_device
```

O `lfq`/`gd` só atua quando a tarefa vai rodar em **CPU**. Para tarefas
marcadas como delegáveis a uma GPU, o runtime pula esses escalonadores e
chama o `parsec_gpu_kernel_scheduler` direto.

## 1. `gd` — Global Dequeue (o "simples")

Arquivo: `parsec/parsec/mca/sched/gd/sched_gd_module.c` (191 linhas).
Cabeçalho (`sched_gd.h:14-16`) diz literalmente *"Global Dequeue
Scheduler"*.

### Estrutura de dados

Existe **uma única fila global**, compartilhada por **todos** os
threads de todos os VPs (virtual processes). É um `parsec_dequeue_t` —
uma deque *lock-free* (sem mutex) que pode ser operada pelas duas
pontas.

A inicialização em `flow_gd_init` (linhas 77-123) é esclarecedora:

- O **primeiro** execution stream de cada VP **cria** a deque
  (`OBJ_NEW(parsec_dequeue_t)`, linhas 90/92).
- Todos os **outros** streams, após uma barreira, simplesmente recebem a
  **mesma referência** (linhas 98-106):
  `es->scheduler_object = LOCAL_SCHED_OBJECT(vp->execution_streams[0])`.
- Ou seja: **1 fila, N consumidores**.

### Inserção: `sched_gd_schedule` (linhas 143-166)

Quando uma tarefa fica pronta, o runtime chama
`sched_gd_schedule(es, task, distance)`:

```c
if( (new_context->task_class->flags & PARSEC_HIGH_PRIORITY_TASK) &&
    (0 == distance) ) {
    parsec_dequeue_chain_front( dq, (parsec_list_item_t*)new_context);
} else {
    parsec_dequeue_chain_back( dq, (parsec_list_item_t*)new_context);
}
```

- Pode inserir uma **cadeia inteira** de tarefas de uma vez (`chain` =
  encadeia lista ligada).
- Alta prioridade + distância 0 → insere pela frente (será consumida
  primeiro).
- Caso contrário → fila FIFO normal (inserção no fim).

### Seleção: `sched_gd_select` (linhas 125-141)

```c
parsec_task_t *context =
    (parsec_task_t*)parsec_dequeue_try_pop_front( sd );
*distance = 0;
return context;
```

Simplesmente tenta **popar da frente** da deque global. `try_pop_front`
significa "não bloqueia; retorna NULL se vazia". Se voltar NULL, aquele
thread fica ocioso.

### Características práticas do `gd`

- **Sem noção de localidade**: nada sabe se a tarefa roda melhor no
  mesmo NUMA/socket que produziu seus dados. É só uma fila FIFO global.
- **Sem prioridade entre tarefas comuns**: só o flag
  `PARSEC_HIGH_PRIORITY_TASK` vai para a frente; o resto é estritamente
  FIFO.
- **Sem hierarquia de filas**: uma única estrutura de dados. Menor
  overhead por tarefa (não há deques locais, nem hbbuffer, nem hwloc),
  **mas** a deque é o ponto quente de contenção quando há muitos
  threads (24 no nosso caso) competindo.
- **Balanceamento de carga perfeito**: como todos puxam da mesma fila,
  nenhum thread fica ocioso enquanto há trabalho.

> Em uma frase: o `gd` é uma **fila FIFO global lock-free**. Cada worker
> pega a próxima tarefa disponível, ponto final. É o escalonador mais
> barato e mais justo, mas ignora localidade de memória.

## 2. `lfq` — Locality First Queue (o "esperto")

Arquivo: `parsec/parsec/mca/sched/lfq/sched_lfq_module.c` (248 linhas).
O cabeçalho chama de *"Local Flat Queue Scheduler"*.

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

- `parsec_hbbuffer_t` = "**hierarchical bounded buffer**". É um buffer
  de cada worker com tamanho fixo (`queue_size = vp->nb_cores * 4`,
  linha 82). Quando enche, transborda (*overflow*) para a fila de
  sistema.
- A fila de sistema é **uma deque global por VP** (criada só pelo
  thread 0, linhas 76-77, e compartilhada, linha 89) — análoga à deque
  do `gd`, mas usada só como "ralo" quando o buffer local enche.

### Inicialização: `flow_lfq_init` (linhas 64-169)

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

### Seleção: `sched_lfq_select` (linhas 171-201) — o coração do "Locality First"

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

1. **Próprio buffer local** (L1 da hierarchy) → máxima localidade: o
   dado provavelmente foi produzido neste mesmo núcleo.
2. **Filas dos vizinhos, em ordem crescente de distância hwloc** →
   rouba trabalho (*work stealing*) do vizinho mais barato primeiro.
3. **Fila de sistema** → quando tudo mais está vazio, busca do ralo
   global compartilhado.

Note o `*distance` retornado: indica "a que distância" a tarefa foi
encontrada (0 = local, crescente = cada vez mais longe). Essa
distância é passada de volta para `sched_lfq_schedule` quando a tarefa
gera sucessoras, para que elas **sejam inseridas na mesma
vizinhança** (afinidade de produtor-consumidor).

### Inserção: `sched_lfq_schedule` (linhas 203-211)

```c
parsec_hbbuffer_push_all(task_queue, new_context, distance);
```

Insere no **próprio buffer local**, passando a `distance` herdada. Se o
buffer encher, o callback
`parsec_mca_sched_push_in_system_queue_wrapper` (em
`sched_local_queues_utils.h:72-82`) faz o overflow cair na **deque de
sistema**.

### Comparação com `pop_best` e prioridade

Tanto no `lfq` (e também nos escalonadores `ltq`, `lhq`, `ap` que
compartilham esse utilitário) a seleção usa
`parsec_execution_context_priority_comparator` — uma função que ordena
tarefas por **prioridade** declarada da task class. Então dentro de
qualquer fila do lfq, quando há múltiplas tarefas, a de maior prioridade
é tirada primeiro (não é FIFO puro dentro do buffer). Mas o `gd` também
insere por prioridade; a diferença é estrutural.

> Em uma frase: o `lfq` é um **work-stealing hierárquico com consciência
> de localidade** — cada worker tem seu buffer local; ao esvaziar, ele
> rouba do vizinho mais próximo (topologicamente) primeiro, e só vai ao
> ralo global em último caso.

## 3. `gd` vs `lfq` lado a lado

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

Para a nossa workload (Cholesky/QR Chameleon, 24 threads, tiles de 480),
o `lfq` tende a vencer o `gd` porque:

- A fatoração tem **forte reuso de tiles** (cada tile é lido/escrito por
  várias tarefas consecutivas);
- Mantê-las no cache do mesmo núcleo (produtor → consumidor no mesmo
  buffer) reduz tráfego de memória;
- 24 threads competindo por uma única deque (`gd`) gera contenção
  mensurável.

## 4. Como fica o caminho quando há uma GPU (RTX 4070)

Aqui mora a parte sutil do experimento. Quando
`PARSEC_MCA_device_cuda_enabled=1` (configurado em
`gpu_doe_sweep.sh:161`), o registro de dispositivos em
`parsec_gpu_init` (`dev_cuda.c:504-609`) cria **uma `gpu_device_t` por
GPU física**, e a registra via `parsec_devices_add` (linha 609). Cada
dispositivo é identificado por um `device_index`.

### 4.1. A decisão CPU-vs-GPU: `parsec_gpu_get_best_device` (`dev_cuda.c:1323-1362`)

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

1. **Afinidade de dado primeiro**: se a tarefa grava em um tile que
   **já foi alocado numa GPU**, ela vai para essa mesma GPU (mantém o
   dado local, evita D2H+H2D).
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

### 4.2. As 6 streams por GPU e o que elas fazem

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
cada operação — nunca blocking (`cudaEventQuery`). É por isso que na
análise (`parsec_tasks_to_parquet.r:195-208`) foi achado o artefato do
"submit duplo": o thread CPU dispara (`prof_event_key_start` no submit)
e a CUDA stream dispara de novo quando o kernel realmente executa — a
solução é guardar só a span da CUDA stream para não inflar em ~16ms.

O `exec_stream` usado para executar é escolhido em **round-robin**
(`dev_cuda.c:2038`):

```c
exec_stream = (exec_stream + 1) % (gpu_device->max_exec_streams - 2);  // max_streams=6 - 2 = 4 streams de execução
```

ou seja, dentre os índices 2..5.

### 4.3. O pipeline em `parsec_gpu_kernel_scheduler` (`dev_cuda.c:1967-2124`)

Esse é o "escalona-escalonador da GPU". Pseudocódigo simplificado:

```
1. Toma o mutex da GPU (atomic). Se != 0, outro thread está mexendo:
   empurra a tarefa em gpu_device->pending (FIFO) e sai (RETURN_ASYNC).   ← lugar onde lfq/gd somem
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
três fases (movein → exec → moveout), e enquanto isso coloca as próximas
tarefas na **`gpu_device->pending`** (uma FIFO, linhas 1993 e 2101).
Ou seja: **a partir do momento em que a tarefa é delegada à GPU, o
escalonador lfq/gd perdeu a jurisdição** — a partir daí é a
`gpu_device->pending` que dita a ordem, e a ordem é **FIFO** (com lapso
de `PARSEC_PUSH_TASK` ordenado por prioridade só na
`exec_stream[0].fifo_pending`, linhas 1427/1366-1371).

### 4.4. A transferência host ↔ device em detalhe

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
  `parsec_gpu_callback_complete_push`, linhas 1377-1426, quando o
  `cudaEventQuery` do h2d-stream retornar success).
- Contabiliza em `gpu_device->super.transferred_data_in += nb_elts`
  (linha 1146) — é o número que aparece em
  `PARSEC_MCA_device_show_statistics=1`.

**D2H — `parsec_gpu_kernel_pop` (`dev_cuda.c:1681-1808`)**:

- Só faz D2H para flows que têm `FLOW_ACCESS_WRITE` **e** o bit
  `gpu_task->pushout & (1<<i)` está setado (linha 1765). Ou seja: o
  runtime tem **controle granular** sobre quais resultados voltam — só
  escreve de volta o que uma tarefa consumidora futura realmente
  precisará na CPU.
- Para dados só-LEITURA, não há cópia de volta: a cópia device vai para
  `gpu_mem_lru` (linha 1744) como **cache**, e fica disponível para a
  próxima tarefa que a ler (afinidade de device — é o cachê LRU
  `gpu_device->gpu_mem_lru`).
- Usa o stream `d2h(1)` (`gpu_stream->cuda_stream` passado é o stream
  1, linhagem 2076-2078).
- `transferred_data_out += nb_elts` (linha 1790).

**A tarefa "W2R" (Write-to-Read)**: quando não há memória livre no
device (`parsec_gpu_data_reserve_device_space` retorna negativo, linhas
1634-1636), o scheduler cria uma tarefa sintética
(`parsec_gpu_create_W2R_task`, linha 2031) do tipo
`GPU_TASK_TYPE_D2HTRANSFER` (constante 111, `dev_cuda.h:45`). Essa
tarefa **força** um D2H de algum bloco sujo da `gpu_mem_owned_lru` para
liberar espaço no device — é a **evicção controlada**. Isso explica o
"thrashing" que `docs/plano_experimentos_gpu.md:17` cita como métrica
causal nos experimentos (*"tráfego H2D/D2H na GPU - thrashing do
pushout do DTD do PaRSEC"*). O `lfq`/`gd` não conseguem reduzir esse
thrashing — só uma política melhor de pushout poderia.

### 4.5. Síntese do caminho completo de UMA tarefa GPU

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

## 5. Implicações para os experimentos

1. **Trocar `lfq` ↔ `gd` muda só a parte CPU.** Para uma fatoração
   Cholesky/QR com a maioria dos gemms na GPU (no nosso caso single GPU
   poti), o que diferencia os escalonadores é **como as tarefas de CPU
   restantes (POTRF, TRSM, pouco paralelas, servidoras de dependência)
   são roteadas**.
2. **A GPU não tem nem lfq nem gd**: é sempre a mesma lógica
   (`get_best_device` + kernel_scheduler). Por isso
   `doe_gpu_tile.r:30` tem razão em dizer que a escolha de GPU é "igual
   p/ lfq/gd".
3. **H2D só acontece quando o `version` difere**. Repetições executam a
   mesma fatoração sobre o mesmo tile com mesma versão → segunda execução
   pode achar cache quente. Isso afeta a repetição média — cuidado ao
   interpretar `gflops` de réplicas sobrepostas na mesma GPU com cache.
4. **O `pushout` seletivo** é o que separa o "modelo DTD do PaRSEC" do
   StarPU `dmda`: no PaRSEC, o runtime decide em tempo de taskpool setup
   quais flows são pushed-out (bit `pushout` é estabelecido quando a
   `parsec_gpu_task_t` é criada pelo Chameleon via `parsec_devices_*`
   registration). O StarPU/dmda decide isso com performance model; o
   PaRSEC decide estaticamente. É daí que vem o "thrashing" que é medido
   com `transfers_summary.csv`.
5. **A stream `h2d(0)` é gargalo único** para todo H2D do device —
   todas as tarefas competem por ela. Igualmente `d2h(1)` para D2H. Só
   os kernels têm 4 streams paralelas. Por isso overlap H2D/D2H é
   natural em GPU sparse. Esse desenho está em `dev_cuda.h:186-191`
   (comentário original: *"4 streams: one for transfers ... to the GPU,
   2 for kernel executions and one for transfers ... into the main
   memory"*).
6. **`PARSEC_MAX_EVENTS_PER_STREAM = 4`** significa que a GPU pode ter
   **até 4 tarefas em voo por stream** sem bloquear o host. Acima disso,
   o thread que está "possuindo" a GPU dorme esperando
   `cudaEventQuery`. Esse é o histórico do artefato de ~16ms no
   `tasks.parquet`.

## 6. Referências de arquivo

- **Enum dos escalonadores**: `parsec/dplasma/testing/common.c:31-41`
  (lfq=1, gd=5).
- **gd**: implementação integral em
  `parsec/parsec/mca/sched/gd/sched_gd_module.c` (init :77, schedule
  :143, select :125). Estrutura: `sched_gd.h:14-16`.
- **lfq**: `parsec/parsec/mca/sched/lfq/sched_lfq_module.c` (init :64
  com hierarquia hwloc, schedule :203, select :171). Estrutura de dados
  em `sched_local_queues_utils.h:24-34`.
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
- **Peso estático do device**: `dev_cuda.c:575-578` (`device_sweight`
  para CPU vs GPU).
- **Reconstrução dos traces**: `scripts/analysis/parsec_tasks_to_parquet.r:53-66`
  (categorias), :161-171 (CPU<N> vs CUDA<dev>_<stream>), :195-226
  (de-duplicação do perfil duplo de kernels GPU).
- **Comentário do plano**: `docs/plano_experimentos_gpu.md:17` (métrica
  de thrashing H2D/D2H).
- **No código do runner**: `scripts/run.sh:77-84` (parsec arm) e
  `scripts/gpu_doe_sweep.sh:159-185` (GPU parsec arm, com
  `device_cuda_enabled` e `mca_pins=task_profiler`).
