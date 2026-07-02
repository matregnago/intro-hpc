# Por que o PaRSEC varia entre repetições na GPU (e o StarPU não)

Documento de referência para o trabalho de Introdução ao HPC. Explica por
que, na figura *"Chameleon GPU - GFLOPS vs N"* (GPU(1), b=500, t=19), o
StarPU (`dmda`/`dmdas`) é muito consistente entre repetições enquanto o
PaRSEC (`gd`/`lfq`) varia bastante — sobretudo o `lfq` em N grande. A
conclusão é que **a variância é inerente ao design do PaRSEC**, não da nossa
implementação GPU/DTD no Chameleon. Ver também
[`escalonadores_parsec.md`](escalonadores_parsec.md),
[`analise_escalonadores_gpu.md`](analise_escalonadores_gpu.md) e
[`cuda_parsec.md`](cuda_parsec.md).

## 1. A observação

Dados de `data/gpu_n_size_poti_20260701_093941/results.csv` (5 repetições
por célula). Coeficiente de variação (DP/média) do GFLOPS, agregado por
config sobre todas as células algo×N:

| Config        | CV médio | CV máx |
|---------------|---------:|-------:|
| `parsec/lfq`  |   ~7,7 % | 13,8 % |
| `parsec/gd`   |   ~2,6 % |  5,7 % |
| `starpu/dmda` |   ~1,1 % |  2,7 % |
| `starpu/dmdas`|   ~0,6 % |  2,4 % |

O `lfq` tem ~7× a variância run-to-run das configs StarPU; o `gd` fica no
meio. O padrão não é ruído gaussiano: em N grande é **bimodal**. Exemplo,
GETRF `lfq` N=50000:

```
rep1: 724.8 GF (114.97 s)   \
rep2: 744.1 GF (111.99 s)    > modo "lento"  ~107-115 s
rep4: 778.0 GF (107.11 s)   /
rep3: 976.4 GF ( 85.35 s)   \  modo "rápido" ~85-90 s
rep5: 929.4 GF ( 89.66 s)   /  (≈ ótimo do StarPU: dmda ~84-85 s aqui)
```

Ou seja: às vezes o PaRSEC acha um bom escalonamento (empata com o StarPU),
às vezes cai ~25-30 % mais lento. Essa assinatura de **dois estados
discretos** é a impressão digital de "device escolhido por corrida de
timing, sem modelo de desempenho para corrigir".

## 2. Causa raiz: o PaRSEC não tem perf-model como o `dmda`

Essa é *a* diferença central. O StarPU `dmda` ("deque model data aware")
mantém um **modelo de desempenho medido** (histórico calibrado) por
codelet/device, estima o tempo de término de cada tarefa em cada worker
(incluindo custo de transferência) e atribui a tarefa ao device que
**minimiza o término previsto**. É um laço de realimentação que se
autocorrige.

O PaRSEC, no caminho GPU que usamos (DTD), **não tem modelo de desempenho
por tarefa algum**. O único "custo" é `ratio * sweight`, onde:

- `ratio` é constante — `1.0` fixo no DTD
  (`parsec/parsec/interfaces/superscalar/insert_function.c:1266`);
- `sweight` é o **pico teórico de flops** do device, calculado uma vez a
  partir de specs do hardware (`parsec/parsec/devices/cuda/dev_cuda.c:575-578`,
  normalizado em `parsec/parsec/devices/device.c:322-359`). Nunca é atualizado
  por medição.

Sem duração medida e sem estimativa de término, a política é determinística
mas cega — não há como aprender que uma dada tarefa terminaria antes na GPU.
São três mecanismos concretos que transformam essa cegueira em variância:

### 2a. Colocação dominada por pinning de residência + corrida first-come

`parsec_gpu_get_best_device` (`dev_cuda.c:1323-1362`) tem dois ramos:

- **Branch A (pinning):** varre as flows de escrita e, se o tile de saída já
  é residente numa GPU (`data_in->original->owner_device > 1`), **fixa a
  tarefa nessa GPU e retorna imediatamente**, sem comparar carga.
- **Branch B (greedy):** só quando o dado ainda não está em GPU, escolhe o
  device que minimiza `parsec_device_load[dev] + ratio*sweight[dev]`.

Numa fatoração, depois que os primeiros tiles ficam residentes, **quase toda
tarefa cai no Branch A** — a colocação segue "onde o tile caiu primeiro". E
"onde caiu primeiro" é decidido por elementos genuinamente não-determinísticos:

- **mutex de dono único da GPU, adquirido first-come** (`dev_cuda.c:1991`):
  o primeiro worker a pegar o mutex vira o gerente do device; os demais só
  empilham em `pending` e retornam. Qual thread vira gerente muda a cada run.
- **fila `pending` drenada em FIFO por ordem de submissão** — a ordenação por
  prioridade é opcional e **desligada por padrão**
  (`device_cuda_sort_pending_tasks=0`, `dev_cuda.c:408-410`).
- **streams de kernel em round-robin puro**, sem considerar ocupação
  (`dev_cuda.c:2038`).
- **ordem de chegada das tarefas ao hook CUDA** definida pelo escalonador de
  CPU (§3).

Como a colocação do grosso das tarefas depende de *qual tile ficou residente
primeiro* — uma corrida de timing — e não de um custo medido, pequenas
diferenças de ordenação entre execuções se propagam em splits CPU/GPU
diferentes e cargas por GPU diferentes.

### 2b. Ponto de virada CPU-vs-GPU sensível a corrida

O único termo dinâmico, `parsec_device_load[dev]`, é um contador grosseiro
de "trabalho atribuído em aberto": incrementa na submissão
(`insert_function.c:1310`) e **decrementa só quando a tarefa termina**
(`dev_cuda.c:2126`). Como `sweight` subestima a vantagem real da GPU, depois
que ~30 GEMMs estão em voo o `weight[GPU]` ultrapassa o `weight[CPU]` no
Branch B e **o próximo GEMM vai para a CPU — onde leva ~100× mais tempo**.
Se e quando esse limiar é cruzado depende das contagens transitórias em voo
no instante da decisão, então ~25 % dos GEMMs podem ser mal-roteados para a
CPU numa execução e menos em outra. É o gatilho direto do modo "lento".

## 3. Por que `lfq` varia mais que `gd`

`gd` e `lfq` são escalonadores **de CPU** — decidem qual thread pega qual
tarefa, e não influenciam a decisão GPU (ver `escalonadores_parsec.md`). Mas
eles definem a **ordem em que as tarefas chegam ao hook CUDA**, e essa ordem
alimenta a corrida de residência do §2a. A diferença:

- **`gd`** (`parsec/parsec/mca/sched/gd/sched_gd_module.c`): um único deque
  global compartilhado por todos os streams. Racy, mas com pouca liberdade de
  reordenação.
- **`lfq`** (`parsec/parsec/mca/sched/lfq/sched_lfq_module.c`): filas locais
  por thread + **work-stealing** de vizinhos. Qual worker ocioso rouba qual
  tarefa varia a cada execução por natureza — muito mais liberdade de
  reordenação, logo mais variação no tile que fica residente primeiro.

Isso casa com os dados: `lfq` (~7,7 %) > `gd` (~2,6 %) > `starpu` (~1 %).

## 4. Não é a nossa camada DTD/GPU

Nossos codelets em `chameleon/runtime/parsec/codelets/codelet_z*.c` são uma
camada fina: registram as variantes CPU+CUDA via
`parsec_dtd_taskpool_insert_task_cuda` e **deixam o PaRSEC decidir**. Eles
**não** setam afinidade de device, **não** setam prioridade específica de GPU
e **não** passam hints de carga — a colocação é 100 % do
`parsec_gpu_get_best_device` com `ratio=1.0` e o `sweight` estático do
próprio PaRSEC. Portanto a variância vem do PaRSEC, não do que implementamos.

A única escolha de política nossa é o **`pushout` incondicional** (D2H eager
em toda flow de escrita, `insert_function.c:1302-1310`): ela prejudica a
**média** (thrashing H2D/D2H, ver `cuda_parsec.md`) e, por forçar um
re-stage H2D em cada tile dependente, **amplifica a sensibilidade a jitter de
escalonamento** — mas não é a raiz da variância. É o ponto mais acionável se
um dia quisermos atacar a variância pelo nosso lado em vez do PaRSEC.

## 5. Consequência na figura

As barras de erro do `plot_n_size.r` eram **min/max** entre repetições, então
uma única rep azarada definia o whisker inteiro e exagerava visualmente o
PaRSEC (ex.: os 30 % de range em GETRF `lfq` N=50000 vinham de um par de reps
lentas). O script foi ajustado para **±1 desvio-padrão** em torno da média
(`scripts/analysis/plot_n_size.r`, com caption explícito). O PaRSEC continua
aparecendo mais variável — porque de fato é —, mas de forma justa, sem um
outlier isolado dominando a barra.

> Nota: para representar a *incerteza da média* (encolhe com o nº de reps)
> em vez da *dispersão run-to-run*, bastaria trocar `gflops_sd` por
> `gflops_sd/sqrt(n_rep)` (±SE). Aqui usamos ±DP de propósito, porque o
> assunto é justamente a variabilidade entre execuções.
