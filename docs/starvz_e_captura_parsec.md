# StarVZ como visualizador único: como capturar do PaRSEC os mesmos tibbles do StarPU

Documento de **design** (sem implementação). Objetivo: entender a fundo como o StarVZ
ingere traces (Phase 1), por que hoje o PaRSEC só vira uma visualização pobre, e projetar
um adaptador que transforme **todas** as informações do PaRSEC nos **mesmos tibbles**
(`application`, `dag`, `tasks`, `variable`, `link`, …) que o StarVZ gera para o StarPU —
de modo a usar o StarVZ como ferramenta **única** de visualização para os dois runtimes.

Traces de referência (job `data/chameleon-gpu-phases_795377/e4_traces/runs/`):

- PaRSEC: `0012_parsec_lfq_dpotrf_n19200_b960_rep1/` (contém `cham_dpotrf-0.prof-wA1QQ1`)
- StarPU: `0013_starpu_dmdas_dpotrf_n19200_b960_rep1/`

---

## 1. Situação atual

Carregando os dois traces com `starvz_read()` e contando linhas dos parquets:

| Tabela          | PaRSEC (0012) | StarPU (0013) | Observação |
|-----------------|:---:|:---:|---|
| `application`   | 3301 | 1540 | PaRSEC inclui 981 `movein` + 771 `moveout` + 9 `cuda` |
| `starpu`        | **0** | 3710 | overhead de runtime/scheduler só no StarPU |
| `variable`      | **0** | 12611 | GFlop/s, banda, fila ao longo do tempo só no StarPU |
| `other_state`   | 70 | 11890 | — |
| `dag`           | **ausente** | 7291 | dependências entre tarefas |
| `tasks`         | **ausente** | 3183 | metadados por tarefa |
| `link`          | **ausente** | 768 | transferências cronometradas (origem/dest/bytes) |
| `memory_state`  | **ausente** | 1960 | DriverCopy/Async/Alloc |
| `data_handles` / `task_handles` | **ausente** | 211 / 5863 | dado por tile / acesso tarefa↔dado |

Ou seja: o StarVZ enxerga o StarPU com riqueza total e o PaRSEC só como um space-time
incompleto. **Este documento mostra que isso NÃO é limitação do PaRSEC**, e sim do
caminho de coleta atual; e projeta como fechar a lacuna.

---

## 2. Como o StarVZ captura informação (Phase 1)

### 2.1 Visão geral do pipeline

A Phase 1 (`starvz/inst/tools/phase1-workflow.sh` + `starvz/R/phase1.R` +
`starvz/R/phase1_parse_csv.R`) faz, em resumo:

```
  Trace bruto                Intermediários (CSV/REC)              Tibbles (parquet)
  ───────────                ───────────────────────              ─────────────────
  FXT (StarPU)  ─fxt2paje─►  paje.sorted.trace.gz ─pj_dump─►  paje.csv.gz
        │                    tasks.rec / data.rec ─rec2csv─►  rec.tasks.csv.gz …
        │                    dag.dot                              │
        ▼                         │ (fatiamento por tipo de linha)│
  paje.trace (PaRSEC)             ▼                               ▼
                          paje.worker_state.csv.gz ──► Application + StarPU + Colors
                          paje.variable.csv.gz     ──► Variable
                          paje.link.csv.gz         ──► Link
                          paje.{memory,comm,other}_state.csv.gz ─► Memory/Comm/Other_state
                          paje.events.csv.gz       ──► Events / Events_data / Events_memory
                          dag.csv.gz               ──► Dag       (full_join com Application)
                          rec.tasks.csv.gz         ──► Tasks / Task_handles
                          rec.data_handles.csv.gz  ──► Data_handles
                          entities.csv             ──► Y (hierarquia) / Zero / Version
```

O coração é o `pj_dump`, que transforma o Paje num CSV onde **cada linha é classificada por
natureza** (`State`, `Variable`, `Link`, `Event`). O `phase1-workflow.sh` (linhas 200-229)
**fatia** esse CSV por natureza/contêiner e gera os intermediários acima; depois o
`phase1-workflow.R` chama as funções de parsing do `phase1_parse_csv.R` para produzir cada
tibble.

### 2.2 Os dois caminhos: StarPU vs PaRSEC

O mesmo script atende os dois runtimes, mas com **desvios importantes** quando
`STARVZ_USE_PAJE_TRACE=1` (caminho PaRSEC):

| Etapa | StarPU | PaRSEC (`STARVZ_USE_PAJE_TRACE=1`) |
|---|---|---|
| Conversão de entrada | `fxt2paje.sh` (FXT → paje + `tasks.rec`/`data.rec`) | usa `paje.trace` pronto (de `dbp2paje`), corrigido por `parsec_paje_fix.sh` |
| `rec2csv` (tasks/data) | **executa** | **PULADO** (linhas 124-126) → sem `Tasks`/`Data_handles` |
| DAG (`dag.dot` → `dag.csv`) | **executa** | **PULADO** (linhas 235-236) → sem `Dag` |
| Fatiamento do Paje | igual | igual |

`parsec_paje_fix.sh` apenas normaliza o Paje do PaRSEC (dedup de headers/definições, escala
de tempo µs→ms via `STARVZ_PARSEC_TIME_SCALE`). Não acrescenta informação.

### 2.3 Tibbles, suas fontes e a chave de ligação

Funções em `starvz/R/phase1_parse_csv.R`:

| Tibble | Função (arquivo:linha aprox.) | Fonte | Colunas-chave de saída |
|---|---|---|---|
| `Application` | `read_worker_csv` (4-312) | `paje.worker_state.csv.gz` | `ResourceId, Node, Resource, ResourceType, Value, Start, End, Duration, JobId, GFlop, Iteration, Outlier, …` |
| `StarPU` | idem (split `Application==FALSE`) | idem | estados de runtime (`Idle`, `Scheduling`, `Executing`, …) |
| `Variable` | `read_vars_set_new_zero` (519-588) | `paje.variable.csv.gz` | `ResourceId, Type, Start, End, Value(double), Node, Resource` |
| `Link` | `read_links` (1057-1112) | `paje.link.csv.gz` | `Container, Type, Start, End, Size, Origin, Dest, Key, …` |
| `Dag` | `read_dag` (984-1054) | `dag.csv.gz` | `JobId, Dependent, Start, End, Cost, Value` |
| `Tasks`/`Task_handles` | `tasks_csv_parser` (822-890) | `rec.tasks.csv.gz` | `JobId, StartTime, EndTime, WorkerId, SubmitOrder, Handles…` |
| `Data_handles` | `data_handles_csv_parser` (763-785) | `rec.data_handles.csv.gz` | `Handle, Name, …` |
| `Y` / `entities` | `hl_y_paje_tree` (phase1.R:260-310) | `entities.csv` | hierarquia + `Height`, `Position` |

**A ligação fundamental — `JobId`.** O `Application` recebe um `JobId` por instância de
tarefa (campo do estado Paje). O `Dag` (`read_dag`) faz literalmente
`full_join(Application, by = "JobId")` para herdar `Start/End/Value/ResourceId`, e o
`Tasks` referencia o **mesmo** `JobId`. Sem um `JobId` único e consistente por tarefa, não
há como reconstruir DAG nem casar metadados.

**Renomeações de `Variable$Type`** (linhas 581-586) que os painéis dependem:
`"Number of Ready Tasks"→"Ready"`, `"…Submitted Uncompleted…"→"Submitted"`,
`"Bandwidth In/Out (MB/s)"→"B. In/Out (MB/s)"`, além de `"GFlops"` e `"Used (MB)"`.
Os painéis filtram por esses nomes: `panel_gflops` → `Type=="GFlops"`; `panel_ready` →
`ResourceId~"sched" & Type~"Ready"`; `panel_gpubandwidth` → `ResourceId~"mpict" & Type~"Out"`;
`panel_usedmemory` → `Type~"Used"`.

---

## 3. Diagnóstico: por que o PaRSEC vira tibbles pobres

Três causas, todas confirmadas nos dados/código:

1. **Os desvios do caminho PaRSEC** (§2.2): `rec2csv` e a conversão do DAG são pulados →
   `Tasks`/`Data_handles`/`Dag` simplesmente não nascem.

2. **O `dbp2paje.c` é "lossy" por construção** (`parsec/tools/profiling/dbp2paje.c:779-860`).
   Ele cria a hierarquia de contêineres (`CT_Appli→CT_P→CT_VP/GPU→CT_T/Stream→CT_ST`) e
   emite **um único tipo de estado** `ST_TS` cujo *valor é o nome da classe da tarefa*
   (`K-<id>` → nome do dicionário: `movein`/`moveout`/`cuda`/`gemm`/`potrf`…). Ele **não
   emite**:
   - `JobId` por instância (o valor do estado é a *classe*, não a tarefa);
   - linhas `Variable` (nada de GFlop/s, banda, fila);
   - linhas `Link` de dependência (o único `LT_TL` é "split event" entre threads, não o DAG).

   Logo, por esse caminho, `variable`/`dag`/`tasks` são **impossíveis**, não importa o
   pós-processamento.

3. **Evidência direta nos parquets atuais** (medida em `0012`):
   - `application$JobId` é **100% `NA`** (1 valor distinto em 3301 linhas). No StarPU são
     1540 ids únicos. → sem chave para DAG/tasks.
   - `ResourceId` do PaRSEC: `M0G0S0-0…M0G0S5-x` (6 *lanes* de stream de GPU) e
     `M0V0T0…M0V0T15` (16 threads de CPU). Os streams de GPU **são** capturados, mas o
     `separate_res`/derivação de `ResourceType` do StarVZ produz lixo (`"MGS-"`, `"MVT-"`)
     porque o esquema de nomes do PaRSEC não bate com o padrão `node_RESOURCE_num` esperado.

---

## 4. O que o PaRSEC consegue capturar de verdade (do código-fonte)

Há **dois pipelines** nativos do PaRSEC, ambos ignorados pelo caminho Paje do StarVZ.

### 4.1 Profiling de eventos (`.prof`)

`parsec/parsec/profiling.{h,c}` (+ formato binário em `parsec_binary_profile.h`). Um
arquivo `.prof-<id>` por rank, com buffers de eventos por thread. Cada evento é um par
início/fim com `key` (tipo no dicionário), `taskpool_id`, `event_id`, `timestamp` e um
`info` opcional.

- **GPU/CUDA** (`parsec/parsec/devices/cuda/dev_cuda.c`): cria 3 tipos de stream por device
  (`GPU d-d`): `movein` (H2D), `moveout` (D2H), `cuda` (kernel). Controlados por flags
  `PARSEC_PROFILE_CUDA_TRACK_DATA_IN/OUT/EXEC`; **default = `EXEC | DATA_OUT`** — por isso a
  trace atual tem `moveout` denso (771 ≈ 1 D2H por kernel, o "eager pushout" de
  [`cuda_parsec.md`](cuda_parsec.md)) e `movein` parcial.
- **PINS** (`parsec/parsec/mca/pins/`): `task_profiler` (instrumenta `EXEC_BEGIN/END`,
  `RELEASE_DEPS_*`, `ACTIVATE_CB_*`, `DATA_FLUSH_*`), `papi` (contadores de hardware),
  `ptg_to_dtd`, etc.
- **MCA params**: `PARSEC_MCA_profile_filename` (liga o `.prof`), `PARSEC_MCA_pins`,
  `PARSEC_MCA_pins_papi_event`, `PARSEC_MCA_profile_rusage`.

### 4.2 DAG / dependências (`.dot`) — caminho separado do profiling

- Compilar o PaRSEC com `PARSEC_PROF_GRAPHER` e rodar com a opção `-. <nome>`
  (registrada em `parsec/parsec.c:354,515,695`) gera `<nome>-<rank>.dot`.
- `parsec/parsec_prof_grapher.c:99-180` emite cada **nó** com
  `tooltip="tpid=…:did=…:tname=…:tid=…"` e cada **aresta** com `flow_src=>flow_dst`.
- `parsec/tools/profiling/python/parsec_dag.py` carrega esses `.dot` num
  `networkx.DiGraph`.

### 4.3 Ferramentas de extração (`parsec/tools/profiling/`)

| Ferramenta | Tipo | Entrada → Saída | Útil para |
|---|---|---|---|
| `dbp2xml` | **C (já compilada)** | `.prof` → XML estruturado de eventos | eventos p/ `Application`/`Tasks`/`Variable` |
| `dbpinfos` | **C (já compilada)** | `.prof` → infos globais | `N, NB, sched, gflops, hostname` |
| `dbp2mem` | **C (já compilada)** | `.prof` → uso de memória | `memory_state` |
| `dbp2paje` | **C (já compilada)** | `.prof` → Paje (lossy, §3) | space-time básico (o atual) |
| `dbp-dot2png` | C | `.profile` + `.dot` → PNG | Gantt **+** DAG no mesmo frame |
| `profile2h5.py` + `pbt2ptt` | Python/Cython | `.prof` → HDF5 (DataFrames) | events/streams/nodes ricos |
| `parsec_trace_tables.py` (PTT) | Python | HDF5 → pandas | análise tabular |
| `parsec_dag.py` / `example-DAG-and-Trace.py` | Python | `.dot` (+HDF5) | DAG + cruzamento por chave |

> **Viabilidade no ambiente atual (checada):** as ferramentas **C** (`dbp2xml`, `dbpinfos`,
> `dbp2mem`, `dbp2paje`) são compiladas e instaladas em `bin` por padrão no
> `nix develop .#parsec` (e o `dbp2paje` já roda — foi ele que gerou o Paje atual). Já o
> caminho **Python/PTT** é frágil aqui: `parsec.nix` constrói com `BUILD_SHARED_LIBS=ON`
> mas **não** provê `pandas`/`numpy`/`tables`/`networkx` no shell, e os scripts
> (`parsec_dag.py`, `example-DAG-and-Trace.py`) são **Python 2** (`print "…"`, `cmp()`),
> dependendo do `2to3` do build.

### 4.4 A chave canônica `tpid_did_tid`

A tripla `(taskpool_id, type, id)` dos eventos HDF5 é **idêntica** à `(tpid, did, tid)` dos
nós do `.dot` (ver `example-DAG-and-Trace.py:42-57`). Essa tripla é o `JobId` natural por
instância de tarefa — a peça que falta hoje (`JobId` NA) e que liga
`Application`↔`Dag`↔`Tasks`.

---

## 5. Design do adaptador PaRSEC → tibbles do StarVZ

### 5.1 Arquitetura recomendada: C-tools → R

Como só as ferramentas do PaRSEC leem `.prof`, e como o caminho Python é frágil neste
ambiente, a rota **robusta** é:

```
  Estágio A (shell .#parsec, ferramentas C já compiladas)
    dbp2xml  cham_dpotrf-0.prof-*   →  events.xml   (todos os eventos: thread/stream,
                                                      key→nome, tpid, did, tid, begin, end, info)
    dbpinfos cham_dpotrf-0.prof-*   →  info.txt     (N, NB, sched, gflops, hostname)
    dbp2mem  cham_dpotrf-0.prof-*   →  mem.txt      (uso de memória → memory_state)
    [Fase B] <nome>-<rank>.dot      →  dag_edges    (src/dst = tpid_did_tid, flow)

  Estágio B (shell default, R + arrow)
    parse events.xml/info.txt/dag_edges  →  monta os tibbles no SCHEMA do StarVZ
                                            grava *.parquet numa pasta lida por starvz_read()
```

- **Por que R no estágio B:** o `arrow` funciona direto no shell default e o destino é o
  schema do StarVZ (R). Onde possível, reutilizar helpers do próprio StarVZ
  (`separate_res`, funções de cor, árvore Y de `hl_y_paje_tree`).
- **Alternativa documentada (Python/PTT):** se um dia provisionarmos
  `pandas`/`numpy`/`tables`/`networkx` no `.#parsec`, o `profile2h5.py` + `parsec_trace_tables`
  + `parsec_dag` dá um caminho mais direto (DataFrames prontos). Fica como opção, não como
  base.
- **Descartado:** estender `dbp2paje.c` em C para emitir Paje mais rico — continuaria sem
  `JobId`/DAG de forma natural e é o ponto mais lossy do fluxo.

### 5.2 Mapeamento por tibble (alvo StarVZ ← fonte PaRSEC ← derivação)

| Tibble StarVZ | Fonte PaRSEC | Como derivar |
|---|---|---|
| `Application` | `events.xml` (tipos de tarefa; heurística "tudo após `PUT_CB`") | `Value`=nome do tipo; `Start/End`=begin/end−`Zero`; `JobId`=`"tpid_did_tid"`; `ResourceId`/`ResourceType` **limpos** a partir de thread/stream (CPU `T#` → `CPU#`; GPU `G#S#` → `CUDA#`); `GFlop` por classe+params |
| `StarPU` (runtime) | `events.xml` de runtime (PINS `task_profiler`) | estados `Scheduling`/`Release deps`/… (requer PINS ligado) |
| `Dag` | `dag_edges` (do `.dot`, **Fase B**) | `(JobId, Dependent)`=`tpid_did_tid`; `full_join(Application)` p/ tempos/Value (igual ao `read_dag`) |
| `Tasks` | `events.xml` | `JobId, StartTime, EndTime, WorkerId`=stream/thread, `SubmitOrder`, `Iteration`, `Name` |
| `Task_handles`/`Data_handles` | flows do `.dot` (`flow_src/flow_dst`, **Fase B**) | qual dado cada tarefa acessa (modo R/W) |
| `Variable` | `events.xml` (exec + `movein`/`moveout`) | **derivado** por janelas de tempo: `GFlops`=Σ flops ativos/janela; `B. In/Out (MB/s)`=Σ bytes/janela; `Ready/Submitted` só com PINS |
| `Link` | eventos `movein`/`moveout` (info = ponteiro/size) | transferências H2D/D2H como links com `Size`, `Origin`, `Dest` |
| `Memory_state` | `dbp2mem` (`mem.txt`) | alocações/cópias ao longo do tempo |
| `Y`/`entities`/`Colors`/`Zero`/`Version` | streams/threads + `Application` | hierarquia MPI→VP/GPU→Thread/Stream; cores por classe; `Zero`=min(begin) |

### 5.3 Fase A (já dá) vs Fase B (precisa re-coleta)

- **Fase A — só com o `.prof` que JÁ existe** (`0012_…/cham_dpotrf-0.prof-wA1QQ1`, sem
  recompilar): reconstrói `Application`, `Tasks`, `Variable` (derivado), `Link`,
  `Memory_state`, com `JobId` real e `ResourceType` limpo (`CPU`/`CUDA`). Já enriquece
  muito a visualização atual do PaRSEC e corrige os rótulos `MGS-/MVT-`.
- **Fase B — exige re-coleta no PCAD**: recompilar o PaRSEC com `PARSEC_PROF_GRAPHER` e
  rodar com `-. <nome>` para gerar os `.dot`. Só então nascem `Dag`/`Task_handles`/
  `Data_handles` → **paridade total** com o StarPU (inclui caminho crítico).

---

## 6. Roteiro de re-coleta (Fase B)

1. **Build**: garantir `PARSEC_PROF_GRAPHER` ligado na compilação do PaRSEC (flag de CMake;
   ver `parsec/CMakeLists.txt` e `parsec.nix`). Confirmar `BUILD_SHARED_LIBS=ON` (já está).
2. **Run** (no `scripts/run.sh`, que hoje só seta `PARSEC_MCA_profile_filename`): acrescentar,
   para PaRSEC, o argumento de aplicação `-- -. cham_<algo>` (gera `cham_<algo>-<rank>.dot`)
   além do `PARSEC_MCA_profile_filename`. Opcional: `PARSEC_MCA_pins=task_profiler` (estados
   de runtime → `StarPU`) e flags de stream CUDA `TRACK_DATA_IN` (para enxergar H2D também).
3. **Pós-processamento** (Estágio A/B do §5): `dbp2xml`/`dbpinfos`/`dbp2mem` + parser do
   `.dot`, depois o builder R → parquets no schema do StarVZ.
4. **Validação**: `starvz_read()` na pasta gerada deve carregar `application$JobId` **não-NA**
   (n_distinct ≈ nº de tarefas) e `ResourceType ∈ {CPU, CUDA}`; `panel_st` desenha lanes
   limpas; `panel_gflops`/`panel_gpubandwidth` mostram curvas (derivadas) não-vazias; e o
   `Dag` permite caminho crítico — comparável lado a lado com o StarPU (`0013`).

> Benefício de diagnóstico: com `Dag` + eventos, o "eager D2H pushout" (771 `moveout`) fica
> **visível** correlacionando cada `moveout` com os sucessores no DAG — exatamente o que
> falta hoje para explicar o gap de desempenho (StarPU ~4150 vs PaRSEC ~1570 GFlops, ver
> [`parsec_cuda_results.md`](parsec_cuda_results.md)).

---

## 7. Apêndice — schemas-alvo dos tibbles principais

Colunas que o builder R deve produzir para casar com o StarVZ (tipos entre parênteses):

- **`Application`**: `ResourceId(factor)`, `Node(factor)`, `Resource(factor)`,
  `ResourceType(factor: CPU|CUDA)`, `Value(factor)`, `Start/End/Duration(double)`,
  `JobId(chr = "tpid_did_tid")`, `GFlop(num)`, `Iteration(int)`, `Outlier(lgl)`,
  `Height/Position(num)`.
- **`Dag`**: `JobId(chr)`, `Dependent(factor)`, `Start/End(double)`, `Cost(num = -(End-Start))`,
  `Value(chr)`. Ligado por `full_join(Application, by="JobId")`.
- **`Tasks`**: `JobId(chr)`, `SubmitOrder(int)`, `StartTime/EndTime(double, −Zero)`,
  `WorkerId(int)`, `Iteration(int)`, `Name(chr)`.
- **`Variable`**: `ResourceId(factor)`, `Type(factor ∈ {GFlops, B. In (MB/s), B. Out (MB/s),
  Used (MB), Ready, Submitted})`, `Start/End(double)`, `Value(double)`,
  `Node/Resource/ResourceType`.
- **`Link`**: `Container(factor)`, `Type(factor)`, `Start/End/Duration(double)`,
  `Size(int)`, `Origin(factor)`, `Dest(factor)`, `Key(factor)`.

---

### Resumo executivo

O StarVZ pode, sim, ser o visualizador único. A pobreza atual do PaRSEC vem de o StarVZ só
consumir um Paje lossy (`dbp2paje`) e pular `rec2csv`/DAG. As informações existem nas fontes
nativas do PaRSEC (`.prof` via `dbp2xml`/`dbpinfos`/`dbp2mem`; DAG via `PARSEC_PROF_GRAPHER`
+ `-.`), e a chave `tpid_did_tid` permite reconstruir `JobId` e ligar tudo. O adaptador
recomendado é **C-tools → R**; **Fase A** já reconstrói a maior parte a partir do `.prof`
existente, e **Fase B** (re-coleta com grapher) entrega paridade total, inclusive o DAG.
