# Transferências H2D/D2H: o que o gráfico mostra e por quê

Explicação do bullet de conclusão: *"StarPU re-encena a matriz por evicção
(H2D ≈ 3× a matriz no poti, onde ela não cabe na VRAM); PaRSEC minimiza H2D
mas paga writeback D2H de todo tile escrito"*.

## TL;DR

Os dois runtimes movem dados demais, mas por razões opostas e com custos
diferentes:

- **StarPU** traz o mesmo tile para a GPU **várias vezes** (H2D alto). No poti
  isso acontece porque a matriz não cabe na VRAM e o runtime evicta tiles que
  vai precisar de novo. E H2D **dói**: o kernel fica esperando o input chegar.
- **PaRSEC** (na nossa camada DTD+CUDA) traz cada tile **uma vez só** (H2D ≈
  1× os dados tocados), mas devolve ao host **todo tile que um kernel de GPU
  escreve**, imediatamente (D2H enorme). Só que D2H **não dói**: sai por um
  stream dedicado enquanto o próximo kernel já roda.

Ou seja: o PaRSEC move mais bytes no total, mas os bytes dele são baratos;
o StarPU move menos, mas os bytes dele estão no caminho crítico.

## Como ler o gráfico (`transfers_d2h.png`)

- **H2D** (host→device): tiles subindo para a GPU antes de um kernel rodar.
- **D2H** (device→host): tiles escritos na GPU voltando para a RAM.
- Linha tracejada = tamanho da matriz (12.8 GB para N=40000 FP64); pontilhada
  = VRAM útil (11.2 GB no poti, 22.8 GB no tupi). Se H2D fica muito acima da
  tracejada, o runtime está re-transferindo dados que já subiram.
- Fonte: **contadores nativos** dos runtimes (`STARPU_BUS_STATS` e o
  `device_show_statistics` do PaRSEC), gravados no `.log` de cada run com
  `TRACE_STATS=1` — não vêm do trace, então não sofrem dos caveats de timing.

## Os números (dpotrf, N=40000)

| | poti (VRAM 11.2 GB) H2D / D2H | tupi (VRAM 22.8 GB) H2D / D2H |
|---|---:|---:|
| starpu:dmda  | 36.4 / 31.1 GB | 22.1 / 16.7 GB |
| starpu:dmdas | 32.5 / 26.7 GB | 20.1 / 14.5 GB |
| parsec:gd    | **6.0** / 61.4 GB | **6.1** / 44.1 GB |
| parsec:lfq   | **6.0** / 68.1 GB | **6.1** / 53.5 GB |

Referências de tamanho: a matriz cheia é 12.8 GB, mas o Cholesky só toca o
triângulo inferior ≈ **6.5 GB**. Em LU (dgetrf_nopiv), que toca a matriz
inteira: StarPU 76.3 GB de H2D no poti (~6× a matriz!), PaRSEC 12.0 GB (1×);
D2H do PaRSEC chega a 123.8 GB (~10× a matriz).

## Lado StarPU: re-encenação por evicção (+ ping-pong de coerência)

O StarPU trata a VRAM como um **cache de tiles com LRU e coerência
MSI**: um tile sobe quando um kernel de GPU precisa dele e fica lá até ser
evictado (falta de espaço) ou invalidado (uma tarefa de CPU escreveu nele).

No **poti**, a matriz (12.8 GB) > VRAM útil (11.2 GB). O working set do
Cholesky right-looking revisita o trailing submatrix a cada iteração k, então
o LRU evicta tiles que serão pedidos de novo na iteração seguinte →
**re-staging**: os mesmos 6.5 GB tocados sobem como 36.4 GB (≈ 5.6× cada
tile em média; "a matriz re-encenada ~3×" do doc do gd usa a matriz cheia
como régua). O D2H do StarPU é o writeback desses tiles sujos ao serem
evictados.

O **tupi** mostra que evicção não é a história toda: lá a matriz cabe
folgada (12.8 < 22.8 GB) e mesmo assim o H2D é 22.1 GB (~3.4× os dados
tocados). Isso é **ping-pong de coerência CPU↔GPU**: ~60% dos flops rodam
nas CPUs, então um tile atualizado por um gemm de CPU invalida a cópia da
GPU e precisa subir de novo quando a GPU o lê na iteração seguinte. No poti
os dois efeitos (capacidade + coerência) se somam.

Por que isso custa makespan: o H2D é **síncrono do ponto de vista do
kernel** — a tarefa só dispara quando os inputs estão no device. Com a GPU
já sendo o gargalo do sistema (eff_gpu ≈ 1.01–1.02 em todos os runs), cada
espera de transferência alonga diretamente a execução. O prefetch do dmda
mitiga, mas não elimina, porque a evicção desfaz o prefetch.

## Lado PaRSEC: LRU eficaz no H2D, pushout incondicional no D2H

- **H2D mínimo**: 6.0–6.1 GB ≈ exatamente 1× o triângulo tocado (0.93×,
  idêntico nas duas máquinas e nos dois schedulers). O log nativo diz que só
  ~3.4% dos bytes *requeridos* foram de fato transferidos — o cache LRU do
  device segura 96.6% dos pedidos. O device layer do PaRSEC também tem
  afinidade de dados (a tarefa tende a ir para onde o dado já está), o que
  reduz o ping-pong que o StarPU sofre.
- **D2H gigante**: é o **pushout eager da nossa camada DTD+CUDA**
  (`cuda_parsec.md`): todo tile **escrito** por um kernel de GPU é copiado de
  volta ao host logo após o kernel, incondicionalmente — mesmo que o tile vá
  ser reescrito na GPU na iteração seguinte. A assinatura está na constância
  do D2H por kernel: **1.86 MB/kernel no poti e 7.45 MB/kernel no tupi ≈
  0.93× o tamanho do tile (b=500 → 2 MB; b=1000 → 8 MB)**, igual em gd, lfq,
  Cholesky e LU. Um tile escrito por kernel, sempre. Por isso o D2H escala
  com o nº de kernels de GPU: 61–68 GB (dpotrf) a 124–126 GB (dgetrf) no poti.
- Por que quase não custa makespan: o pushout sai num **stream de D2H
  dedicado**, sobreposto à computação; o kernel seguinte não espera. O custo
  existe (ocupa PCIe e pode atrasar H2D concorrente), mas fica majoritariamente
  **fora do caminho crítico**.

## A assimetria que fecha o argumento

| | StarPU | PaRSEC (nossa camada DTD) |
|---|---|---|
| H2D | alto (evicção + coerência) | mínimo (LRU + afinidade) |
| D2H | moderado (writeback ao evictar) | enorme (pushout de todo tile escrito) |
| custo do H2D | **no caminho crítico** (kernel espera) | idem, mas quase não há |
| custo do D2H | fora do caminho crítico | fora do caminho crítico |

No poti, onde a GPU já é fraca em FP64 e ainda é o gargalo, o H2D extra do
StarPU vira makespan — uma das duas pernas da vantagem do `gd` (a outra é a
CPU ~26% mais lenta por tarefa; ver `por_que_gd_vence.md` §4a/§4b). No tupi a
pressão de VRAM some, o H2D do StarPU cai e o `dmda/dmdas` dispara na frente.

## Caveats

- 1 repetição por config nos jobs de trace (mas os contadores são volumes
  determinísticos pela estrutura do algoritmo — variam pouco entre runs).
- O pushout eager é comportamento da **nossa** implementação DTD+CUDA, não
  uma fatalidade do PaRSEC: um writeback lazy (só ao evictar/finalizar)
  cortaria o D2H em ~5–10× sem mudar o H2D. É um trabalho futuro honesto.
- "~3× a matriz" no slide usa a matriz cheia (12.8 GB) como régua para o
  dpotrf no poti (36.4/12.8 ≈ 2.8). Se alguém perguntar pela régua "dados
  tocados" (6.5 GB), a razão é ~5.6× — o argumento só fica mais forte.

## Frase pronta para a apresentação

> "Os contadores nativos de PCIe mostram os dois modelos de memória: o StarPU
> usa a VRAM como cache com evicção — quando a matriz não cabe, ele re-encena
> os mesmos tiles várias vezes, e essa espera fica no caminho crítico do
> kernel. Nossa camada CUDA do PaRSEC faz o oposto: sobe cada tile uma vez,
> mas devolve ao host todo tile escrito, um por kernel — só que isso sai num
> stream assíncrono e quase não custa makespan."
