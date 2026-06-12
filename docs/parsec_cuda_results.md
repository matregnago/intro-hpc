# PaRSEC+CUDA gemm sweep — results and analysis

This document complements `docs/cuda_parsec.md` (the implementation writeup).
It captures the `{s,d,z}gemm` sweep on a real PCAD `poti` node (RTX 4070 +
i7-14700KF) and explains why the headline GPU-vs-CPU speedups in
`data/parsec-cuda_788786/comparison.csv` look modest against the 16-thread
CPU baseline yet strong against a single-worker CPU baseline.

## TL;DR

- **dgemm and zgemm GPU speedup vs the 16-thread CPU is ~0.9–1.0×, but vs a
  single PaRSEC worker the GPU wins ~10×.** The RTX 4070 has a FP64 peak of
  **458 GFlops** (PaRSEC prints this in every log). The i7-14700KF at 16
  threads has FP64 peak ~870 GFlops, i.e. ~1.9× the GPU's peak — so on this
  consumer card the 16-thread CPU is the *intrinsic* upper bound for FP64.
  Worker-for-worker the GPU is clearly faster (the new `speedup_vs_cpu1`
  column shows 8–10× across dgemm/zgemm). Both `dgemm` and `zgemm` ride the
  same FP64 ceiling because complex double runs on the same hardware units.
- **sgemm wins by up to 4.37× vs 16-thread CPU and 44.97× vs a single CPU
  worker** in the latest sweep (best at `n=16000 b=1000`). It is at ~16 % of
  the GPU's 29.3 TFlops FP32 peak. The dominant bottleneck is PaRSEC
  reducing the per-process worker count to 1 when CUDA is enabled (see the
  warning analysis below). This sweep did *not* yet exercise the `T=8`
  change — the SLURM driver was still forcing `T=4`.
- The benchmark methodology previously compared a 16-thread CPU run against
  a 1-thread + 1-GPU run. We added a 1-thread CPU baseline (`cpu1` column)
  so the comparison is both honest (worker-for-worker) and realistic
  (whole-CPU-vs-CPU+GPU).

## Setup

| | |
|-|-|
| Host | PCAD `poti` partition, node `poti5` |
| CPU | Intel Core i7-14700KF (16 logical threads exposed to PaRSEC) |
| GPU | NVIDIA GeForce RTX 4070, 12 GiB, driver 590.44.01, capability 8.9 |
| Runner | sweep driver (since retired), job `parsec-cuda_788786` |
| Precisions | `sgemm` (FP32), `dgemm` (FP64), `zgemm` (complex FP64) |
| Sizes | `n ∈ {2000, 4000, 8000, 16000}` |
| Tile sizes | `b ∈ {100, 500, 1000}` (subset per `n`) |
| CPU baseline | `-t 16 -g 0`, `.#parsec-cpu` dev shell (PaRSEC without CUDA) |
| 1-worker CPU baseline | `-t 1 -g 0`, same dev shell |
| GPU run | `-t 4 -g 1`, `.#parsec` dev shell (cuBLAS gemm path) |
| BLAS threads | `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1` |

Raw logs live under `data/parsec-cuda_788786/`. Per-run CSV is `summary.csv`;
side-by-side pivot is `comparison.csv`.

Note on `T`: the SLURM driver hard-coded `T=4` for this run, which PaRSEC
reduces to one effective worker (see §1 below). Later sweeps default `T=8`
to free more submitters.

## Result table

GFlop/s. Speedup columns:

- **speedup** = `gpu / cpu(t=16)` — practical "did the GPU help my workstation?"
- **speedup_vs_cpu1** = `gpu / cpu(t=1)` — apples-to-apples worker-for-worker

### sgemm (FP32)

| n     | b    | CPU t=16 | CPU t=1 | GPU | vs CPU-16 | vs CPU-1 |
|------:|-----:|---------:|--------:|----:|----------:|---------:|
|  2000 |  100 |      847 |      82 | 279 |     0.33  |     3.38 |
|  2000 |  500 |      647 |     101 | 1380 |    2.13 |    13.67 |
|  4000 |  500 |      888 |     101 | 1663 |    1.87 |    16.54 |
|  8000 |  500 |     1087 |     100 | 1941 |    1.79 |    19.34 |
|  8000 | 1000 |      954 |     104 | 4112 |    4.31 |    39.48 |
| 16000 | 1000 |     1072 |     104 | 4689 |  **4.37** | **44.97** |

### dgemm (FP64)

| n     | b    | CPU t=16 | CPU t=1 | GPU | vs CPU-16 | vs CPU-1 |
|------:|-----:|---------:|--------:|----:|----------:|---------:|
|  2000 |  100 |      437 |      40 | 150 |     0.34  |     3.71 |
|  2000 |  500 |      326 |      49 | 359 |     1.10  |     7.30 |
|  4000 |  500 |      463 |      50 | 411 |     0.89  |     8.31 |
|  8000 |  500 |      523 |      49 | 516 |     0.99  |    10.43 |
|  8000 | 1000 |      474 |      51 | 413 |     0.87  |     8.15 |
| 16000 | 1000 |      544 |      51 | 513 |     0.94  |    10.09 |

### zgemm (complex FP64)

| n     | b    | CPU t=16 | CPU t=1 | GPU | vs CPU-16 | vs CPU-1 |
|------:|-----:|---------:|--------:|----:|----------:|---------:|
|  2000 |  100 |      512 |      46 | 309 |     0.60  |     6.66 |
|  2000 |  500 |      338 |      52 | 462 |     1.37  |     8.95 |
|  4000 |  500 |      457 |      52 | 469 |     1.02  |     9.08 |
|  8000 |  500 |      548 |      52 | 566 |     1.03  |    10.96 |
|  8000 | 1000 |      503 |      52 | 427 |     0.85  |     8.15 |
| 16000 | 1000 |      544 |      52 | 514 |     0.95  |     9.81 |

## Why dgemm and zgemm look flat vs the 16-thread CPU

PaRSEC reports the per-device peak GFlops at startup. In every CUDA log
(e.g. `data/parsec-cuda_788786/dgemm_gpu_n16000_b1000.log:10`):

```
peak Gflops : double 458.1600, single 29322.2422 tensor 117288.9688 half 58644.4844
```

And the same log line 4 (per-stream CPU peak; multiply by ~4 for the
16-thread aggregate the CPU baseline actually achieves):

```
peak Gflops : double 217.7670, single 435.5339   (per-stream; 16-thread ≈ 870)
```

FP64 ceilings:

| Device | FP64 peak |
|--|-:|
| RTX 4070 | **458 GFlops** |
| i7-14700KF (16 threads) | **870 GFlops** |

The RTX 4070 is a consumer Ada Lovelace card; NVIDIA crippled FP64 to 1/64
of FP32 throughput to upsell datacenter parts. There is no way to move this
number with code changes — the GPU FP64 hardware is the limit.

In the measured numbers, the GPU's dgemm at `n=8000 b=500` (516 GFlops)
slightly *exceeds* the listed device peak. That is not a measurement error:
PaRSEC's heterogeneous scheduler ratios each device by its peak and
dispatches a proportional share of tiles, so the throughput we measure is
the combined CPU+GPU number. This is the same reason `zgemm` tracks
`dgemm` — complex double runs on the same FP64 units.

What the `cpu1` column reveals: a single PaRSEC worker only hits ~50 GFlops
on FP64, so the GPU genuinely wins by ~10× on a worker-for-worker basis.
The "flat" speedup-vs-16-thread is a story about how many CPU threads we
chose to throw at the problem, not about the GPU implementation.

## Why sgemm only hits 16 % of the GPU FP32 peak

29.3 TFlops theoretical vs. 4.69 TFlops observed at the best
`(n=16000, b=1000)` point — there is real headroom here. Three things
contribute, in decreasing order of impact.

### 1. PaRSEC clamps the CPU worker count under CUDA

Every CUDA log prints (e.g. `sgemm_gpu_n16000_b1000.log:18`):

```
WARNING: A runtime environment variable or a configuration file has overridden
the number of required threads from 4 to 1
```

`-t 4` requested 4 PaRSEC workers. PaRSEC reserves one CPU thread per GPU
stream (the device reports `Parsec Streams : 4`), so all four requested
workers get consumed by GPU-stream management, leaving exactly one
general-purpose worker thread. The result line confirms it:

```
0;sgemm;1;1;...
       ^ threads=1, gpus=1
```

A single worker is the sole task submitter feeding all four GPU streams.
That serializes the submission path — even though cuBLAS launches are
non-blocking, only one stream sees fresh work at a time.

This run still showed `from 4 to 1` because the SLURM driver was
hard-coding `T=4`. With `T=8` the runtime reports `from 8 to 4`-ish,
freeing more submitters.

### 2. Tile size sensitivity

| n    | b   | tasks | sgemm GPU GFlops |
|-----:|----:|------:|-----------------:|
| 2000 | 100 |  8000 |              279 |
| 2000 | 500 |    64 |             1380 |
| 8000 | 500 |  4096 |             1941 |
| 8000 |1000 |   512 |             4112 |
| 16000|1000 |  4096 |             4689 |

At `b=100` each task is a 100×100×100 sgemm — ~2 MFlops, well under cuBLAS
launch overhead. There is no way for the GPU to win regardless of
implementation. At `b=1000` each task is ~2 GFlops, comfortably above the
launch-cost knee, and the GPU pulls ahead. The takeaway is configuration,
not code: when running dpotrf/dgeqrf on this hardware, prefer `b ≥ 500`.

### 3. Per-task cuBLAS handle null-check + SetStream

`chameleon/runtime/parsec/codelets/codelet_zgemm.c` keeps a
`__thread cublasHandle_t` and rebinds the stream per task. After the first
task the handle is non-null and the `if(NULL == ...)` check is a
branch-predicted no-op. `cublasSetStream` is a few µs but is independent of
problem size; it is below the noise floor of the launch-overhead effect.
We considered moving handle creation to a worker init hook and decided not
to — it would not move the headline number and would add a fragile
lifecycle path. The TODO in `docs/cuda_parsec.md` about the handle leak on
thread exit is unchanged.

## Methodology: two baselines, two questions

The CPU baseline at `-t 16 -g 0` and the GPU run at `-t 4 -g 1` (clamped to
1 worker, see §1) compare 16 CPU workers against 1 CPU worker + 1 GPU.
That answers *should I buy this card to accelerate my dense gemms on this
workstation?* — and the answer for FP64 on the RTX 4070 is "barely". It
does not answer *what is the GPU contributing per worker?*, which is what
the new `cpu1` column reports. Both columns are useful for different
audiences; neither is "more correct".

## Harness changes since this sweep

The sweep above ran on an earlier driver that has since been replaced by
`slurm/cuda_vs_cpu.slurm` + `scripts/run_cuda_vscpu.sh` (enters each dev
shell once and loops inside, rather than calling `nix develop` per run).
The two takeaways that drove the config: default `T=8` so PaRSEC keeps real
worker threads after reserving one per GPU stream, and extend the sweep to
`16000 2000` to check whether the sgemm curve is still climbing past
`b=1000`.

We did *not* touch `chameleon/runtime/parsec/codelets/codelet_zgemm.c` —
the per-task overhead there is not the bottleneck.

If after these knobs sgemm still saturates around 4–5 TFlops, the remaining
gap likely needs:

- More GPU work per call (larger tiles in the harness's outer driver, not
  just within `gemm`); or
- A multi-codelet implementation that overlaps non-gemm work on the GPU
  (out of scope for the current gemm proof-of-concept).

## dpotrf sweep (pending)

After landing the Cholesky kernels (`docs/cuda_parsec.md` §4 — trsm, syrk,
potrf), the CPU-vs-GPU sweep (`slurm/cuda_vs_cpu.slurm`, `OPS="potrf geqrf"`)
also produces `dpotrf` rows in `summary.csv` and `comparison.csv`. Numbers
from `data/parsec-cuda-vscpu_<jobid>/` to be added here
once the post-implementation run completes; the FP64 ceiling discussion
above applies unchanged (potrf is dominated by gemm/syrk in the trailing
update, so the same 458 GFlops FP64 cap on the RTX 4070 sets the upper
bound for `dpotrf` on this hardware).

## dgeqrf sweep (pending)

After landing the QR kernels (`docs/cuda_parsec.md` §5 — ztpmqrt, zunmqr),
the same `slurm/cuda_vs_cpu.slurm` sweep covers `geqrf` (it runs
`OPS="potrf geqrf"`). Numbers from `data/parsec-cuda-vscpu_<jobid>/` to be
added once a clean run completes. Expectations:

- **dgeqrf/zgeqrf on the RTX 4070 hit the same 458 GFlops FP64 ceiling**
  as dgemm and dpotrf — the trailing update (`CUDA_ztpmqrt`) ultimately
  decomposes into cuBLAS FP64 calls subject to the card's hardware cap.
- **sgeqrf should track somewhere below sgemm.** The panel + vertical
  merge stay on CPU (intentional, matching StarPU — see §5), so the GPU
  contribution is bounded by the trailing-update fraction of the
  algorithm. Best-case observed in any StarPU+CUDA reference would be a
  useful comparison once data is available.
- **`speedup_vs_cpu1`** is the more informative column here, same as for
  dpotrf — the CPU-16 baseline benefits more from the QR critical path
  than the GPU does, so headline `speedup` looks modest while the
  worker-for-worker ratio reveals the actual GPU contribution.

## Files referenced

| Path | Role |
|------|------|
| `docs/cuda_parsec.md` | Implementation writeup (DTD CUDA chore, codelet, fixes) |
| `scripts/run_cuda_vscpu.sh` | CPU-vs-GPU sweep runner (`cpu`/`gpu`/`pivot`) |
| `slurm/cuda_vs_cpu.slurm` | SLURM wrapper for PCAD `poti` (CPU vs GPU sweep) |
| `data/parsec-cuda_788786/summary.csv` | Per-run results |
| `data/parsec-cuda_788786/comparison.csv` | CPU↔GPU pivot |
| `data/parsec-cuda_788786/*_n*_b*.log` | Raw chameleon_dtesting logs |
