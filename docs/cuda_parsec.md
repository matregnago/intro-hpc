# CUDA backend for PaRSEC inside Chameleon

## The gap

Chameleon's StarPU runtime ships cuBLAS-backed codelets out of the box
(`chameleon/runtime/starpu/codelets/codelet_zgemm.c`, `cl_zgemm_cuda_func`).
The PaRSEC runtime did not: every codelet under `chameleon/runtime/parsec/codelets/`
only registered a CPU function via `parsec_dtd_taskpool_insert_task`, so even
on a machine with an NVIDIA GPU, Chameleon-on-PaRSEC was CPU-only.

The reason this is not a one-codelet fix: PaRSEC's DTD (Dynamic Task
Discovery) interface — which Chameleon's PaRSEC backend uses everywhere — was
hard-wired to a single shared CPU chore table
(`parsec/parsec/interfaces/superscalar/insert_function.c:dtd_chore[]`). DTD
tasks always picked `chore_id = 0`, so there was no mechanism to attach a
CUDA implementation alongside a CPU one. JDF files have `BODY [type=CUDA]`
and the PTG compiler emits dispatch logic for that, but DTD has no
equivalent.

## What changed

### 1. PaRSEC DTD now supports per-task-class CUDA chores

`parsec/parsec/interfaces/superscalar/insert_function.c` and the related
headers gained:

- Per-task-class `incarnations[]` arrays (3 slots: CPU, optional CUDA, END).
  Previously every task class shared `dtd_chore[]`; now each task class owns
  its own table.
- A new public API:
  `parsec_dtd_taskpool_insert_task_cuda(tp, cpu_fp, cuda_fp, priority, name, …)`.
  It works exactly like `parsec_dtd_taskpool_insert_task` but registers a
  CUDA chore on first call for the task class.
- The CUDA chore hook (`hook_of_dtd_task_cuda`) follows the JDF-generated
  pattern (`parsec/parsec/interfaces/ptg/ptg-compiler/jdf2c.c:5042`+): pick the
  best CUDA device, build a `parsec_gpu_task_t`, set its `submit` callback to
  `gpu_kernel_submit_dtd_task_cuda` (which forwards to the user's CUDA function
  with the active `cudaStream_t`), and hand the task to
  `parsec_gpu_kernel_scheduler`.
- If no GPU is available at dispatch time (`PARSEC_MCA_device_cuda_enabled=0`,
  or `parsec_gpu_get_best_device()` returns < 2), the CUDA hook returns
  `PARSEC_HOOK_RETURN_NEXT` and the runtime falls back to the CPU chore.

### 1a. Bug fix: registration ordering for DTD taskpools

Once 1 was in place the CUDA hook still never won at dispatch time: the
verbose log printed `Device 2 (NVIDIA GeForce RTX 4070) disabled for taskpool
<ptr> @parsec_cuda_taskpool_register:362` for every Chameleon taskpool and
dgemm ran on the CPU. Root cause:

1. Chameleon registers the DTD taskpool **immediately** after creating it
   (`chameleon/runtime/parsec/control/runtime_async.c:30-32` —
   `parsec_dtd_taskpool_new()` then `parsec_context_add_taskpool()`).
2. DTD task classes are created **lazily** on the first
   `parsec_dtd_taskpool_insert_task[_cuda]` call, so at register time the
   taskpool has zero task classes.
3. `parsec_cuda_taskpool_register` (`parsec/parsec/devices/cuda/dev_cuda.c:324-365`)
   walks `tp->task_classes_array` looking for a chore with `type ==
   PARSEC_DEV_CUDA`. With zero task classes the loop body never runs, `rc`
   stays `PARSEC_ERR_NOT_FOUND`, and the CUDA bit is cleared from
   `tp->devices_index_mask` (line 361).
4. `INSERT_TASK_zgemm` later attaches the CUDA chore to the dgemm task
   class, but the taskpool's mask is already missing the CUDA bit.
5. At dispatch the CUDA hook calls `parsec_gpu_get_best_device`
   (`dev_cuda.c:1323`), which at line 1350 does `if(!(tp->devices_index_mask
   & (1 << dev_index))) continue;` — every CUDA device is skipped, best
   device falls back to CPU (index 0), the hook returns
   `PARSEC_HOOK_RETURN_NEXT`, and PaRSEC runs the CPU chore.

**Fix, part 1 (mask).** `parsec_dtd_taskpool_insert_task_cuda` now, on the
first CUDA insertion for a task class (the `NULL == dtd_tc->cuda_fpointer`
branch), walks all registered devices and for each CUDA device re-sets the
bit on `tp->devices_index_mask` and re-invokes the device's
`device_taskpool_register` hook so dyld resolution is also re-run. This runs
exactly once per task class (gated by `cuda_fpointer`), so it adds no
per-task overhead.

**Fix, part 2 (incarnation order).** With the mask fix alone, dgemm still
ran on the CPU. The scheduler at `parsec/parsec/scheduling.c:152-189`
iterates `incarnations[]` in order and only advances past a chore when its
hook returns `PARSEC_HOOK_RETURN_NEXT`. The DTD CPU hook
(`hook_of_dtd_task`) always returns `DONE`, so the CUDA chore at
`incarnations[1]` was never reached — CPU at slot 0 won every dispatch.

Mirroring the JDF convention (a `BODY [type=CUDA]` listed before the default
CPU body is preferred), `parsec_dtd_taskpool_insert_task_cuda` now installs
the CUDA chore at slot 0 and shifts the CPU chore to slot 1.
`hook_of_dtd_task_cuda` already returns `PARSEC_HOOK_RETURN_NEXT` when no
GPU is available (`parsec_gpu_get_best_device(...) < 2`), so the CPU
fallback still works.

**Fix, part 3 (data pointers).** With the CUDA hook reached, cuBLAS crashed
with `cudaEventQuery an illegal memory access was encountered`. Root cause:
`parsec_dtd_unpack_args` returns `this_task->data[i].data_in` — the host
copy. When PaRSEC's GPU scheduler stages a tile to the device it stores the
device-side `parsec_data_copy_t` in `this_task->data[i].data_out`, but the
DTD unpack never looks there. The CUDA codelet (`codelet_zgemm.c`) now
overrides A/B/C with `parsec_data_copy_get_ptr(this_task->data[i].data_out)`
after calling `parsec_dtd_unpack_args` — VALUE args still come from the
generic unpack, only the data pointers are replaced.

A cleaner long-term fix would be to add a `parsec_dtd_unpack_args_cuda`
inside PaRSEC that reads `data_out` for the PASSED_BY_REF flows; we chose
the in-codelet override because Chameleon currently has only one CUDA
codelet (dgemm) and we wanted to keep the PaRSEC patch minimal.

**Verified on the RTX 4070 (PCAD `poti` partition).** With all three fixes
in place the cuBLAS path is functionally correct: `chameleon_dtesting`
completes without `cudaEventQuery ... illegal memory access`, the GPU is
exercised across `{s,d,z}gemm`, and PaRSEC's heterogeneous scheduler
dispatches tiles to both CPU and GPU. The full sweep and the
interpretation of its GFlop/s and CPU↔GPU speedup numbers — including why
dgemm/zgemm look essentially flat on this consumer card — live in
`docs/parsec_cuda_results.md`. To re-verify after future changes, submit
`slurm/cuda_vs_cpu.slurm` (full CPU-vs-GPU sweep) or run
`bash scripts/run_cuda_vscpu.sh <cpu|gpu>` inside the matching dev shell
and check the resulting log for the absence of CUDA errors and a sensible
GFlop/s number.

### 2. Chameleon's dgemm has a cuBLAS path

`chameleon/runtime/parsec/codelets/codelet_zgemm.c` now declares both
`CORE_zgemm_parsec` (CPU) and `CORE_zgemm_parsec_cuda` (cuBLAS-backed via the
existing `CUDA_zgemm` wrapper in `chameleon/gpucublas/`). `INSERT_TASK_zgemm`
calls `parsec_dtd_taskpool_insert_task_cuda` when `CHAMELEON_USE_CUDA` is set,
falling back to the CPU-only insert otherwise.

Each worker thread keeps a `__thread cublasHandle_t` and rebinds the CUDA
stream per task — adequate for a proof of concept; handles leak on thread exit.

### 3. Nix build enables CUDA

`flake.nix`, `parsec.nix`, and `chameleon.nix` now thread `cudaPackages`
through and set `-DPARSEC_GPU_WITH_CUDA=ON` / `-DCHAMELEON_USE_CUDA=ON`.
Because the local `parsec/` and `chameleon/` source trees are untracked,
builds require `nix build .#chameleon-parsec --impure` (the flake reads them
through `builtins.getEnv "PWD"`).

### 4. Cholesky kernels on the GPU (trsm, syrk, potrf)

With the gemm path validated end-to-end (see
`data/parsec-cuda_788841/comparison.csv`), the three remaining kernels in
the Cholesky factorization got the same treatment:

- `chameleon/runtime/parsec/codelets/codelet_ztrsm.c` —
  `CORE_ztrsm_parsec_cuda` forwards to `CUDA_ztrsm` (pure cuBLAS, no
  workspace). Two PASSED_BY_REF flows (`T` INPUT, `B` INOUT) overridden from
  `data_out`, per-translation-unit `__thread cublasHandle_t` lazily created
  on first call.
- `chameleon/runtime/parsec/codelets/codelet_zsyrk.c` —
  `CORE_zsyrk_parsec_cuda` forwards to `CUDA_zsyrk` (also pure cuBLAS).
  Same shape as trsm: data_out override for `A` INPUT and `C` INOUT,
  `__thread cublasHandle_t`.
- `chameleon/runtime/parsec/codelets/codelet_zpotrf.c` —
  `CORE_zpotrf_parsec_cuda` forwards to `CUDA_zpotrf` (cuSOLVER). Diverges
  from the gemm pattern in three ways, all per-thread state:
  1. **cuSOLVER handle.** `__thread cusolverDnHandle_t` created on first
     call via `cusolverDnCreate`, stream rebound per task with
     `cusolverDnSetStream`.
  2. **Workspace.** `cusolverDnZpotrf_bufferSize` returns `lwork` (in
     `cuDoubleComplex` elements) for the current tile size; a
     `__thread cuDoubleComplex *` buffer is grown lazily via
     `cudaFree`/`cudaMalloc`. Tile size is constant across a run, so the
     realloc fires exactly once per worker thread.
  3. **Device info.** `__thread int *d_info` (4 bytes) holds the cuSOLVER
     return code. After the kernel call, an async D→H `cudaMemcpyAsync` +
     `cudaStreamSynchronize` mirrors StarPU's
     `chameleon/runtime/starpu/codelets/codelet_zpotrf.c` cuda body so that
     non-positive-definite tiles still surface through
     `RUNTIME_sequence_flush`. This sync is unavoidable for potrf — it is
     the only kernel in the Cholesky family that needs a host-visible
     status code per call.

No PaRSEC source changes were needed. The DTD CUDA chore plumbing landed
for gemm (mask fix, incarnation order, `data_out` pointers — §1 and §1a)
is the entire machinery these new codelets rely on; `INSERT_TASK_z<kernel>`
just switches to `parsec_dtd_taskpool_insert_task_cuda` under
`#ifdef CHAMELEON_USE_CUDA` and passes both function pointers.

### 5. QR kernels on the GPU (ztpmqrt, zunmqr)

`dgeqrf`'s task graph (see `chameleon/compute/pzgeqrf.c`) emits four
kernels: `zgeqrt` (panel), `ztpqrt` (vertical merge), `ztpmqrt` (trailing
update), and `zunmqr` (apply Q to RHS, used only when `D` is given). Of
these, only the last two are reasonable candidates for GPU dispatch — and
that is exactly the split upstream Chameleon uses for StarPU:

| StarPU codelet | CODELETS macro |
|---|---|
| `codelet_zgeqrt.c:62` | `CODELETS_CPU(zgeqrt, ...)` — CPU-only |
| `codelet_ztpqrt.c:53` | `CODELETS_CPU(ztpqrt, ...)` — CPU-only |
| `codelet_ztpmqrt.c:87` | `CODELETS(ztpmqrt, ..., cl_ztpmqrt_cuda_func, STARPU_CUDA_ASYNC)` |
| `codelet_zunmqr.c:91` | `CODELETS(zunmqr, ..., cl_zunmqr_cuda_func, STARPU_CUDA_ASYNC)` |

This is the standard hybrid QR pattern: the panel + vertical merge are on
the critical path but small per call and would lose to H↔D staging cost
on the GPU; the trailing update is embarrassingly parallel across columns
and dominates flops. Putting only those two on the GPU is "full coverage"
by upstream's definition.

Two new PaRSEC codelets, both following the gemm/trsm/syrk shape with one
extension (per-thread GPU scratch, sized lazily — same idea as potrf's
workspace):

- `chameleon/runtime/parsec/codelets/codelet_ztpmqrt.c` —
  `CORE_ztpmqrt_parsec_cuda` forwards to `CUDA_ztpmqrt` (pure cuBLAS,
  internally a `CUDA_zparfb` loop). Four PASSED_BY_REF flows (V, T, A,
  B → `data[0..3]`) overridden from `data_out`. A `__thread cuDoubleComplex *`
  scratch is grown to `3*ib*nb` elements (StarPU's documented size for
  this codelet) on first call; constant across a run, so the realloc
  fires once per worker.
- `chameleon/runtime/parsec/codelets/codelet_zunmqr.c` —
  `CORE_zunmqr_parsec_cuda` forwards to `CUDA_zunmqrt` (pure cuBLAS,
  internally a `CUDA_zlarfb` loop). Three PASSED_BY_REF flows (A, T,
  C → `data[0..2]`), scratch sized `ib*ldwork`.

`zgeqrt` and `ztpqrt` deliberately keep their CPU-only `INSERT_TASK`
calls. PaRSEC's heterogeneous scheduler will dispatch them to CPU
workers in parallel with ztpmqrt/zunmqr running on the GPU — exactly
what StarPU does. Note that this means dgeqrf's critical path (panel
factorization + vertical merges) executes on CPU; the GPU contribution
is the trailing update on every panel.

No PaRSEC source changes, no new wrappers in
`chameleon/gpucublas/compute/`, and no changes to `CUDA_ztpmqrt` /
`CUDA_zunmqrt` (both already used by StarPU).

## What's still on CPU

All three task graphs (`dgemm`, `dpotrf`, `dgeqrf`) now exercise the GPU.
Within `dgeqrf` the panel factorization (`zgeqrt`) and the vertical
merge (`ztpqrt`) intentionally remain CPU-only — that matches upstream
Chameleon's StarPU coverage (see the table in §5). The remaining gap
is not a PaRSEC port question but a wrapper question: no
`cuda_ztpqrt.c` exists in `chameleon/gpucublas/compute/`, and
`cuda_zgeqrt.c` is gated by `#if defined(CHAMELEON_USE_MAGMA)`
(MAGMA is not in our Nix build). Closing those gaps would require new
GPU implementations of those kernels, not just runtime plumbing.

## Running it

```bash
sbatch slurm/cuda_vs_cpu.slurm
```

The SLURM job enters the `.#parsec-cpu` and `.#parsec` dev shells once each
and runs the full CPU-vs-GPU sweep inside them via
`scripts/run_cuda_vscpu.sh`. To run a single side interactively, enter the
dev shell and call the runner directly, e.g.
`nix develop .#parsec --command bash -c 'RESULTS_DIR=out bash scripts/run_cuda_vscpu.sh gpu'`.
Override geometry/precision with `SWEEP_PAIRS`, `PRECISIONS`, `OPS`, `REPS`.

## Files

| File | Role |
| ---- | ---- |
| `parsec/parsec/interfaces/superscalar/insert_function.c` | DTD CUDA chore registration + dispatch |
| `parsec/parsec/interfaces/superscalar/insert_function.h` | `parsec_dtd_taskpool_insert_task_cuda` declaration |
| `parsec/parsec/interfaces/superscalar/insert_function_internal.h` | `cuda_fpointer` field on `parsec_dtd_task_class_s` |
| `chameleon/runtime/parsec/codelets/codelet_zgemm.c` | cuBLAS-backed gemm codelet |
| `chameleon/runtime/parsec/codelets/codelet_ztrsm.c` | cuBLAS-backed trsm codelet |
| `chameleon/runtime/parsec/codelets/codelet_zsyrk.c` | cuBLAS-backed syrk codelet |
| `chameleon/runtime/parsec/codelets/codelet_zpotrf.c` | cuSOLVER-backed potrf codelet + workspace/devInfo lifecycle |
| `chameleon/runtime/parsec/codelets/codelet_ztpmqrt.c` | cuBLAS-backed QR trailing-update codelet + 3*ib*nb scratch |
| `chameleon/runtime/parsec/codelets/codelet_zunmqr.c` | cuBLAS-backed apply-Q codelet + ib*ldwork scratch |
| `chameleon/gpucublas/compute/cuda_z{gemm,trsm,syrk,potrf,tpmqrt,unmqrt}.c` | (reused, unchanged) cuBLAS/cuSOLVER calls |
| `flake.nix`, `parsec.nix`, `chameleon.nix` | CUDA-aware Nix builds |
| `scripts/run_cuda_vscpu.sh` | in-shell CPU-vs-GPU sweep runner (`cpu`/`gpu`/`pivot`) |
| `slurm/cuda_vs_cpu.slurm` | SLURM job: PaRSEC CPU vs GPU sweep on PCAD `poti` |

## Open items

- **Eager D2H pushout defeats GPU data residency (biggest perf gap, not yet
  fixed).** `hook_of_dtd_task_cuda` in
  `parsec/parsec/interfaces/superscalar/insert_function.c` sets
  `gpu_task->pushout` on *every* write flow, so each GPU task copies its output
  back to host immediately. Always correct, but in a factorization a produced
  tile feeds the next GPU task, so this round-trips every tile H↔D. Measured on
  a local RTX 4060 Ti with the GPU-stream trace (spotrf, FP32, n=19200, b=960):
  919 GPU kernels → **919 D2H + 1129 H2D**, i.e. exactly one writeback per
  kernel. This is the main reason the PaRSEC spotrf GPU path trails StarPU
  (StarPU ~4150 GFlops vs PaRSEC ~1570 at n=19200; the FP32 gemm-heavy QR path,
  where tiles are reused less, is competitive). The JDF/PTG backend keeps tiles
  resident and only pushes out when a flow's *consumer* is a different task
  class (`jdf2c.c:5127`), co-designed with a scheduler that keeps same-class
  chains on the device.
  *Attempted fix (reverted):* a `PARSEC_DTD_PUSHOUT=0` "keep resident + stage
  back on demand" path — GPU hook leaves `pushout=0`, the DTD CPU hook
  (`hook_of_dtd_task`) and the data flush pull tiles to host before a host read.
  Two iterations both failed validation: (1) routing the pull through
  `parsec_data_transfer_ownership_to_copy()` **deadlocks** (it bumps the source
  copy's reader count, never balanced on the DTD CPU path, and the next GPU
  writer waits forever for readers to drain); (2) a hand-rolled version-compare
  + `cudaMemcpy` D2H + manual invalidate/version-bump stopped the deadlock but
  **returns wrong results** (`spotrf -c` FAILS) — the home-grown versioning does
  not interoperate with the GPU stage-in's version protocol
  (`dev_cuda.c:1096/1122/1150`). Conclusion: a correct fix must drive the demand
  D2H through PaRSEC's existing data-copy version/coherence machinery (mirroring
  the PTG CPU body at `jdf2c.c:5306` which calls
  `parsec_data_transfer_ownership_to_copy(...,0,...)` *plus* the real staging),
  not bypass it. The upside is real: the (incorrect) resident run hit ~4300
  GFlops on sgeqrf vs ~1150 with eager pushout.
- **dgemm/zgemm are FP64-hardware-bound on the RTX 4070.** Not actionable in
  software. The card's listed FP64 peak (458 GFlops, see
  `docs/parsec_cuda_results.md`) is below the 16-thread CPU's FP64 peak
  (~870 GFlops), so the speedup ceiling for these precisions on this
  hardware is ~1×. Moving the number requires a different GPU.
- **The `threads=1` in result lines is a reporting stub, not a real worker
  reduction.** `chameleon_dtesting`/`stesting` print `threads = CHAMELEON_GetThreadNbr()`,
  which on the PaRSEC backend calls `RUNTIME_thread_size()` —
  hardcoded to `return 1` with a `//TODO fixme` in
  `chameleon/runtime/parsec/control/runtime_control.c`. `RUNTIME_init` actually
  passes the requested core count to `parsec_init(ncpus, ...)`, and a trace
  confirms PaRSEC runs with all `-t` worker threads (e.g. 16 `M0V0T*` threads
  for `-t 16`). So the "overridden the number of required threads from N to 1"
  warning and the `threads=1` column are cosmetic; earlier notes that read this
  as PaRSEC throttling to one submitter were mistaken. (Worth fixing the stub so
  the column is honest.)
- **cuBLAS handle lifecycle.** Per-thread handles are still leaked at thread
  exit. Move to a registered fini hook if this becomes load-bearing. We
  considered an init hook to pre-create handles and decided against it: the
  per-task null-check is not on the critical path and an init hook would
  add a fragile lifecycle without moving the headline number.
- **zgeqrt / ztpqrt remain CPU-only — by design.** This matches StarPU's
  upstream coverage (both `codelet_zgeqrt.c` and `codelet_ztpqrt.c` are
  `CODELETS_CPU`). Putting them on the GPU would require new wrappers in
  `chameleon/gpucublas/compute/` (or enabling MAGMA in the Nix build for
  the existing `CUDA_zgeqrt`). Not a runtime gap — a kernel-library gap.
- **Validation against StarPU.** Compare PaRSEC+CUDA GFlops at a given
  problem size against a StarPU+CUDA run (`chameleon_dtesting` inside
  `nix develop .#starpu` with `-g 1`). StarPU is expected to be subject to
  the same FP64 ceiling for dgemm/zgemm; the comparison is most informative
  for sgemm.
