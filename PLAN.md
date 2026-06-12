# Problem
Today, we are facing a problem where Chameleon executed with PaRSEC does not use CUDA accelerated kernels to execute. With the StarPU runtime, Chameleon can execute with GPU.

We need to modify the source code of PaRSEC and Chameleon to use CUDA accelerated kernels to match StarPU runtime implementation inside Chameleon.

Currently, PaRSEC has support for CUDA but it is not used by Chameleon.

## Plan
1. Compile PaRSEC with CUDA support inside Nix.
2. Compile Chameleon with PaRSEC with CUDA support inside Nix.
3. Implement a CUDA backend for PaRSEC inside Chameleon. You dont need to implement all the kernels. But a sufficient amount to test the new CUDA implementation inside Chameleon.

## What i know
Inside the Chameleon repository, there are the runtime/ directory where all the supported runtime backends are implemented inside Chameleon.

## This machine
The current machine that we are talking to and making the development doesnt have a NVIDIA GPU, so we are not testing here. I have access to a machine with a NVIDIA GPU, but i cannot give you the access to it. I can test it myself manually.

I would like a test script in bash so i can test the CUDA implementation in this machine manually. It has a RTX 4070 on it.

For testing i'm executing this command with this variables

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

PARSEC_MCA_device_cuda_enabled=1 PARSEC_MCA_device_show_capabilities=1 PARSEC_MCA_device_cuda_verbose=20 PARSEC_MCA_device_verbose=20 chameleon_dtesting -o dgemm -n 2000 -b 100 -g 1 -t 4

So i need the script to run this command with the variables set.

## Goal
The goal is to implement a CUDA backend for PaRSEC inside Chameleon and compile it with Nix. You can write .md files to document your work and thought process for me and you remember all the issues and have more knowlage about the repositories.

## Current State
The CUDA backend for PaRSEC inside Chameleon is now functional end-to-end on
PCAD `poti` (RTX 4070). The gemm CUDA path is verified across `{s,d,z}gemm`
at multiple `(n, b)` geometries; PaRSEC's heterogeneous scheduler dispatches
tiles to the GPU, and `chameleon_dtesting` completes without
`cudaEventQuery ... illegal memory access`.

- Implementation details and the DTD CUDA chore patches:
  `docs/cuda_parsec.md`.
- Full sweep results and analysis (including why dgemm/zgemm speedups are
  ~1× on this consumer card and where sgemm has remaining headroom):
  `docs/parsec_cuda_results.md`.
- Raw logs from the verification run: `data/parsec-cuda_788764/`.

## Docs
THe development docs is located in /docs/. It is your responsability to read and maintain this docs up to date with the changes in the codebase.
