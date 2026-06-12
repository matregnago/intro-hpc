{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    # Pinned snapshot whose cudaPackages.backendStdenv is gcc13-compatible
    # (current nixos-unstable moved to gcc14 which nvcc 12.9 cannot parse).
    cudaNixpkgs.url = "github:nixos/nixpkgs/1da52dd49a127ad74486b135898da2cef8c62665";
    starpu.url = "github:Sacolle/nix-starpu";
    flake-utils.url = "github:numtide/flake-utils";
    starvz.url = "github:matregnago/starvz";
  };

  outputs =
    {
      self,
      nixpkgs,
      cudaNixpkgs,
      starpu,
      starvz,
      flake-utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs { inherit system; };

        cudapkgs = import cudaNixpkgs {
          inherit system;
          config = {
            allowUnfree = true;
            cudaSupport = true;
          };
        };

        starpuPkg = starpu.packages.${system}.default.override {
          enableMPI = false;
          enableTrace = true;
          compileAsRelease = true;
        };

        starpuPkgCuda = starpu.packages.${system}.default.override {
          enableMPI = false;
          enableTrace = true;
          compileAsRelease = true;
          enableCUDA = true;
        };
        parsecPkg = cudapkgs.callPackage ./parsec.nix {
          stdenv = cudapkgs.gcc13Stdenv;
          cudaPackages = cudapkgs.cudaPackages;
        };

        chameleon_starpu = pkgs.callPackage ./chameleon.nix {
          starpu = starpuPkg;
          parsec = parsecPkg;
          runtime = "starpu";
          enableCuda = false;
          cudaPackages = null;
          autoAddDriverRunpath = null;
        };

        chameleon_starpu_cuda = cudapkgs.callPackage ./chameleon.nix {
          stdenv = cudapkgs.gcc13Stdenv;
          starpu = starpuPkgCuda;
          parsec = parsecPkg;
          runtime = "starpu";
          cudaPackages = cudapkgs.cudaPackages;
          enableCuda = true;
        };
        chameleon_parsec = cudapkgs.callPackage ./chameleon.nix {
          stdenv = cudapkgs.gcc13Stdenv;
          starpu = starpuPkg;
          parsec = parsecPkg;
          runtime = "parsec";
          cudaPackages = cudapkgs.cudaPackages;
        };

        parsecPkgCpu = pkgs.callPackage ./parsec.nix {
          enableCuda = false;
          cudaPackages = null;
          autoAddDriverRunpath = null;
        };

        chameleon_parsec_cpu = pkgs.callPackage ./chameleon.nix {
          starpu = starpuPkg;
          parsec = parsecPkgCpu;
          runtime = "parsec";
          enableCuda = false;
          cudaPackages = null;
          autoAddDriverRunpath = null;
        };

        starvz_tools = starvz.packages.${system}.default;
        starvz_r = starvz.packages.${system}.starvz;
        rEnv = pkgs.rWrapper.override {
          packages = [
            starvz_r
          ]
          ++ (with pkgs.rPackages; [
            languageserver
            lintr
            here
            DoE_base
            FrF2
            tidyverse
            janitor
          ]);
        };
      in
      {
        packages = {
          parsec = parsecPkg;
          parsec-cpu = parsecPkgCpu;
          starpu-cuda = starpuPkgCuda;
          chameleon-starpu = chameleon_starpu;
          chameleon-starpu-cuda = chameleon_starpu_cuda;
          chameleon-parsec = chameleon_parsec;
          chameleon-parsec-cpu = chameleon_parsec_cpu;
          default = chameleon_starpu;
        };
        devShells = {
          parsec = cudapkgs.mkShell {
            buildInputs = [
              chameleon_parsec
              parsecPkg
            ];
            shellHook = ''
              # Host's libcuda.so + libnvidia* aren't in the Nix store; copy
              # them into a tempdir and prepend to LD_LIBRARY_PATH so PaRSEC's
              # cudaGetDeviceCount finds a driver matching the cudart it was
              # linked against. Mirrors what teste_cuda does.
              DRIVER_DIR="$(mktemp -d)"
              for d in /usr/lib /usr/lib/x86_64-linux-gnu /usr/lib64; do
                [ -d "$d" ] || continue
                for p in "$d"/libcuda.so* "$d"/libnvidia*.so* "$d"/libnvcuvid.so*; do
                  [ -e "$p" ] || [ -L "$p" ] && cp -Lf "$p" "$DRIVER_DIR/" 2>/dev/null || true
                done
              done
              export LD_LIBRARY_PATH="$DRIVER_DIR:''${LD_LIBRARY_PATH:-}"
              unset CUDA_VISIBLE_DEVICES
            '';
          };
          starpu-cuda = cudapkgs.mkShell {
            buildInputs = [
              chameleon_starpu_cuda
              starpuPkgCuda
            ];
            shellHook = ''
              # Host's libcuda.so + libnvidia* aren't in the Nix store; copy
              # them into a tempdir and prepend to LD_LIBRARY_PATH so PaRSEC's
              # cudaGetDeviceCount finds a driver matching the cudart it was
              # linked against. Mirrors what teste_cuda does.
              DRIVER_DIR="$(mktemp -d)"
              for d in /usr/lib /usr/lib/x86_64-linux-gnu /usr/lib64; do
                [ -d "$d" ] || continue
                for p in "$d"/libcuda.so* "$d"/libnvidia*.so* "$d"/libnvcuvid.so*; do
                  [ -e "$p" ] || [ -L "$p" ] && cp -Lf "$p" "$DRIVER_DIR/" 2>/dev/null || true
                done
              done
              export LD_LIBRARY_PATH="$DRIVER_DIR:''${LD_LIBRARY_PATH:-}"
              unset CUDA_VISIBLE_DEVICES
            '';
          };
          parsec-cpu = pkgs.mkShell {
            buildInputs = [
              chameleon_parsec_cpu
              parsecPkgCpu
            ];
          };
          starpu = pkgs.mkShell {
            buildInputs = [
              chameleon_starpu
            ];
          };
          default = pkgs.mkShell {
            buildInputs = [
              pkgs.marp-cli
              pkgs.chromium
              starvz_tools
              rEnv
            ];
          };
        };
      }
    );
}
