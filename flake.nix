{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
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
          enableCUDA = true;
        };
        parsecPkg = cudapkgs.callPackage ./parsec.nix {
          stdenv = cudapkgs.gcc13Stdenv;
          cudaPackages = cudapkgs.cudaPackages;
        };

        chameleon_starpu = cudapkgs.callPackage ./chameleon.nix {
          stdenv = cudapkgs.gcc13Stdenv;
          starpu = starpuPkg;
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
          chameleon-starpu = chameleon_starpu;
          chameleon-parsec = chameleon_parsec;
          starvz-tools = starvz_tools;
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
              # On WSL2 the genuine GPU driver lives in /usr/lib/wsl/lib (an
              # isolated dir holding only libcuda/libnvidia*). Point at it
              # directly: copying libcuda.so out of it breaks its dlopen of the
              # sibling WSL driver shims, so cudaGetDeviceCount would still fail
              # with "driver version is insufficient" and fall back to the stub.
              if [ -e /usr/lib/wsl/lib/libcuda.so.1 ]; then
                export LD_LIBRARY_PATH="/usr/lib/wsl/lib:''${LD_LIBRARY_PATH}"
              fi
              unset CUDA_VISIBLE_DEVICES
            '';
          };
          starpu = cudapkgs.mkShell {
            buildInputs = [
              chameleon_starpu
              starpuPkg
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
              # On WSL2 the genuine GPU driver lives in /usr/lib/wsl/lib (an
              # isolated dir holding only libcuda/libnvidia*). Point at it
              # directly: copying libcuda.so out of it breaks its dlopen of the
              # sibling WSL driver shims, so cudaGetDeviceCount would still fail
              # with "driver version is insufficient" and fall back to the stub.
              if [ -e /usr/lib/wsl/lib/libcuda.so.1 ]; then
                export LD_LIBRARY_PATH="/usr/lib/wsl/lib:''${LD_LIBRARY_PATH}"
              fi
              unset CUDA_VISIBLE_DEVICES
            '';
          };
          default = pkgs.mkShell {
            buildInputs = [
              pkgs.marp-cli
              pkgs.chromium
              pkgs.graphviz
              starvz_tools
              rEnv
            ];
          };
        };
      }
    );
}
