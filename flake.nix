{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    # pinned commit used by nix-starpu for CUDA 13 (glibc 2.40-36 workaround)
    cudaNixpkgs.url = "github:nixos/nixpkgs/1da52dd49a127ad74486b135898da2cef8c62665";
    starpu.url = "github:Sacolle/nix-starpu";
    flake-utils.url = "github:numtide/flake-utils";
    starvz.url = "github:schnorr/starvz";
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
            cudaVersion = "13";
          };
        };

        starpuPkg = starpu.packages.${system}.default.override {
          enableMPI = false;
          enableTrace = true;
        };
        starpuCudaPkg = starpu.packages.${system}.default.override {
          enableMPI = false;
          enableTrace = true;
          enableCUDA = true;
        };
        parsecPkg = pkgs.callPackage ./parsec.nix {
          stdenv = pkgs.gcc13Stdenv;
        };

        parsecCudaPkg = cudapkgs.callPackage ./parsec_cuda.nix {
          stdenv = cudapkgs.gcc13Stdenv;
          cudaPackages = cudapkgs.cudaPackages;
        };

        chameleon_starpu = pkgs.callPackage ./chameleon.nix {
          starpu = starpuPkg;
          parsec = parsecPkg;
          runtime = "starpu";
        };

        chameleon_starpu_cuda = pkgs.callPackage ./chameleon.nix {
          starpu = starpuCudaPkg;
          parsec = parsecPkg;
          runtime = "starpu";
        };

        chameleon_parsec = pkgs.callPackage ./chameleon.nix {
          starpu = starpuPkg;
          parsec = parsecPkg;
          runtime = "parsec";
        };

        chameleon_parsec_cuda = cudapkgs.callPackage ./chameleon_parsec_cuda.nix {
          stdenv = cudapkgs.gcc13Stdenv;
          parsec = parsecCudaPkg;
          cudaPackages = cudapkgs.cudaPackages;
        };

        starvz_tools = starvz.packages.${system}.default;
        starvz_r = starvz.packages.${system}.starvz;
        pyEnv = pkgs.python3.withPackages (ps: with ps; [
          matplotlib
          pandas
          numpy
        ]);
      in
      {
        packages = {
          parsec = parsecPkg;
          parsec-cuda = parsecCudaPkg;
          chameleon-starpu = chameleon_starpu;
          chameleon-parsec-cuda = chameleon_parsec_cuda;
          chameleon-starpu-cuda = chameleon_starpu_cuda;
          chameleon-parsec = chameleon_parsec;
          default = chameleon_starpu;
        };
        devShells = {
          parsec = pkgs.mkShell {
            buildInputs = [
              chameleon_parsec
              pkgs.eztrace
            ];
          };
          starpu = pkgs.mkShell {
            buildInputs = [
              chameleon_starpu
              starvz_tools
              starvz_r
            ];
          };
          starpu_cuda = pkgs.mkShell {
            buildInputs = [
              chameleon_starpu_cuda
              starvz_tools
              starvz_r
            ];
          };
          parsec_cuda = cudapkgs.mkShell {
            buildInputs = [
              chameleon_parsec_cuda
              pkgs.eztrace
            ];
            shellHook = ''
              DRIVER_DIR="$(mktemp -d)"
              # Arch Linux: /usr/lib — Ubuntu/PCAD: /usr/lib/x86_64-linux-gnu and /usr/lib64
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
          default = pkgs.mkShell {
            buildInputs = [
              pkgs.marp-cli
              pkgs.chromium
              pyEnv
            ];
          };
        };
      }
    );
}
