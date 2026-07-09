{
  description = "chameleon with parsec and starpu flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    cudaNixpkgs.url = "github:nixos/nixpkgs/1da52dd49a127ad74486b135898da2cef8c62665";
    starpu.url = "github:Sacolle/nix-starpu";
    starvz.url = "github:schnorr/starvz";
    nix-gl-host.url = "github:numtide/nix-gl-host";
  };

  outputs =
    {
      self,
      nixpkgs,
      cudaNixpkgs,
      starpu,
      starvz,
      nix-gl-host,
    }:
    let
      system = "x86_64-linux";
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
      nixglhost = nix-gl-host.defaultPackage.${system};
    in
    {

      devShells.${system} = {
        parsec = cudapkgs.mkShell {
          buildInputs = [
            chameleon_parsec
            parsecPkg
            nixglhost
          ];
        };
        starpu = cudapkgs.mkShell {
          buildInputs = [
            chameleon_starpu
            starpuPkg
            nixglhost
          ];
        };
        default = pkgs.mkShell {
          buildInputs = [
            pkgs.marp-cli
            pkgs.chromium
            pkgs.graphviz
            pkgs.just
            starvz_tools
            rEnv
            pkgs.texliveFull
          ];
        };
      };
    };
}
