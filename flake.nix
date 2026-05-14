{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    starpu.url = "github:Sacolle/nix-starpu";
    flake-utils.url = "github:numtide/flake-utils";
    starvz.url = "github:matregnago/starvz";
  };

  outputs =
    {
      self,
      nixpkgs,
      starpu,
      starvz,
      flake-utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs { inherit system; };

        starpuPkg = starpu.packages.${system}.default.override {
          enableMPI = false;
          enableTrace = true;
          compileAsRelease = true;
        };

        parsecPkg = pkgs.callPackage ./parsec.nix {
          stdenv = pkgs.gcc13Stdenv;
        };

        chameleon_starpu = pkgs.callPackage ./chameleon.nix {
          starpu = starpuPkg;
          parsec = parsecPkg;
          runtime = "starpu";
        };

        chameleon_parsec = pkgs.callPackage ./chameleon.nix {
          starpu = starpuPkg;
          parsec = parsecPkg;
          runtime = "parsec";
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
          default = chameleon_starpu;
        };
        devShells = {
          parsec = pkgs.mkShell {
            buildInputs = [
              chameleon_parsec
              parsecPkg
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
