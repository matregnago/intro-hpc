{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    starpu.url = "github:Sacolle/nix-starpu";
    flake-utils.url = "github:numtide/flake-utils";
    starvz.url = "github:schnorr/starvz";
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
        chameleon = pkgs.callPackage ./chameleon.nix {
          starpu = starpu.packages.${system}.default.override {
            enableMPI = false;
          };
          parsec = pkgs.callPackage ./parsec.nix { };
        };
        starvz_pkg = starvz.packages.${system}.default;
      in
      {
        packages.parsec = pkgs.callPackage ./parsec.nix { };
        packages.chameleon = chameleon;
        devShells.default = pkgs.mkShell {
          buildInputs = [
            chameleon
            starvz_pkg
            pkgs.marp-cli
            pkgs.chromium
            pkgs.eztrace
          ];
        };
      }
    );
}
