{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    starpu.url = "github:Sacolle/nix-starpu";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, starpu, flake-utils }: 
   flake-utils.lib.eachDefaultSystem (system: 
   let 
    		pkgs = import nixpkgs { inherit system; };
        chameleon = pkgs.callPackage ./chameleon.nix {
          starpu = starpu.packages.${system}.default.override {
            enableMPI = false;
            };
            parsec = pkgs.callPackage ./parsec.nix {};
        };
      in 
      {
        packages.parsec = pkgs.callPackage ./parsec.nix {};
        packages.chameleon = chameleon;
        devShells.default = pkgs.mkShell {
         buildInputs = [
          chameleon
         ];
        };
      }
   );
}
