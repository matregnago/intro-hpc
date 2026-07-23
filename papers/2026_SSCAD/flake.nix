{
  description = "Paper SBC escrito em Org-mode, exportado e compilado para PDF";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        paper = import ./default.nix { inherit pkgs; };
      in
      {
        packages.default = paper.pdf;
      }
    );
}
