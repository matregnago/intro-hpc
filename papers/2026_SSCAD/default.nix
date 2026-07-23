{ pkgs ? import <nixpkgs> { } }:

rec {
  emacs = pkgs.emacs-nox.pkgs.withPackages (epkgs: [ epkgs.org-ref ]);

  pdf = pkgs.stdenv.mkDerivation {
    name = "sbc-paper";
    srcs = [ ./. ../plots ];
    sourceRoot = "papers/2026_SSCAD/";
    nativeBuildInputs = [ emacs pkgs.texliveFull ];
    buildPhase = ''
      emacs --batch --load export.el
      latexmk -pdf SSCAD.tex
    '';
    installPhase = ''
      mkdir -p $out
      cp sbc.pdf $out/
    '';
  };
}
