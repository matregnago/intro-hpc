{ pkgs ? import <nixpkgs> { } }:

rec {
  emacs = pkgs.emacs-nox.pkgs.withPackages (epkgs: [ epkgs.org-ref ]);

  pdf = pkgs.stdenv.mkDerivation {
    name = "sbc-paper";
    # SSCAD.org references figures as ../../plots/final/..., mirroring the
    # real repo layout (papers/2026_SSCAD/ two levels under the repo root
    # that holds plots/). Recreate that nesting in the sandbox instead of
    # relying on the default flat unpack (which drops both srcs as unrelated
    # top-level dirs and breaks those relative links).
    unpackPhase = ''
      mkdir -p papers/2026_SSCAD
      cp -r ${./.}/. papers/2026_SSCAD/
      cp -r ${../../plots} plots
      chmod -R u+w papers plots
    '';
    sourceRoot = "papers/2026_SSCAD";
    nativeBuildInputs = [ emacs pkgs.texliveFull ];
    buildPhase = ''
      emacs --batch --load export.el
      latexmk -pdf SSCAD.tex
    '';
    installPhase = ''
      mkdir -p $out
      cp SSCAD.pdf $out/
    '';
  };
}
