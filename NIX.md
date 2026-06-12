# Nix
The Nix package manager is used in this repo to reproduce the environment. Today, the derivations for PaRSEC and Chameleon are based from their respective remote repository. 

## Changes to make
We'll need to make some changes in both software. So, first of all, we need to point the nix derivations to the local source code instead of the remote repositories. They are located in parsec/ and chameleon/ respectively.
