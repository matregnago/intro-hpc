# Passos para rodar o projeto

1. Instalar o `Nix`:

```bash
sh <(curl --proto '=https' --tlsv1.2 -L https://nixos.org/nix/install) --daemon
```

2. Habilitar os flakes:

Adicione a linha a seguir em `~/.config/nix/nix.conf` ou `/etc/nix/nix.conf`:

```bash
experimental-features = nix-command flakes
```

3. Clonar o repositório:

```bash
git clone https://github.com/matregnago/intro-hpc.git
cd intro-hpc
```

4. Entrar no ambiente do `nix` de desenvolvimento:

```bash
nix develop
```

## Slides

Os slides estão localizados na pasta `slides`. Eles são escritos em `markdown` e compilados para pdf com o [Marp](https://marp.app) a partir dos seguintes comandos:

```bash
cd slides/part1
marp --pdf part1.md --allow-local-files
```
Isso gera um arquivo chamado `part1.pdf`.

## DOE e Análise dos Dados
Os scripts estão disponíveis na pasta `scripts` e podem ser rodados com o `RScript`. Exemplo:

```bash
Rscript scripts/analysis/plot_traces.r
```