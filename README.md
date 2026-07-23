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
Os scripts estão na pasta `scripts`:

- `scripts/run.sh` — executa o design (`DESIGN_FILE`) para um runtime (chamado pelos jobs em `slurm/`);
- `scripts/process_data/` — converte os traces brutos em parquets StarVZ (`trace_phase1.sh`, `parsec_phase1.sh` e os conversores `parsec_*_to_parquet.r`);
- `scripts/analysis/` — gera as figuras a partir dos resultados/parquets. Exemplo:

```bash
Rscript scripts/analysis/plot_n_size_compare.r
scripts/process_data/parsec_phase1.sh data/<job>/runs/00*_parsec_*
```