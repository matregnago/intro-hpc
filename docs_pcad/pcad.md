# Como rodar os experimentos no PCAD

## Configuração SSH

Configure um alias no ssh para facilitar o acesso ao PCAD em `.ssh/config`:

```bash
Host pcad
	User usuario
	Hostname gppd-hpc.inf.ufrgs.br
	ForwardAgent no
```

Agora você pode usar `ssh pcad` ao invés de `ssh usuario@gppd-hpc.inf.ufrgs.br`.

## Configuração do Nix no PCAD

No PCAD, o Nix deve ser instalado com o [`nixw`](nixw). Coloque `nixw` em `$HOME/bin` e adicione `$HOME/bin` ao seu `PATH`. Todo comando nix deve ser precedido por `nixw`.

## Setup inicial

Clone o repositório para a sua `$HOME`:

```bash
git clone https://github.com/matregnago/intro-hpc.git
cd intro-hpc
```

## Lançar os experimentos

Submeta os jobs ao SLURM:

```bash
bash slurm/block_size.slurm
bash slurm/n_size.slurm
bash slurm/traces.slurm
```