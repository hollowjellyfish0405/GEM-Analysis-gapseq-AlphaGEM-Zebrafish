"""Launch the upstream rebuild tools in separate, configured environments."""
import argparse
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PROTEOME = ROOT / "input/UP000001570.faa"


def make_alpha_fasta(source, target):
    seen = set()
    with source.open() as reader, target.open("w") as writer:
        for line in reader:
            if line.startswith(">"):
                token = line[1:].split()[0]
                fields = token.split("|")
                accession = fields[1] if len(fields) == 3 else token
                if accession in seen:
                    raise ValueError(f"Repeated FASTA accession {accession}")
                seen.add(accession)
                line = f">{accession}\n"
            writer.write(line)
    if not seen:
        raise ValueError(f"No FASTA entries in {source}")


def execute(command, directory):
    print(shlex.join(map(str, command)), flush=True)
    subprocess.run(command, cwd=directory, check=True)


def gapseq(args):
    executable = shutil.which(args.executable)
    if executable is None:
        raise FileNotFoundError(f"gapseq executable unavailable: {args.executable}")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    fasta = out / "Bsubtilis168.faa"
    shutil.copyfile(args.fasta.resolve(), fasta)
    command = [Path(executable).absolute(), "doall", fasta.name]
    if args.medium:
        command.append(args.medium.resolve(strict=True))
    command.append("Bacteria")
    execute(command, out)
    model = out / "Bsubtilis168.xml"
    if not model.exists():
        raise FileNotFoundError(f"gapseq produced no final SBML at {model}; check SBML support")
    print(model)


def alpha(args):
    repo = args.repo.resolve(strict=True)
    cli = repo / "AlphaGEM.py"
    if not cli.is_file():
        raise FileNotFoundError(cli)
    work = repo / "working" / args.name
    if work.exists():
        raise FileExistsError(f"Choose an unused job name: {work}")
    table = args.table.resolve(strict=True)
    structures = args.structures.resolve(strict=True)
    if not structures.is_dir():
        raise NotADirectoryError(structures)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    fasta = out / "Bsubtilis168_accessions.fasta"
    make_alpha_fasta(args.fasta.resolve(strict=True), fasta)
    command = [sys.executable, cli, "--mode", "structure_alignment", "--refname", "ecoli",
               "--name", args.name, "--fasta", fasta, "--structure", structures,
               "--list", table, "--grothmedium", "min", "--solver", args.solver]
    execute(command, repo)
    print(f"AlphaGEM working directory: {work}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="Use 'reconstruct.py gapseq --help' or 'reconstruct.py alphagem --help' "
               "for the arguments to each workflow.",
    )
    runners = parser.add_subparsers(dest="tool", required=True)
    gap = runners.add_parser("gapseq", help="Rebuild using gapseq")
    gap.add_argument("--executable", default="gapseq")
    gap.add_argument("--fasta", type=Path, default=PROTEOME)
    gap.add_argument("--out", type=Path, default=ROOT / "rebuild/gapseq_168")
    gap.add_argument("--medium", type=Path)
    gap.set_defaults(run=gapseq)
    ag = runners.add_parser("alphagem", help="Rebuild using AlphaGEM structure alignment")
    ag.add_argument("--repo", type=Path, required=True)
    ag.add_argument("--structures", type=Path, required=True)
    ag.add_argument("--table", type=Path, required=True)
    ag.add_argument("--fasta", type=Path, default=PROTEOME)
    ag.add_argument("--name", default="Bsubtilis168")
    ag.add_argument("--out", type=Path, default=ROOT / "rebuild/alphagem_168")
    ag.add_argument("--solver", choices=("auto", "gurobi", "copt"), default="auto")
    ag.set_defaults(run=alpha)
    args = parser.parse_args(sys.argv[1:] or ["--help"])
    args.run(args)


if __name__ == "__main__":
    main()
