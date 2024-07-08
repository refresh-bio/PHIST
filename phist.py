#!/usr/bin/env python3
"""A tool to predict prokaryotic hosts for phage (meta)genomic sequences. 
PHIST links viruses to hosts based on the number of k-mers shared between 
their sequences.

Copyright (C) 2021 A. Zielezinski, S. Deorowicz, and A. Gudys
https://github.com/refresh-bio/PHIST
"""

from __future__ import annotations
import argparse
import multiprocessing
import platform
from pathlib import Path
import subprocess
import sys

__version__ = '1.2.1'


def get_parser() -> argparse.ArgumentParser:
    desc = f'PHIST predicts hosts from phage (meta)genomic data'
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('virus_path',
                   help='Input FASTA file or directory with FASTA files (plain or gzip)')
    p.add_argument('host_dir', metavar='host_dir',
                   help='Input directory w/ host FASTA files (plain or gzip)')
    p.add_argument('out_dir', metavar='out_dir', nargs='+',
                   help='Output directory (will be created if it does not exist)')
    p.add_argument('-k', dest='k', type=int,
                   default=25, help='k-mer length [default =  %(default)s]')
    p.add_argument('-t', dest='num_threads', type=int,
                   default=multiprocessing.cpu_count(),
                   help='Number of threads [default = %(default)s]')
    p.add_argument('--keep_temp', action="store_true",
                   help='Keep temporary kmer-db files [%(default)s]')
    p.add_argument('--version', action='version',
                   version=__version__,
                   help="Show tool's version number and exit")

    # Display help if the script is run without arguments.
    if len(sys.argv[1:]) == 0:
        p.print_help()
        p.exit()
    return p


def validate_args(parser: argparse.ArgumentParser) -> argparse.Namespace:
    """Validates arguments provided by the user.

    Returns:
        Arguments provided by the users.
    Raises:
        argparse.ArgumentParser.error if arguments are invalid.
    """
    args = parser.parse_args()
    
    # Validate k-mer length
    if args.k < 3 or args.k > 30:
        parser.error(f'K-mer length should be in range 3-30.')

    # Validate virus input
    v_path = Path(args.virus_path)
    if not v_path.exists():
        parser.error(f'Virus input does not exist: {v_path}')

    # Validate host input
    hdir_path = Path(args.host_dir)
    if not hdir_path.exists() or not hdir_path.is_dir():
        parser.error(f'Input host directory does not exist: {hdir_path}')

    args.v_path = v_path
    args.hdir_path = hdir_path

    # Validate output files
    if len(args.out_dir) > 1:
        args.outtable_path = Path(args.out_dir[0])
        args.outpred_path = Path(args.out_dir[1])
        args.out_dir = args.outtable_path.parent
    else:
        args.out_dir = Path(args.out_dir[0])
        args.outtable_path = args.out_dir / 'common_kmers.csv'
        args.outpred_path = args.out_dir / 'predictions.csv'
    return args


if __name__ == '__main__':
    
    PHIST_DIR = Path(__file__).resolve().parent

    if platform.system() == "Windows":
        kmer_exec = PHIST_DIR.joinpath('kmer-db', 'src', 'x64', 'Release', 'kmer-db.exe')
        util_exec = PHIST_DIR.joinpath('utils', 'x64', 'Release', 'phist.exe')
    else:
        kmer_exec = PHIST_DIR.joinpath('kmer-db', 'kmer-db')
        util_exec = PHIST_DIR.joinpath('utils', 'phist')

    parser = get_parser()
    args = validate_args(parser)

    print(
        f'PHIST  v{__version__}\n',
        'A. Zielezinski, S. Deorowicz, A. Gudys (c) 2021\n\n')

    v_path = args.v_path
    hdir_path = args.hdir_path
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Paths to temp files
    vlst_path = out_dir / 'virus.list'
    hlst_path = out_dir / 'host.list'
    db_path = out_dir / 'virus.kdb'

    # Create virus.lst
    with open(vlst_path, 'w') as oh:
        if v_path.is_dir():
            for f in sorted(v_path.rglob('*')):
                if f.is_file():
                    oh.write(f'{f}\n')
        else:
            oh.write(f'{v_path}')

    # Create host.list.
    with open(hlst_path, 'w') as oh:
        for f in sorted(hdir_path.rglob('*')):
            if f.is_file():
                oh.write(f"{f}\n")
        oh.close()

    # Kmer-db build
    cmd = [
        f'{kmer_exec}',
        'build',
        '-k',
        f'{args.k}',
        '-t',
        f'{args.num_threads}',
        f'{vlst_path}',
        f'{db_path}',
    ]
    if v_path.is_file():
        cmd.insert(6, '-multisample-fasta')
    subprocess.run(cmd)

    # Kmer-db new2all
    cmd = [
        f'{kmer_exec}',
        'new2all',
        '-sparse',
        '-t',
        f'{args.num_threads}',
        f'{db_path}',
        f'{hlst_path}',
        f'{args.outtable_path}',
    ]
    subprocess.run(cmd)

    # Remove temp files.
    if not args.keep_temp:
        vlst_path.unlink()
        hlst_path.unlink()
        db_path.unlink()  
    
    # Postprocessing
    cmd = [
        f'{util_exec}',
        f'{args.outtable_path}',
        f'{args.outpred_path}',
    ]
    subprocess.run(cmd)
