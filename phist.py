# 
# PHIST
# Copyright (C) 2021, A. Zielezinski, S. Deorowicz, and A. Gudys
# https://github.com/refresh-bio/PHIST
# 
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this
# program. If not, see https://www.gnu.org/licenses/.
# 
import argparse
import multiprocessing
import subprocess
from pathlib import Path
import platform


__version__ = '1.0.0'

def get_arguments():
    desc = f'PHIST predicts hosts from phage (meta)genomic data'
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('virus_dir', metavar='virus_dir',
                   help='Input directory with virus FASTA files (gzipped or not)')
    p.add_argument('host_dir', metavar='host_dir',
                   help='Input directory with host FASTA files (gzipped or not)')
    p.add_argument('out_table', metavar='out_table',
                   help='Output CSV file with common k-mers table')
    p.add_argument('out_predictions', metavar='out_predictions',
                   help='Output CSV file with hosts predictions')
    p.add_argument('-k', '--k', dest='k', type=int,
                   default=25, help='k-mer length [default =  %(default)s]')
    p.add_argument('-t', '--t', dest='num_threads', type=int,
                   default=multiprocessing.cpu_count(),
                   help='Number of threads [default = %(default)s]')
    p.add_argument('--version', action='version',
                   version='PHIST v' + __version__,
                   help="Show tool's version number and exit")
    args = p.parse_args()
    return args


if __name__ == '__main__':   
    print(
        f'PHIST  v{__version__}\n',
		'A. Zielezinski, S. Deorowicz, A. Gudys (c) 2021\n\n')
    
    PHIST_DIR = Path(__file__).resolve().parent

    if platform.system() == "Windows":
        kmer_exec = PHIST_DIR.joinpath('kmer-db', 'src', 'x64', 'Release', 'kmer-db.exe')
        util_exec = PHIST_DIR.joinpath('utils', 'x64', 'Release', 'phist.exe')
    else:
        kmer_exec = PHIST_DIR.joinpath('kmer-db', 'kmer-db')
        util_exec = PHIST_DIR.joinpath('utils', 'phist')

    args = get_arguments()

    vdir_path = Path(args.virus_dir)
    hdir_path = Path(args.host_dir)
    outtable_path = Path(args.out_table)
    outpred_path = Path(args.out_predictions)

    # Paths to temp files.
    temp_path = outtable_path.resolve().parent
    vlst_path = temp_path.joinpath('vir.list')
    hlst_path = temp_path.joinpath('host.list')
    db_path = temp_path.joinpath('vir.db')

    # Create vir.list and host.list.
    for lst_path, dir_path in [(vlst_path, vdir_path), (hlst_path, hdir_path)]:
        oh = open(lst_path, 'w')
        for f in sorted(dir_path.iterdir()):
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
        f'{outtable_path}',
    ]
    subprocess.run(cmd)

    # Remove temp files.
    vlst_path.unlink()
    hlst_path.unlink()
    db_path.unlink()  
    
    # Postprocessing
    cmd = [
        f'{util_exec}',
        f'{outtable_path}',
        f'{outpred_path}',
    ]
    subprocess.run(cmd)
