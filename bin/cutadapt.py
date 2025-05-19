#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module provides a portion of the pipeline for trimming adapter contaminated reads with the feature set generated
in feature_extraction.py. It cross-references the read IDs from a cutadapt info file and tags the reads as poor quality based
on the feature thresholds - start pos of adapter and mean quality score.

Author: Nick Wlodychak
Version: 0.1
Date: 2025-05-19
"""

import subprocess
from pathlib import Path
import argparse
import glob
import os
from utils import run_command

def get_args():
    parser = argparse.ArgumentParser(description="""

    DESCRIPTION

        """, formatter_class= argparse.RawTextHelpFormatter)
    parser.add_argument('--platform',
                    type=str,
                    default=None,
                    help='Type of sequencing method. [illumina, nanopore]',
                    required=True)
    parser.add_argument('--target_dir',
                    type=str,
                    default='.',
                    help='Path to dir containing sample FASTQs',
                    required=True) 
    parser.add_argument('--sample_id',
                    type=str,
                    default=None,
                    help='Sample basename e.g. Sample01 - No R1 or R2 required',
                    required=True)
    parser.add_argument('--min_len',
                    type=int,
                    default=100,
                    help='Minimum length to retain - Default 100bp',
                    required=True)
    parser.add_argument('--max_len',
                    type=int,
                    default=None,
                    help='Maximum length to retain.',
                    required=True)
    parser.add_argument('--min_avg_qual',
                    type=int,
                    default=None,
                    help='Minimum Q30 score to retain.',
                    required=True)
    parser.add_argument("--trim_adapters",
                    help='Trim Illumina Adapters',
                    type=bool,
                    required=False,
                    action='store_true')
    parser.add_argument("--no_porechop",
                    help='Skip porechop for Nanopore',
                    type=bool,
                    required=False,
                    action='store_true')
    parser.add_argument('--adapter_fa',
                    type=str,
                    default=None,
                    help='file containing the adatpers in a fasta.',
                    required=False)
    parser.add_argument('--threads',
                    type=int,
                    default=1,
                    help='threads to process',
                    required=False)

    args = parser.parse_args()
    return args

args = get_args()

platform = args.platform.strip().lower()
args.sample_id = args.sample_id.replace(".fastq.gz", "").replace(".fq.gz", "").replace(".fastq", "").replace(".fq", "")

# potential for other platforms
platforms = ['illumina', 'nanopore']


if platform in platforms:
    print(f'Trimming {platform} reads...')
else:
    print(f'Unrecognized platform! Check support {platforms}')
    exit(1)

fq = glob.glob(f"{args.target_dir}/{args.sample_id}.f*q.gz")
fq = fq.sort()

if len(fq) == 2:
    paired = True
elif len(fq) == 1:
    paired = False
else:
    print("No fastqs found!")
    exit(1)

trimmed = []

if platform == 'illumina':
    if args.trim_adapters:
        # paried / illumina trim
        if paired:
            call = (f"cutadapt "
                    f"--cores {args.threads} "
                    f"-a {args.adapter_fa} "
                    f"-A {args.adapter_fa} "
                    f"-o {args.sample_id}_R1.trimmed.fastq.gz "
                    f"-p {args.sample_id}_R2.trimmed.fastq.gz "
                    f"{fq[0]} "
                    f"{fq[1]}")
            run_command(call)
            trimmed.append(f'{args.sample_id}_R1.trimmed.fastq.gz')
            trimmed.append(f'{args.sample_id}_R2.trimmed.fastq.gz')
        
        # single / illumina trim
        else:
            call = (f"cutadapt "
                    f"--cores {args.threads} "
                    f"-a {args.adapter_fa} "
                    f"-A {args.adapter_fa} "
                    f"-o {args.sample_id}.trimmed.fastq.gz "
                    f"{fq[0]}")
            run_command(call)
            trimmed.append(f'{args.sample_id}.trimmed.fastq.gz')
    else:
        # paried / polya trim
        if paired:
            call = (f"cutadapt "
                    f"--cores {args.threads} "
                    f"--poly-a "
                    f"-o {args.sample_id}_R1.trimmed.fastq.gz "
                    f"-p {args.sample_id}_R2.trimmed.fastq.gz "
                    f"{fq[0]} "
                    f"{fq[1]}")
            run_command(call)
            trimmed.append(f'{args.sample_id}_R1.trimmed.fastq.gz')
            trimmed.append(f'{args.sample_id}_R2.trimmed.fastq.gz')
        
        # single / polya trim
        else:
            call = (f"cutadapt "
                    f"--cores {args.threads} "
                    f"--poly-a "
                    f"-o {args.sample_id}.trimmed.fastq.gz "
                    f"{fq[0]}")
            run_command(call)
            trimmed.append(f'{args.sample_id}.trimmed.fastq.gz')

elif platform == 'nanopore' and not args.no_porechop:
    call1 = f"porechop -i {fq[0]} -o {args.sample_id}.chopped.fastq.gz"
    run_command(call1)

    call2 = f"""filtlong --min_length {args.min_len} --min_mean_q {args.min_avg_qual} \
                --target_bases {args.target_bases} {f'{args.sample_id}.chopped.fastq.gz'} | \
                cutadapt -j {args.threads} -m {args.min_len} -M {args.max_len} -o {args.sample_id}.trimmed.fastq.gz """
    run_command(call2)

    trimmed.append(f'{args.sample_id}.trimmed.fastq.gz')


elif platform == 'nanopore' and args.no_porechop:
    call = f"""filtlong --min_length {args.min_len} --min_mean_q {args.min_avg_qual} \
                --target_bases {args.target_bases} {f'{args.sample_id}.chopped.fastq.gz'} | \
                cutadapt -j {args.threads} -m {args.min_len} -M {args.max_len} -o {args.sample_id}.trimmed.fastq.gz """
    run_command(call)
    trimmed.append(f'{args.sample_id}.trimmed.fastq.gz')

print(f"Complete - {fq}")
exit(0)