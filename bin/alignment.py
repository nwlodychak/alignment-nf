#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module provides the portion of the pipeline to perform alignment against a reference genome. We can provide different types of aligners
to mostly universialize the process.

Author: Nick Wlodychak
Version: 0.1
Date: 2025-05-19
"""

import subprocess
import polars as pl
from pathlib import Path
import argparse
import glob
import os

def get_args():
    parser = argparse.ArgumentParser(description="""

    DESCRIPTION

        """, formatter_class= argparse.RawTextHelpFormatter)
    parser.add_argument('--platform',
                    type=str,
                    default=None,
                    help='Type of sequencing method. [ILLUMINA, NANOPORE, PACBIO]',
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
    parser.add_argument('--threads',
                    type=int,
                    default=1,
                    help='threads to process',
                    required=False)
    parser.add_argument('--aligner',
                    type=str,
                    default=None,
                    help='Type of aligner to use. [BWA, BOWTIE, STAR]',
                    required=True)

    args = parser.parse_args()
    return args


args = get_args()

fq = glob.glob(f"{args.trimming_directory}/{args.sample_id}*.f*q.gz")
fq = fq.sort()

aligner = args.aligner.strip().upper()

if aligner in ['BWA', 'BOWTIE', 'STAR']:
    print(f'Aligning reads by {aligner}...')
else:
    print('Unrecognized aligner!')
    exit(1)

if len(fq) == 2:
    paired = True
elif len(fq) == 1:
    paired = False
else:
    print("No fastqs found!")
    exit(1)

if aligner == "BWA":
    if paired:
        call = f'bwa mem {args.ref_genome} {args.sample_id}_R1.trimmed.fastq.gz {args.sample_id}_R2.trimmed.fastq.gz > {args.sample_id}.sam -t {args.threads} samtools view -bS {args.sample_id}.sam > {args.sample_id}.bam samtools flagstat {args.sample_id}.bam'
        subprocess.call(call, shell = True)
    else:
        call = f'bwa mem {args.ref_genome} {args.sample_id}.trimmed.fastq.gz > {args.sample_id}.sam samtools view -bS {args.sample_id}.sam t {args.threads} > {args.sample_id}.bam samtools flagstat {args.sample_id}.bam'
        subprocess.call(call, shell = True)

if aligner == "BOWTIE":
    if paired:
        call = f'bwa mem {args.ref_genome} {args.sample_id}_R1.trimmed.fastq.gz {args.sample_id}_R2.trimmed.fastq.gz > {args.sample_id}.sam -t {args.threads} samtools view -bS {args.sample_id}.sam > {args.sample_id}.bam samtools flagstat {args.sample_id}.bam'
        subprocess.call(call, shell = True)
    else:
        call = f'bwa mem {args.ref_genome} {args.sample_id}.trimmed.fastq.gz > {args.sample_id}.sam samtools view -bS {args.sample_id}.sam t {args.threads} > {args.sample_id}.bam samtools flagstat {args.sample_id}.bam'
        subprocess.call(call, shell = True)