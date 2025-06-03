import argparse
import os
import logging
from utils import run_command


def get_args():
    parser = argparse.ArgumentParser(description="""

    DESCRIPTION

        """, formatter_class= argparse.RawTextHelpFormatter)
    parser.add_argument('--sample_id',
                    type=str,
                    default=None,
                    help='Sample file - e.g SAMPLE001_R1.fastq.gz',
                    required=True)
    parser.add_argument('--fastqc_dir',
                    type=str,
                    default='.',
                    help='Path to dir containing sample FASTQs',
                    required=True)
    parser.add_argument('--multi_dir',
                    type=str,
                    default='.',
                    help='Path to dir containing sample FASTQs',
                    required=True)
    args = parser.parse_args()
    return args

def fastqc(fq, outdir):
    """
    Generate fastqc files
    :param fq: full path for fastq file
    :param outdir: where you want the fastqc files to go
    :return: .html and .zip of fastq
    """  
    try:
        run_command(f'fastqc {fq} --outdir {outdir}')
    except FileNotFoundError as e:
        logging.error(f"An error occurred: {e}")
        raise


def multiqc():
    """
    Generate multiqc files
    :param outdir: results location
    :return: .html and .zip of fastq
    """
    try:
        run_command(f'multiqc -O {multi_dir} {fastqc_indir}')
    except FileNotFoundError as e:
        logging.error(f"An error occurred: {e}")
        raise
