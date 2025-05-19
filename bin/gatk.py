#!/usr/bin/env python3

import os
import argparse
import subprocess
import sys

def get_args():
    parser = argparse.ArgumentParser(description="""

    DESCRIPTION

        """, formatter_class= argparse.RawTextHelpFormatter)
    parser.add_argument('--threads',
                    type=int,
                    default=1,
                    help='threads to process',
                    required=False)

    args = parser.parse_args()
    return args

args = get_args()


class GATKWrapper:
    def __init__(self, gatk_path=None, java_options=None):
        """Initialize the GATK wrapper with paths and default options."""
        self.gatk_path = gatk_path if gatk_path else "gatk"
        self.java_options = java_options if java_options else "-Xmx8G"
        
    def run_command(self, command):
        """Execute a shell command and return the process."""
        print(f"Running: {command}")
        process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, text=True)
        
        if process.returncode != 0:
            print(f"Error executing command: {command}")
            print(f"STDERR: {process.stderr}")
            return False
        return True
        
    def haplotype_caller(self, reference, input_bam, output_vcf, extra_args=None):
        """Run the HaplotypeCaller tool on a BAM file."""
        cmd = f"{self.gatk_path} --java-options '{self.java_options}' HaplotypeCaller "
        cmd += f"-R {reference} -I {input_bam} -O {output_vcf} "
        
        if extra_args:
            cmd += extra_args
            
        return self.run_command(cmd)
    
    def base_recalibrator(self, reference, input_bam, known_sites, output_table, extra_args=None):
        """Run the BaseRecalibrator tool on a BAM file."""
        cmd = f"{self.gatk_path} --java-options '{self.java_options}' BaseRecalibrator "
        cmd += f"-R {reference} -I {input_bam} --known-sites {known_sites} -O {output_table} "
        
        if extra_args:
            cmd += extra_args
            
        return self.run_command(cmd)
    
    def apply_bqsr(self, reference, input_bam, recal_table, output_bam, extra_args=None):
        """Apply base quality score recalibration to a BAM file."""
        cmd = f"{self.gatk_path} --java-options '{self.java_options}' ApplyBQSR "
        cmd += f"-R {reference} -I {input_bam} --bqsr-recal-file {recal_table} -O {output_bam} "
        
        if extra_args:
            cmd += extra_args
            
        return self.run_command(cmd)
    
    def variant_filtration(self, reference, input_vcf, output_vcf, filter_expression, filter_name, extra_args=None):
        """Filter variants according to specified expressions."""
        cmd = f"{self.gatk_path} --java-options '{self.java_options}' VariantFiltration "
        cmd += f"-R {reference} -V {input_vcf} -O {output_vcf} "
        cmd += f"--filter-expression '{filter_expression}' --filter-name '{filter_name}' "
        
        if extra_args:
            cmd += extra_args
            
        return self.run_command(cmd)

def main():
    parser = argparse.ArgumentParser(description="GATK Wrapper for BAM processing")
    parser.add_argument("--gatk-path", help="Path to GATK executable", default="gatk")
    parser.add_argument("--java-options", help="Java options for GATK", default="-Xmx8G")
    
    subparsers = parser.add_subparsers(dest="command", help="GATK command to run")
    
    # HaplotypeCaller command
    hc_parser = subparsers.add_parser("HaplotypeCaller", help="Run HaplotypeCaller")
    hc_parser.add_argument("-R", "--reference", required=True, help="Reference genome")
    hc_parser.add_argument("-I", "--input", required=True, help="Input BAM file")
    hc_parser.add_argument("-O", "--output", required=True, help="Output VCF file")
    hc_parser.add_argument("--extra-args", help="Additional arguments for HaplotypeCaller")
    
    # BaseRecalibrator command
    br_parser = subparsers.add_parser("BaseRecalibrator", help="Run BaseRecalibrator")
    br_parser.add_argument("-R", "--reference", required=True, help="Reference genome")
    br_parser.add_argument("-I", "--input", required=True, help="Input BAM file")
    br_parser.add_argument("--known-sites", required=True, help="Known sites VCF")
    br_parser.add_argument("-O", "--output", required=True, help="Output recalibration table")
    br_parser.add_argument("--extra-args", help="Additional arguments for BaseRecalibrator")
    
    # ApplyBQSR command
    bqsr_parser = subparsers.add_parser("ApplyBQSR", help="Run ApplyBQSR")
    bqsr_parser.add_argument("-R", "--reference", required=True, help="Reference genome")
    bqsr_parser.add_argument("-I", "--input", required=True, help="Input BAM file")
    bqsr_parser.add_argument("--bqsr-recal-file", required=True, help="Recalibration table")
    bqsr_parser.add_argument("-O", "--output", required=True, help="Output BAM file")
    bqsr_parser.add_argument("--extra-args", help="Additional arguments for ApplyBQSR")
    
    # VariantFiltration command
    vf_parser = subparsers.add_parser("VariantFiltration", help="Run VariantFiltration")
    vf_parser.add_argument("-R", "--reference", required=True, help="Reference genome")
    vf_parser.add_argument("-V", "--input", required=True, help="Input VCF file")
    vf_parser.add_argument("-O", "--output", required=True, help="Output filtered VCF file")
    vf_parser.add_argument("--filter-expression", required=True, help="Filter expression")
    vf_parser.add_argument("--filter-name", required=True, help="Filter name")
    vf_parser.add_argument("--extra-args", help="Additional arguments for VariantFiltration")
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    gatk = GATKWrapper(args.gatk_path, args.java_options)
    
    if args.command == "HaplotypeCaller":
        gatk.haplotype_caller(args.reference, args.input, args.output, args.extra_args)
    elif args.command == "BaseRecalibrator":
        gatk.base_recalibrator(args.reference, args.input, args.known_sites, args.output, args.extra_args)
    elif args.command == "ApplyBQSR":
        gatk.apply_bqsr(args.reference, args.input, args.bqsr_recal_file, args.output, args.extra_args)
    elif args.command == "VariantFiltration":
        gatk.variant_filtration(args.reference, args.input, args.output, args.filter_expression, args.filter_name, args.extra_args)

if __name__ == "__main__":
    main()
