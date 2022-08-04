#!/usr/bin/env python3
__version__ = "0.0.1"
__author__ = "Sangjin Lee"

# modules
import basic.util
import basic.caller
import basic.vcflib 
import basic.mutlib 
from basic.parse_args import parse_args

def main():
    parser, options = parse_args(program_version=__version__)
    if options.sub == "call":  # call somatic substitutions
        basic.util.check_num_threads(options.threads)
        print("himut is calling somatic single base substitutions")
        basic.caller.call_somatic_substitutions(
            options.bam, # input # bam_file
            options.vcf, # germline mutations
            options.region, # target contigs/scaffolds/chromosomes
            options.region_list, # target contigs/scaffolds/chromosomes fofn
            options.min_mapq, # int: 0-60
            options.min_trim, # float: 0.01 - 0.1
            options.min_sequence_identity, # float: blast sequence identity
            options.min_hq_base_proportion, # float: proportion of BQ=93 bases
            options.min_alignment_proportion, # float: query alignment length / query read length
            options.common_snps, # common snps
            options.panel_of_normals, # panel of normals
            options.min_bq, # minimum base quality score: int
            options.mismatch_window, # mismatch window size
            options.max_mismatch_count, # maximum number of mismatches within a window
            options.min_ref_count, # number of reads supporting the reference base
            options.min_alt_count, # number of reads supporting the alterantive base
            options.min_hap_count, # number of reads supporting h0 and h1 haplotype
            options.threads, # maxminum number of threads
            options.phase, # bool
            __version__, # str
            options.out, # output # himut vcf file
        )
    else:
        print("The subcommand does not exist!\n")
        parser.print_help()
        parser.exit()

if __name__ == "__main__":
    main()
    basic.util.exit()
