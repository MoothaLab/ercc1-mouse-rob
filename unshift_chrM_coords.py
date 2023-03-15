#! /usr/bin/env python

import argparse
import os
from plumbum import local
import pysam

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_file')
    parser.add_argument('--out_bam_file')
    parser.add_argument('--chrM_shift', type=int, default=8000)
    parser.add_argument('--sort_nthreads', type=int, default=1)
    parser.add_argument('--sort_mem_gb_per_thread', type=int, default=2)
    args = parser.parse_args()

    if args.out_bam_file is None:
        args.out_bam_file = os.path.splitext(bam_file)[0] + '.chrM_unshifted.unsorted.bam'
    else:
        args.out_bam_file = os.path.splitext(args.out_bam_file)[0] + '.unsorted.bam'

    with pysam.AlignmentFile(args.bam_file) as bam_in:
        chrM_len = bam_in.get_reference_length('chrM')
        with pysam.AlignmentFile(args.out_bam_file, 'wb', template=bam_in) as bam_out:
            #iterate over the bam file
            for record in bam_in:
                #unshift the chrM records
                if record.reference_name == 'chrM':
                    shifted_map_coord = record.reference_start
                    record.reference_start = (shifted_map_coord + args.chrM_shift) % chrM_len
                    if record.is_paired:
                        shifted_mate_coord = record.next_reference_start
                        record.next_reference_start = (shifted_mate_coord + args.chrM_shift) % chrM_len
                    #for now, just skip any read pairs with a mate that maps directly over the
                    #reference genome breakpoint
                    if (record.reference_end > chrM_len
                        or (record.is_paired and (record.next_reference_start + record.query_length) >= chrM_len)):
                        continue
                bam_out.write(record)

    #sort and index unshifted bam file
    sorted_bam = args.out_bam_file.replace('.unsorted.bam', '.bam')
    local['samtools']('sort',
                      '-@', str(args.sort_nthreads),
                      '-m', str(args.sort_mem_gb_per_thread) + 'G',
                      '-o', sorted_bam,
                      args.out_bam_file)
    local['samtools']('index', sorted_bam)

    #clean up intermediate unsorted bam
    local['rm'](args.out_bam_file)
