#! /usr/bin/env python

import argparse
import os
import pandas
from plumbum import local
import pysam

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script will look at reads overlapping NUMTs and, for any reads with a secondary mtDNA alignment that is equally good (same MAPQ score), will change the mtDNA alignment to be the primary alignment.")
    parser.add_argument('input_bam')
    parser.add_argument('numt_bed')
    parser.add_argument('--output_bam', help='Path to bam file for the output. If not specified, the script will be named based on the input_bam.')
    args = parser.parse_args()

    #identify reads that map to both a NUMT and chrM and make the chrM reads primary
    numt_regions = pandas.read_csv(args.numt_bed, sep='\t', names=['chrom', 'start', 'stop'])
    with pysam.AlignmentFile(args.input_bam, 'rb') as bam_in:
        mito_secondaries = {}
        #catalog chrM alignments that are secondary
        for record in bam_in.fetch(contig='chrM'):
            if record.is_secondary:
                #index them by their name, their reference, and whether they are read1 or read2
                dict_key = (record.query_name, record.reference_name, record.is_read1)
                #in case there are multiple secondary chrM alignments for a given read, keep the
                #  one with the higher aligner-assigned alignment score (SAM tag "AS")
                if (dict_key in mito_secondaries and
                    record.get_tag('AS') <= mito_secondaries[dict_key].get_tag('AS')):
                        continue
                mito_secondaries[dict_key] = record

        #iterate over the NUMT regions and check whether any primary alignments there are for
        # any of the previously-identified reads with secondary chrM alignments
        amended_reads = {}
        for numt_idx in range(numt_regions.shape[0]):
            numt = numt_regions.iloc[numt_idx]
            numt_reads = {(elt.query_name, elt.reference_name, elt.is_read1):elt
                          for elt in bam_in.fetch(contig=numt['chrom'], start=numt['start'], stop=numt['stop'])}
            for key, record in numt_reads.items():
                #if this is a primary alignment and both mate pairs are overlapping the NUMT region, check for
                # a corresponding chrM secondary alignment
                if not record.is_secondary and ((key[0], key[1], not key[2]) in numt_reads):
                    dict_key = (record.query_name, 'chrM', record.is_read1)
                    try:
                        chrM_read = mito_secondaries[dict_key]
                    except KeyError:
                        continue
                    else:
                        #if a primary alignment for a read with a secondary alignment to chrM is discovered
                        # and the alignment score is the same for both alignments, prioritize the chrM alignment
                        if (record.mapping_quality == chrM_read.mapping_quality
                            and record.get_tag('AS') == chrM_read.get_tag('AS')):
                            #modify the chrM_read from mito_secondaries in-place
                            chrM_read.is_secondary = False
                            chrM_read.mapping_quality = 31
                            #change this NUMT record and add to amended reads
                            record.is_secondary = True
                            amended_reads[(record.query_name, record.reference_name, record.is_read1)] = record

    #write out a new BAM file with the NUMT records corrected
    tmp_bam = os.path.splitext(args.input_bam)[0] + '.chrM_prioritize.bam'
    with pysam.AlignmentFile(args.input_bam, 'rb') as bam_in:
        with pysam.AlignmentFile(tmp_bam, 'wb', template=bam_in) as bam_out:
            for record in bam_in:
                record_key = (record.query_name, record.reference_name, record.is_read1)
                #if this read is an amended read, write the corrected alignment record
                if (record_key in amended_reads and
                    record.reference_start == amended_reads[record_key].reference_start):
                    bam_out.write(amended_reads[record_key])
                #if this read is a mito-mapping read, write the (possibly corrected) alignment record
                elif (record_key in mito_secondaries and
                      record.reference_start == mito_secondaries[record_key].reference_start):
                    bam_out.write(mito_secondaries[record_key])
                #write all other reads as they are
                else:
                    bam_out.write(record)

    #rename the BAM according to the user options
    if args.output_bam:
        local['mv'](tmp_bam, args.output_bam)
