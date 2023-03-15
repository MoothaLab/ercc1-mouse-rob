#! /usr/bin/env python

import argparse
import os
from plumbum import local
import pysam

def pick_best_alignment(bam_records):
    ''' Takes a collection of pysam BAM records that represent multiple alignments of the
        same read and chooses the one with the best MAPQ score.
    '''
    if bam_records[0].is_paired:
        if len(bam_records) == 2:
            return bam_records
        else:
            #sort the read 1 alignments by mapping quality
            mate1_list = sorted([elt for elt in bam_records if elt.is_read1 and elt.is_proper_pair], key=lambda x:x.mapping_quality, reverse=True)
            if len(mate1_list) == 0:
                raise Exception('No reads mapped in proper pairs.')
            #find the matching mate pair alignment for the best read 1 alignment
            for mate1 in mate1_list:
                mate2_list = [elt for elt in bam_records if elt.is_read2 and mate1.next_reference_start == elt.reference_start]
                if len(mate2_list) == 0:
                    continue
                else:
                    mate2 = mate2_list[0]
                    return [mate1, mate2]
            else:
                raise Exception('Could not find mapped pair.')
    else:
        if len(bam_records) == 1:
            return bam_records
        else:
            #for single end reads, just sort by mapping quality and return the alignment with the best one
            best_mapq = sorted([elt for elt in bam_records], key=lambda x:x.mapping_quality)[-1]
            return [best_mapq]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_file1')
    parser.add_argument('bam_file2')
    parser.add_argument('output_bam')
    parser.add_argument('--java_maxmem_gb', type=int, default=8)
    args = parser.parse_args()

    #merge the bam files and sort by read name
    tmp_namesorted = os.path.splitext(args.output_bam)[0] + '.tmp_merged_namesorted.bam'
    local['gatk']('--java-options', '-Xmx{!s}G -Djava.io.tmpdir={!s}'.format(args.java_maxmem_gb, os.path.dirname(args.output_bam)),
                  'MergeSamFiles',
                  '--INPUT', args.bam_file1, '--INPUT', args.bam_file2,
                  '--OUTPUT', tmp_namesorted,
                  '--SORT_ORDER', 'queryname',
                  '--VALIDATION_STRINGENCY=LENIENT',
                  '--TMP_DIR', os.path.dirname(args.output_bam))
    #remove any identical alignments from the merged set of reads
    tmp_namesorted_dedup = tmp_namesorted.replace('.bam', '.uniq.bam')
    cmd = local['samtools']['view', '-h', tmp_namesorted] | local['uniq'] | local['samtools']['view', '-Sbh', '-'] > tmp_namesorted_dedup
    cmd()

    #iterate over the merged, sorted, and deduplicated bam file and choose the best primary alignment
    # for each read
    out_bam_namesorted = tmp_namesorted_dedup.replace('.bam', '.reconciled.bam')
    with pysam.AlignmentFile(tmp_namesorted_dedup) as bam_in:
        with pysam.AlignmentFile(out_bam_namesorted, 'wb', template=bam_in) as bam_out:
            to_compare = []
            prev_qname = None
            for record in bam_in:
                #only compare primary alignments
                if record.is_secondary or record.is_supplementary or record.is_unmapped:
                    continue
                #if we have collected all alignments for the current read, choose the best one and write it out
                if prev_qname and prev_qname != record.query_name:
                    if len(to_compare) > 1:
                        try:
                            for aln in pick_best_alignment(to_compare):
                                bam_out.write(aln)
                        except Exception as err:
                            print('Warning: problem reconciling reads. Exception: {!s} Reads: {!s}'
                                  .format(err, '\n'.join([str(elt) for elt in to_compare])))
                    else:
                        print('Only one record for the following read. Ignoring. {!s}'.format(to_compare[0]))
                    to_compare = [record]
                    prev_qname = record.query_name
                #otherwise, keep collecting alignments for the current read
                else:
                    to_compare.append(record)
                    prev_qname = record.query_name
            #handle the last read in the file
            else:
                if len(to_compare) > 0:
                    for aln in pick_best_alignment(to_compare):
                        bam_out.write(aln)

    #sort and index the new bam file
    local['samtools']('sort', '-o', args.output_bam, out_bam_namesorted)
    local['samtools']('index', args.output_bam)
    #clean up intermediate files
    local['rm'](tmp_namesorted, tmp_namesorted_dedup, out_bam_namesorted)
