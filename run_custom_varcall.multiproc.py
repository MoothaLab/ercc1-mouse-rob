import argparse
from collections import OrderedDict, Counter
import csv
import functools
import glob
import itertools
import multiprocessing as mp
import numpy
import os
import pandas
import plumbum
from plumbum import local
import pysam
import random
from scipy import stats
import sys
import time

def index_bam_func(path):
    return pysam.index(path)

def split_bam_by_cell(bam_path, tmp_dir, cb_tag='read_name', no_umi=False, chrom=None, ncores=1):
    bam_data = pysam.AlignmentFile(bam_path)

    if isinstance(chrom, str):
        chrom = [chrom]
    elif chrom is None:
        chrom = bam_data.references
    elif not isinstance(chrom, (list, tuple)):
        raise Exception('chrom param must be a string, list, or tuple')

    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    out_bams = {}
    tmp_bam_paths = []
    try:
        #split bam file based on the cell barcode
        for chrom_elt in sorted(chrom):
            chrom_bam_data = bam_data.fetch(region=chrom_elt)
            for record in chrom_bam_data:
                if cb_tag == 'read_name':
                    cell_bc = record.query_name.split('_')[-1 if no_umi is True else -2]
                else:
                    cell_bc = record.get_tag(cb_tag)
                try:
                    out_bams[cell_bc].write(record)
                except KeyError:
                    tmp_bam_path = os.path.join(tmp_dir, cell_bc+'.bam')
                    tmp_bam_paths.append(tmp_bam_path)
                    out_bams[cell_bc] = pysam.AlignmentFile(tmp_bam_path, mode='wb',
                                                            template=bam_data)
                    out_bams[cell_bc].write(record)
    #make sure the AlignmentFile objects are properly closed and cleaned up
    finally:
        for bam_obj in out_bams.values():
            bam_obj.close()

    #index the temporary bam files
    if ncores > 1:
        with mp.Pool(processes=ncores) as pool:
            pool.map(index_bam_func, tmp_bam_paths)
    else:
        for bam_path in tmp_bam_paths:
            pysam.index(bam_path)

    #return the list of temporary bams
    return tmp_bam_paths

def call_consensus_sequence(bam_path, ref_path, fetch_region='chrM', proper_pairs=True,
                            mapq_filter=20, min_pos_cov=20, min_mean_bq=23.8):
    ''' Takes a bam file and the reference sequence the reads are aligned to, iterates 
         over the reference positions, and returns a dict of dicts that has chromosome
         name as the top level key, and then the second level keys are reference 
         positions where the consensus bp from the reads differs from the reference. 
         The values in the interior dicts are the consensus base calls.
    '''
    bam_data = pysam.AlignmentFile(bam_path)
    ref_data = pysam.FastaFile(ref_path)
    #always get the full reference chromosome sequence regardless of whether a subregion was requested
    ref_sequence = ref_data.fetch(region=fetch_region.split(':')[0])

    #set alignment filters
    BAM_FUNMAP = 4
    BAM_FSECONDARY = 256
    BAM_FQCFAIL = 512
    BAM_FDUP = 1024
    flag_filter = (BAM_FUNMAP
                   | BAM_FSECONDARY
                   | BAM_FQCFAIL
                   | BAM_FDUP)

    #get the samtools pileup data
    pile_obj = bam_data.pileup(region=fetch_region, truncate=True,
                               stepper='samtools', fasta_file=ref_data,
                               ignore_overlaps=True, flag_filter=flag_filter,
                               ignore_orphans=proper_pairs,
                               min_mapping_quality=mapq_filter, compute_baq=True,
                               redo_baq=False, adjust_capq_threshold=0,
                               max_depth=100000000)

    consensus_ref = {}
    chrom = fetch_region.split(':')[0]
    for pilecol in pile_obj:
        ref_pos = pilecol.reference_pos
        try:
            ref_allele = ref_sequence[ref_pos].upper()
        except IndexError:
            #this shouldn't happen because we retrieved the whole reference chr sequence
            print('Really??', chrom, ref_pos)
            raise

        pos_seq = numpy.array([elt.upper() for elt in pilecol.get_query_sequences()], dtype=object)
        #reads with spliced alignments will be returned for all bases within the spliced region, but these
        # reads have no base information, so get the indices to ignore them
        non_skip_idx = [idx for idx, elt in enumerate(pos_seq) if elt in ['A', 'C', 'G', 'T', 'N']]

        #don't attempt to call a consensus reference allele if there is insufficient coverage
        if len(non_skip_idx) < min_pos_cov:
            continue

        #skip low-quality positions
        pos_mean_bq = numpy.array(pilecol.get_query_qualities())[non_skip_idx].mean()
        if pos_mean_bq < min_mean_bq:
            continue

        #otherwise, what is the consensus base here?
        try:
            pos_consensus_ref_allele = sorted(Counter(pos_seq[non_skip_idx]).items(), key=lambda x:x[1])[-1][0]
        except IndexError:
            if len(pos_seq) > 0:
                print(ref_pos, pos_seq)
                raise
            else:
                continue

        #if the consensus allele is the reference allele, just move on
        if (pos_consensus_ref_allele == ref_allele
            or (ref_allele == 'N' and pos_consensus_ref_allele == '')):
            continue

        #if everything checks out, then map this as a new consensus base
        consensus_ref[ref_pos] = pos_consensus_ref_allele
    return {chrom:consensus_ref}

def compile_bp_dataframe(bam_path, ref_path, fetch_region='chrM',
                         cb_tag='read_name', no_umi=False,
                         proper_pairs=True, mapq_filter=20, min_pos_cov=20,
                         min_mean_bq=23.8, min_allele_strand_cov=3,
                         max_strand_frac=0.7, skip_ref_allele=True,
                         homopolymer_filt=-1, het_threshold=0.2,
                         het_threshold_chrM=0.01, out_file=None, ignore_overlaps=False,
                         all_pos_w_suff_cov=False, bulk_consensus_ref=None):
    ''' The main variant calling function. Iterates over the pysam pileup output to call
         variants that pass a set of heuristic filters and returns the results as a
         pandas DataFrame with a row for each detected variant.
    '''
    bp_data = OrderedDict([('chrom', []),
                           ('pos', []),
                           ('cell', []),
                           ('ref_allele', []),
                           ('allele', []),
                           ('fwd_count', []),
                           ('rev_count', []),
                           ('pos_cov', []),
                           ('fwd_bq_mean', []),
                           ('rev_bq_mean', []),
                           ('homo_poly_len', []),
                           ('het', []),
                           ('mut_type', [])])
    if bulk_consensus_ref is not None:
        bp_data['consensus_ref_allele'] = []

    bam_data = pysam.AlignmentFile(bam_path)
    ref_data = pysam.FastaFile(ref_path)
    ref_sequence = ref_data.fetch(region=fetch_region)

    #set flag filters
    BAM_FUNMAP = 4
    BAM_FSECONDARY = 256
    BAM_FQCFAIL = 512
    BAM_FDUP = 1024
    flag_filter = (BAM_FUNMAP
                   | BAM_FSECONDARY
                   | BAM_FQCFAIL
                   | BAM_FDUP)

    ##NOTE: we don't want samtools to do the ignoring of overlaps because that messes with the strand balance
    # calculation! I manually find overlapping mate pairs and correctly adjust the coverage instead.
    pile_obj = bam_data.pileup(region=fetch_region, truncate=True,
                               stepper='samtools', fasta_file=ref_data,
                               ignore_overlaps=False, flag_filter=flag_filter,
                               ignore_orphans=proper_pairs,
                               min_mapping_quality=mapq_filter, compute_baq=True,
                               redo_baq=False, adjust_capq_threshold=0,
                               max_depth=100000000)

    chrom = fetch_region.split(':')[0]
    base_list = ['A', 'C', 'G', 'T']
    het_filt = het_threshold_chrM if chrom == 'chrM' else het_threshold

    #keep track of how many variants are rejected by the heuristic filters
    filter_stats = {}
    def increment_filter_stats(key):
        try:
            filter_stats[key] += 1
        except KeyError:
            filter_stats[key] = 1

    #main loop -- iterate over the pileup output
    for pilecol in pile_obj:
        ref_pos = pilecol.reference_pos
        try:
            ref_allele = ref_sequence[ref_pos].upper()
        except IndexError:
            #we should never encounter an out of bounds reference coordinate
            print('Really??', chrom, ref_pos)
            raise
        #get the base calls from the reads piled up at this location
        pos_seq = numpy.array([elt.upper() for elt in pilecol.get_query_sequences()], dtype=object)
        #reads with spliced alignments will be returned for all bases within the spliced region, but these
        # reads have no base information, so get the indices to ignore them for calculating coverage
        non_skip_idx = [idx for idx, elt in enumerate(pos_seq) if elt in ['A', 'C', 'G', 'T', 'N']]
        pos_seq = pos_seq[non_skip_idx]
        #get their corresponding base quality scores
        pos_bq = numpy.array(pilecol.get_query_qualities())[non_skip_idx]
        #get the original BAM records
        pos_reads = numpy.array([elt for elt in pilecol.pileups])[non_skip_idx]
        #list the read names
        pos_qnames = numpy.array([elt for elt in pilecol.get_query_names()], dtype=object)[non_skip_idx]
        #check which reads have indels at this location
        pos_indel = numpy.array([elt.indel for elt in pos_reads])
        #get the strand of each read in the pileup
        pos_str = numpy.array([elt.alignment.is_reverse for elt in pos_reads])

        #check whether we are calling variants in single cells or in bulk data
        if cb_tag == 'bulk':
            pos_cells = ['bulk']
        #single cell data and the cell barcode is in the read name
        elif cb_tag == 'read_name':
            pos_cells = numpy.array([elt.split('_')[-1 if no_umi is True else -2] for elt in pos_qnames], dtype=object)
        #single cell data and the cell barcode is in the provided BAM tag name
        else:
            pos_cells = numpy.array([elt.alignment.get_tag(cb_tag) for elt in pos_reads], dtype=object)
        
        pos_cells_set = set(pos_cells)
        #some read pairs have insert lengths that are so short that either mate pair could support
        #a variant call, but we only want to count them once for the coverage
        pos_qnames_counter = Counter(pos_qnames)
        pos_reads_w_mates = numpy.array([1 if pos_qnames_counter[elt] > 1 else 0 for elt in pos_qnames])
        #record the bulk consensus allele here if requested
        # (such alleles differ from the reference, but will be present in all cells because they are
        #  most likely inherited and not somatic)
        if bulk_consensus_ref is None:
            pos_consensus_ref_allele = '-'
        else:
            pos_consensus_ref_allele = bulk_consensus_ref.get(ref_pos, '-')

        #iterate over the reads for each cell and call variants
        for cell_name in pos_cells_set:
            if len(pos_cells_set) == 1:
                cell_bq = pos_bq
                cell_seq = pos_seq
                cell_reads = pos_reads
                cell_indel = pos_indel
                cell_str = pos_str
                cell_reads_w_mates = pos_reads_w_mates
            else:
                cell_idx = numpy.where(pos_cells == cell_name)[0]
                cell_bq = pos_bq[cell_idx]
                cell_seq = pos_seq[cell_idx]
                cell_reads = pos_reads[cell_idx]
                cell_indel = pos_indel[cell_idx]
                cell_str = pos_str[cell_idx]
                cell_reads_w_mates = pos_reads_w_mates[cell_idx]
            #compute the coverage at this position, correcting for mate pairs that overlap
            mates_correction = numpy.sum(cell_reads_w_mates)/2
            cell_cov = cell_reads.shape[0] - mates_correction
            #skip this position in this cell if coverage is too low
            if cell_cov < min_pos_cov:
                increment_filter_stats('min_pos_cov')
                continue

            #handle indel variants
            if numpy.any(cell_indel != 0):
                #check for insertions at this position
                ins_coords = numpy.where(cell_indel > 0)[0]
                if len(ins_coords) > 0:
                    #there could be multiple insertion alleles of different lengths
                    ins_lens = sorted(set(cell_indel[ins_coords]))
                    for ins_len in ins_lens:
                        ins_len_coords = numpy.where(cell_indel == ins_len)[0]
                        #check that the insertion alleles are supported by reads on both strands
                        ins_fwd_count = numpy.sum(~cell_str[ins_len_coords])
                        ins_rev_count = numpy.sum(cell_str[ins_len_coords])
                        strand_frac = max(ins_fwd_count, ins_rev_count)/(ins_fwd_count + ins_rev_count)
                        if (ins_fwd_count < min_allele_strand_cov or ins_rev_count < min_allele_strand_cov) or strand_frac > max_strand_frac:
                            increment_filter_stats('strand_balance')
                            continue
                        #avoid double-counting insertion alleles for mate pairs that overlap
                        ins_mates_correction = numpy.sum(cell_reads_w_mates[ins_len_coords])/2
                        #calculate the variant heteroplasmy
                        het_val = ((ins_fwd_count + ins_rev_count) - ins_mates_correction)/cell_cov
                        if het_val < het_filt:
                            increment_filter_stats('min_het')
                            continue
                        #get the insertion sequence
                        try:
                            ins_types = []
                            for cell_read in cell_reads[ins_len_coords]:
                                if cell_read.query_position is None:
                                    query_pos_idx = cell_read.query_position_or_next
                                else:
                                    query_pos_idx = cell_read.query_position + 1
                                ins_types.append(list(cell_read.alignment.query_sequence[query_pos_idx:query_pos_idx + ins_len]))
                        except TypeError:
                            print(chrom, ref_pos)
                            print(ins_len_coords)
                            print([cell_read.alignment.query_position_or_next for cell_read in cell_reads[ins_len_coords]])
                            raise
                        #make a consensus sequence for insertions of this length, assuming that
                        # any variants are sequencing errors
                        ins_consensus = ''.join(stats.mode(numpy.vstack(ins_types), axis=0).mode.flatten())
                        #record the information about the detected insertion
                        bp_data['chrom'].append(chrom)
                        bp_data['pos'].append(ref_pos)
                        bp_data['cell'].append(cell_name)
                        bp_data['allele'].append(ins_consensus)
                        bp_data['ref_allele'].append('-')
                        bp_data['fwd_count'].append(ins_fwd_count)
                        bp_data['rev_count'].append(ins_rev_count)
                        bp_data['pos_cov'].append(cell_cov)
                        bp_data['fwd_bq_mean'].append(numpy.nan)
                        bp_data['rev_bq_mean'].append(numpy.nan)
                        bp_data['homo_poly_len'].append(numpy.nan)
                        bp_data['het'].append(het_val)
                        bp_data['mut_type'].append('ins')
                        try:
                            bp_data['consensus_ref_allele'].append(pos_consensus_ref_allele)
                        except KeyError:
                            pass

                #handle any deletion variants in similar way to insertions
                del_coords = numpy.where(cell_indel < 0)[0]
                if len(del_coords) > 0:
                    del_lens = sorted(set(cell_indel[del_coords]))
                    for del_len in del_lens:
                        del_len_coords = numpy.where(cell_indel == del_len)[0]
                        #check that deletion alleles are supported by reads on both strands
                        del_fwd_count = numpy.sum(~cell_str[del_len_coords])
                        del_rev_count = numpy.sum(cell_str[del_len_coords])
                        strand_frac = max(del_fwd_count, del_rev_count)/(del_fwd_count + del_rev_count)
                        if (del_fwd_count < min_allele_strand_cov or del_rev_count < min_allele_strand_cov) or strand_frac > max_strand_frac:
                            increment_filter_stats('strand_balance')
                            continue
                        #avoid double-counting deletion alleles for mate pairs that overlap
                        del_mates_correction = numpy.sum(cell_reads_w_mates[del_len_coords])/2
                        #calculate the variant heteroplasmy
                        het_val = ((del_fwd_count + del_rev_count) - del_mates_correction)/cell_cov
                        if het_val < het_filt:
                            increment_filter_stats('min_het')
                            continue
                        #record the data for this deletion allele
                        bp_data['chrom'].append(chrom)
                        #deletions are reported in the pileup output at the base before the first deleted base
                        bp_data['pos'].append(ref_pos + 1)
                        bp_data['cell'].append(cell_name)
                        bp_data['allele'].append('-')
                        bp_data['ref_allele'].append(ref_sequence[(ref_pos + 1):(ref_pos + 1) + (0-del_len)])
                        bp_data['fwd_count'].append(del_fwd_count)
                        bp_data['rev_count'].append(del_rev_count)
                        bp_data['pos_cov'].append(cell_cov)
                        bp_data['fwd_bq_mean'].append(numpy.nan)
                        bp_data['rev_bq_mean'].append(numpy.nan)
                        bp_data['homo_poly_len'].append(numpy.nan)
                        bp_data['het'].append(het_val)
                        bp_data['mut_type'].append('del')
                        try:
                            bp_data['consensus_ref_allele'].append(pos_consensus_ref_allele)
                        except KeyError:
                            pass

            #unless otherwise requested, skip all positions with no difference from the reference
            if not all_pos_w_suff_cov and (len(set(cell_seq)) == 1 and cell_seq[0] == ref_allele):
                continue

            #iterate over each base to detect SNVs of each type
            for allele_base in base_list:
                #ignore skip_ref_allele option if all_pos_w_suff_cov is on
                if not (all_pos_w_suff_cov is True):
                    if skip_ref_allele is True and allele_base == ref_allele:
                        continue

                #if an allele doesn't exist, just continue
                allele_idx = numpy.where(cell_seq == allele_base)[0]
                if len(allele_idx) == 0:
                    continue

                #otherwise, check that the allele is supported by reads from both strands
                allele_fwd_count = numpy.sum(~cell_str[allele_idx])
                allele_rev_count = numpy.sum(cell_str[allele_idx])
                strand_frac = max(allele_fwd_count, allele_rev_count)/(allele_fwd_count + allele_rev_count)
                if (allele_fwd_count < min_allele_strand_cov or allele_rev_count < min_allele_strand_cov) or strand_frac > max_strand_frac:
                    increment_filter_stats('strand_balance')
                    continue
                #make sure not to double-count the allele from overlapping mate pairs
                allele_mates_correction = numpy.sum(cell_reads_w_mates[allele_idx])/2
                #and calculate the heteroplasmy (the cell_cov value was already separately
                #  corrected for overlapping mate pairs)
                het_val = ((allele_fwd_count + allele_rev_count) - allele_mates_correction)/cell_cov
                if het_val < het_filt:
                    increment_filter_stats('min_het')
                    continue

                #check the base calling quality of the variant alleles
                allele_fwd_mean_bq = numpy.mean(cell_bq[allele_idx][~cell_str[allele_idx]])
                allele_rev_mean_bq = numpy.mean(cell_bq[allele_idx][cell_str[allele_idx]])
                if ((allele_fwd_mean_bq < min_mean_bq or allele_rev_mean_bq < min_mean_bq)
                    and het_val > (1/500)):
                    increment_filter_stats('min_mean_bq')
                    break

                #detect whether this allele is at the end of a homopolymer stretch of nucleotides
                dist = 0
                rev_homopoly = True
                fwd_homopoly = True
                while rev_homopoly or fwd_homopoly:
                    dist += 1
                    if rev_homopoly:
                        rev_homopoly = ref_sequence[ref_pos - dist] == allele_base
                    if fwd_homopoly:
                        fwd_homopoly = ref_sequence[ref_pos + dist] == allele_base
                allele_homopoly_len = dist
                if homopolymer_filt > 0 and allele_homopoly_len > homopolymer_filt:
                    increment_filter_stats('homopolymer')
                    continue

                #record the data for this allele
                bp_data['chrom'].append(chrom)
                bp_data['pos'].append(ref_pos)
                bp_data['cell'].append(cell_name)
                bp_data['allele'].append(allele_base)
                bp_data['ref_allele'].append(ref_allele)
                bp_data['fwd_count'].append(allele_fwd_count)
                bp_data['rev_count'].append(allele_rev_count)
                bp_data['pos_cov'].append(cell_cov)
                bp_data['fwd_bq_mean'].append(allele_fwd_mean_bq)
                bp_data['rev_bq_mean'].append(allele_rev_mean_bq)
                bp_data['homo_poly_len'].append(allele_homopoly_len)
                bp_data['het'].append(het_val)
                if ((bulk_consensus_ref is not None)
                    and (allele_base == pos_consensus_ref_allele)
                    and (pos_consensus_ref_allele != ref_allele)):
                    bp_data['mut_type'].append('cons')
                else:
                    bp_data['mut_type'].append('snv' if allele_base != ref_allele else 'ref')
                try:
                    bp_data['consensus_ref_allele'].append(pos_consensus_ref_allele)
                except KeyError:
                    pass
    #compile the results into a dataframe
    var_df = pandas.DataFrame(bp_data)
    #write it to a file
    if out_file is not None:
        var_df.to_csv(out_file, sep='\t', quoting=csv.QUOTE_NONE)
    #write the numbers of alleles filtered out to stdout
    print(chrom)
    print('\n'.join(sorted(['{!s}\t{!s}'.format(*elt) for elt in filter_stats.items()])))
    #return the dataframe
    return var_df

def run_local_realignment(bam_file, ref_fasta, gtf_file=None, nthreads=1, tmpdir=os.getcwd()):
    ''' Local realignment of reads can improve the accuracy of the read alignments, so
         can be useful to do before variant calling. This function takes the BAM file 
         passed in for variant calling and runs ABRA2 for local realignment.
    '''
    out_bam = os.path.splitext(bam_file)[0] + '.abra2.bam'
    abra_log = os.path.join(os.path.dirname(bam_file), 'abra2.log')
    #if a GTF file is provided, analyze with RNA-seq mode
    if gtf_file is not None:
        abra_cmd = local['abra2']['--in', bam_file,
                                  '--out', out_bam,
                                  '--ref', ref_fasta,
                                  '--threads', str(nthreads),
                                  '--tmpdir', tmpdir,
                                  '--junctions', 'bam',
                                  '--gtf', gtf_file,
                                  '--dist', '500000',
                                  '--sua'] > abra_log
    else:
        abra_cmd = local['abra2']['--in', bam_file,
                                  '--out', out_bam,
                                  '--ref', ref_fasta,
                                  '--threads', str(nthreads),
                                  '--tmpdir', tmpdir] > abra_log
    abra_cmd()
    local['samtools']('index', out_bam)
    return out_bam

#This maps chromosome names from UCSC reference format to non-UCSC reference format
CHROM_MAP = {'chrM':'MT', 'chrX':'X', 'chrY':'Y', 'chr1':'1', 'chr2':'2',
             'chr3':'3', 'chr4':'4', 'chr5':'5', 'chr6':'6', 'chr7':'7',
             'chr8':'8', 'chr9':'9', 'chr10':'10', 'chr11':'11',
             'chr12':'12', 'chr13':'13', 'chr14':'14', 'chr15':'15',
             'chr16':'16', 'chr17':'17', 'chr18':'18', 'chr19':'19',
             'chr20':'20', 'chr21':'21', 'chr22':'22'}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_list', help='BAM paths to call variants on. This can be a single BAM path or the path to a file containing a list of BAM files separated by newlines.') 
    parser.add_argument('--ref_fasta', help='Path to the reference FASTA file.')
    parser.add_argument('--ref_gtf', help='Only for running ABRA2 for local realignment of reads')
    parser.add_argument('--vep_species', help='Species name for running VEP. Most likely either \'mus_musculus\' or \'homo_sapiens\'')
    parser.add_argument('--vep_dir_cache', default=os.path.expandvars('$HOME/.vep/'), help='This is the directory that contains the VEP cache.')
    parser.add_argument('--out_root')
    parser.add_argument('--out_base_prefix', default='var_call_result.all')
    parser.add_argument('--chrom', nargs='+', default=['chrM'])
    parser.add_argument('--min_pos_cov', type=int, default=20)
    parser.add_argument('--min_mean_bq', type=int, default=23.8)
    parser.add_argument('--min_allele_strand_cov', type=int, default=3)
    parser.add_argument('--max_strand_frac', type=float, default=0.7)
    parser.add_argument('--het_threshold', type=float, default=0.2)
    parser.add_argument('--het_threshold_chrM', type=float, default=0.01)
    parser.add_argument('--homopolymer_filt', type=int, default=-1, help='Ignore variants occurring in homopolymer runs of this length. Negative values turn this filter off. [default: %(default)s]')
    parser.add_argument('--min_mapq_read_filter', type=int, default=20)
    parser.add_argument('--report_ref_allele', action='store_true', default=False)
    parser.add_argument('--ncores', type=int, default=1)
    parser.add_argument('--gb_mem_per_core', type=int, default=8)
    parser.add_argument('--hrs_per_job', type=int, default=24)
    parser.add_argument('--max_concurrent_jobs', type=int, default=50)
    parser.add_argument('--bams_per_job', type=int)
    parser.add_argument('--start_idx', type=int)
    parser.add_argument('--split_by_cell', action='store_true', default=False)
    parser.add_argument('--skip_completed', action='store_true', default=False)
    parser.add_argument('--min_num_chrM_reads', type=int, default=3000)
    parser.add_argument('--cb_tag', default='read_name', help='Specify a BAM record tag name at which the cell barcode is stored. If the cell barcode is attached to the end of the read name with an underscore, then specify "read_name" for this option (this is the default).')
    parser.add_argument('--no_umi', action='store_true', default=False, help='If there is not UMI sequence, then take the last element of the read name as the cell barcode, not the second to last.')
    parser.add_argument('--local_realignment', action='store_true', default=False, help='If this option is provided, do a local realignment with ABRA2 before calling variants.')
    parser.add_argument('--pipeline_completion_flag', help='Provide a completion flag file name to check for before calling variants. If provided, the variant calling will skip over any BAM files that are not in a directory with a completed flag file.')
    parser.add_argument('--double_count_mate_overlaps', action='store_true', default=False, help='Note: this option is deprecated and currently does nothing.')
    parser.add_argument('--output_all_pos_w_sufficient_cov', action='store_true', default=False)
    parser.add_argument('--call_bulk_consensus_ref', action='store_true', default=False)
    parser.add_argument('--just_call_bulk_consensus_ref', action='store_true', default=False)
    parser.add_argument('--ignore_strand_balance', action='store_true', default=False)
    args = parser.parse_args()

    #get the input bam files, either a single path or a file containing a list of bam paths
    if args.input_list.endswith('.bam'):
        bam_list = [os.path.abspath(args.input_list)]
        args.start_idx = 1
        args.bams_per_job = 1
    else:
        bam_list = [os.path.abspath(elt) for elt in numpy.loadtxt(args.input_list, dtype=object)]

    #split the input list of bams into chunks and start UGER jobs to call variants on each chunk
    if args.start_idx is None:
        num_per_job = int(numpy.ceil(len(bam_list)/args.max_concurrent_jobs))

        conda = os.path.basename(os.environ['CONDA_PREFIX'])
        var_call_cmd = 'conda activate {conda!s}; python {script!s} --input_list={input_list!s} --ref_fasta={ref_fasta!s} --ref_gtf={ref_gtf!s} --vep_species={vep_species!s} --vep_dir_cache={vep_dir_cache!s} --out_root={out_root!s} --out_base_prefix={out_prefix!s} --chrom {chrom!s} --min_pos_cov={min_pos_cov!s} --min_mean_bq={min_mean_bq!s} --min_allele_strand_cov={min_allele_strand_cov!s} --max_strand_frac={max_strand_frac!s} --het_threshold={het_threshold!s} --het_threshold_chrM={het_chrM!s} --homopolymer_filt={hom_filt!s} --min_num_chrM_reads={min_chrM_reads!s} --min_mapq_read_filter={min_mapq!s} --cb_tag={cb_tag!s} --ncores={ncores!s} --gb_mem_per_core={mem_per_core!s} --bams_per_job={num_per_job!s} --start_idx=$SGE_TASK_ID'.format(conda=conda, script=__file__, input_list=args.input_list, ref_fasta=args.ref_fasta, ref_gtf=args.ref_gtf, vep_species=args.vep_species, vep_dir_cache=args.vep_dir_cache, out_root=args.out_root, out_prefix=args.out_base_prefix, chrom=' '.join(args.chrom), min_pos_cov=args.min_pos_cov, min_mean_bq=args.min_mean_bq, min_allele_strand_cov=args.min_allele_strand_cov, max_strand_frac=args.max_strand_frac, het_threshold=args.het_threshold, het_chrM=args.het_threshold_chrM, hom_filt=args.homopolymer_filt, min_chrM_reads=args.min_num_chrM_reads, min_mapq=args.min_mapq_read_filter, cb_tag=args.cb_tag, ncores=args.ncores, mem_per_core=args.gb_mem_per_core, num_per_job=num_per_job)
        if args.pipeline_completion_flag is not None:
            var_call_cmd += ' --pipeline_completion_flag={!s}'.format(args.pipeline_completion_flag)
        if args.skip_completed is True:
            var_call_cmd += ' --skip_completed'
        if args.no_umi is True:
            var_call_cmd += ' --no_umi'
        if args.split_by_cell is True:
            var_call_cmd += ' --split_by_cell'
        if args.local_realignment is True:
            var_call_cmd += ' --local_realignment'
        if args.report_ref_allele is True:
            var_call_cmd += ' --report_ref_allele'
        if args.double_count_mate_overlaps is True:
            var_call_cmd += ' --double_count_mate_overlaps'
        if args.output_all_pos_w_sufficient_cov is True:
            var_call_cmd += ' --output_all_pos_w_sufficient_cov'
        if args.call_bulk_consensus_ref is True:
            var_call_cmd += ' --call_bulk_consensus_ref'
        if args.just_call_bulk_consensus_ref is True:
            var_call_cmd += ' --just_call_bulk_consensus_ref'
        if args.ignore_strand_balance is True:
            var_call_cmd += '--ignore_strand_balance'
        #run multi-core jobs
        if args.ncores > 1:
            qsub_cmd = local['qsub']['-cwd', '-j', 'y', '-b', 'y', '-N', 'run_var_call', '-o', args.out_root, '-l', 'h_rt={!s}:00:00'.format(args.hrs_per_job), '-l', 'h_vmem={!s}G'.format(args.gb_mem_per_core), '-pe', 'smp', str(args.ncores), '-binding', 'linear:{!s}'.format(args.ncores), '-t', '{!s}-{!s}:{!s}'.format(1, len(bam_list), num_per_job), var_call_cmd]
        #run single-core jobs
        else:
            qsub_cmd = local['qsub']['-cwd', '-j', 'y', '-b', 'y', '-N', 'run_var_call', '-o', args.out_root, '-l', 'h_rt={!s}:00:00'.format(args.hrs_per_job), '-l', 'h_vmem={!s}G'.format(args.gb_mem_per_core), '-t', '{!s}-{!s}:{!s}'.format(1, len(bam_list), num_per_job), var_call_cmd]
        qsub_cmd()
    #otherwise, this input should be processed locally, so iterate over the bam input(s) and call variants
    else:
        start_idx = args.start_idx - 1
        stop_idx = min(start_idx + args.bams_per_job, len(bam_list))
        for dir_idx in range(start_idx, stop_idx):
            bam_path = bam_list[dir_idx]
            bam_dir = os.path.dirname(bam_path)

            #skip this directory if the alignment pipeline hasn't completed successfully
            if args.pipeline_completion_flag and not os.path.isfile(os.path.join(bam_dir, args.pipeline_completion_flag)):
                continue

            #skip this directory if there are not enough chrM reads
            if 'all' in args.chrom:
                bam_data = pysam.AlignmentFile(bam_path)
                args.chrom = [elt for elt in bam_data.references if elt in CHROM_MAP]
            if 'chrM' in args.chrom:
                bam_read_count = int(local['samtools']('view', '-c', bam_path, 'chrM'))
                if bam_read_count < args.min_num_chrM_reads:
                    continue

            #print the bam directory that is currently being processed
            print(bam_dir)

            #do local realignment if requested
            if args.local_realignment is True:
                bam_path = run_local_realignment(bam_path, args.ref_fasta, gtf_file=args.ref_gtf, nthreads=args.ncores, tmpdir=bam_dir)

            #set output file paths
            out_prefix = '{!s}.{!s}'.format(os.path.splitext(os.path.basename(bam_path))[0], args.out_base_prefix)
            vep_output = os.path.join(bam_dir, '{!s}.vep_output.txt'.format(out_prefix))
            #if requested, skip any files that already exist
            if args.skip_completed is True and os.path.isfile(vep_output):
                continue

            #set some variant calling filter parameters
            if args.ignore_strand_balance is True:
                min_allele_strand_cov = 0
                max_strand_frac = 1.0
            else:
                min_allele_strand_cov = args.min_allele_strand_cov
                max_strand_frac = args.max_strand_frac

            #call the bulk consensus bases if requested
            if args.call_bulk_consensus_ref is True or args.just_call_bulk_consensus_ref is True:
                #use a process pool if we have access to multiple cores
                if args.ncores > 1:
                    pool_args = [((bam_path, args.ref_fasta), {'fetch_region':chrom_elt,
                                                               'proper_pairs':True,
                                                               'mapq_filter':args.min_mapq_read_filter,
                                                               'min_pos_cov':args.min_pos_cov,
                                                               'min_mean_bq':args.min_mean_bq,
                                                               'min_allele_strand_cov':min_allele_strand_cov,
                                                               'max_strand_frac':max_strand_frac,
                                                               'ignore_overlaps':not args.double_count_mate_overlaps})
                                  for chrom_elt in args.chrom]
                    def apply_pool_args(args, kwargs):
                        return call_consensus_sequence(*args, **kwargs)
                    with mp.Pool(processes=args.ncores) as pool:
                        consensus_ref_list = pool.starmap(apply_pool_args, pool_args)
                #otherwise, use a single process
                else:
                    consensus_ref_list = [call_consensus_sequence(bam_path, args.ref_fasta,
                                                                  fetch_region=chrom_elt, proper_pairs=True,
                                                                  mapq_filter=args.min_mapq_read_filter,
                                                                  min_pos_cov=args.min_pos_cov,
                                                                  min_mean_bq=args.min_mean_bq,
                                                                  min_allele_strand_cov=min_allele_strand_cov,
                                                                  max_strand_frac=max_strand_frac,
                                                                  ignore_overlaps=not args.double_count_mate_overlaps)
                                          for chrom_elt in args.chrom]
                #merge the consensus-calling results for each chromosome
                consensus_ref = {}
                for elt in consensus_ref_list:
                    consensus_ref.update(elt)
                #output the changed coordinates and a modified reference FASTA if requested
                if args.just_call_bulk_consensus_ref is True:
                    out_fasta = os.path.join(bam_dir, '{!s}.fasta'.format(out_prefix))
                    out_meta = os.path.splitext(out_fasta)[0] + '.consensus_ref_variants.tab'
                    ref_data = pysam.FastaFile(args.ref_fasta)
                    with open(out_fasta, 'w') as out, open(out_meta, 'w') as out2:
                        out2.write('chrom\tpos\tref_allele\tvar_allele\n')
                        for chrom, chrom_vars in consensus_ref.items():
                            #always get the full reference chromosome sequence regardless
                            #of whether a subregion was requested
                            chrom_seq = list(ref_data.fetch(region=chrom))
                            for ref_idx, var_allele in chrom_vars.items():
                                out2.write('{!s}\t{!s}\t{!s}\t{!s}\n'.format(chrom, ref_idx, chrom_seq[ref_idx], var_allele))
                                chrom_seq[ref_idx] = var_allele
                            out.write('>{!s}\n{!s}\n'.format(chrom, ''.join(chrom_seq)))
                    sys.exit()
            else:
                consensus_ref = None

            #call variants
            var_call_out = os.path.join(bam_dir, '{!s}.tab'.format(out_prefix))
            #use a process pool if we have access to more than one core
            if args.ncores > 1:
                if args.split_by_cell is True:
                    tmp_dir = os.path.join(bam_dir, '{!s}.tmp_cell_split'.format(out_prefix))
                    tmp_bam_paths = split_bam_by_cell(bam_path, tmp_dir, cb_tag=args.cb_tag, no_umi=args.no_umi, chrom=args.chrom, ncores=args.ncores)
                    pool_args = itertools.product(tmp_bam_paths, args.chrom)
                else:
                    pool_args = itertools.product([bam_path], args.chrom)
                pool_args = [((bpath, args.ref_fasta), {'fetch_region':chrom_elt,
                                                        'cb_tag':args.cb_tag,
                                                        'no_umi':args.no_umi,
                                                        'proper_pairs':True,
                                                        'mapq_filter':args.min_mapq_read_filter,
                                                        'min_pos_cov':args.min_pos_cov,
                                                        'min_mean_bq':args.min_mean_bq,
                                                        'min_allele_strand_cov':min_allele_strand_cov,
                                                        'max_strand_frac':max_strand_frac,
                                                        'skip_ref_allele':not args.report_ref_allele,
                                                        'homopolymer_filt':args.homopolymer_filt,
                                                        'het_threshold':args.het_threshold,
                                                        'het_threshold_chrM':args.het_threshold_chrM,
                                                        'out_file':None,
                                                        'ignore_overlaps':not args.double_count_mate_overlaps,
                                                        'all_pos_w_suff_cov':args.output_all_pos_w_sufficient_cov,
                                                        'bulk_consensus_ref':None if consensus_ref is None else consensus_ref.get(chrom_elt, {})})
                             for bpath, chrom_elt in pool_args]
                def apply_pool_args(args, kwargs):
                    return compile_bp_dataframe(*args, **kwargs)
                with mp.Pool(processes=args.ncores) as pool:
                    var_info = pool.starmap(apply_pool_args, pool_args)
            #otherwise call variants in a single process
            else:
                var_info = []
                for chrom in args.chrom:
                    var_df = compile_bp_dataframe(bam_path, args.ref_fasta,
                                                  fetch_region=chrom,
                                                  cb_tag=args.cb_tag,
                                                  no_umi=args.no_umi,
                                                  proper_pairs=True,
                                                  mapq_filter=args.min_mapq_read_filter,
                                                  min_pos_cov=args.min_pos_cov,
                                                  min_mean_bq=args.min_mean_bq,
                                                  min_allele_strand_cov=min_allele_strand_cov,
                                                  max_strand_frac=max_strand_frac,
                                                  skip_ref_allele=not args.report_ref_allele,
                                                  homopolymer_filt=args.homopolymer_filt,
                                                  het_threshold=args.het_threshold,
                                                  het_threshold_chrM=args.het_threshold_chrM,
                                                  out_file=None,
                                                  ignore_overlaps=not args.double_count_mate_overlaps,
                                                  all_pos_w_suff_cov=args.output_all_pos_w_sufficient_cov,
                                                  bulk_consensus_ref=None if consensus_ref is None else consensus_ref.get(chrom, {}))
                    var_info.append(var_df)
            #compile all of the variant calls for each chromosome and all cells into a single dataframe
            # and save it to a file
            var_df = pandas.concat(var_info).reset_index()
            var_df.to_csv(var_call_out, sep='\t', index=False, quoting=csv.QUOTE_NONE)

            #skip VEP if there are no variants that make it through filtering
            if var_df.shape[0] == 0:
                continue

            #otherwise, write the variants in VEP input format
            vep_input = os.path.join(bam_dir, '{!s}.vep_input.txt'.format(out_prefix))
            unique_vars = var_df[var_df['mut_type'].isin(['ins', 'del', 'snv', 'cons'])][['chrom', 'pos', 'ref_allele', 'allele']].drop_duplicates().sort_values(by=['chrom', 'pos', 'allele']).reset_index()
            with open(vep_input, 'w') as out:
                for idx in range(unique_vars.shape[0]):
                    cur_pos = int(unique_vars['pos'][idx] + 1)
                    #insertion
                    if unique_vars['ref_allele'][idx] == '-':
                        out.write('{chrom!s}\t{start!s}\t{stop!s}\t{ref!s}/{var!s}\t1\n'.
                              format(chrom=CHROM_MAP[unique_vars['chrom'][idx]],
                                     start=cur_pos + 1,
                                     stop=cur_pos,
                                     ref=unique_vars['ref_allele'][idx],
                                     var=unique_vars['allele'][idx]))
                    #deletion
                    elif unique_vars['allele'][idx] == '-':
                        out.write('{chrom!s}\t{start!s}\t{stop!s}\t{ref!s}/{var!s}\t1\n'.
                              format(chrom=CHROM_MAP[unique_vars['chrom'][idx]],
                                     start=cur_pos,
                                     stop=cur_pos + len(unique_vars['ref_allele'][idx]) - 1,
                                     ref=unique_vars['ref_allele'][idx],
                                     var=unique_vars['allele'][idx]))
                    #SNV
                    else:
                        out.write('{chrom!s}\t{start!s}\t{stop!s}\t{ref!s}/{var!s}\t1\n'.
                                  format(chrom=CHROM_MAP[unique_vars['chrom'][idx]],
                                         start=cur_pos,
                                         stop=cur_pos,
                                         ref=unique_vars['ref_allele'][idx],
                                         var=unique_vars['allele'][idx]))

            #run VEP on the variant calls
            vep_cmd = local['vep']['--appris',
                                   '--biotype',
                                   '--buffer_size', '500',
                                   '--check_existing',
                                   '--distance', '0',
                                   '--mane',
                                   '--regulatory',
                                   '--sift', 'b',
                                   '--species', args.vep_species,
                                   '--symbol',
                                   '--transcript_version',
                                   '--tsl',
                                   '--offline',
                                   '--dir_cache', args.vep_cache_dir,
                                   '--input_file', vep_input,
                                   '--output_file', vep_output,
                                   '--force_overwrite']
            try:
                vep_cmd()
            except plumbum.commands.processes.ProcessExecutionError as err:
                print(repr(err))
