#! /usr/bin/env python

import argparse
import importlib
import os
from plumbum import local
import pypiper
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('read1_fastq')
    parser.add_argument('read2_fastq')
    parser.add_argument('--out_dir', default=os.getcwd())
    parser.add_argument('--sample_name')
    parser.add_argument('--aln_reference')
    parser.add_argument('--aln_reference_sizes')
    parser.add_argument('--ref_fasta', help='FASTA file of reference for the local realignment to use.')
    parser.add_argument('--aln_nthreads', type=int, default=1)
    parser.add_argument('--pipeline_jvmmemavail', type=int, default=8)
    parser.add_argument('--vcf_for_bqsr')
    parser.add_argument('--max_concurrent_subjobs', type=int, default=50)
    parser.add_argument('--pipeline_recover', action='store_true', default=False)
    args = parser.parse_args()

    code_dir = os.path.dirname(__file__)
    if args.sample_name is None:
        input_basename = os.path.splitext(os.path.basename(args.read1_fastq).rstrip('.gz'))[0]
        args.sample_name = input_basename
    out_prefix = os.path.join(args.out_dir, args.sample_name)
    pm = pypiper.PipelineManager(name="bulk_DNA", outfolder=args.out_dir,
                                 recover=args.pipeline_recover)


##########
# TRIM ILLUMINA ADAPTERS
##########
    pm.timestamp('### Trim Illumina adapters')
    adapters_loc = os.path.join(os.environ['CONDA_PREFIX'], 'share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa')
    trimmed_fq1 = args.read1_fastq.replace('.fastq', '.trimmed.fastq')
    trimmed_fq2 = args.read2_fastq.replace('.fastq', '.trimmed.fastq')
    trim_args = ['PE', args.read1_fastq,
                 args.read2_fastq,
                 trimmed_fq1,
                 trimmed_fq1.replace('.fastq', '.unpaired.fastq'),
                 trimmed_fq2,
                 trimmed_fq2.replace('.fastq', '.unpaired.fastq'),
                 'ILLUMINACLIP:{!s}:2:30:10:1:true'.format(adapters_loc),
                 'TRAILING:3', 'SLIDINGWINDOW:4:10', 'MINLEN:20']
    trim_cmd = local['trimmomatic'].__getitem__(trim_args)
    pm.run(str(trim_cmd), target=trimmed_fq1)


##########
# ALIGN READS WITH BOWTIE2
##########
    pm.timestamp('### Align reads')
    out_bam = out_prefix + '.filt.q10.bam'
    out_bam_norg = out_prefix + '.filt.q10.no_rg.bam'
    aln_args = ['-p {!s}'.format(args.aln_nthreads),
                '-X 2000',
                '-x {!s}'.format(args.aln_reference),
                '-1', trimmed_fq1,
                '-2', trimmed_fq2]
    aln_cmd = local['bowtie2'].__getitem__(aln_args) | local['samtools']['view', '-bhSu', '-f 3', '-F 12', '-q 10', '-'] | local['samtools']['sort', '-T', args.out_dir, '-o', out_bam_norg, '-@ {!s}'.format(args.aln_nthreads)]
    addrg_cmd = local['gatk']['--java-options', '-Xmx{!s}G'.format(args.pipeline_jvmmemavail),
                              'AddOrReplaceReadGroups',
                              '-I', out_bam_norg,
                              '-O', out_bam,
                              '--RGLB={!s}'.format(args.sample_name),
                              '--RGPL=ILLUMINA',
                              '--RGPU={!s}'.format(args.sample_name),
                              '--RGSM={!s}'.format(args.sample_name),
                              '--VALIDATION_STRINGENCY=LENIENT'] # Necessary because for some reason bowtie2 doesn't put its version number in the @PG header line...
    index_cmd = local['samtools']['index', out_bam]
    pm.run([str(aln_cmd), str(addrg_cmd), str(index_cmd)], target=out_bam+'.bai')
    #clean up the non-deduped bam at pipeline completion
    pm.clean_add(out_bam_norg)
    pm.clean_add(out_bam)
    pm.clean_add(out_bam + '.bai')


##########
# LOCAL REALIGNMENT WITH ABRA2
##########
    pm.timestamp('### Run local realignment with ABRA2')

    #NOTE: ABRA2 will re-sort the chromosomes in the BAM file in lexicographic order. If your sequence
    # lines in the header are in another order, this will cause problems because ABRA2 copies the
    # header from the input bam but then writes the BAM records for the chromosomes in lexicographic
    # order. This violates the BAM specification because the reads from each chromosome must appear
    # in the same order as the chromosomes appear in the header...

    local_realign_bam = os.path.splitext(out_bam)[0] + '.abra2.bam'
    abra_log = os.path.splitext(out_bam)[0] + '.abra2.log'
    abra2_args = ['--in', out_bam,
                  '--out', local_realign_bam,
                  '--ref', args.ref_fasta,
                  '--threads', str(args.aln_nthreads),
                  '--tmpdir', os.path.dirname(out_prefix)]

    abra_cmd = local['abra2'].__getitem__(abra2_args) > abra_log

    index_cmd = local['samtools']['index', local_realign_bam]
    pm.run([str(abra_cmd), str(index_cmd)], target=local_realign_bam+'.bai')


##########
# BASE QUALITY SCORE RECALIBRATION -- if requested
##########
    if args.vcf_for_bqsr is not None:
        pm.timestamp('### Run GATK to recalibrate base quality scores')

        #NOTE: BQSR requires read groups to run, so check for that and add them if they are not found

        bqsr1_table = out_prefix + '.bqsr1_recal.table'
        bqsr1_cmd = local['gatk']['--java-options',
                                  '-Xmx{!s}G'.format(args.pipeline_jvmmemavail),
                                  'BaseRecalibrator',
                                  '-I', local_realign_bam,
                                  '-R', args.ref_fasta,
                                  '--known-sites', args.vcf_for_bqsr,
                                  '-O', bqsr1_table]

        recal_bam = local_realign_bam.replace('.bam', '.bqsr.bam')
        apply_bqsr_cmd = local['gatk']['--java-options',
                                       '-Xmx{!s}G'.format(args.pipeline_jvmmemavail),
                                       'ApplyBQSR',
                                       '-I', local_realign_bam,
                                       '-R', args.ref_fasta,
                                       '--bqsr-recal-file', bqsr1_table,
                                       '-O', recal_bam]

        bqsr2_table = out_prefix + '.bqsr2_recal.table'
        bqsr2_cmd = local['gatk']['--java-options',
                                  '-Xmx{!s}G'.format(args.pipeline_jvmmemavail),
                                  'BaseRecalibrator',
                                  '-I', recal_bam,
                                  '-R', args.ref_fasta,
                                  '--known-sites', args.vcf_for_bqsr,
                                  '-O', bqsr2_table]

        plots_pdf = out_prefix + '.bqsr_plots.pdf'
        plot_cmd = local['gatk']['--java-options',
                                 '-Xmx{!s}G'.format(args.pipeline_jvmmemavail),
                                 'AnalyzeCovariates',
                                 '-before', bqsr1_table,
                                 '-after', bqsr2_table,
                                 '-plots', plots_pdf]
        pm.run([str(bqsr1_cmd), str(apply_bqsr_cmd),
                str(bqsr2_cmd), str(plot_cmd)], target=plots_pdf)
    else:
        recal_bam = local_realign_bam


##########
# MARK DUPLICATE READS
##########
    pm.timestamp('### Remove duplicates')

    dedup_bam = os.path.splitext(recal_bam)[0] + '.dedup.bam'
    dedup_bam_index = dedup_bam + '.bai'
    dedup_metrics = dedup_bam.replace('.bam', '.metrics.txt')

    dedup_cmd = local['gatk']['--java-options',
                              '-Xmx{!s}G'.format(args.pipeline_jvmmemavail),
                              'MarkDuplicates',
                              '-I', local_realign_bam,
                              '-O', dedup_bam,
                              '--METRICS_FILE={!s}'.format(dedup_metrics)]
    index_cmd = local['samtools']['index', dedup_bam]

    pm.run([str(dedup_cmd), str(index_cmd)], target=dedup_bam_index)


##########
# COMPUTE COVERAGE
##########
    pm.timestamp('### Compute coverage')

    #read coverage
    out_bdg = out_prefix + '.coverage.bdg'
    out_bw = out_prefix + '.coverage.bw'
    cov_cmd = (local['samtools']['view', '-h', '-f', '3', '-F', '3852', '-q', '30', dedup_bam]
               | local['bedtools']['genomecov',
                                   '-ibam', 'stdin',
                                   '-bg',
                                   '-split']
               | local['grep']['-v', '_']
               | local['sort']['-k1,1', '-k2,2n'] > out_bdg)
    bw_cmd = local['bedGraphToBigWig'][out_bdg, args.aln_reference_sizes, out_bw]
    gzip_cmd = local['gzip'][out_bdg]
    pm.run([str(cov_cmd), str(bw_cmd), str(gzip_cmd)], target=out_bdg + '.gz')
    pm.clean_add(out_bdg + '.gz')

    #paired end fragment coverage
    out_bdg = out_prefix + '.coverage.pc.bdg'
    out_bw = out_prefix + '.coverage.pc.bw'
    cov_cmd = (local['samtools']['view', '-h', '-f', '3', '-F', '3852', '-q', '30', dedup_bam]
               | local['bedtools']['genomecov',
                                   '-ibam', 'stdin',
                                   '-bg',
                                   '-pc',
                                   '-split']
               | local['grep']['-v', '_']
               | local['sort']['-k1,1', '-k2,2n'] > out_bdg)
    bw_cmd = local['bedGraphToBigWig'][out_bdg, args.aln_reference_sizes, out_bw]
    gzip_cmd = local['gzip'][out_bdg]
    pm.run([str(cov_cmd), str(bw_cmd), str(gzip_cmd)], target=out_bdg + '.gz')
    pm.clean_add(out_bdg + '.gz')


##########
# END OF PIPELINE
##########
    pm.stop_pipeline()
