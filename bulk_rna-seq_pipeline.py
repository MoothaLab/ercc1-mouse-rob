#! /usr/bin/env python

import argparse
import os
from plumbum import local
import pypiper

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq1')
    parser.add_argument('fastq2', nargs='?', default=None)
    parser.add_argument('--out_dir', default=os.getcwd())
    parser.add_argument('--aln_rdna_reference')
    parser.add_argument('--aln_reference')
    parser.add_argument('--aln_reference_fasta')
    parser.add_argument('--aln_reference_sizes')
    parser.add_argument('--aln_nthreads', type=int, default=1)
    parser.add_argument('--aln_annot_gtf')
    parser.add_argument('--local_realignment', action='store_true', default=False, help='If specified, run a local ralignment of mapped reads using ABRA2 before doing deduplication and base quality recalibration.')
    parser.add_argument('--numts_bed')
    parser.add_argument('--vcf_for_bqsr')
    parser.add_argument('--featureCounts_t', default='exon')
    parser.add_argument('--max_concurrent_subjobs', type=int, default=50)
    parser.add_argument('--sample_name')
    parser.add_argument('--pipeline_recover', action='store_true', default=False)
    parser.add_argument('--pipeline_no_cleanup', action='store_true', default=False)
    parser.add_argument('--pipeline_run_all', action='store_true', default=False)
    parser.add_argument('--pipeline_jvmmemavail', default='30g')
    args = parser.parse_args()

    if args.vcf_for_bqsr is not None and args.aln_reference_fasta is None:
        raise Exception('If --vcf_for_bqsr is supplied for running Base Quality Score Recalibration with GATK, the user must also supply a path to an indexed fasta file for the reference genome with --aln_reference_fasta.')

    code_dir = os.path.dirname(__file__)
    if args.sample_name is None:
        input_basename = os.path.splitext(os.path.basename(args.fastq1).rstrip('.gz'))[0]
        args.sample_name = input_basename
    out_prefix = os.path.join(args.out_dir, args.sample_name)
    pm = pypiper.PipelineManager(name="bulk_RNA-seq", outfolder=args.out_dir,
                                 recover=args.pipeline_recover,
                                 dirty=args.pipeline_no_cleanup,
                                 new_start=args.pipeline_run_all)

##########
# DEPLETE rDNA IF REQUESTED
##########
    if args.aln_rdna_reference is not None:
        pm.timestamp('### Deplete rDNA reads')

        out_prefix_rdna = out_prefix + '_rDNA'
        out_bam = out_prefix_rdna + '_Aligned.out.bam'
        nordna_fq1 = out_prefix + '_R1.rDNA_depleted.fastq.gz'
        #single end run
        if args.fastq2 is None:
            cmd_args = ['--runThreadN {!s}'.format(args.aln_nthreads),
                        '--genomeDir {!s}'.format(args.aln_rdna_reference),
                        '--readFilesIn', args.fastq1,
                        '--readFilesCommand {!s}'.format('zcat' if args.fastq1.endswith('.gz') else 'cat'),
                        '--outFileNamePrefix {!s}_'.format(out_prefix_rdna),
                        '--outSAMattrRGline ID:{0!s}\tPL:ILLUMINA\tSM:{0!s}'.format(args.sample_name),
                        '--runRNGseed 777',
                        '--outSAMtype BAM Unsorted', #note: with unmapped reads included, STAR has sorting issues
                        '--outSAMunmapped Within']
            aln_cmd = local['STAR'].__getitem__(cmd_args)
            extract_cmd = local['samtools']['fastq', '-n', '-1', nordna_fq1, '-F', '2', out_bam]
            cleanup_cmd = local['rm']['-r', out_prefix + '__STARtmp']
        #paired end run
        else:
            cmd_args = ['--runThreadN {!s}'.format(args.aln_nthreads),
                        '--genomeDir {!s}'.format(args.aln_rdna_reference),
                        '--readFilesIn', args.fastq1, args.fastq2,
                        '--readFilesCommand {!s}'.format('zcat' if args.fastq1.endswith('.gz') else 'cat'),
                        '--outFileNamePrefix {!s}_'.format(out_prefix_rdna),
                        '--outSAMattrRGline ID:{0!s}\tPL:ILLUMINA\tSM:{0!s}'.format(args.sample_name),
                        '--runRNGseed 777',
                        '--outSAMtype BAM Unsorted', #note: with unmapped reads included, STAR has sorting issues
                        '--outSAMunmapped Within']
            nordna_fq2 = out_prefix + '_R2.rDNA_depleted.fastq.gz'
            aln_cmd = local['STAR'].__getitem__(cmd_args)
            extract_cmd = local['samtools']['fastq', '-n', '-1', nordna_fq1, '-2', nordna_fq2, '-F', '2', out_bam]
            cleanup_cmd = local['rm']['-r', out_prefix + '__STARtmp']
        #align the reads to the rDNA reference, index the bam file, convert the unmapped reads to fastq,
        # and clean up tmp files.
        pm.run([str(aln_cmd), str(index_cmd), str(extract_cmd), str(cleanup_cmd)], target=nordna_fq1)
        args.fastq1 = nordna_fq1
        args.fastq2 = nordna_fq2 if args.fastq2 is not None else None
        pm.clean_add(out_bam)
        pm.clean_add(out_bam + '.bai')

##########
# MAIN ALIGNMENT STEP
##########
    pm.timestamp('### Align reads')
    #single end run
    if args.fastq2 is None:
        cmd_args = ['--runThreadN {!s}'.format(args.aln_nthreads),
                    '--genomeDir {!s}'.format(args.aln_reference),
                    '--readFilesIn', args.fastq1,
                    '--readFilesCommand {!s}'.format('zcat' if args.fastq1.endswith('.gz') else 'cat'),
                    '--outFileNamePrefix {!s}_'.format(out_prefix),
                    '--outSAMattrRGline ID:{0!s}\tPL:ILLUMINA\tSM:{0!s}'.format(args.sample_name),
                    '--runRNGseed 777',
                    '--outSAMtype BAM SortedByCoordinate',
                    '--outSAMattributes NH HI NM MD AS XS nM',
                    '--twopassMode', 'Basic']
    #paired end run
    else:
        cmd_args = ['--runThreadN {!s}'.format(args.aln_nthreads),
                    '--genomeDir {!s}'.format(args.aln_reference),
                    '--readFilesIn', args.fastq1, args.fastq2,
                    '--readFilesCommand {!s}'.format('zcat' if args.fastq1.endswith('.gz') else 'cat'),
                    '--outFileNamePrefix {!s}_'.format(out_prefix),
                    '--outSAMattrRGline ID:{0!s}\tPL:ILLUMINA\tSM:{0!s}'.format(args.sample_name),
                    '--runRNGseed 777',
                    '--outSAMtype BAM SortedByCoordinate',
                    '--outSAMattributes NH HI NM MD AS XS nM',
                    '--twopassMode', 'Basic']
    out_bam = out_prefix + '_Aligned.sortedByCoord.out.bam'
    aln_cmd = local['STAR'].__getitem__(cmd_args)
    index_cmd = local['samtools']['index', out_bam]
    cleanup_cmd = local['rm']['-r', out_prefix + '__STARtmp']
    # align the reads, index the resulting BAM, and clean up intermediate files
    pm.run([str(aln_cmd), str(index_cmd), str(cleanup_cmd)], target=out_bam+'.bai')
    pm.clean_add(out_bam)
    pm.clean_add(out_bam + '.bai')


##########
# PRIORITIZE chrM MAPPING OVER NUMTS -- for any reads mapping equally well to the nuclear genome
#   and to the mitochondrial genome, set the primary alignment to be the one to chrM. The assumption
#   here is that most NUMTs are not highly expressed but chrM is.
##########
    pm.timestamp('### Prioritize chrM read over NUMTs')

    prioritize_script = os.path.join(os.path.dirname(__file__), 'prioritize_chrM_over_numts.py')
    prioritized_bam = os.path.splitext(out_bam)[0] + '.prioritized.bam'
    prioritize_cmd = local['python'][prioritize_script,
                                     out_bam,
                                     args.numts_bed,
                                     '--output_bam={!s}'.format(prioritized_bam)]
    index_cmd = local['samtools']['index', prioritized_bam]
    #run the priortizing script and index the resulting BAM file
    pm.run([str(prioritize_cmd), str(index_cmd)], target=prioritized_bam+'.bai')
    pm.clean_add(prioritized_bam)
    pm.clean_add(prioritized_bam + '.bai')


##########
# LOCAL REALIGNMENT WITH ABRA2 IF REQUESTED
##########
    if args.local_realignment is True:
        pm.timestamp('### Run local realignment with ABRA2')

        local_realign_bam = os.path.splitext(prioritized_bam)[0] + '.abra2.bam'
        abra_log = os.path.splitext(prioritized_bam)[0] + '.abra2.log'
        abra_cmd = local['abra2']['--in', prioritized_bam,
                                  '--out', local_realign_bam,
                                  '--ref', args.aln_reference_fasta,
                                  '--threads', str(args.aln_nthreads),
                                  '--tmpdir', os.path.dirname(out_prefix),
                                  '--junctions', 'bam',
                                  '--gtf', args.aln_annot_gtf,
                                  '--dist', '500000',
                                  '--sua'] > abra_log
        index_cmd = local['samtools']['index', local_realign_bam]
        #run the local realignment and index the resulting BAM file
        pm.run([str(abra_cmd), str(index_cmd)], target=local_realign_bam+'.bai')
        pm.clean_add(local_realign_bam)
        pm.clean_add(local_realign_bam + '.bai')
        prioritized_bam = local_realign_bam


##########
# REMOVE DUPLICATE READS
##########
    pm.timestamp('### Remove duplicates')

    dedup_bam = os.path.splitext(prioritized_bam)[0] + '.dedup.bam'
    dedup_bam_index = dedup_bam + '.bai'
    dedup_metrics = dedup_bam.replace('.bam', '.metrics.txt')

    dedup_cmd = local['gatk']['--java-options',
                              '-Xmx{!s}'.format(args.pipeline_jvmmemavail),
                              'MarkDuplicates',
                              '-I', prioritized_bam,
                              '-O', dedup_bam,
                              '--METRICS_FILE={!s}'.format(dedup_metrics),
                              '--REMOVE_DUPLICATES=true']
    index_cmd = local['samtools']['index', dedup_bam]
    #run the deduplication and index the resulting BAM file
    pm.run([str(dedup_cmd), str(index_cmd)], target=dedup_bam_index)
    pm.clean_add(dedup_bam)
    pm.clean_add(dedup_bam_index)


##########
# RUN BASE QUALITY SCORE RECALIBRATION IF REQUESTED
##########
    if args.vcf_for_bqsr is not None:
        pm.timestamp('### Run GATK to recalibrate base quality scores')

        #train the recalibration model
        bqsr1_table = out_prefix + '.bqsr1_recal.table'
        bqsr1_cmd = local['gatk']['--java-options',
                                  '-Xmx{!s}'.format(args.pipeline_jvmmemavail),
                                  'BaseRecalibrator',
                                  '-I', dedup_bam,
                                  '-R', args.aln_reference_fasta,
                                  '--known-sites', args.vcf_for_bqsr,
                                  '-O', bqsr1_table]
        #apply it to our aligned reads
        dedup_recal = dedup_bam.replace('.bam', '.bqsr.bam')
        apply_bqsr_cmd = local['gatk']['--java-options',
                                       '-Xmx{!s}'.format(args.pipeline_jvmmemavail),
                                       'ApplyBQSR',
                                       '-I', dedup_bam,
                                       '-R', args.aln_reference_fasta,
                                       '--bqsr-recal-file', bqsr1_table,
                                       '-O', dedup_recal]

        mv_cmd1 = local['mv'][dedup_recal, dedup_bam]
        mv_cmd2 = local['mv'][dedup_recal.replace('.bam', '.bai'), dedup_bam_index]

        #train another recalibration model on the corrected reads to see how we did
        bqsr2_table = out_prefix + '.bqsr2_recal.table'
        bqsr2_cmd = local['gatk']['--java-options',
                                  '-Xmx{!s}'.format(args.pipeline_jvmmemavail),
                                  'BaseRecalibrator',
                                  '-I', dedup_bam,
                                  '-R', args.aln_reference_fasta,
                                  '--known-sites', args.vcf_for_bqsr,
                                  '-O', bqsr2_table]

        #generate a report about the base quality score recalibration
        plots_pdf = out_prefix + '.bqsr_plots.pdf'
        plot_cmd = local['gatk']['--java-options',
                                 '-Xmx{!s}'.format(args.pipeline_jvmmemavail),
                                 'AnalyzeCovariates',
                                 '-before', bqsr1_table,
                                 '-after', bqsr2_table,
                                 '-plots', plots_pdf]
        pm.run([str(bqsr1_cmd), str(apply_bqsr_cmd),
                str(mv_cmd1), str(mv_cmd2),
                str(bqsr2_cmd), str(plot_cmd)], target=plots_pdf)


##########
# ANNOTATE READS WITH GENES
##########
    pm.timestamp('### Annotate reads with genes')

    annot_bam = dedup_bam + '.featureCounts.bam'
    annot_bam_index = annot_bam + '.bai'
    annot_cmd = local['featureCounts']['-a {!s}'.format(args.aln_annot_gtf),
                                       '-t {!s}'.format(args.featureCounts_t),
                                       '-o {!s}.featureCounts.csv'.format(out_prefix),
                                       '-R BAM',
                                       '-T {!s}'.format(args.aln_nthreads),
                                       '-p', #force it to use paired end setting and count fragments
                                       '-M', #count multimapping reads
                                       '--primary', #but only the primary ones
                                       dedup_bam]
    sort_cmd = local['samtools']['sort', '-@ {!s}'.format(args.aln_nthreads),
                                 '-o', annot_bam + '.sorted', annot_bam]
    mv_cmd = local['mv'][annot_bam + '.sorted', annot_bam]
    mv_cmd2 = local['mv'][out_prefix, out_prefix + '.feature_counts.counts.tsv']
    index_cmd = local['samtools']['index', annot_bam]
    #run featureCounts to annotate reads and count reads per gene, then sort and index the resulting BAM
    pm.run([str(annot_cmd), str(sort_cmd), str(mv_cmd), str(index_cmd)],
           target=annot_bam_index)


##########
# COMPUTE THE COVERAGE ACROSS THE GENOME
##########
    pm.timestamp('### Compute coverage')

    #read coverage
    out_bdg = out_prefix + '.coverage.bdg'
    out_bw = out_prefix + '.coverage.bw'
    #single end run
    if args.fastq2 is None:
        cov_cmd = (local['samtools']['view', '-h', '-F', '3852', '-q', '30', annot_bam]
                   | local['bedtools']['genomecov',
                                       '-ibam', 'stdin',
                                       '-bg',
                                       '-split']
                   | local['grep']['-v', '_']
                   | local['sort']['-k1,1', '-k2,2n'] > out_bdg)
    #paired end run
    else:
        cov_cmd = (local['samtools']['view', '-h', '-f', '3', '-F', '3852', '-q', '30', annot_bam]
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

    #paired end fragment coverage (not just read coverage)
    if args.fastq2 is not None:
        out_bdg = out_prefix + '.coverage.pc.bdg'
        out_bw = out_prefix + '.coverage.pc.bw'
        cov_cmd = (local['samtools']['view', '-h', '-f', '3', '-F', '3852', '-q', '30', annot_bam]
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
