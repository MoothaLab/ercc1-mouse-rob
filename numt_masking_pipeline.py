#! /usr/bin/env python

'''
This NUMT-masking pipeline is inspired by the NUMT-masking strategy in Lareau, et al. Nat Biotech. 2023.
Here, we add steps to align reads to both the standard chrM reference and a shifted chrM reference to 
account for the circular nature of the chrM genome, and also a step to call sample-specific homoplasmies
that differ from the chrM reference so that we can generate a sample-specific reference.
'''

import argparse
import itertools
import os
from plumbum import local
import pypiper
import random
import stat
import sys

sys.path.append(os.path.dirname(__file__))
from wait_for_qsub_jobs import wait_for_jobs, wait_for_jobs_check_return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq1')
    parser.add_argument('--fastq2')
    parser.add_argument('--out_dir', default=os.getcwd())
    parser.add_argument('--index_type', help='The aligner index to create. Either bowtie2 or star.')
    parser.add_argument('--reference_fasta')
    parser.add_argument('--star_aln_annot_gtf', help='Only needed if we are aligning RNA-seq data with STAR')
    parser.add_argument('--no_chrM_reference_fasta')
    parser.add_argument('--indexed_no_chrM_reference')
    parser.add_argument('--aln_nthreads', type=int, default=1)
    parser.add_argument('--aln_gb_per_thread', type=int, default=16)
    parser.add_argument('--sample_name')
    parser.add_argument('--art_error_model_r1', help='Specify a custom error model for read1 built with empirical data using the art_profiler_illumina tool.')
    parser.add_argument('--art_error_model_r2', help='Specify a custom error model for read1 built with empirical data using the art_profiler_illumina tool.')
    parser.add_argument('--art_seqSys', default='NS50', help='If no custom error profile is supplied through the --art_error_model_r1/2 option(s), use this built-in one (NextSeq NS50 by default).')
    parser.add_argument('--pipeline_recover', action='store_true', default=False)
    parser.add_argument('--pipeline_no_cleanup', action='store_true', default=False)
    parser.add_argument('--pipeline_run_all', action='store_true', default=False)
    parser.add_argument('--pipeline_jvmmemavail', type=int, help='Available memory in GB', default=30)
    args = parser.parse_args()

    code_dir = os.path.dirname(__file__)
    if args.sample_name is None:
        input_basename = os.path.splitext(os.path.basename(args.fastq1).rstrip('.gz'))[0]
        args.sample_name = input_basename
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
    out_prefix = os.path.join(args.out_dir, args.sample_name)
    pm = pypiper.PipelineManager(name="NUMT_masking", outfolder=args.out_dir,
                                 recover=args.pipeline_recover,
                                 dirty=args.pipeline_no_cleanup,
                                 new_start=args.pipeline_run_all)


    ##########
    # MAKE SHIFTED CHRM REFERENCE
    ##########
    pm.timestamp('### Make shifted chrM')

    script_path = os.path.join(code_dir, 'shift_chrM.py')
    nonshifted_chrM = out_prefix + '.chrM.fasta'
    nonshift_cmd = local['python'][script_path,
                                   args.reference_fasta,
                                   '--shifted_fasta={!s}'.format(nonshifted_chrM),
                                   '--split_idx=0']

    shifted_chrM = out_prefix + '.chrM_shifted.fasta'
    shift_cmd = local['python'][script_path,
                                args.reference_fasta,
                                '--shifted_fasta={!s}'.format(shifted_chrM),
                                '--split_idx=8000']
    pm.run([str(nonshift_cmd), str(shift_cmd)], target=shifted_chrM)


    ##########
    # INDEX THE CHRM REFERENCES FOR ALIGNMENT
    ##########
    pm.timestamp('### Make aligner indices')

    if args.index_type == 'bowtie2':
        nonshifted_ref = os.path.join(args.out_dir, 'chrM')
        index_cmd1 = local['bowtie2-build']['-f', '--seed=50', '--threads=1', nonshifted_chrM, nonshifted_ref]

        shifted_ref = os.path.join(args.out_dir, 'chrM_shifted')
        index_cmd2 = local['bowtie2-build']['-f', '--seed=50', '--threads=1', shifted_chrM, shifted_ref]

        pm.run([str(index_cmd1), str(index_cmd2)], target=shifted_ref + '.1.bt2')
    elif args.index_type == 'star':
        nonshifted_ref = os.path.join(args.out_dir, 'chrM')
        index_cmd1 = local['STAR']['--runThreadN=1',
                                   '--runMode', 'genomeGenerate',
                                   '--genomeDir={!s}'.format(nonshifted_ref),
                                   '--genomeFastaFiles', nonshifted_chrM,
                                   '--sjdbGTFfile', args.star_aln_annot_gtf,
                                   '--genomeSAindexNbases', '6'] #tuning this parameter because of the following warning when using the default parameter: WARNING: --genomeSAindexNbases 14 is too large for the genome size=262144, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 8. Also, note that the manual states that this value should be chosen based on this formula min(14, log2(GenomeLength)/2 - 1), which actually works out to about 6 for the mitochondrial genome. I am seeing SegFaults in one of the samples so I may try values lower than 8.
        
        shifted_ref = os.path.join(args.out_dir, 'chrM_shifted')
        index_cmd2 = local['STAR']['--runThreadN=1',
                                   '--runMode', 'genomeGenerate',
                                   '--genomeDir={!s}'.format(shifted_ref),
                                   '--genomeFastaFiles', shifted_chrM,
                                   '--sjdbGTFfile', args.star_aln_annot_gtf,
                                   '--genomeSAindexNbases', '6'] #tuning this parameter because of the following warning when using the default parameter: WARNING: --genomeSAindexNbases 14 is too large for the genome size=262144, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 8
        
        pm.run([str(index_cmd1), str(index_cmd2)], target=os.path.join(shifted_ref, 'SA'))
    else:
        raise Exception('Currently only bowtie2 and star indexing is implemented.')


    ##########
    # ALIGN TO THE CHRM REFERENCES
    ##########
    pm.timestamp('### Align to shifted and non-shifted chrM')

    job_names = []
    qsub_jobs = []
    for ref_path in [nonshifted_ref, shifted_ref]:
        if args.index_type == 'bowtie2':
            out_bam = ref_path + '.bam'
            aln_args = ['-p', '{!s}'.format(args.aln_nthreads),
                        '-X', '2000',
                        '-x', '{!s}'.format(ref_path),
                        '-1', args.fastq1,
                        '-2', args.fastq2]
            aln_cmd = local['bowtie2'].__getitem__(aln_args) | local['samtools']['view', '-bhSu', '-f', '3', '-F', '12', '-q', '10', '-'] | local['samtools']['sort', '-T', args.out_dir, '-o', out_bam, '-@', '{!s}'.format(args.aln_nthreads)]
            index_cmd = local['samtools']['index', out_bam]
            cmd_str = '{!s}; {!s};'.format(str(aln_cmd), str(index_cmd))

        elif args.index_type == 'star':
            out_bam_prefix = ref_path
            fastq_list = [elt for elt in [args.fastq1, args.fastq2] if elt is not None]
            cmd_args = ['--runThreadN', '{!s}'.format(args.aln_nthreads),
                        '--genomeDir', '{!s}'.format(ref_path),
                        '--readFilesIn'] + fastq_list + [
                        '--readFilesCommand', '{!s}'.format('zcat' if fastq_list[0].endswith('.gz') else 'cat'),
                        '--outFileNamePrefix', '{!s}_'.format(out_bam_prefix),
                        '--runRNGseed', '777',
                        '--outSAMtype', 'BAM', 'Unsorted',
                        '--outSAMattributes', 'NH', 'HI', 'NM', 'MD', 'AS', 'XS', 'nM',
                        '--twopassMode', 'Basic']
            aln_cmd = local['STAR'].__getitem__(cmd_args)
            out_unsorted = out_bam_prefix + '_Aligned.out.bam'
            out_bam = out_bam_prefix + '.bam'
            sort_cmd = (local['samtools']['view', '-bhu', '-q', '10', out_unsorted]
                        | local['samtools']['sort', '-@', str(args.aln_nthreads - 1),
                                            '-m', '{!s}G'.format(args.aln_gb_per_thread),
                                            '-o', out_bam, '-'])
            index_cmd = local['samtools']['index', out_bam]
            cmd_str = '{!s}; {!s}; {!s};'.format(str(aln_cmd), str(sort_cmd), str(index_cmd))

        #if the pipeline failed while aligning to the second (i.e. the shifted) reference, don't
        # also re-align to the first (non-shifted) reference
        if (not args.pipeline_run_all) and os.path.exists(out_bam) and os.path.exists(out_bam + '.bai'):
            continue
        conda_env = os.path.basename(os.environ['CONDA_PREFIX'])
        job_out = ref_path + '_{!s}_aln.{:04}.out'.format(args.sample_name, random.randint(0,9999))
        job_name = os.path.basename(job_out)
        qsub_path = os.path.join(args.out_dir, job_name + '.qsub.sh')
        with open(qsub_path, 'w') as out:
            out.write('#! /usr/bin/env bash\n\n')
            out.write('# run in conda env: {!s}\n\n'.format(conda_env))
            out.write(cmd_str + '\n')
        os.chmod(qsub_path, 0o774)
        cmd_to_submit = 'conda run -n {!s} {!s}'.format(conda_env, qsub_path)
        job_array_args = ['-j', 'y',
                          '-b', 'y',
                          '-cwd',
                          '-N', job_name,
                          '-o', job_out,
                          '-l', 'h_vmem={!s}G'.format(args.aln_gb_per_thread),
                          '-l', 'h_rt=12:00:00']
        if args.aln_nthreads > 1:
            job_array_args.extend(['-pe', 'smp', str(args.aln_nthreads), '-binding', 'linear:{!s}'.format(args.aln_nthreads)])
        job_array_args.append(cmd_to_submit)
        qsub_job = local['qsub'].__getitem__(job_array_args)
        print(str(qsub_job))
        qsub_jobs.append(qsub_job)
        job_names.append(job_name)

    #run the alignments in parallel as UGER jobs, then wait for them to complete
    # via the wait_for_qsub_jobs.py before continuing.
    if len(job_names) > 0:
        script_path = os.path.join(code_dir, 'wait_for_qsub_jobs.py')
        wait_cmds = [local['python'][script_path, elt, '--uge_proj=broad', '--num_jobs=1'] for elt in job_names]

        pm.run([str(elt) for elt in qsub_jobs] + [str(elt) for elt in wait_cmds], target=out_bam + '.bai')


    ###########
    # UNSHIFT THE COORDINATES FROM THE SHIFTED REFERENCE
    ###########
    pm.timestamp('### Unshift chrM coords')

    script_path = os.path.join(code_dir, 'unshift_chrM_coords.py')
    shifted_bam = shifted_ref + '.bam'
    unshifted_bam = shifted_ref + '.unshifted.bam'
    unshift_cmd = local['python'][script_path,
                                  shifted_bam,
                                  '--out_bam_file={!s}'.format(unshifted_bam),
                                  '--chrM_shift=8000']
    pm.run(str(unshift_cmd), target=unshifted_bam + '.bai')


    ##########
    # SELECT THE BEST MAPPING LOCATION FOR READS THAT MAP TO DIFFERENT COORDINATES IN THE
    #   SHIFTED AND UNSHIFTED ALIGNMENTS
    ##########
    pm.timestamp('### Reconcile shifted and unshifted alignments.')

    script_path = os.path.join(code_dir, 'reconcile_alignments.py')
    nonshifted_bam = nonshifted_ref + '.bam'
    reconciled_bam = out_prefix + '.chrM_reconciled.bam'
    reconcile_cmd = local['python'][script_path,
                                    nonshifted_bam,
                                    unshifted_bam,
                                    reconciled_bam,
                                    '--java_maxmem_gb={!s}'.format(args.pipeline_jvmmemavail)]
    pm.run(str(reconcile_cmd), target=reconciled_bam + '.bai')


    ##########
    # ITERATE OVER THE CHRM PILEUP AND CALL SAMPLE HOMOPLASMIES THAT DIFFER FROM REFERENCE
    ##########
    pm.timestamp('### Call consensus chrM reference')

    script_path = os.path.join(CODE_DIR, 'run_custom_varcall.multiproc.py')
    consensus_ref_fasta = os.path.splitext(reconciled_bam)[0] + '.consensus_ref.fasta'
    consensus_ref_cmd = local['python'][script_path,
                                        '--input_list={!s}'.format(reconciled_bam),
                                        '--ref_fasta={!s}'.format(args.reference_fasta),
                                        '--out_root={!s}'.format(args.out_dir),
                                        '--out_base_prefix=consensus_ref',
                                        '--chrom=chrM',
                                        '--just_call_bulk_consensus_ref']
    pm.run(str(consensus_ref_cmd), target=consensus_ref_fasta)

    ##########
    # COMPILE SOME MAPPING STATISTICS NEEDED FOR SIMULATING READS
    ##########
    pm.timestamp('### Get chrM mapping stats')
    out_stats = os.path.splitext(reconciled_bam)[0] + '.stats.txt'
    samtools_cmd = local['samtools']['stats', reconciled_bam, 'chrM'] > out_stats
    pm.run(str(samtools_cmd), target=out_stats)
    stats_txt = (local['grep']['^SN', out_stats] | local['cut']['-f', '2-'])()
    stats_dict = dict([elt.strip().split('\t')[:2] for elt in stats_txt.strip().split('\n')])
    read_len = stats_dict['average length:']
    frag_len = stats_dict['insert size average:']
    frag_len_sdev = stats_dict['insert size standard deviation:']

    ##########
    # SIMULATE CHRM READS USING THE ART PROGRAM
    ##########
    pm.timestamp('### Make synthetic sequencing result with ART')

    synthetic_fq_base = os.path.splitext(consensus_ref_fasta)[0]
    art_args = ['--noALN',
                '--len', read_len,
                '--rcount', '10000000',
                '--in', consensus_ref_fasta,
                '--out', synthetic_fq_base]
    if ((args.fastq2 is not None)
        and (args.art_error_model_r1 is not None)):
        if args.art_error_model_r2 is None:
            raise Exception('If supplying a custom ART error model for paired end data, you must supply an error model for read1 and one for read2.')
        art_args.extend(['--qprof1', args.art_error_model_r1,
                         '--qprof2', args.art_error_model_r2])
    elif args.art_error_model_r1 is not None:
        art_args.extend(['--qprof1', args.art_error_model_r1])
    else:
        art_args.extend(['--seqSys', args.art_seqSys])

    if args.fastq2 is not None:
        target_loc = synthetic_fq_base + '2.fq'
        art_args.extend(['--paired',
                         '--mflen', frag_len,
                         '--sdev', frag_len_sdev])
    else:
        target_loc = synthetic_fq_base + '.fq'

    art_cmd = local['art_illumina'].__getitem__(art_args)
    pm.run(str(art_cmd), target=target_loc)


    ##########
    # ALIGN THE SIMULATED READS TO THE NUCLEAR GENOME REFERENCE
    ##########
    pm.timestamp('### Align synthetic data to nuclear genome reference')

    if args.index_type == 'bowtie2':
        out_bam = out_prefix + '.synthetic_reads_aligned_to_nuc.bam'
        aln_args = ['-p {!s}'.format(args.aln_nthreads),
                    '-X 2000',
                    '-x {!s}'.format(args.indexed_no_chrM_reference),
                    '-1', synthetic_fq_base + '1.fq',
                    '-2', synthetic_fq_base + '2.fq']
        aln_cmd = local['bowtie2'].__getitem__(aln_args) | local['samtools']['view', '-bhSu', '-f 3', '-F 12', '-q 10', '-'] | local['samtools']['sort', '-T', args.out_dir, '-o', out_bam, '-@ {!s}'.format(args.aln_nthreads)]
        index_cmd = local['samtools']['index', out_bam]
        cmd_str = '{!s}; {!s};'.format(str(aln_cmd), str(index_cmd))

    elif args.index_type == 'star':
        out_bam_prefix = out_prefix + '.synthetic_reads_aligned_to_nuc'
        fastq_list = [elt for elt in [args.fastq1, args.fastq2] if elt is not None]
        cmd_args = ['--runThreadN {!s}'.format(args.aln_nthreads),
                    '--genomeDir {!s}'.format(args.indexed_no_chrM_reference),
                    '--readFilesIn'] + fastq_list + [
                    '--readFilesCommand {!s}'.format('zcat' if fastq_list[0].endswith('.gz') else 'cat'),
                    '--outFileNamePrefix {!s}_'.format(out_bam_prefix),
                    '--runRNGseed 777',
                    '--outSAMtype BAM Unsorted',
                    '--outSAMattributes NH HI NM MD AS XS nM',
                    '--twopassMode Basic']
        aln_cmd = local['STAR'].__getitem__(cmd_args)
        out_unsorted = out_bam_prefix + '_Aligned.out.bam'
        out_bam = out_bam_prefix + '_Aligned.sortedByCoord.out.bam'
        sort_cmd = (local['samtools']['view', '-bhu', '-q 10', out_unsorted]
                    | local['samtools']['sort', '-@', str(args.aln_nthreads - 1),
                                        '-m', '{!s}G'.format(args.aln_gb_per_thread),
                                        '-o', out_bam, '-'])
        index_cmd = local['samtools']['index', out_bam]

        cmd_str = '{!s}; {!s}; {!s};'.format(str(aln_cmd), str(sort_cmd), str(index_cmd))


    conda_env = os.path.basename(os.environ['CONDA_PREFIX'])
    #write the qsub script
    qsub_path = out_bam.replace('.bam', '.qsub.sh')
    with open(qsub_path, 'w') as out:
        out.write('#! /usr/bin/env bash\n\n')
        out.write('# run in conda env: {!s}\n\n'.format(conda_env))
        out.write(cmd_str + '\n')
    os.chmod(qsub_path, 0o774)
    cmd_to_submit = 'conda run -n {!s} {!s}'.format(conda_env, qsub_path)
    job_out = os.path.splitext(out_bam)[0] + '_aln.{:04}.out'.format(random.randint(0,9999))
    job_name = os.path.basename(job_out)
    job_array_args = ['-j', 'y',
                      '-b', 'y',
                      '-cwd',
                      '-N', job_name,
                      '-o', job_out,
                      '-l', 'h_vmem={!s}G'.format(args.aln_gb_per_thread),
                      '-l', 'h_rt=12:00:00']
    if args.aln_nthreads > 1:
        job_array_args.extend(['-pe', 'smp', str(args.aln_nthreads), '-binding', 'linear:{!s}'.format(args.aln_nthreads)])
    job_array_args.append(cmd_to_submit)
    #submit the alignment as a UGER job, then wait for it to complete using the wait_for_qsub_jobs.py script
    if args.pipeline_run_all or not (os.path.exists(out_bam) and os.path.exists(out_bam + '.bai')):
        qsub_job = local['qsub'].__getitem__(job_array_args)

        script_path = os.path.join(code_dir, 'wait_for_qsub_jobs.py')
        wait_cmd = local['python'][script_path, job_name, '--uge_proj=broad', '--num_jobs=1']

        pm.run([str(qsub_job), str(wait_cmd)], target=out_bam + '.bai')


    ##########
    # USE MACS2 to call NUMTs as regions where the simulated chrM reads align
    ##########
    pm.timestamp('### Call NUMT loci')

    if 'hg' in os.path.basename(args.indexed_no_chrM_reference):
        gsize = '2.7e9'
    elif 'mm' in os.path.basename(args.indexed_no_chrM_reference):
        gsize = '1.87e9'
    else:
        raise Exception('Could not determine whether {!s} is a human or mouse reference.'.format(args.indexed_no_chrM_reference))
    narrowPeaks = os.path.splitext(out_bam)[0] + '_peaks.narrowPeak'
    macs2_cmd = local['macs2']['callpeak', '-t', out_bam, '--nomodel', '--nolambda', '--keep-dup', 'all', '-g', gsize, '-n', os.path.splitext(out_bam)[0]]

    pm.run(str(macs2_cmd), target=narrowPeaks)


    ##########
    # CREATE A NEW REFERENCE WITH THE NUMT REGIONS MASKED
    ##########
    pm.timestamp('### Mask NUMT loci')

    merged_fasta = out_prefix + '.merged.fasta'
    masked_fasta = os.path.splitext(merged_fasta)[0] + '.numt_masked.fasta'
    merge_fasta_cmd = local['cat'][args.no_chrM_reference_fasta, consensus_ref_fasta] > merged_fasta
    mask_cmd = local['bedtools']['maskfasta', '-fi', merged_fasta, '-bed', narrowPeaks, '-fo', masked_fasta]
    faidx_cmd = local['samtools']['faidx', masked_fasta]
    seqdict_cmd = local['gatk']['--java-options', '-Xmx{!s}G'.format(args.pipeline_jvmmemavail),
                                'CreateSequenceDictionary',
                                '-R', masked_fasta]
    pm.run([str(merge_fasta_cmd), str(mask_cmd), str(faidx_cmd), str(seqdict_cmd)], target=masked_fasta)


    ##########
    # INDEX THE NEW REFERENCE
    ##########
    pm.timestamp('### Index new reference consisting of masked nuclear genome and consensus chrM')

    indexed_ref = os.path.splitext(masked_fasta)[0]
    if args.index_type == 'bowtie2':
        index_cmd = local['bowtie2-build']['-f', '--seed=51', '--threads={!s}'.format(args.aln_nthreads), masked_fasta, indexed_ref]
        target_loc = indexed_ref + '.1.bt2'
    elif args.index_type == 'star':
        index_cmd = local['STAR']['--runThreadN={!s}'.format(args.aln_nthreads),
                                  '--runMode', 'genomeGenerate',
                                  '--genomeDir={!s}'.format(indexed_ref),
                                  '--genomeFastaFiles', masked_fasta,
                                  '--sjdbGTFfile', args.star_aln_annot_gtf]
        target_loc = os.path.join(indexed_ref, 'SA')
    else:
        raise Exception('Only bowtie2 and star index types are implemented.')

    conda_env = os.path.basename(os.environ['CONDA_PREFIX'])
    cmd_to_submit = 'conda activate {!s}; {!s};'.format(conda_env, str(index_cmd))
    job_out = indexed_ref + '.index_job.{:04}.out'.format(random.randint(0,9999))
    job_name = os.path.basename(job_out) + '2'
    job_array_args = ['-j', 'y', '-b', 'y', '-cwd', '-N', job_name, '-o', job_out, '-l', 'h_vmem={!s}G'.format(args.aln_gb_per_thread), '-l', 'h_rt=6:00:00']
    if args.aln_nthreads > 1:
        job_array_args.extend(['-pe', 'smp', str(args.aln_nthreads), '-binding', 'linear:{!s}'.format(args.aln_nthreads)])
    job_array_args.append(cmd_to_submit)
    if args.pipeline_run_all or not os.path.exists(target_loc):
        qsub_cmd = local['qsub'].__getitem__(job_array_args)
        qsub_cmd()
        script_path = os.path.join(code_dir, 'wait_for_qsub_jobs.py')
        wait_cmd = local['python'][script_path, job_name, '--uge_proj=broad', '--num_jobs=1']

        pm.run([str(qsub_cmd), str(wait_cmd)], target=target_loc)

    if args.index_type == 'bowtie2':
        sizes_cmd = (local['grep']['^@SQ', indexed_ref + '.dict'] | local['sed']['-E', 's/.*SN:([a-zA-Z0-9]+)\tLN:([0-9]+).*/\1\t\2/']) > indexed_ref + '.sizes'
        pm.run(str(sizes_cmd), target=indexed_ref + '.sizes'

    ##########
    # END OF PIPELINE
    ##########
    pm.stop_pipeline()
