#! /usr/bin/env python

import argparse
import os
from plumbum import local
import pwd
import time

def num_jobs(jobname, uge_proj='broad'):
    ''' Get the number of jobs currently running by parsing qstat output.
    '''
    qstat_cmd = local['qstat']['-u', pwd.getpwuid(os.getuid())[0], '-j', jobname] | local['grep']['-P', 'project:\s*{!s}'.format(uge_proj)]
    return len([elt for elt in qstat_cmd(retcode=None).strip().split('\n') if elt])

def wait_for_jobs(jobname, uge_proj='broad', num_jobs=1):
    ''' Block in this function until a certain number of jobs finish.
    '''
    qstat_cmd = local['qstat']['-u', pwd.getpwuid(os.getuid())[0], '-j', jobname] | local['grep']['-P', 'project:\s*{!s}'.format(uge_proj)]
    while len([elt for elt in qstat_cmd(retcode=None).strip().split('\n') if elt]) >= num_jobs:
        time.sleep(10)
    time.sleep(10)
    #sometimes qstat errors can make it seem like there are no jobs, so double-check before returning
    while len([elt for elt in qstat_cmd(retcode=None).strip().split('\n') if elt]) >= num_jobs:
        time.sleep(10)
    return

def wait_for_jobs_check_return(job_name, uge_proj='broad', num_jobs=1):
    #wait for jobs
    wait_for_jobs(job_name, uge_proj=uge_proj, num_jobs=num_jobs)

    #check that all tasks exited successfully
    exit_codes_cmd = local['qacct']['-j', job_name] | local['grep']['exit_status']
    for attempt in range(10):
        try:
            exit_codes_cmd_res = exit_codes_cmd()
        except:
            time.sleep(2)
        else:
            break
    else:
        raise Exception('Could not identify error codes of subjobs. Failed command: {!s}'.format(str(exit_codes_cmd)))
    exit_codes = [int(elt.strip().split()[1]) for elt in exit_codes_cmd_res.strip().split('\n')]
    if sum(exit_codes) > 0:
        raise Exception('Job array did not complete successfully. Exiting.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('job_name', help='Job identifier, can be a job id, job name, or pattern with the same syntax as the -j option for qstat and qacct.')
    parser.add_argument('--uge_proj', default='broad')
    parser.add_argument('--num_jobs', type=int, default=1)
    parser.add_argument('--no_exit_check', action='store_true', default=False)
    args = parser.parse_args()

    if args.no_exit_check is True:
        wait_for_jobs(args.job_name, uge_proj=args.uge_proj, num_jobs=args.num_jobs)
    else:
        wait_for_jobs_check_return(args.job_name, uge_proj=args.uge_proj, num_jobs=args.num_jobs)
