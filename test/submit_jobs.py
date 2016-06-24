'''
:Date: 27 Apr 2016
:Author: Public Health England
'''

import argparse
import glob
import os
from posix import wait
from subprocess import CalledProcessError
import subprocess
from time import sleep


def get_my_jobs(jobs, job_name=None):
    cmd = ["qstat"]

    current_jobs = set()

    output = subprocess.check_output(cmd)

    c = 0
    for line in output.split("\n"):
        c += 1
        if c <= 2:
            continue

        if job_name and job_name in line:
            current_jobs.add(line.split()[0])
        else:
            current_jobs.add(line.split(" ")[0])

    return jobs & current_jobs

def get_job_id(line):
    line = line.replace("Your job ", "")

    return line[:line.index("(") - 1]

def get_args():
    args = argparse.ArgumentParser()

    args.add_argument("config_dir")

    return args

def main():
    args = vars(get_args().parse_args())

    base_dir = os.path.join("test", "scripts")

    if not os.path.exists("results"):
        os.makedirs("results")

    jobs = set()

    for config in glob.glob(os.path.join(args["config_dir"], "*.yml")):
        config_name = os.path.basename(config)
        params = {"R1": os.path.join(base_dir, "data", "R1.fastq.gz"),
                    "R2": os.path.join(base_dir, "data", "R2.fastq.gz"),
                    "REF": os.path.join(base_dir, "data", "reference.fa"),
                    "outdir": os.path.join("results", config_name.replace(".yml", "")),
                    "CONFIG": config,
                    "sample_name": "test_sample"
                    }

        bash = r"""
#!/bin/bash

. /etc/profile.d/modules.sh

pwd

module load anaconda samtools/1.3 bwa/0.7.13 jdk picard-tools gatk bowtie2
    
PYTHONPATH=.:$PYTHONPATH scripts/phenix.py run_snp_pipeline --keep-temp -r1 %(R1)s -r2 %(R2)s -r %(REF)s -c %(CONFIG)s --sample-name %(sample_name)s -o %(outdir)s
    
""" % params

        job_file = config_name.replace("yml", "sh")
        with open(job_file, "wb") as fp:
            fp.write(bash)

        cmd = ["qsub", "-cwd", "-N", "test_job", "-q", "test.q", job_file]

        try:
            output = subprocess.check_output(cmd)
            jobs.add(get_job_id(output))

        except CalledProcessError:
            continue

    while jobs:

        sleep(5)

        jobs = get_my_jobs(jobs)
        if jobs:
            print "Waiting for [%s] jobs to finish" % ",".join(jobs)

    return 0

if __name__ == '__main__':
    exit(main())
