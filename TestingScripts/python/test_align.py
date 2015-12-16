#!/usr/bin/env python
import argparse
import re
import subprocess as sub
import sys

BASE_CMD = ['/home/shane/github/RapMap/build/src/rapmap',
            'quasimap',
            '-i /home/shane/github/RapMap/output/',
            '-1 /home/shane/github/RapMap/data/1M.1.fastq',
            '-2 /home/shane/github/RapMap/data/1M.2.fastq',
            '-n']
ALIGN_CMD = BASE_CMD + ['-a']

regex = r"Elapsed time: (\d*\.\d+|\d+)s"


def run(cmd):
    proc = sub.Popen(cmd, stdout=sub.PIPE, stderr=sub.PIPE)
    _, stderr = proc.communicate()
    proc.wait()
    time_str = stderr.split()[-1][:-1]
    times = re.findall(regex, stderr)
    total = sum(map(float, times))
    return total


def main(runs, threads):
    print('Taking average of %s runs with %s threads.' % (runs, threads))
    config = ['-t %s' % threads, '-z %s' % runs]
    baseTime = run(BASE_CMD + config)/runs
    print('Base time average:  %s' % baseTime)
    alignTime = run(ALIGN_CMD + config)/runs
    print('Align time average: %s' % alignTime)
    print('Align/Base: %s' % (alignTime/baseTime))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='RapMap run time averager.')
    parser.add_argument('--runs',
                        type=int,
                        default=10,
                        help='Number of times to run RapMap.')
    parser.add_argument('--threads',
                        type=int,
                        default=8,
                        help='Number of threads to use.')
    args = parser.parse_args()
    main(args.runs, args.threads)
