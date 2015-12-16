#!/usr/bin/env python
import argparse
import re
import subprocess as sub
import sys

BASE_CMD = ['/home/shane/github/RapMap/build/src/rapmap',
			'quasimap',
			'-i /home/shane/github/RapMap/output/',
			'-r /home/shane/github/RapMap/data/1M.1.fastq',
			'-t 8',
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


def main(num_runs):
    print('Taking average of %s runs.' % num_runs)
    baseTime = run(BASE_CMD + ['-z %s' % num_runs])
    print('Base time average:  %s' % (baseTime/num_runs))
    alignTime = run(ALIGN_CMD + ['-z %s' % num_runs])
    print('Align time average: %s' % (alignTime/num_runs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='RapMap run time averager.')
    parser.add_argument('--runs',
                        type=int,
                        default=10,
                        help='Number of times to run RapMap.')
    args = parser.parse_args()
    main(args.runs)
