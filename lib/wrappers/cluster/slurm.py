#!/usr/bin/env python3

# submits slurm jobs

import os
import sys
import math
#from itertools import izip # for py2.7 and old method
from snakemake.utils import read_job_properties

# below is a snakemake way to do this
jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)


cluster = job_properties.get('cluster', {})

group = cluster.get('group', 'snakemake') # job group
queue = cluster.get('queue', None) # job queue
account = cluster.get('account', None) # job account
time = cluster.get('time', '02:00:00') # time limit
job_name = cluster.get('name', 'snakemake') # job name

# Check threads in cluster config first then get job properties
threads = cluster.get('threads', None)
if not threads:
    threads = job_properties.get('threads', 1)

memory = cluster.get('memory', 10) # default to 10Gb
if memory:
    # mem is specified as total memory in Gb, but we
    # need to give it to SGE as memory per thread in Mb
    memory = int(math.floor(float(memory * 1e9) / (threads * 1e6)))

output = cluster.get('output', None)
if output:
    output = os.path.realpath(output)
    tdir = os.path.dirname(output)
    if not os.path.exists(tdir): os.makedirs(tdir)
error = cluster.get('error', None)
if error:
    error = os.path.realpath(error)
    tdir = os.path.dirname(error)
    if not os.path.exists(tdir): os.makedirs(tdir)


# build the command
cmd = 'sbatch'
if queue:
    cmd = '{} --partition={}'.format(cmd, queue)
if account:
    cmd = '{} -A {}'.format(cmd, account)
cmd = '{} --time={}'.format(cmd, time)
cmd = '{} --job-name="{}"'.format(cmd, job_name)
cmd = '{} --mem={}'.format(cmd, memory)

if threads > 1:
    cmd = '{} --cpus-per-task={}'.format(cmd, threads)
if output:
    output = output.replace(",", "__")
    cmd = '{} --output={}'.format(cmd, output)
if error:
    error = error.replace(",", "__")
    cmd = '{} --error={}'.format(cmd, error)
cmd = '{} {}'.format(cmd, jobscript)

# run the command
os.system(cmd)
