#!/usr/bin/env python3

import os
import sys
import math
import subprocess
from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

cluster = job_properties.get('cluster', {})

queue = cluster.get('queue', None)
account = cluster.get('account', None)
time = cluster.get('time', '02:00:00')
job_name = cluster.get('name', 'snakemake')

threads = cluster.get('threads', None)
if not threads:
    threads = job_properties.get('threads', 1)

memory = cluster.get('memory', 10)
if memory:
    memory = int(math.floor(float(memory * 1e9) / (threads * 1e6)))

output = cluster.get('output', None)
if output:
    output = os.path.realpath(output)
    os.makedirs(os.path.dirname(output), exist_ok=True)

error = cluster.get('error', None)
if error:
    error = os.path.realpath(error)
    os.makedirs(os.path.dirname(error), exist_ok=True)

cmd = ['sbatch']
if queue:
    cmd += ['--partition', queue]
if account:
    cmd += ['-A', account]
cmd += ['--time', time]
cmd += ['--job-name', job_name]
cmd += ['--mem', str(memory)]

if threads > 1:
    cmd += ['--cpus-per-task', str(threads)]
if output:
    output = output.replace(",", "__")
    cmd += ['--output', output]
if error:
    error = error.replace(",", "__")
    cmd += ['--error', error]

# helpful: make sbatch return only the job id
cmd += ['--parsable', jobscript]

print("SUBMIT:", " ".join(cmd), file=sys.stderr, flush=True)

result = subprocess.run(cmd, text=True, capture_output=True)

print("SBATCH STDOUT:", result.stdout.strip(), file=sys.stderr, flush=True)
print("SBATCH STDERR:", result.stderr.strip(), file=sys.stderr, flush=True)

if result.returncode != 0:
    print(f"SBATCH FAILED with code {result.returncode}", file=sys.stderr, flush=True)
    sys.exit(result.returncode)

# optional: also echo job id to stdout for Snakemake/status tools
print(result.stdout.strip())