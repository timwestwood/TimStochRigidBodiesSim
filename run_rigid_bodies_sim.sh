#!/bin/bash
#
# This script should be used to run the code compiled through the makefile
# rather than calling ./rigid_bodies directly. This is to ensure that we always use
# the correct environment variables.
#
# If this script doesn't run, ensure that it has the correct permissions by
# running "chmod u+x run_rigid_bodies_sim.sh" or run it without permissions using
# "bash run_rigid_bodies_sim.sh" ("bash" can be replaced by "sh" etc. if needed).
#
# NOTE: Don't edit this file on Windows as the different line return will break
# the script on Unix systems. Having the Unix returns doesn't matter on Windows
# as I run through cygwin anyway.
#
export OPENBLAS_NUM_THREADS=1
read -p 'Number of processes: ' num_in
mpiexec -n $num_in ./rigid_bodies
