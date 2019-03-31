#!/bin/bash -l

# Suggested Batch script to run the Facet-based Radar Altimeter Echo Model for Sea Ice
# on Myriad.
#
# Based on openmp.sh by:
# Owain Kenway, Research Computing, 20/Sept/2010

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request 1 hours of wallclock time for testing (max is 48 hours)
# (format hours:minutes:seconds).
#$ -l h_rt=30:00:00

# 3. Request 2 gigabyte of RAM per core - may need to adjust this!
#$ -l mem=2G

# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=15G

# 5. Select 36 threads. On Myriad this is the maximium you can set (one whole node)
#$ -pe smp 36

# 6. Reserve one Matlab licence - this stops your job starting and failing when no
#    licences are available.
#$ -l matlab=1

# 7. The way Matlab threads work requires Matlab to not share nodes with other
# jobs.
#$ -ac exclusive

# 8. Set the name of the job.
#$ -N JobNL
# 8.5 Tells me when stuff messes up in .e file
#$ -notify

# 9. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
#
# This needs to be the dirtectory where the module Matlab scripts are.
#
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/zcapcde/Scratch/


# 12. Run Matlab job

module load xorg-utils/X11R7.7
module load matlab/full/r2018b/9.5
module list

Matlab_infile=SHELL_LP1.m

echo ""
echo "Running matlab -nosplash -nodisplay < $Matlab_infile ..."
echo ""
/usr/bin/time --verbose matlab -nosplash -nodesktop -nodisplay < $Matlab_infile

# End
