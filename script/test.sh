#!/bin/bash
#$ -S /bin/bash
#$ -N DecError
#$ -cwd
#$ -e ./log
#$ -o ./out
#$ -V
#$ -m e
#$ -t 1:2
date
sleep 10 &
wait
date
