#!/bin/bash -ue
sleep 1
printf 'test21.fastq.gz '
gunzip -c test21.fastq.gz | wc -l
