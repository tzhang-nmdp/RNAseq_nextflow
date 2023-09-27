#!/bin/bash -ue
sleep 1
printf 'test_2.fastq.gz '
gunzip -c test_2.fastq.gz | wc -l
