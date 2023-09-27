#!/bin/bash -ue
sleep 2
printf 'test_1.fastq.gz '
gunzip -c test_1.fastq.gz | wc -l
