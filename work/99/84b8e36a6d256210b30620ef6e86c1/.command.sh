#!/bin/bash -ue
println 
This is a multi-line string
using triple quotes.

printf 'test21.fastq.gz '
gunzip -c test21.fastq.gz | wc -l
