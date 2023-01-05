To be checked with 

```
gzip test_R1.fastq.gz ## kept it uncompressed for the sake of visualization/gitlab
Rscript 01_harmonize_fastqs.R -r1 test_R1.fastq.gz -r2 test_R1.fastq.gz -o test
# and compare with the expected_output folder
```
