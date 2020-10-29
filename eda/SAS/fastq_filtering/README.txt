the shell script to filter reads for one sample is:
  scripts/run_filter.sh
which will, among other things, invoke
  scripts/ffq.TS_PE_SE_RNA.pl

The only argument it requires is the identifier for the Sample
For example, assuming the necessary raw fastq files are in place as expected, you might execute:

cd fastq_filtering
nohup ./run_filter.sh C100_S3 > nohup.C100_S3.out 2> nohup.C100_S3.log &

This would produce a C100_S3 subdirectory in fastq_filtering, containing several output files
  of various sorts and purposes.

See scripts/run_filter.sh for more details.
