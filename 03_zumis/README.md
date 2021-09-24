Steps:
0. Install the required soft, including `libhdf5-dev`
1. Run the indexing for mouse/genome/alien
2. Edit the yaml file, i.e. as `test.yaml` (this example has no BC whitelisting)
3. Uncompress the relevant GTF
3. Run, e.g. with the proper binaries within the $PATH (no conda), with `./zUMIs.sh  -y test.yaml -d /home/ubuntu/sandbox/zUMIs/ 2>&1 | tee -a logs/test_counting.log`
