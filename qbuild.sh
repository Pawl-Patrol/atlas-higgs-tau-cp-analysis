cd ./build
make

cd ../run
ATestRun_eljob.py -c /samples/cp-even-hadhad -s cp-even-hadhad $@

cd ../stuff_in_the_wrong_dir/
python script.py