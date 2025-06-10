cd ./build
make

cd ../run
ATestRun_eljob.py -c ../samples/cp-even -s cp-even $@
ATestRun_eljob.py -c ../samples/cp-odd -s cp-odd $@

cd ../stuff_in_the_wrong_dir/
python script.py