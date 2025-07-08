cd ./build
make

cd ../run
ATestRun_eljob.py -c /samples/cp-even-hadhad -s cp-even-hadhad $@
ATestRun_eljob.py -c /samples/cp-odd-hadhad -s cp-odd-hadhad $@
# ATestRun_eljob.py -c /samples/cp-even-hadlep -s cp-even-hadlep $@
# ATestRun_eljob.py -c /samples/cp-odd-hadlep -s cp-odd-hadlep $@

cd ../stuff_in_the_wrong_dir/
python script.py