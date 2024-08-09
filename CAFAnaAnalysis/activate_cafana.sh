source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh

dirName=v09_75_03
cd ${dirName}
source localProducts_larsoft_v09_75_03_e20_prof/setup

mrbsetenv

mrb i -j4 # recompile

mrbslp
cd ..