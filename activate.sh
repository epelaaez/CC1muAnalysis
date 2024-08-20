source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh

dirName=v09_75_03
cd ${dirName}
source localProducts_larsoft_v09_75_03_e20_prof/setup

mrbsetenv
cd ..

kx509
voms-proxy-init -noregen -rfc -voms 'fermilab:/fermilab/sbnd/Role=Analysis'
export TERM=screen