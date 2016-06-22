# source /grid/fermiapp/products/common/etc/setups.sh
# setup mu2e
# setup ifdhc
echo 'importing setups for smp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh
setup perl 5.18.1+6
setup Y2Nstack 1.0.6+18
setup diffimg gwdevel8
setup ftools v6.17
setup autoscan v3.1+0
setup easyaccess
setup extralibs 1.0

echo 'sourcing fermiapp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
source /grid/fermiapp/products/common/etc/setups.sh
setup mu2e
setup ifdhc
echo 'copying zip to worker node!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
ifdh cp /pnfs/des/persistent/desdm/djbrout/pysmp.zip .
echo 'unzippppppppppppp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
unzip pysmp.zip
cd pysmp
pwd
echo 'running smp.py!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

python smp.py --index=0 -f g --nozpt --fermigrid
