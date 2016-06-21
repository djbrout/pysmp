source /grid/fermiapp/products/common/etc/setups.sh
setup mu2e
setup ifdhc
touch test.txt
ifdh mkdir /pnfs/des/persistent/desdm/djbrout
ifdh cp -D test.txt /pnfs/des/persistent/desdm/djbrout
