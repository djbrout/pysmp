cd /data/des41.a/data/djbrout/pysmp/
 
source /data/des41.a/data/djbrout/pysmp/setup_scripts/setup_fermilab_des41.a.sh

echo python smp.py --index=$1 -f $2
python smp.py --index=$1 -f $2


