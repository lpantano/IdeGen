DIRNAME=`dirname $BASH_SOURCE`
source $DIRNAME"/include.sh"

REF=$1
MASK=$2
REFDIR=$3

if [ -e "$REF" -a -f "$MASK" -a -e "$REFDIR" ] ; then 
 	printf "everything seems correct\n"
else
	printf "ERROR I/O\n"
 	ls "$REFDIR" "$REF"  "$MASK"
 	exit 1
fi


bedtools  maskfasta -fi $REF -bed $MASK -fo $REFDIR/genome.mask.fa

