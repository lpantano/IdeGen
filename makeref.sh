set -e 
FILE=$2
FILEINV=$3
NAME=$1
REFDIR=$5
TEMPDIR=$6
REF=$4


if [ -e "$REF" -a -e "$REFDIR" -a -e  "$CHAIN" -a -e "$FILE" -a -e "$FILEINV" ] ; then 
	printf "ERROR I/O"
	`ls "$REF" "$REFDIR" "$CHAIN" "$FILE" "$FILEINV" `
	exit 1
fi

mkdir -p "$TEMPDIR"
mkdir -p "$REFDIR"

echo GET REF.FA
#get position from file
grep $NAME $FILE | awk '{start=$2-2000000;if ($2-2000000<1){start=1}; print $1"\t"start"\t"$3+2000000"\t.\t.\t+"}' > $REFDIR/ref.bed

cat $REFDIR/ref.bed >>mask.bed

chr=`awk '{print $1}' $REFDIR/ref.bed`
bedtools getfasta -s -fi $REF -bed $REFDIR/ref.bed -fo $REFDIR/ref.fa
sed -i 's/:.*//' $REFDIR/ref.fa

echo GET POSITION AT $FILEINV
POS1=`grep Inversion $FILEINV | head -1 | sed 's/\.\./ /' | awk '{print $2}'`
POS1=`echo $[POS1-1000]`
POS2=`grep Inversion $FILEINV | tail -1 | sed 's/\.\./ /' | awk '{print $2}'`
POS2=`echo $[POS2+1000]`
echo $POS1 $POS2

SEQ=`sed 's/ //g' $FILEINV | grep "^[0-9]" | sed 's/[0-9]//g' |awk 1 ORS=''`
printf ">ORIINV\n"$SEQ > $TEMPDIR/oriseq.fa

INV=${SEQ:$POS1:$POS2-$POS1}

LEFT=${SEQ:$[POS1-1000]:1000}

RIGHT=${SEQ:$POS2-1:1000}

########################################LEFT
echo "GET LEFT SEQ"
printf ">INVL\n"$LEFT > $TEMPDIR/left.fa

/home/lpantano/invregions/blat  $REFDIR/ref.fa $TEMPDIR/left.fa $TEMPDIR/outl

CHRREF=`grep "^[0-9]" $TEMPDIR/outl |  awk '{print $14}'| head -1`
LEFTREF=`grep "^[0-9]" $TEMPDIR/outl |  awk '{print $17+($11-$13)}'| head -1`

STARTLEFTREF=`echo $[LEFTREF-1000000]`
printf $CHRREF"\t"$STARTLEFTREF"\t"$LEFTREF | awk '{print $0}' >$TEMPDIR/leftref.bed
bedtools getfasta  -fi $REFDIR/ref.fa -bed $TEMPDIR/leftref.bed -fo $TEMPDIR/leftref.fa
LEFTINV=`tail -1 $TEMPDIR/leftref.fa`

########################################RIGHT
echo GET RIGHT SEQ
printf ">INVR\n"$RIGHT > $TEMPDIR/right.fa

/home/lpantano/invregions/blat  $REFDIR/ref.fa $TEMPDIR/right.fa $TEMPDIR/outr

CHRREF=`grep "^[0-9]" $TEMPDIR/outr |  awk '{print $14}'| head -1`
RIGHTREF=`grep "^[0-9]" $TEMPDIR/outr |  awk '{print $16-$12}'| head -1`

ENDRIGHTREF=`echo $[RIGHTREF+1000000]`
printf $CHRREF"\t"$RIGHTREF"\t"$ENDRIGHTREF | awk '{print $0}'>$TEMPDIR/rightref.bed
bedtools getfasta  -fi $REFDIR/ref.fa -bed $TEMPDIR/rightref.bed -fo $TEMPDIR/rightref.fa
RIGHTINV=`tail -1 $TEMPDIR/rightref.fa`

########################################
echo GET ALL SEQ
ALL=$LEFTINV$INV$RIGHTINV
printf ">"$NAME"\n"$ALL > $NAME.fa
#echo DO CHECK
#/home/lpantano/invregions/blat -out=blast $TEMPDIR/oriseq.fa all.fa  check.blast
#all.fa is reference to detect new transcripts
######################################
