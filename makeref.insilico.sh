DIRNAME=`dirname $BASH_SOURCE`
source $DIRNAME"/include.sh"

NAME=$1
FILE=$2
REF=$3
SIZE=$4
TEMPDIR=$5
REFDIR=$6


# if [ -e "$REF" -a -f "$FILE" -a -e "$REF" ] ; then 
#  	printf "evreything seems correct\n"
# else
# 	printf "ERROR I/O\n"
#  	ls "$REFDIR" "$REF"  "$FILE"
#  	exit 1
# fi

echo "USING INSILICO SEQUENCE"


mkdir -p "$TEMPDIR"
mkdir -p "$REFDIR"

echo "GET REF.FA"
#get position from file
chr=`grep -w $NAME $FILE | cut -f 1`
max=`grep -w $chr $SIZE | cut -f 2`
printf "chr $chr has size $max\n"

grep $NAME $FILE | awk '{start=$2-1000000;if (start<1){start=1}; print $1"\t"start"\t"$2"\t.\t.\t+"}' > $TEMPDIR/left.bed
grep $NAME $FILE | awk -v max=$max '{end=$3+1000000;if (end>max){end=max-100}; print $1"\t"$3-1"\t"end"\t.\t.\t+"}' > "$TEMPDIR/rigth.bed"
grep $NAME $FILE | awk -v name=$NAME  '{print $1"\t"$2"\t"$3-1"\t"name"\t.\t-"}' > "$TEMPDIR/inv.bed"
grep $NAME $FILE | awk '{print $1"\t"$2"\t"$3-1}' > "$$REFDIR/$NAME.bed"


cat $TEMPDIR/*.bed >>"$REFDIR"/mask.bed

bedtools getfasta  -fi $REF -bed "$TEMPDIR/left.bed" -fo "$TEMPDIR/left.fa"
bedtools getfasta  -fi $REF -bed "$TEMPDIR/rigth.bed" -fo "$TEMPDIR/rigth.fa"
bedtools getfasta  -s -fi $REF -bed "$TEMPDIR/inv.bed" -fo "$TEMPDIR/inv.fa"


SEQ1=`grep -v ">" $TEMPDIR/left.fa |awk 1 ORS=''`
SEQ2=`grep -v ">" $TEMPDIR/inv.fa |awk 1 ORS=''`
SEQ3=`grep -v ">" $TEMPDIR/rigth.fa |awk 1 ORS=''`
printf ">$NAME\n"$SEQ1$SEQ2$SEQ3"\n" > "$REFDIR/$NAME.fa"

rm "$TEMPDIR"/*fa
rm "$TEMPDIR"/*bed
rm "$TEMPDIR"/out*
#all.fa is reference to detect new transcripts
######################################
