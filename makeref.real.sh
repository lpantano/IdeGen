DIRNAME=`dirname $BASH_SOURCE`
source $DIRNAME"/include.sh"

FILE=$2
FILEINV=$3
NAME=$1
REF=$4
SIZE=$5
REFDIR=$7
TEMPDIR=$6


if [ "$#" -ne 7 ] ; then
	printf "need 7 arguments\n"
	exit 1
fi

# if [ -f "$REF" -a -e "$REFDIR" -a -f "$FILE" -a -f "$FILEINV" -a -f "$SIZE" ] ; then 
# 	printf "ERROR I/O"
# 	`ls "$REF" "$REFDIR" "$SIZE" "$FILE" "$FILEINV" `
# 	exit 1
# fi

echo "USING REAL SEQUENCE"

mkdir -p "$TEMPDIR"
mkdir -p "$REFDIR"

SEQ=`sed 's/ //g' $FILEINV | grep "^[0-9]" | sed 's/[0-9]//g' |awk 1 ORS=''`
printf ">ORIINV\n"$SEQ > $TEMPDIR/oriseq.fa


echo "GET REF.FA"
#get position from file
chr=`grep -w $NAME $FILE | cut -f 1` #get chr where inv is
max=`grep -w $chr $SIZE | cut -f 2` #get size of the chr
printf "chr $chr has size $max\n"

grep $NAME $FILE | awk -v max=$max '{start=$2-2000000;end=$3+2000000;if ($2-2000000<1){start=1};if (end>max){end=max-1000}; print $1"\t"start"\t"end"\t.\t.\t+"}' > $TEMPDIR/ref.bed #region containing the inv with 2 MB flanking regions

cat $TEMPDIR/ref.bed
maxref=`cat $TEMPDIR/ref.bed | awk '{print  $3-$2}'` #get max size of the previous region
printf "ref region has size $maxref\n"

cat $TEMPDIR/ref.bed >>"$REFDIR/mask.bed" #add region to mask.bed, used to mask the genome

bedtools getfasta -s -fi $REF -bed "$TEMPDIR/ref.bed" -fo "$TEMPDIR/ref.fa" #get region containing the inv
sed -i 's/:.*//' $TEMPDIR/ref.fa #get the chr as the name of the fasta

echo "GET POSITION AT $FILEINV"
POS1=`grep -B 1 "label=BP" $FILEINV | head -1 | sed 's/\.\./ /' |sed 's/\^/ /' | sed 's/[A-Za-z()_]/ /g'| awk '{print $1}'` #find first BP
POS2=`grep -B 1  "label=BP" $FILEINV | tail -2 | head -1  | sed 's/\.\./ /' |sed 's/\^/ /' |  sed 's/[A-Za-z()_]/ /g'| awk '{print $NF}'` #find last BP
echo "BP $POS1 $POS2"
INV=${SEQ:$POS1:$POS2-$POS1} #get sequence among BPs
printf ">INVREAL\n"$INV > $TEMPDIR/realinv.fa #create FASTA with previous sequence


POS1=`echo $POS1 | awk '{start=$1-1000;if (start<0){start=$1};print start}'` #define flanking of BP1 region to be used to find the absolute position in the GENOME
POS2=`echo $[POS2+1000]` #define flanking of BP2 region to be used to find the absolute position in the GENOME
echo "estimated $POS1 $POS2"

LENS=`echo "$SEQ" | awk '{print length($1)}'` #determine length of the sequence in gbk file
POS2=`echo $LENS $POS2 | awk '{if ($1>$2){print $2}else{print $2-1000}}'` #change BP2 if downstream 1Kb flanking region is longer than the length of the sequence in the file

echo "final $POS1 $POS2"
echo "length origianl $LENS"

LEN1=`echo $POS1 | awk '{len=1000;if ($1-1000<0){len=$1-1};print len}'` #truncate the size of the sequence that it will be used to map the left part of the inversion to the GENOME. 1kb by default, but shorter if there is no more sequence in file.
LEN2=`echo $POS2 $LENS | awk '{len=1000;if ($1+1000>$2){len=$2-$1-1};print len}'` #same than prevous but affecting the right part

printf "len1 $LEN1; len2 $LEN2\n"

INV=${SEQ:$POS1:$POS2-$POS1} #get inversion
LEFT=${SEQ:$[POS1-$LEN1]:$LEN1} #get left
RIGHT=${SEQ:$POS2-1:$LEN2} #get right

########################################LEFT
echo "GET LEFT SEQ"
printf ">INVL\n"$LEFT > $TEMPDIR/left.fa #create file

/home/lpantano/invregions/blat -tileSize=18 -stepSize=50 -maxIntron=1000 $TEMPDIR/ref.fa $TEMPDIR/left.fa $TEMPDIR/outl #determine the position in the GENOME

CHRREF=`grep "^[0-9]" $TEMPDIR/outl |  awk '{print $14}'| head -1` #get chr
LEFTREF=`grep "^[0-9]" $TEMPDIR/outl |  awk '{print $17+($11-$13)}'| head -1` #get end position of the flanking region

STARTLEFTREF=`echo $LEFTREF | awk '{start=$1-1000000;if (start<0){start=1};print start}'` #add as much as 1Mb to the left side
printf $CHRREF"\t"$STARTLEFTREF"\t"$LEFTREF | awk '{print $0}' >$TEMPDIR/leftref.bed
bedtools getfasta  -fi $TEMPDIR/ref.fa -bed $TEMPDIR/leftref.bed -fo $TEMPDIR/leftref.fa #get the sequence form the GENOME
LEFTINV=`tail -1 $TEMPDIR/leftref.fa`

########################################RIGHT
echo "GET RIGHT SEQ"
printf ">INVR\n"$RIGHT > $TEMPDIR/right.fa #create file

/home/lpantano/invregions/blat  -tileSize=18 -stepSize=50 -maxIntron=1000 $TEMPDIR/ref.fa $TEMPDIR/right.fa $TEMPDIR/outr #determine the position in the GENOME

CHRREF=`grep "^[0-9]" $TEMPDIR/outr |  awk '{print $14}'| head -1` #get chr
RIGHTREF=`grep "^[0-9]" $TEMPDIR/outr |  awk '{print $16-$12}'| head -1` #get start position of the flanking region

ENDRIGHTREF=`echo $RIGHTREF | awk -v maxref=$maxref '{end=$RIGHTREF+1000000;if (end>maxref){end=maxref };print end}'` #add as much as 1Mb to the tight side
printf $CHRREF"\t"$RIGHTREF"\t"$ENDRIGHTREF | awk '{print $0}'>$TEMPDIR/rightref.bed
bedtools getfasta  -fi $TEMPDIR/ref.fa -bed $TEMPDIR/rightref.bed -fo $TEMPDIR/rightref.fa #get the sequence form the GENOME
RIGHTINV=`tail -1 $TEMPDIR/rightref.fa`

########################################CREATE REGION
echo "GET ALL SEQ"
ALL=$LEFTINV$INV$RIGHTINV
printf ">"$NAME"\n"$ALL > "$REFDIR/$NAME.fa" #join all fragments to create the virtual inversion

echo "GET POSITION INV AT NEW CHR"
/home/lpantano/invregions/blat  -tileSize=18 -stepSize=50 -maxIntron=1000 "$TEMPDIR/realinv.fa" "$REFDIR/$NAME.fa"   "$TEMPDIR"/outinv #determine the position of the inversion in the new FASTA
grep "^[0-9]" "$TEMPDIR"/outinv |  awk '{print $10"\t"$12"\t"$13}'| head -1 >"$REFDIR"/"$NAME".bed #get positions


rm "$TEMPDIR"/*fa*
rm "$TEMPDIR"/*bed
rm "$TEMPDIR"/out*

#echo DO CHECK
#/home/lpantano/invregions/blat -out=blastq $TEMPDIR/oriseq.fa all.fa  check.blast
#all.fa is reference to detect new transcripts
######################################
