DIRNAME=`dirname $BASH_SOURCE`
source $DIRNAME"/include.sh"

NAME=$1
POS=$2
REFDIR=$3
ANNFILE=$4
TEMPDIR=$5
OUTDIR=$6
INVFILE=$REFDIR/$NAME.bed

#masking positions
CHR=`grep -w $NAME $POS | cut -f 1`
START=`grep -w $NAME $POS | cut -f 2`
END=`grep -w $NAME $POS | cut -f 3`

#change inverted region coordinates to genomic coordinates
awk -v start=$START -v end=$END -v chr=$CHR '{s=$2+start;e=$3+start;print chr"\t"s"\t"e"\t"}' $INVFILE >  $TEMPDIR/TMPinvregion.bed
#sed -i 's/$NAME/$CHR/'  $TEMPDIR/TMPinvregion.bed
#get outside regions 
bedtools subtract -a $REFDIR/mask.bed -b $TEMPDIR/TMPinvregion.bed > $TEMPDIR/TMPousides.bed
bedtools intersect -f 1 -a $ANNFILE -b $TEMPDIR/TMPousides.bed | awk '$3=="exon"' > $TEMPDIR/TMPgenes.flanko.gtf
#remove start position to genes.flank
awk -v start=$START 'FS="\t" {s=$4-start+1;e=$5-start+1;print  $1"\t"$2"\t"$3"\t"s"\t"e"\t"$6"\t"$7"\t"$8"\t"$9}' $TEMPDIR/TMPgenes.flanko.gtf >$TEMPDIR/TMPgenes.shift.flanko.gtf
#change chr name
sed -i s/$CHR/$NAME/ $TEMPDIR/TMPgenes.shift.flanko.gtf

#get inner exons
bedtools intersect -f 1 -a $ANNFILE -b $TEMPDIR/TMPinvregion.bed | awk '$3=="exon"' > $TEMPDIR/TMPgenes.flanki.gtf
#inverted region
STARTI=`cut -f 2 $INVFILE`
ENDI=`cut -f 3 $INVFILE`
#add end_inv-exon_pos + start invregion
awk -v start=$START -v end=$END -v starti=$STARTI -v endi=$ENDI 'FS="\t" {s=((start+endi)-$4)+starti+1;e=((start+endi)-$5)+starti+1;c="+";if ($7=="+"){c="-"};print  $1"\t"$2"\t"$3"\t"e"\t"s"\t"$6"\t"c"\t"$8"\t"$9}' $TEMPDIR/TMPgenes.flanki.gtf >$TEMPDIR/TMPgenes.shift.flanki.gtf
#change chromosome
sed -i s/$CHR/$NAME/ $TEMPDIR/TMPgenes.shift.flanki.gtf

cat $TEMPDIR/TMPgenes.shift.flanki.gtf $TEMPDIR/TMPgenes.shift.flanko.gtf > $OUTDIR/$NAME.knowngenes.gtf

rm $TEMPDIR/TMP*