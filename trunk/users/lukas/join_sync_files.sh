#!/bin/bash
# usage: join_sync_files.sh infa.sync infb.sync > outfile.sync
# creates temp files with CHRBPS as key column, and joins on all common fields 
# and also only keeping all lines exclusively in the second file with the missing ones from the first replaced by 0:0:0:0:0:0
# padds empty lines from  file 1 with 0s
# to only retrieve common lines remove: -a $F -e "0:0:0:0:0:0" from the join line
# 
INFA=$1
INFB=$2
# use first or second file to join onto
F=2
awk '{key=sprintf("%s%09d",$1,$2); print key,"\t",$0}' $INFA > ${INFA}".temp"
awk '{key=sprintf("%s%09d",$1,$2); print key,"\t",$0}' $INFB > ${INFB}".temp"
wA=`head -1 ${INFA}.temp | wc -w`
wB=`head -1 ${INFB}.temp | wc -w`
#CMD='$a='$wA'; $b='$wB'; $c=""; for $i (2..$a) {$c .=" 1.".$i}; for $i (5..$b){$c.=" 2.".$i}; print $c;'
CMD='$a='$wA'; $b='$wB'; $c=""; for $i (2..4){$c .=" '$F'.".$i} ; for $i (5..$a) {$c .=" 1.".$i}; for $i (5..$b){$c.=" 2.".$i}; print $c;'
#echo $CMD
OUTFORM=$(perl -e "$CMD") 
#echo $OUTFORM
join -j 1 -a $F -e "0:0:0:0:0:0" -o "$OUTFORM"  ${INFA}".temp"  ${INFB}".temp"
# to not pad and only print common lines comment the last and uncomment the next line
#join -j 1 -o "$OUTFORM"  ${INFA}".temp"  ${INFB}".temp"
rm -f ${INFA}".temp"  ${INFB}".temp"

 
