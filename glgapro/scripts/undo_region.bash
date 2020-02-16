#!/bin/bash


export fname=$1


# getting the number of all lines that start woth #
hash_no=`awk '($1 == "#"){print(NR)}' ${fname}`

# find the number of # lines
N=0
for i in  $hash_no ; do
((N++))
done

# echo $N
# echo $hash_no
# if there is only two #, then remove the file, only one region
if [ "$N" = 2 ]; then
  
   rm ${fname}
 
# else, replicate the entire file except the last region, after two last # lines
else

  ((i--))  # since there are two # begning lines for each header
#   echo $i
  rand=$(( ( RANDOM % 100 )  + 1 ))  # a random number between 1-100
  awk -v n="$i" '(NR<n){print($0)}' ${fname}  > ${fname}_${rand}_tmp_undo_region # ${fname}
  mv ${fname}_${rand}_tmp_undo_region ${fname}
 
fi



