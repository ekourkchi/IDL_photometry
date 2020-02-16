#!/bin/bash


export fname=$1
export index=$2


ind_n=`awk -v ind="$index" '($1 == ind){print(NR)}'  ${fname}`



if [ "$ind_n" = '' ]; then

echo 
echo '**********************************'
echo "Warning from: remove_region.bash"
echo "Could not find such region in "${fname}
echo "Ignored your command !!!!!"
echo '**********************************'
echo 

else

echo 
echo '**********************************'
echo "Message from: remove_region.bash"
echo "Found such region in "${fname}
echo "Please wait .... "
echo '**********************************'
echo

ind_n=`echo $ind_n | awk '{print($1)}'`

# getting the number of all lines that start woth #
hash_no=`awk '($1 == "#"){print(NR)}' ${fname}`

# find the number of # lines
N=0
bol=1
b=0
for i in  $hash_no ; do
((N++))

if [ "$i" -lt  "$ind_n" ]; then
a=$i
fi

if [ "$bol" = 1 ]; then
if [ "$i" -gt  "$ind_n" ]; then
b=$i
bol=0
fi
fi

done


# 
#  echo $N
#  echo $hash_no
#  echo $ind_n, $a, $b
# if there is only two #, then remove the file, only one region
if [ "$N" = 2 ]; then
  
   rm ${fname}
#    echo "N=2"
 
# else, replicate the entire file except the last region, after two last # lines
elif [ "$b" = 0 ]; then

  ((i--))  # since there are two # begning lines for each header
#    echo 'b=0' $i
  rand=$(( ( RANDOM % 100 )  + 1 ))  # a random number between 1-100
  awk -v n="$i" '(NR<n){print($0)}' ${fname}  > ${fname}_${rand}_tmp_remove_region 
  mv ${fname}_${rand}_tmp_remove_region ${fname}

else 

   rand=$(( ( RANDOM % 100 )  + 1 ))  # a random number between 1-100
   awk -v n="$a" '(NR<n-1){print($0)}' ${fname}  > ${fname}_${rand}_tmp_remove_region 
   awk -v n="$b" '(NR>=n){print($0)}' ${fname}  >> ${fname}_${rand}_tmp_remove_region 
   mv ${fname}_${rand}_tmp_remove_region ${fname}
fi


fi


