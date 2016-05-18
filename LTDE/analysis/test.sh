#!/usr/bin/env bash

## declare an array variable
declare -a arr1=("element1" "element2" "element3")
declare -a arr2=("R1" "R2")

## now loop through the above array

for i in "${arr1[@]}"
do
    for j in "${arr2[@]}"
    do
    #  echo $i
      #echo $j
      echo "$(pwd)/$i/$j.txt"
    done
   # or do whatever with individual element of the array
  echo ""
done

# You can access them using echo "${arr[0]}", "${arr[1]}" also
