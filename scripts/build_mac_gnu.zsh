#!/bin/zsh

if [[ $1 == "" ]]
then
    echo "\nNo filename. Correct usage: build_mac_gnu.zsh filename.f90\n"
    exit
fi

file="${1##*/}"
filename="${file%.f90}"

echo "\nLoading LAPACK libraries\n"
LIST2=""
for file in lapack/lib/*.o
do
  echo "Loading" $file
  LIST2=$LIST2" "$file
done

echo "\n============================\n"
echo "Loading BEAM App libraries\n"
LIST1=""
for file in lib/*.o
do
  echo "Loading" $file
  LIST1=$LIST1" "$file
done

echo "\n============================\n"
echo "Building program " $1
comm="gfortran -Ilib -Ilapack/lib "$1$LIST2$LIST1" -o "$filename

eval $comm

echo "Finished building program\n"
