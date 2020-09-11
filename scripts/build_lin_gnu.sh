if [ "$1" = "" ]
then
  echo "No filename. Correct usage: bash scripts/build_lin_gnu.sh filename.f90"
  exit
fi

file="${1##*/}"
filename="${file%.f90}"

LIST1=""
for file in lib/*.o
do
  LIST1=$LIST1" "$file
done

LIST2=""
for file in lapack/lib/*.o
do
  LIST2=$LIST2" "$file
done

gfortran -Ilib -Ilapack/lib $1 $LIST2 $LIST1 -o $filename".exe"

echo "Finished building program"
