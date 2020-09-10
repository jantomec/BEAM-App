#!/bin/zsh

if [[ $1 == "" ]]
then
    echo "No filename. Correct usage: build_mac_gnu.zsh filename.f90"
    exit
fi
echo $1