#!/bin/bash

subc_nums='10 30 50 70 90 110 130 150 170 190 210 230 242'

for subc_num in $subc_nums ; do

fileorig="VHTMUMIMO_emulation_loop_subcarriers_NCG.m"

filenew="VHTMUMIMO_emulation_loop_subcarriers_${subc_num}_NCG.m"

cp $fileorig $filenew
echo $filenew

orig_string="numSTSelected = 150;"
new_string="numSTSelected = ${subc_num};"
echo $new_string

sed -i -e "s/$orig_string/$new_string/g" $filenew

done