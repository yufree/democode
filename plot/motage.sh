#!/bin/bash
mkdir final
_now=$(date +"%m_%d_%Y_%m_%s")
i=0
for f in *.tiff; do
  convert "$f" miff:-
  ((i++))
done | montage -       \
   -tile 2x2           \
   -geometry 800x800   \
   final/huge_$_now.tiff
 # system('motage.sh')