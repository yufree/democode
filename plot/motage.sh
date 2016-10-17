#!/bin/bash
mkdir final
_now=$(date +"%m_%d_%Y_%m_%s")
i=0
for f in *.png; do
  convert "$f" miff:-
  ((i++))
done | montage -       \
   -tile 2x2           \
   -geometry 800x800   \
   final/huge_$_now.png
 # system('motage.sh')
 
# 用convert命令将图片进行标注，然后用montage组合成图片。其中 -gravity 表示注释位置，与matlab一致，-pointsize表示字体像素， -fill 表示字体颜色。

convert file1.png \
	-fill black \
	-pointsize 28 \
	label:' (A) ' \
	-gravity northwest -geometry +10+12 \
	-composite file1.png
convert file2.png \
	-fill black \
	-pointsize 28 \
	label:' (B) ' \
	-gravity northwest -geometry +10+12 \
	-composite file2.png
convert file3.png \
	-fill black \
	-pointsize 28 \
	label:' (C) ' \
	-gravity northwest -geometry +10+12 \
	-composite file3.png
convert file4.png \
	-fill black \
	-pointsize 28 \
	label:' (D) ' \
	-gravity northwest -geometry +10+12 \
	-composite file4.png

montage file1.tiff file2.tiff file3.tiff file4.tiff\
		-tile 2x2\  
		-geometry 800x800 +1+1 \
		huge.tiff

# 也可以用append去在图片外写标题

convert file1.tiff \
        -fill black \
        -pointsize 80 \
        - label:'A) this is a long title' \
        +swap  -gravity NorthWest -append \
        file1.tiff
convert file2.tiff \
        -fill black \
        -pointsize 80 \
        - label:'B) this is a long title' \
        +swap  -gravity NorthWest -append \
        file2.tiff
convert file3.tiff \
        -fill black \
        -pointsize 80 \
        - label:'C) this is a long title' \
        +swap  -gravity NorthWest -append \
        file3.tiff
convert file4.tiff \
        -fill black \
        -pointsize 80 \
        - label:'D) this is a long title' \
        +swap  -gravity NorthWest -append \
        file4.tiff

montage file1.tiff file2.tiff file3.tiff file4.tiff \
        -tile 2x2  \
        -geometry 800x800 +1+1 \
        huge.tiff

# 这是用 montage 标签的方法加label，不过只能加到图片下方

montage -label 'a) plot1' file1.tiff\
	-label 'b) plot2' file2.tiff\
	-label 'c) plot3(mm)' file3.tiff\
	-label 'd) plot4' file4.tiff\
	-tile 2x2\  
	-geometry 800x800 +1+1 \
	huge.png
