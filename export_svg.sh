for i 
do
 inkscape --export-area-drawing     \
 --export-png=${i%.*}.png  \
 --export-width=500                \
 --export-height=500 ${i} 
done
