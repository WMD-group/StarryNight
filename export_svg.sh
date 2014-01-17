echo "WARNING: hard coded to size of lattice!"

for i 
do
 inkscape --export-area=0:0:100:100    \
 --export-png=${i%.*}.png  \
 --export-width=1000                \
 --export-height=1000 ${i} 
done
