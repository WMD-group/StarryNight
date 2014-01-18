echo "WARNING: hard coded to size of lattice!"

LATTICE=50
PIXELS=500

for i 
do
 inkscape --export-area=0:0:${LATTICE}:${LATTICE}    \
 --export-png=${i%.*}.png  \
 --export-width=${PIXELS}                \
 --export-height=${PIXELS} ${i} 
done
