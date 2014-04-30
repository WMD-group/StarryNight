echo "WARNING: hard coded to size of lattice!"

# Sensible Linux
INKSCAPE="inkscape"
# MAC OS X
INKSCAPE="/Applications/Inkscape.app/Contents/Resources/script"
PWD=` pwd `


LATTICE=50
PIXELS=500

for i 
do
 ${INKSCAPE} --export-area=0:0:${LATTICE}:${LATTICE}    \
 --export-png=${PWD}/${i%.*}.png  \
 --export-width=${PIXELS}                \
 --export-height=${PIXELS} ${PWD}/${i} 
done
