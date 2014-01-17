#echo "WARNING: hard coded to size of lattice!"

for i 
do
    echo -n "."
# turn on background transparent for coalescing into combined movie with colour plots... otherwise white for direct movie viewing
# -background transparent
    convert -antialias -background transparent -density 1200 -scale 500x "${i}" ${i%.*}.png  
done
