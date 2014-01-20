#echo "WARNING: hard coded to size of lattice!"

for i 
do
    echo -n "."
#Get rid of triangle head to line so convert doesn't bork on bender
    cp -a "${i}" "${i}".tmp #Unique names per process to avoid race condition :^)
    grep -v "<marker" "${i}".tmp | sed 's/marker-end="url(#triangle)"//g' > "${i}"
    rm  "${i}".tmp
#turn on background transparent for coalescing into combined movie with colour plots... otherwise white for direct movie viewing
# -background transparent
    convert -antialias -background transparent -density 1200 -scale 500x "${i}" ${i%.*}.png  
done
