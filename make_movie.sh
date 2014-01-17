for file in MC_step*.pnm
do
    echo -n "."
    convert -scale 500x500 "${file}" "${file%.*}.jpg" #Nb: scale = no interpolation between pixels
done

mencoder -ovc lavc mf://MC_step*.jpg -o movie.avi
