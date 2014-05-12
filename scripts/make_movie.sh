for file in MC-PNG_step*.png
do
    echo -n "."
    convert -scale 500x500 "${file}" "${file%.*}.jpg" #Nb: scale = no interpolation between pixels
done

mplayer -fps 10 mf://MC-PNG_step*.jpg
#mencoder -ovc lavc mf://MC-PNG_step*.jpg -o movie.avi
