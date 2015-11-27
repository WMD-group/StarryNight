# Files look like: 
# T_0001_Strain_1.000000.log               T_0001_Strain_1.000000.png               T_0001_Strain_1.000000_MC-PNG_final.png  T_0001_Strain_1.000000_MC-SVG_final.svg

mkdir movieframes

for T in ` seq -w 0 1000 `
do
    echo -n "."
    echo $T
    lattice="T_${T}_Strain_?.000000_MC-PNG_final.png"
    potential="T_${T}_Strain_?.000000.png"
    convert -scale 480x480 $lattice movieframes/L_${T}.png
    convert -scale 480x480 $potential movieframes/P_${T}.png
    convert +append movieframes/L_${T}.png movieframes/P_${T}.png movieframes/M_${T}.png
    #montage movieframes/L_${T}.png movieframes/P_${T}.png -tile 2x1 -geometry +0+0 movieframes/M_${T}.png # Works, but very slow
    #    convert -scale 500x500 "${file}" "movieframes/${file%.*}.png" #Nb: scale = no interpolation between pixels
done

ffmpeg -framerate 10 -i movieframes/M_%04d.png -s:v 960x480 -c:v libx264 -profile:v high -crf 23 -pix_fmt yuv420p -r 10 TRNM_movie_lattice.mp4
#mplayer -fps 10 mf://movieframes/*.png 
#mencoder -ovc lavc mf://MC-PNG_step*.jpg -o movie.avi
