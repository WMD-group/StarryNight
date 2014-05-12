# coalesce / merge Inkscape SVG renders + PNG checkerboards into a super duper visual

for i in Dipole_pot_xy_T*MC-SVG_final.png
do

 convert -scale 500x500 -coalesce "${i/SVG/PNG}" "${i}" out.png
 mv out-1.png "${i%.png}_COMBINED.png"
 echo -n "."
done

mplayer mf://*COMBINED.png -fps 10
