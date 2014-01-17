# coalesce / merge Inkscape SVG renders + PNG checkerboards into a super duper visual

for i in MC-SVG_step_????.png
do

 convert -scale 1000x1000 -coalesce "${i/SVG/PNG}" "${i}" out.png
 mv out-1.png "${i%.png}_COMBINED.png"
 echo -n "."
done

mplayer mf://*COMBINED.png -fps 10
