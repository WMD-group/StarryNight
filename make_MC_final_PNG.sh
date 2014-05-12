. export_svg.sh MC-SVG_final.svg
convert -scale 1000x1000 -coalesce MC-PNG_final.png MC-SVG_final.png out.png
mv out-1.png MC-SVG_final_COMBINED.png
