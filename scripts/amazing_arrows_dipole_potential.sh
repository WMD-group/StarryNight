./starrynight
. export_svg.sh MC-SVG_final.svg
python fourier_transform_potential.py
convert -scale 1000x1000 -coalesce fourier_transform_potential_data.png MC-SVG_final.png out.png
mv out-1.png amazing_arrows_dipole_potential.png

open amazing_arrows_dipole_potential.png
. coalesce_this.sh MC-SVG_final.png
open MC-SVG_final_COMBINED.png

