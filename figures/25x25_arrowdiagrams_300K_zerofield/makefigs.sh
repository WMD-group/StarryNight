
./starrynight
. export_svg.sh MC-SVG_step_0000*.svg MC-SVG_step_0049*.svg 
. coalesce.sh
cp -a MC-SVG_step_0000_COMBINED.png MC-SVG_step_0049_COMBINED.png ~/REPOSITORY/jarv_writings/StarryNight-Paper/figures/25x25_arrowdiagrams_300K_zerofield/

gnuplot plot_potential.gpt
cp -a initial_pot_xy.png final_pot_xy.png ~/REPOSITORY/jarv_writings/StarryNight-Paper/figures/25x25_arrowdiagrams_300K_zerofield/
