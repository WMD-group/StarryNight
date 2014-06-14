# Load files from Starrynight
load dipoles.xyz
run dipoles.py

#Load bounding box code & draw a box
run drawBoundingBox.py
drawBoundingBox battenberg

#Setup for view / raytracing
set bg_rgb, white
set antialias, 2
set ray_shadows, off

# I hate perspective .'. orthoscopic view
set orthoscopic, on

# From QuteMol example on Pymol Wiki
set light_count,8
set spec_count,1
set shininess, 10
set specular, 0.25
set ambient,0
set direct,0
set reflect,1.5
set ray_shadow_decay_factor, 0.1
set ray_shadow_decay_range, 2
unset depth_cue

