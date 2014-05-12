import matplotlib.pyplot as plt

import numpy
import math
from IPython import embed #iPython magic for interactive session...

# General setup of figure...
fig=plt.figure()
ax=fig.add_subplot(111)

# Awk one liner to generate test data (checkerboard):
#  awk 'BEGIN{for (i=0;i<25;i=i+1) { for (j=0;j<25;j=j+1) printf ("%d %d %f\n",i,j,(i+j)%2); }}' > final_pot_xy.dat
# (funky sin strips - use of sin is a bit of an abuse, but gives nice off-diag 2D FFT)
#  awk 'BEGIN{for (i=0;i<25;i=i+1) { for (j=0;j<25;j=j+1) printf ("%d %d %f\n",i,j,sin((i+j)%10)); }}' > final_pot_xy.dat
#data = numpy.genfromtxt("initial_pot_xy.dat") # Typically randomly oriented dipoles - useful to check power spectrum / code
data = numpy.genfromtxt("final_pot_xy.dat")

nrows, ncols = 25,25
grid=data[:,2].reshape((nrows,ncols))

# 6 FOOT HIGH NOTE:
#   The AXES are screwy. Y is reversed; but this is also the case for the PNG output + the SVG output (reversed at the generation stage).
#   Presently everything agrees. BE CAREFUL!

print data

plt.imshow(grid,extent=(data[:,0].min(), data[:,0].max(), data[:,1].max(), data[:,1].min()),
                   interpolation='nearest', cmap=plt.cm.RdBu) #RdBu) #, cmap=cm.gist_rainbow)
#plt.colorbar()
#plt.tight_layout(pad=0.3) #, w_pad=0.5, h_pad=1.0) # Magic incantation for non-terrible plots

# Make plot with no border for overlaying with direct outputs...
plt.axis('off')
plt.gca().set_frame_on(False)
plt.gca().set_xticks([])
plt.gca().set_yticks([])

plt.show()
fig.savefig("fourier_transform_potential_data.png",bbox_inches='tight', pad_inches=0)

# _____ _____ _____
# |  ___|  ___|_   _|
# | |_  | |_    | |
# |  _| |  _|   | |
# |_|   |_|     |_|
#

fig=plt.figure()
ax=fig.add_subplot(111)

fftdata=numpy.fft.fft2(grid)

fftdata=numpy.fft.fftshift(fftdata)
fftdata=fftdata.real
print fftdata 
#grid=fftdata[:,2].reshape((nrows,ncols))

plt.imshow(abs(fftdata),extent=(data[:,0].min(), data[:,0].max(), data[:,1].max(), data[:,1].min()),
                   interpolation='nearest', cmap=plt.cm.PuBuGn) #, cmap=cm.gist_rainbow)
#plt.colorbar()
plt.axis('off')
plt.gca().set_frame_on(False)
plt.gca().set_xticks([])
plt.gca().set_yticks([])

plt.show()
fig.savefig("fourier_transform_potential_transform.png",bbox_inches='tight', pad_inches=0)

  #   ######        ####### ####### #######
 ##   #     #       #       #          #
# #   #     #       #       #          #
  #   #     # ##### #####   #####      #
  #   #     #       #       #          #
  #   #     #       #       #          #
##### ######        #       #          #

fftdata=numpy.fft.rfft2(grid)

#fftdata=numpy.fft.fftshift(fftdata)
fftdata=abs(fftdata.real)
print fftdata 

trace = numpy.sum(fftdata,axis=0)
trace = abs(trace)

print trace

plt.plot(trace)
plt.show()

fig.savefig("fourier_transform_potential_1D_FFT_final.png",bbox_inches='tight', pad_inches=0)

data = numpy.genfromtxt("initial_pot_xy.dat")
nrows, ncols = 100,100
grid=data[:,2].reshape((nrows,ncols))
fftdata=numpy.fft.rfft2(grid)
fftdata=abs(fftdata.real)
trace = numpy.sum(fftdata,axis=0)
trace = abs(trace)
plt.plot(trace)

plt.show()

fig.savefig("fourier_transform_potential_1D_FFT_initial.png",bbox_inches='tight', pad_inches=0)


plt.plot(fftdata)
plt.show()

end

H,xedges,yedges = numpy.histogram2d(data[:,1],data[:,2],bins=36)
H.shape, xedges.shape, yedges.shape
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
#extent = [0,1,0,1] # Fractional coordinates .'. plot 0 to 1
    
#Contours - via http://micropore.wordpress.com/2011/10/01/2d-density-plot-or-2d-histogram/
# - Data are too noisy for this to be useful. Also, they're upside down?! Weird!
#fig.subplots_adjust(bottom=0.15,left=0.15)
#levels = (5.0e1, 4.0e1, 3.0e1, 2.0e1)
#cset = plt.contour(H, levels, origin='lower',colors=['black','green','blue','red'],linewidths=(1.9, 1.6, 1.5, 1.4),extent=extent)
#plt.clabel(cset, inline=1, fontsize=10, fmt='%1.0i')
#for c in cset.collections:
#    c.set_linestyle('solid')

plt.imshow(H,extent=extent,interpolation='nearest')
plt.colorbar()
plt.show()

    # Pb_I_2dhistogram.py
#fig.savefig(datafile.split(".")[0]+'_2dhistogram.png')
plt.close()
