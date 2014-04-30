import matplotlib.pyplot as plt

import numpy
import math
from IPython import embed #iPython magic for interactive session...

# General setup of figure...
fig=plt.figure()
ax=fig.add_subplot(121)

data = numpy.genfromtxt("final_pot_xy.dat")

nrows, ncols = 100,100
grid=data[:,2].reshape((nrows,ncols))
print data

plt.imshow(grid,extent=(data[:,0].min(), data[:,0].max(), data[:,1].max(), data[:,1].min()),
                   interpolation='nearest', cmap=plt.cm.RdBu) #RdBu) #, cmap=cm.gist_rainbow)
#plt.colorbar()
#plt.show()

ax=fig.add_subplot(122)

fftdata=numpy.fft.fft2(grid)

fftdata=numpy.fft.fftshift(fftdata)
fftdata=fftdata.real
print fftdata 
#grid=fftdata[:,2].reshape((nrows,ncols))

plt.imshow(abs(fftdata),extent=(data[:,0].min(), data[:,0].max(), data[:,1].max(), data[:,1].min()),
                   interpolation='nearest', cmap=plt.cm.PuBuGn) #, cmap=cm.gist_rainbow)
#plt.colorbar()
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
