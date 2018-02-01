#import recombination
include("recombination.jl")

T=300           # Temperature; Kelvin

kb=8.6173324E-5 # in units of eV
β=1/(kb*T)      # Thermodynamic Beta; units

#N=1000
#pot=0.13*randn(N) # Normal distribution
#pot=0.13*sin(range(0,N)) # Sinusoidal variations

N,pot = recombination.starrynight_read_potential(ARGS[1])

function potential_statistics(N,pot)
    println("Boltzmann statistics...")
    recombination.calc_recombination(N, pot, En -> exp.(-β*En) )
    println("Fermi-Dirac statistics...")
    recombination.calc_recombination(N, pot, En -> 1./(exp.(β*En) + 1.0) )
    println("Fermi-Dirac, chem pot = 0.8 eV...")
    recombination.calc_recombination(N, pot, En -> 1./(exp.(β*(En+0.8)) + 1.0) )
    recombination.calc_mobility(N,pot)
end

potential_statistics(N,pot)

# Apply knowledge of polaron Gaussian wavefunction by Gaussian blurring the underlying potential

using Images
n=50# Terrible hard-wired hack
pot=reshape(pot,n,n,n)
#print(pot)

## Make up a pot to test Gaussian filter
#n=10
#pot=zeros(n,n,n)
#pot[n/2,n/2,n/2]=1.0


const l=49 # size of Gaussian Kernel; make this larger to avoid edge effects
# defaults to: l = 4*ceil(Int,σ)+1) ; which is why it had a discretisation error
# c.f.
# https://github.com/JuliaImages/ImageFiltering.jl/blob/f9195cc14f24dc5bcb3d45adb1608a7d7d103018/src/kernelfactors.jl#L267

#print(pot)
for i=0:0.05:10
    println(STDERR,"Gaussian filter with Sigma=$i")
    R=recombination.calc_recombination(N,
                                       imfilter(pot,reflect(Kernel.gaussian([i,i,i], [l,l,l])), "circular"), 
                                       En -> 1./(exp.(β*(En+0.8)) + 1.0) )
    # circular gives periodic boundaries - actually not sure of this, may have broken.
    # reflect() makes it operation a convolution
    # (i,i,i) define the standard deviation of the Gaussian
    # (l,l,l) are the size of the Kernel (stencil); set larger than default to avoid discretisation error
    println("$i $R") # STDOUT, for later printing.
end
