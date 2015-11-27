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
    recombination.calc_recombination(N, pot, En -> exp(-β*En) )
    println("Fermi-Dirac statistics...")
    recombination.calc_recombination(N, pot, En -> 1./(exp(β*En) + 1.0) )
    println("Fermi-Dirac, chem pot = 1.0 eV...")
    recombination.calc_recombination(N, pot, En -> 1./(exp(β*(En+1.0)) + 1.0) )
    recombination.calc_mobility(N,pot)
end

potential_statistics(N,pot)

using Images
n=20 # Terrible hard-wired hack
pot=reshape(pot,n,n,n)
#print(pot)

## Make up a pot to test Gaussian filter
#n=10
#pot=zeros(n,n,n)
#pot[n/2,n/2,n/2]=1.0


#print(pot)
for i=0:0.05:10
    println(STDERR,"Gaussian filter with Sigma=$i")
    R=recombination.calc_recombination(N,imfilter_gaussian(pot,[i,i,i]), En -> 1./(exp(β*(En+1.0)) + 1.0) )
    println("$i $R")
end
