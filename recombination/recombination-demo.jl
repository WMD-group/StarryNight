import recombination

T=300           # Temperature; Kelvin

kb=8.6173324E-5 # in units of eV
β=1/(kb*T)      # Thermodynamic Beta; units

#N=1000
#pot=0.13*randn(N) # Normal distribution
#pot=0.13*sin(range(0,N)) # Sinusoidal variations

N,pot = recombination.starrynight_read_potential(ARGS[1])

println("Boltzmann statistics...")
recombination.calc_recombination(N, pot, En -> exp(-β*En) )
println("Fermi-Dirac statistics...")
recombination.calc_recombination(N, pot, En -> 1./(exp(β*En) + 1.0) )
println("Fermi-Dirac, chem pot = 1.0 eV...")
recombination.calc_recombination(N, pot, En -> 1./(exp(β*(En+1.0)) + 1.0) )

recombination.calc_mobility(N,pot)

