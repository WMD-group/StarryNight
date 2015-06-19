#Read in potential...

lattice=readdlm(ARGS[1])
pot=lattice[:,4] # just potential, in internal units of Starrynight

T=300           # Temperature; Kelvin

kb=8.6173324E-5 # in units of eV
β=1/(kb*T)      # Thermodynamic Beta; units

Å=1E-10 # Angstrom
r=6.29Å # Sensible figure for MA for consistency; AW
D=3.336E-30 #Debye in SI
ε0=8.854E-12 #Units: C2N−1m−2
D1=2.18D #Methyl-Ammonia; b3lyp/6-31g*

# OK; let's calc Dipole Potential with simple form, contribution from one dipole at distance 1
factor=1/(4*pi*ε0) * (D1*r)/(r^3)
@printf "Conversion factor from Starrynight potential internal units --> V : %f\n" factor
# Scale from Starrynight potential --> volts
pot*=factor

# Apply effective dielectric screening...
ɛr=4.5
pot/=ɛr

N=length(pot) # number of elements

@printf "OK; data loaded and converted to Volts. Variance: %e SD: %e\n\n" var(pot) std(pot)
print(hist(pot,50))

N=1000
#pot=0.13*randn(N) # Normal distribution
#pot=0.13*sin(range(0,N)) # Sinusoidal variations

function recombination(kernel)
    # Calculate partition function by sum over state energies
    Zh=sum(kernel(-pot))
    Ze=sum(kernel(pot))
    # Charge densities from partition fn.
    ρh=(1/Zh) * kernel(-pot)
    ρe=(1/Ze) * kernel(pot)

    @printf "Zh: %e Ze: %e e: %e h: %e\n" Zh Ze sum(ρe) sum(ρh)
    # Factors of N make these values independent of number of elements in the
    # potential, self-consistent, and equal to R=1 for a flat potential.
    @printf "R=N*N/(Zh*Ze): %e\t" N*N/(Zh*Ze)
    @printf "R=N*sum(ρe.*ρh): %e" N*sum(ρh.*ρe)
    println()
    println()
end

println("Boltzmann statistics...")
recombination( En -> exp(-β*En) )
println("Fermi-Dirac statistics...")
recombination( En -> 1./(exp(β*En) + 1.0) )
println("Fermi-Dirac, chem pot = 1.0 eV...")
recombination( En -> 1./(exp(β*(En+1.0)) + 1.0) )

