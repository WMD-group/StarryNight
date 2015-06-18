#Read in potential...

lattice=readdlm(ARGS[1])
pot=lattice[:,4] # just potential, in eV

kb=8.6173324E-5 # in units of eV
β=1/(kb*300)    # Thermodynamic Beta; units

Å=1E-10
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

function partition(f)
    Zh=sum(f(-pot))
    Ze=sum(f(pot))
    ρh=(1/Zh) * f(-pot)
    ρe=(1/Ze) * f(pot)

    @printf "Zh: %e Ze: %e e: %e h: %e\n" Zh Ze sum(ρe) sum(ρh)
    # Factors of N make these values independent of number of elements in the
    # potential, self-consistent, and equal to R=1 for a flat potential.
    @printf "R=N*N/(Zh*Ze): %e\t" N*N/(Zh*Ze)
    @printf "R=N*sum(ρe.*ρh): %e" N*sum(ρh.*ρe)
    println()
    println()
end

println("Boltzmann statistics...")
partition( x -> exp(-β*x) )
println("Fermi-Dirac statistics...")
partition( x -> 1./(exp(β*x) + 1.0) )
println("Fermi-Dirac, chem pot = 1.0 eV...")
partition( x -> 1./(exp(β*(x+1.0)) + 1.0) )

