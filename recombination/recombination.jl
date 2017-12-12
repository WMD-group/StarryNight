module recombination

using StatsBase

# Physical constants / units
const kb=8.6173324E-5 # in units of eV
const ε0=8.854E-12 #Units: C2N−1m−2
const D=3.336E-30 #Debye in SI

const Å=1E-10 # Angstrom

# Specific to MAPI
const r=6.29Å # Sensible figure for MA for consistency; AW
const D1=2.18D #Methyl-Ammonia; b3lyp/6-31g*

#Read in potential...
function starrynight_read_potential(filename)
    lattice=readdlm(filename)
    pot=lattice[:,4] # just potential, in internal units of Starrynight

    T=300           # Temperature; Kelvin
    β=1/(kb*T)      # Thermodynamic Beta; units

    # OK; let's calc Dipole Potential with simple form, contribution from one dipole at distance 1
    factor=1/(4*pi*ε0) * (D1*r)/(r^3)
    @printf STDERR "Conversion factor from Starrynight potential internal units --> V : %f\n" factor
    # Scale from Starrynight potential --> volts
    pot*=factor

    # Apply effective dielectric screening...
    # Calculated value DFT: 5.6 --> 6.5 (dep. on MA alignment)
    # See: APL Mat. 1, 042111 (2013); http://dx.doi.org/10.1063/1.4824147
    # Measured value, Ellipso in NIR: 5.0
    ɛr=4.1
    pot/=ɛr

    N=length(pot) # number of elements

    @printf STDERR "OK; data loaded and converted to Volts. Variance: %e SD: %e\n\n" var(pot) std(pot)
#    print(hist(pot,50)) # I think this used to make a nice little histogram; API has changed

    describe(pot)

    return N, pot
end

#N=1000
#pot=0.13*randn(N) # Normal distribution
#pot=0.13*sin(range(0,N)) # Sinusoidal variations

function calc_recombination(N,pot,kernel)
    # Calculate partition function by sum over state energies
    Zh=sum(kernel(-pot))
    Ze=sum(kernel(pot))
    # Charge densities from partition fn.
    ρh=(1/Zh) * kernel(-pot)
    ρe=(1/Ze) * kernel(pot)

    @printf STDERR "Zh: %e Ze: %e e: %e h: %e\n" Zh Ze sum(ρe) sum(ρh)
    # Factors of N make these values independent of number of elements in the
    # potential, self-consistent, and equal to R=1 for a flat potential.
    @printf STDERR "R=N*N/(Zh*Ze): %e\t" N*N/(Zh*Ze)
    @printf STDERR "R=N*sum(ρe.*ρh): %e\n" N*sum(ρh.*ρe)
    return N*sum(ρh.*ρe)
end

function calc_mobility(N::Int,pot)
    potsorted=sort(reshape(pot,length(pot)))
    trap=potsorted[div(N,4)]-potsorted[1] # Diff in energy between deepest point and first quartile, 25% CDF
    @printf "Trap Depth: %e Bottom(eV): %e 25\%%(eV): %e\n" trap potsorted[1] potsorted[div(N,4)]
    # OK; now need some FD distribution to get trapped vs. free charges.

    T=300           # Temperature; Kelvin
    β=1/(kb*T)      # Thermodynamic Beta; units

    @printf "Boltzmann ratio given trap 0->25\%% trap depth @ 300 K %e \n" exp(trap*β)
end

end 
