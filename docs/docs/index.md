# StarryNight Documentation 

StarryNight are Monte Carlo codes to simulate dipole-dipole interactions and
(anti-)ferroelectric domains in a hybrid organic-inorganic perovskite solar cell.

This is (currently) achieved with the Metropolish algorithm. The simulation
volume is represented as a 3D lattice with individual sites posessing dipoles. 
The Hamiltonian for the system is then the dipole-dipole interaction (summed to
a cut-off); a local 'cage strain' evaluated with the nearest neighbours
(motivated by ab-initio DFT calculations); and an applied field. 
Random mixtures of dipoles of use specification can be used. 

Though motivated by a mechanistic model of fixed methylammonium dipoles
rotating in 3D space, the same physics is produced by interacting dipoles of
e.g. Jahn-Teller distorted octahedra. 

From the finite temperature Metropolis simulation, an electrostatic potential
can be reconstructed. This can then be used to infer device behaviour
(particularly, recombination and mobility).

The codes are written in C. They currently use a Mersenne-Twister
Pseudo-Random-Number-Generator to power the Metropolis algorithm.

## Publications

- [Molecular ferroelectric contributions to anomalous hysteresis in hybrid perovskite solar cells](http://scitation.aip.org/content/aip/journal/aplmater/2/8/10.1063/1.4890246) APL Materials (2014)
- [The dynamics of methylammonium ions in hybrid organic–inorganic perovskite solar cells](http://www.nature.com/ncomms/2015/150529/ncomms8124/abs/ncomms8124.html) Nature Communications (2015)
- [Role of microstructure in the electron–hole interaction of hybrid lead halide perovskites](http://www.nature.com/nphoton/journal/v9/n10/abs/nphoton.2015.151.html) Nature Photonics (2015)

## Development Notes

2016 - Master has now been moved to the Icarius branch, the 2014-2015 extension
to 3D, solid solutions, many further analysis tools, electrostatic potentials,
Fermi-Dirac/Boltzmann hole/electron populations.

2014-05-31 - Started work on Icarius

2014-01-29
Added dependency on libconfig for lightweight config file parsing.
Mac OSX: `brew install libconfig`
Debian: `sudo apt-get install libconfig-dev`
