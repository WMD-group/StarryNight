[![DOI](https://zenodo.org/badge/16991744.svg)](https://zenodo.org/badge/latestdoi/16991744)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.org/WMD-group/StarryNight.svg?branch=master)](https://travis-ci.org/WMD-group/StarryNight)

# Starry Night

Monte Carlo codes to simulate dipole-dipole interactions and ferroelectric domains in a hybrid organic-inorganic perovskite solar cell. 

We have started to work on a more complete description and documentation,
available online: [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://wmd-group.github.io/StarryNight/)

![Dipole Domains](dipole-domains.png)

Requirements 
-----------------
A C compiler such as [gcc](https://gcc.gnu.org) is required for the main code, while various scripts and post-processing tools use a combination of [python](https://www.python.org), [julia](https://julialang.org), and [gnuplot](http://www.gnuplot.info).

[libconfig](https://github.com/hyperrealm/libconfig) is used for lightweight config file parsing.

Mac OSX: `brew install libconfig`

Debian: `sudo apt-get install libconfig-dev`

Installation
-----------------
`make` compiles the binary file `starrynight`

`make test` will run a series of test calculations based on the configuration file `starrynight.cfg` in the main directory

Publications
------------
- [Molecular ferroelectric contributions to anomalous hysteresis in hybrid perovskite solar cells](http://scitation.aip.org/content/aip/journal/aplmater/2/8/10.1063/1.4890246) APL Materials (2014)
- [The dynamics of methylammonium ions in hybrid organic–inorganic perovskite solar cells](http://www.nature.com/ncomms/2015/150529/ncomms8124/abs/ncomms8124.html) Nature Communications (2015)
- [Role of microstructure in the electron–hole interaction of hybrid lead halide perovskites](http://www.nature.com/nphoton/journal/v9/n10/abs/nphoton.2015.151.html) Nature Photonics (2015)
- [Dielectric and ferroic properties of metal halide perovskites](https://aip.scitation.org/doi/10.1063/1.5079633) APL Materials (2019)

Development Notes
-----------------
2016 - Extended to 3D, solid solutions, many further analysis tools, electrostatic
potentials, Fermi-Dirac/Boltzmann hole/electron populations

2014-05-31 - Started work on Icarius

2014-01-29 - Added dependency on libconfig
