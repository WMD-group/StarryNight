/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

#include <math.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libconfig.h>

#include "mt19937ar-cok.c" //Code _included_ to allow more global optimisation

// Prototypes...
static int rand_int(int SPAN);

#include "starrynight-config.c" //Global variables & config file reader function  
#include "starrynight-lattice.c" //Lattice initialisation / zeroing / sphere picker fn; dot product
#include "starrynight-analysis.c" //Analysis functions, and output routines
#include "starrynight-montecarlo-core.c" // Core simulation

// this analysis function run before MC moves start.
void analysis_initial()
{
    if(CalculateEfield) lattice_Efield_XYZ("initial_lattice_efield.xyz");
    if(CalculateEfield) lattice_Efieldoffset_XYZ("initial_lattice_efieldoffset.xyz");
    if(CalculatePotential) lattice_potential_XYZ("initial_lattice_potential.xyz"); // potential distro
    if(SavePotentialCube) lattice_potential_cube("initial_lattice_potential.cube");
    if(SaveDipolesSVG) outputlattice_svg("initial-SVG.svg");
    if(CalculatePotential) outputpotential_png("initial_pot.png"); //"final_pot.png");
    if(SaveDipolesXYZ) outputlattice_xyz("initial_dipoles.xyz");

    if(CalculateRadialOrderParameter) radial_order_parameter("rdf.dat");
    // output initialised lattice - mainly for debugging
    if(SaveDipolesPNG) outputlattice_ppm_hsv("initial.png");
    if(CalculatePotential) outputpotential_png("initial_pot.png"); //"final_pot.png");
    //outputlattice_xyz("initial_dipoles.xyz");
    //outputlattice_xyz_overprint("initial_overprint.xyz");

    if (DisplayDumbTerminal) outputlattice_dumb_terminal(); //Party like it's 1980
    if (CalculateRecombination) recombination_calculator(stderr);
}

// this analysis function run every MEGASTEPS, for intermediate data collection
void analysis_midpoint(int MCstep, FILE *log)
{
    char name[100],prefix[100]; 

    // FIXME: JMF 2017-10 - this commented out code should either be made .cfg
    // options, or deleted.

    // Log some data... Nb: Slow as does a NxN summation of lattice energy
    // contributions!
    //        lattice_potential_log(log);
    //fprintf(log,"%lu %f %f %f\n",ACCEPT+REJECT,lattice_energy(),Efield,Eangle); //FIXME: lattice_energy all broken, data worthless presently.
    // TODO: some kind of dipole distribution? Would I have to bin it
    // myself? (boring.)
    // TODO: Split Energy into different contributions... would be nice to
    // see polarisation delta.E spike when the field flips

    // Update the (interactive) user what we're up to
    //fprintf(stderr,".");
    //fprintf(stderr,"\n");
    if(DisplayDumbTerminal) outputlattice_dumb_terminal(); //Party like it's 1980
    if(CalculateRecombination) recombination_calculator(log);
    sprintf(name,"RDF-%.4d.dat",T);
    if(CalculateRadialOrderParameter) radial_order_parameter(name); // appends to file 

    //fprintf(stderr,"Efield: x %f y %f z %f | Dipole %f CageStrain %f K %f\n",Efield.x,Efield.y,Efield.z,Dipole,CageStrain,K);
    //            fprintf(stderr,"dipole_fraction: %f T: %d Landau: %f\n",dipole_fraction,T,landau_order());
    //            fprintf(stdout,"Moves: %d CageStrain: %f T: %d Landau: %f\n",i*(MCMinorSteps/(X*Y*Z)),CageStrain,T,landau_order());
    //fprintf(stderr,"\n");
    //fprintf(stdout, "T: %d Efield: x %f Polar: %f\n",T,Efield.x,polarisation());
    //fprintf(stderr,"\n");

    sprintf(prefix,"T_%04d_i_%03d",T,MCstep); // ie. T_0300_t_002 - for batch runs

    sprintf(name,"%s_efield.xyz",prefix);
    if(CalculateEfield) lattice_Efield_XYZ(name);

    sprintf(name,"%s_potential.xyz",prefix);
    if(CalculatePotential) lattice_potential_XYZ(name); // potential distro; text format (3D)

    sprintf(name,"%s_potential.cube",prefix);
    if(SavePotentialCube) lattice_potential_cube(name);
 
    sprintf(name,"%s_potential.png",prefix);
    if(CalculatePotential) outputpotential_png(name); // 'surface of the moon' picture (2D)

    sprintf(name,"%s_MC-PNG_final.png",prefix); 
    if(SaveDipolesPNG) outputlattice_ppm_hsv(name); // coloured squares (2D)

    sprintf(name,"%s_MC-SVG_final.svg",prefix);
    if(SaveDipolesSVG) outputlattice_svg(name); // arrows (2D)

}

// this analysis function runs at simulation end
void analysis_final()
{
    // FIXME: JMF 2017-10 - Either delete commented code or make a config option!

    // Final data output / summaries.
    //    outputlattice_ppm_hsv("MC-PNG_final.png");
    //    outputlattice_svg("MC-SVG_final.svg");

    //lattice_potential_log(log);
    //lattice_angle_log(log);

    //    lattice_potential_XY("final_pot_xy.dat");

    //outputlattice_xyz("dipoles.xyz");
    //outputlattice_xyz_overprint("overprint.xyz");
    //outputlattice_pymol_cgo("dipoles.py");
}


int main(int argc, char *argv[])
{
    int i,j,k, x,y; // for loop iterators
    int tic,toc,tac;    // keep track of time for user interface; how many MC moves per second

    char name[100],prefix[100]; 
    char const *LOGFILE = NULL; //for output filenames
    // Yes, I know, 100 chars are enough for any segfault ^_^

    fprintf(stderr,"Starry Night - Monte Carlo brushstrokes.\n");

    fprintf(stderr,"Loading config...\n");
    load_config();

    fprintf(stderr,"Checking %d command line options for config file overrides...\n",argc-1);
    // Now override with command line options if supplied...
    //  This enables simple parallel via the Makefile and 'GNU Parallel'
    if (argc>1)
    {
        sscanf(argv[1],"%d",&T);
        fprintf(stderr,"Command line temperature: T = %d\n",T);
    }
    if (argc>2)
    {
        sscanf(argv[2],"%lf",&CageStrain);
        fprintf(stderr,"Command Line CageStrain: CageStrain = %lf\n",CageStrain);
    }

    // Allocate lattice on the heap; which is made up of 'dipole' structs
    fprintf(stderr,"Memory allocation for lattice with X=%d Y=%d Z=%d\n",X,Y,Z);
    lattice = (struct dipole ***)malloc(sizeof(struct dipole **)*X);
    for (x=0;x<X;x++) 
    {
        lattice[x]=(struct dipole **)malloc(sizeof(struct dipole *)*Y);
        for (y=0;y<Y;y++)
            lattice[x][y]=(struct dipole *)malloc(sizeof(struct dipole)*Z);
    } 
    fprintf(stderr,"Lattice allocated");

    // LOGFILE -- If we're going to do some actual science, we better have one...
    sprintf(name,"Recombination_T_%04d.log",T);
    FILE *log;
    LOGFILE=name;
    log=fopen(LOGFILE,"w");
    fprintf(stderr,"Log file '%s' opened. ",LOGFILE);

    //Fire up the twister!
    int SEED=0xDEADBEEF + T; 
    // By adding T to the SEED, the different temperature ensembles have
    // a different starting config, while still being reproducible 
    // Consider time(NULL), for independent runs.
    init_genrand(SEED);  // reproducible data :)
    fprintf(stderr,"Mersenne Twister initialised... seed: %X\t",SEED);
    fprintf(log,"# Starrynight - simulation run on time(NULL)= %ld\n# Mersenne Twister Seed: %X\n",time(NULL),SEED);

    gen_neighbour(); //generate neighbour list (out to dipole cut-off) for fast iteration in energy calculator

//    { random, ferroelectric, buckled, antiferro_wall, ferro_wall, antiferro_slip, spectrum, slab_delete}; 
    // C-function pointer to chosen initial lattice
    void (*initialise_lattice)();
    initialise_lattice= & initialise_lattice_random; // default
    if (strcmp(InitialLattice,"random")==0)
        {initialise_lattice =  & initialise_lattice_random;}
    if (strcmp(InitialLattice,"buckled")==0)
        {initialise_lattice =  & initialise_lattice_buckled;}
    if (strcmp(InitialLattice,"ferroelectric")==0)
        {initialise_lattice =  & initialise_lattice_ferroelectric;} 
    if (strcmp(InitialLattice,"antiferro_wall")==0)
        {initialise_lattice =  & initialise_lattice_antiferro_wall;} 
    if (strcmp(InitialLattice,"ferro_wall")==0)
        {initialise_lattice =  & initialise_lattice_ferro_wall;} 
    if (strcmp(InitialLattice,"antiferro_slip")==0)
        {initialise_lattice =  & initialise_lattice_antiferro_slip;}
    if (strcmp(InitialLattice,"spectrum")==0)
        {initialise_lattice =  & initialise_lattice_spectrum;} 
    if (strcmp(InitialLattice,"slab_delete")==0)
        {initialise_lattice =  & initialise_lattice_slab_delete;} 

    initialise_lattice(); //populate with random dipoles
    fprintf(stderr,"Lattice initialised...");
    solid_solution(); //populate dipole strengths on top of this
    fprintf(stderr,"Solid solution formed...\n");

    if(DisplayDumbTerminal) outputlattice_dumb_terminal(); 
    analysis_initial(); // output initial lattice analysis

    fprintf(stderr,"\n\tMC startup. 'Do I dare disturb the universe?'\n");

    fprintf(stderr,"'.' is %e MC moves attempted.\n",(double)MCMinorSteps);

    beta=1/((float)T/300.0);

    // Equilibriated before Hysterisis scan
    fprintf(stderr,"Equilibriation MC moves... %e\n",(double)MCMinorSteps*(double)MCEqmSteps);
    for (i=0;i<MCEqmSteps;i++)
    {
        fprintf(stderr,",");
        MC_moves(MCMinorSteps);
    }

    if(CalculateEfield) lattice_Efield_XYZ("equilib_lattice_efield.xyz");
    if(SaveDipolesSVG)   outputlattice_svg("equilib-SVG.svg");
    if(CalculatePotential) outputpotential_png("equilib_pot.png"); //"final_pot.png");

    // FIXME: This commented code applies a time-varying field. Either delete
    // or make a config option.

    //    double AMP; double PHASE;
    //    for (AMP=0.01; AMP<=0.05; AMP+=0.01)
    //        for (PHASE=0; PHASE<=2*M_PI; PHASE+=M_PI/16) // DOESN'T SAW TOOTH CURRENTLY!
    //    for (T=0;T<500;T+=1) //I know, I know... shouldn't hard code this.
    {
        //        Efield.x=AMP*sin(PHASE);

        beta=1/((float)T/300.0); 
        // recalculate thermodynamic Beta (used internally), in case T has changed (for loop)

        // Do some MC moves!
        // THIS IS THE SIMULATION CORE / HOT LOOP.
        for (i=0;i<MCMegaSteps;i++)
        {
            //            initialise_lattice(); // RESET LATTICE!
            tic=clock(); // measured in CLOCKS_PER_SECs of a second
            MC_moves(MCMinorSteps);
            toc=clock();

            analysis_midpoint(i,log);
            fflush(stdout); // flush the output buffer, so we can live-graph / it's saved if we interupt

            tac=clock();
            fprintf(stderr,"MC Moves (per second): %f MHz\n",
                    1e-6*(double)(MCMinorSteps)/(double)(toc-tic)*(double)CLOCKS_PER_SEC);
            fprintf(stderr,"Output routines: %f s ; Efficiency of MC moves vs. analysis %.2f\%%\n",
                    (double)(tac-toc)/(double)CLOCKS_PER_SEC,
                    100.0*(double)(toc-tic)/(double)(tac-tic));

            // Manipulate the run conditions depending on simulation time
            //        if (i==100) { DIM=3;}  // ESCAPE FROM FLATLAND
            //        if (i==200) { Efield.z=1.0;}      // relax back to nothing
            //        if (i==300) {Efield.z=0.0; Efield.x=1.0;}
        }
    } 
    // OK; we're finished...

    fprintf(stderr,"\n");

    analysis_final();

    fprintf(stderr,"Monte Carlo moves - ACCEPT: %lu REJECT: %lu ratio: %f\n",ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\n\n");

    return 0;
}

