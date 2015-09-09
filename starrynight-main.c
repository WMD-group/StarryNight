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
#include <stdlib.h>
#include <libconfig.h>

#include "mt19937ar-cok.c" //Code _included_ to allow more global optimisation

// Prototypes...
static int rand_int(int SPAN);

#include "starrynight-config.c" //Global variables & config file reader function  
#include "starrynight-lattice.c" //Lattice initialisation / zeroing / sphere picker fn; dot product
#include "starrynight-analysis.c" //Analysis functions, and output routines

static void gen_neighbour();

static double site_energy(int x, int y, int z, struct dipole *newdipole, struct dipole *olddipole);

static void MC_moves(int moves);
static void MC_move();
static void MC_move_openmp();

int main(int argc, char *argv[])
{
    int i,j,k, x,y; // for loop iterators
    int tic,toc;    // keep track of time for user interface; how many MC moves per second

    char name[100],prefix[100]; 
    char const *LOGFILE = NULL; //for output filenames
    // Yes, I know, 100 chars are enough for any segfault ^_^

    fprintf(stderr,"Starry Night - Monte Carlo brushstrokes.\n");

    fprintf(stderr,"Loading config...\n");
    load_config();

    // Now override with command line options if supplied...
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
    
    // LOGFILE -- If we're going to do some actual science, we better have one...
    sprintf(name,"Recombination_T_%04d_CageStrain_%f.log",T,CageStrain);
    FILE *log;
    LOGFILE=name;
    log=fopen(LOGFILE,"w");
    fprintf(stderr,"Log file '%s' opened. ",LOGFILE);

    //Fire up the twister!
    init_genrand(0xDEADBEEF); //314159265);  // reproducible data :)
    //init_genrand(time(NULL)); // seeded with current time
    fprintf(stderr,"Mersenne Twister initialised...\t");

    gen_neighbour(); //generate neighbour list for fast iteration in energy calculator
    fprintf(stderr,"Neighbour list generated...\t");

    void (*initialise_lattice)() =  & initialise_lattice_antiferro_wall ; // C-function pointer to chosen initial lattice
    // FIXME: C Foo might confuse people? Explain more? Turn into a config
    // option?

    initialise_lattice(); //populate with random dipoles
    fprintf(stderr,"Lattice initialised...");

    solid_solution(); //populate dipole strengths on top of this
    fprintf(stderr,"Solid solution formed...");

    if(CalculateEfield) lattice_Efield_XYZ("initial_lattice_efield.xyz");
    if(CalculateEfield) lattice_Efieldoffset_XYZ("initial_lattice_efieldoffset.xyz");
    if(CalculatePotential) lattice_potential_XYZ("initial_lattice_potential.xyz"); // potential distro
    if(SaveDipolesSVG) outputlattice_svg("initial-SVG.svg");
    if(CalculatePotential) outputpotential_png("initial_pot.png"); //"final_pot.png");
    if(SaveDipolesXYZ) outputlattice_xyz("initial_dipoles.xyz");
 
    if(CalculateRadialOrderParameter) radial_order_parameter();
    // output initialised lattice - mainly for debugging
    if(SaveDipolesPNG) outputlattice_ppm_hsv("initial.png");
    if(CalculatePotential) outputpotential_png("initial_pot.png"); //"final_pot.png");
    //outputlattice_xyz("initial_dipoles.xyz");
    //outputlattice_xyz_overprint("initial_overprint.xyz");
    
    if (DisplayDumbTerminal) outputlattice_dumb_terminal(); //Party like it's 1980
    if (CalculateRecombination) recombination_calculator(stderr);

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
 
//    double AMP; double PHASE;
//    for (AMP=0.01; AMP<=0.05; AMP+=0.01)
//        for (PHASE=0; PHASE<=2*M_PI; PHASE+=M_PI/16) // DOESN'T SAW TOOTH CURRENTLY!
    //    for (T=0;T<500;T+=1) //I know, I know... shouldn't hard code this.
    {
//        Efield.x=AMP*sin(PHASE);

        beta=1/((float)T/300.0); // recalculate beta (used internally) based
//        on T-dep forloop

        for (i=0;i<MCMegaSteps;i++)
        {
/*
// Crazy code to iterate through temperatures; dep on i, with good coverage of
//  range
            // Alright, this is the plan
            // First we take our variable
            // Then we bit reverse it as binary
            // Pretty confusing, but means it will fill in the temperature
            // range with maximum coverage, rather than linear ramping
            unsigned char r=i;
            r=(r&0xF0)>>4 | (r&0x0F)<<4;
            r=(r&0xCC)>>2 | (r&0x33)<<2;
            r=(r&0xAA)>>1 | (r&0x55)<<1;

            T=r*2;
            beta=1/((float)T/300.0);  
*/  

            // Do some MC moves!

//            initialise_lattice(); // RESET LATTICE!
            tic=clock();
            MC_moves(MCMinorSteps);
            toc=clock();

//            fprintf(stderr,"Clocks: tic: %d toc: %d\n",tic,toc);

            // Log some data... Nb: Slow as does a NxN summation of lattice energy
            // contributions!
            //        lattice_potential_log(log);
            //fprintf(log,"%lu %f %f %f\n",ACCEPT+REJECT,lattice_energy(),Efield,Eangle); //FIXME: lattice_energy all broken, data worthless presently.
            // TODO: some kind of dipole distribution? Would I have to bin it
            // myself? (boring.)
            // TODO: Split Energy into different contributions... would be nice to
            // see polarisation delta.E spike when the field flips

            // Log some pretty pictures...
            //        sprintf(name,"MC-PNG_step_%.4d.png",i);
            //        outputlattice_ppm_hsv(name);

            //        sprintf(name,"MC-SVG_step_%.4d.svg",i);
            //        outputlattice_svg(name);


            // Update the (interactive) user what we're up to
            //fprintf(stderr,".");
            //fprintf(stderr,"\n");
            if(DisplayDumbTerminal) outputlattice_dumb_terminal(); //Party like it's 1980
            if(CalculateRecombination) recombination_calculator(log);
            if(CalculateRadialOrderParameter) radial_order_parameter(); // outputs directly to Terminal

            //fprintf(stderr,"Efield: x %f y %f z %f | Dipole %f CageStrain %f K %f\n",Efield.x,Efield.y,Efield.z,Dipole,CageStrain,K);
//            fprintf(stderr,"dipole_fraction: %f T: %d Landau: %f\n",dipole_fraction,T,landau_order());
//            fprintf(stdout,"Moves: %d CageStrain: %f T: %d Landau: %f\n",i*(MCMinorSteps/(X*Y*Z)),CageStrain,T,landau_order());
            //fprintf(stderr,"\n");
            //fprintf(stdout, "T: %d Efield: x %f Polar: %f\n",T,Efield.x,polarisation());
            //fprintf(stderr,"\n");
            fflush(stdout); // flush the output buffer, so we can live-graph / it's saved if we interupt
            
            fprintf(stderr,"MC Moves (per second): %f MHz\n",1e-6*(double)(MCMinorSteps)/(double)(toc-tic)*(double)CLOCKS_PER_SEC);

            sprintf(prefix,"T_%04d_i_%03d_CageStrain_%f",T,i,CageStrain);

            sprintf(name,"%s_lattice_efield.xyz",prefix);
            if(CalculateEfield) lattice_Efield_XYZ(name);

            sprintf(name,"%s_lattice_potential.xyz",prefix);
            if(CalculatePotential) lattice_potential_XYZ(name); // potential distro
    
            sprintf(name,"%s_MC-PNG_final.png",prefix);
            if(SaveDipolesPNG) outputlattice_ppm_hsv(name);

            sprintf(name,"%s_MC-SVG_final.svg",prefix);
            if(SaveDipolesSVG) outputlattice_svg(name);

            sprintf(name,"%s_potential.png",prefix);
            if(CalculatePotential) outputpotential_png(name); //"final_pot.png");

            // Manipulate the run conditions depending on simulation time
            //        if (i==100) { DIM=3;}  // ESCAPE FROM FLATLAND
            //        if (i==200) { Efield.z=1.0;}      // relax back to nothing
            //        if (i==300) {Efield.z=0.0; Efield.x=1.0;}
        }
    } 
    // OK; we're finished...

    fprintf(stderr,"\n");

    // Final data output / summaries.
//    outputlattice_ppm_hsv("MC-PNG_final.png");
//    outputlattice_svg("MC-SVG_final.svg");

    //lattice_potential_log(log);
    //lattice_angle_log(log);

    //    lattice_potential_XY("final_pot_xy.dat");

    //outputlattice_xyz("dipoles.xyz");
    //outputlattice_xyz_overprint("overprint.xyz");
    //outputlattice_pymol_cgo("dipoles.py");

    fprintf(stderr,"Monte Carlo moves - ACCEPT: %lu REJECT: %lu ratio: %f\n",ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\n\n");

    return 0;
}

static int rand_int(int SPAN) // TODO: profile this to make sure it runs at an OK speed.
{
    return((int)( (unsigned long) genrand_int32() % (unsigned long)SPAN));
}

// The following code builds a neighbour list for evaluation of energy;
//   this should result in a speedup as it avoids the for loops + 'ifs' during
//   MC

struct {
    int dx;
    int dy;
    int dz;
    float d;
} neighbours[10000];
int neighbour=0; //count of neighbours

static void gen_neighbour()
{

    int dx,dy,dz=0;
    float d;
 
    for (dx=-DipoleCutOff;dx<=DipoleCutOff;dx++)
        for (dy=-DipoleCutOff;dy<=DipoleCutOff;dy++)
#if(Z>1) //i.e. 3D in Z
            for (dz=-DipoleCutOff;dz<=DipoleCutOff;dz++) //NB: conditional zDipoleCutOff to allow for 2D version
#endif
            {
                if (dx==0 && dy==0 && dz==0)
                    continue; //no infinities / self interactions please!

                d=sqrt((float) dx*dx + dy*dy + dz*dz); //that old chestnut; distance in Euler space

                if (d>(float)DipoleCutOff) continue; // Cutoff in d

                // store precomputed list of neighbours
                neighbours[neighbour].dx=dx; neighbours[neighbour].dy=dy; neighbours[neighbour].dz=dz;
                neighbours[neighbour].d=d;
                neighbour++;
            }
    fprintf(stderr,"Neighbour list generated: %d neighbours with %d DipoleCutOff.\n",neighbour,DipoleCutOff);
}

static double site_energy(int x, int y, int z, struct dipole *newdipole, struct dipole *olddipole)
{
    int dx,dy,dz=0;
    float d;
    double dE=0.0;
    struct dipole *testdipole, n;

    // This now iterates over the neighbour list of neighbours[0..neighbour]
    // Which contains all the precomputed dx,dy,dz for a spherical cutoff, and
    // the sqrt distances etc.

    // Sum over near neighbours for dipole-dipole interaction
    int i;
    #pragma omp parallel for private(dx,dy,dz,d,n) reduction(+:dE) schedule(static,1)
// NB: Works, but only modest speed gains!
    for (i=0;i<neighbour;i++)
    {
        // read in dirn to neighbours + precomputed values
        dx=neighbours[i].dx; dy=neighbours[i].dy; dz=neighbours[i].dz;
        d=neighbours[i].d;

                testdipole=& lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z];

                n.x=(float)dx/d; n.y=(float)dy/d; n.z=(float)dz/d; //normalised diff. vector

    //True dipole like
dE+= (olddipole->length * testdipole->length) * 
    (
        ( dot(newdipole,testdipole) - 3*dot(&n,newdipole)*dot(&n,testdipole) ) -
        ( dot(olddipole,testdipole) - 3*dot(&n,olddipole)*dot(&n,testdipole) ) 
    ) / (d*d*d); 

                // Ferroelectric / Potts model - vector form
                //            dE+= - Dipole * dot(newdipole,testdipole) / (d*d*d)
                //                + Dipole * dot(olddipole,testdipole) / (d*d*d);

                // Now reborn as our cage-strain term!
                if ((dx*dx+dy*dy+dz*dz)==1) //only nearest neighbour
                    dE+= - CageStrain* dot(newdipole,testdipole)
                        + CageStrain * dot(olddipole,testdipole); // signs to energetic drive alignment of vectors (dot product = more +ve, dE = -ve)

            }

    // Interaction of dipole with (unshielded) E-field
    dE+= + dot(newdipole, & Efield)
        - dot(olddipole, & Efield);
    //fprintf(stderr,"%f\n",dot(newdipole, & Efield));

    if (K>0.0) // rarely used anymore; 2D lattice epitaxial strain term
    {
    // along .x projection, squared
    n.x=1.0; n.y=0.0; n.z=0.0;
    dE +=   - K*fabs(dot(newdipole,&n))
        + K*fabs(dot(olddipole,&n));
    // along .y projection, squared
    n.x=0.0; n.y=1.0; n.z=0.0;
    dE +=   - K*fabs(dot(newdipole,&n))
        + K*fabs(dot(olddipole,&n));
    }

    // point charge at centre of space
//    n.x=x-(X/2); n.y=y-(Y/2); n.z=z-(Z/2);
//    dE += 1.0 * (dot(newdipole,&n) - dot(olddipole,&n) ) / ((x-X/2)^2 - (y-Y/2)^2 - (z-Z/2)^2);
    
    return(dE); 
}

static void MC_moves(int moves)
{
    int i;
    //moves/=8; //hard coded domain decomp.
    for (i=0;i<moves;i++)
        MC_move();
}

static void MC_move()
{
    int x, y, z;
    float dE=0.0;
    struct dipole newdipole, *olddipole;

    // Choose random dipole / lattice location

    x=rand_int(X);
    y=rand_int(Y);
    z=rand_int(Z);

    if (lattice[x][y][z].length==0.0) return; //dipole zero length .'. not present

    // random new orientation. 
    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well
    if (ConstrainToX)
        random_X_point(& newdipole); //consider any <100> vector
    else
        random_sphere_point(& newdipole);    

    newdipole.length = lattice[x][y][z].length; // preserve length / i.d. of dipole
    olddipole=& lattice[x][y][z];

    //calc site energy
    dE=site_energy(x,y,z, & newdipole,olddipole);

    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
    {
        lattice[x][y][z].x=newdipole.x;
        lattice[x][y][z].y=newdipole.y;
        lattice[x][y][z].z=newdipole.z;
//      lattice[x][y][z].length=newdipole.length; // never changes with current
//      algorithms.

        ACCEPT++;
    }
    else
        REJECT++;
}

static void MC_move_openmp()
{
    int rx, ry, rz;
    int DDx=2, DDy=2, DDz=2; //domain decomposition in X,Y,Z
    struct dipole newdipole, *olddipole;

    // Choose random dipole / lattice location

    // within segmented domains
    rx=rand_int(X/DDx);
    ry=rand_int(Y/DDy);
    rz=rand_int(Z/DDz);

    //if (lattice[x][y][z].length==0.0) return; //dipole zero length .'. not present
    // random new orientation. 
    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well
    if (ConstrainToX)
        random_X_point(& newdipole); //consider any <100> vector
    else
        random_sphere_point(& newdipole);    

    int Dx,Dy,Dz; //iterature over domains
    #pragma omp parallel for firstprivate(newdipole,olddipole) collapse(3) schedule(static,1)
    for (Dx=0;Dx<DDx;Dx++)
        for (Dy=0;Dy<DDy;Dy++)
            for (Dz=0;Dz<DDz;Dz++)
{
    float dE=0.0;
    int x,y,z;
    // move into our constrained domain
    x=rx+Dx*(X/DDx);
    y=ry+Dy*(Y/DDy);
    z=rz+Dz*(Z/DDz);

    newdipole.length = lattice[x][y][z].length; // preserve length / i.d. of dipole
    olddipole=& lattice[x][y][z];

    //calc site energy
    dE=site_energy(x,y,z, & newdipole,olddipole);

#pragma omp critical
    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
    {
        lattice[x][y][z].x=newdipole.x;
        lattice[x][y][z].y=newdipole.y;
        lattice[x][y][z].z=newdipole.z;
//      lattice[x][y][z].length=newdipole.length; // never changes with current
//      algorithms.

   //     ACCEPT++;
    }
   // else
   //     REJECT++;
}
}
