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

static double site_energy(int x, int y, int z, struct dipole *newdipole, struct dipole *olddipole);
static void MC_move();

int main(int argc, char *argv[])
{
    int i,j,k, x,y; //for loop iterators

    double P=0.0;

    int tic,toc;

    char name[100];
    char const *LOGFILE = NULL; //for output filenames
    // Yes, I know, 50 chars are enough for any segfault ^_^

    fprintf(stderr,"Starry Night - Monte Carlo brushstrokes.\n");

    fprintf(stderr,"Loading config...\n");
    load_config();

    // Now override with command line options if supplied...
    if (argc>1)
    {
        sscanf(argv[1],"%d",&T);
        fprintf(stderr,"Command line Temperature = %d\n",T); 
    }
/*
    if (argc>1)
    {
        sscanf(argv[1],"%d",&T);
        fprintf(stderr,"Command line temperature: T = %d\n",T);
    }
    if (argc>2)
    {
        sscanf(argv[2],"%lf",&Dipole);
        fprintf(stderr,"Command Line Dipole: Dipole = %lf\n",Dipole);
    }
    */
    sprintf(name,"T_%d_DipoleFraction_%f.log",T,dipole_fraction);

    // If we're going to do some actual science, we better have a logfile...
    FILE *log;
    LOGFILE=name;
    log=fopen(LOGFILE,"w");
    fprintf(stderr,"Log file '%s' opened. ",LOGFILE);

    //Fire up the twister!
    init_genrand(0xDEADBEEF); //314159265);  // reproducible data :)
    //init_genrand(time(NULL)); // seeded with current time
    fprintf(stderr,"Twister initialised. ");

    //initialise_lattice(); //populate wiht random dipoles
    //initialise_lattice_spectrum(); //dipoles to test output routines
    
    //initialise_lattice_antiferro_wall(); //already-paired to test simulator
    initialise_lattice_ferro_wall(); // ferroelectric bi-partition domain; for domain wall creep
    //initialise_lattice_antiferro_slip();

    //initialise_lattice_slab_delete(); 

    fprintf(stderr,"Lattice initialised.");

    fprintf (stderr,"LOCAL ORDER: \n");
    radial_order_parameter();
    // output initialised lattice - mainly for debugging
/*    outputlattice_ppm_hsv("initial.png");
    outputlattice_svg("initial-SVG.svg");
    outputpotential_png("initial_pot.png"); //"final_pot.png");
    outputlattice_xyz("initial_dipoles.xyz");
    outputlattice_xyz_overprint("initial_overprint.xyz");
*/
    outputlattice_dumb_terminal(); //Party like it's 1980

    //lattice_potential_XY("initial_pot_xy.dat"); // potential distro

    fprintf(stderr,"\n\tMC startup. 'Do I dare disturb the universe?'\n");

    fprintf(stderr,"'.' is %d MC moves attempted.\n",MCMinorSteps);

    fprintf(log,"# ACCEPT+REJECT, Efield, Eangle, E_dipole, E_strain, E_field, (E_dipole+E_strain+E_field)\n");

    beta=1/((float)T/300.0);

    //old code - now read in option, so I can parallise externally
    //    for (Efield.x=0.1; Efield.x<3.0; Efield.x+=0.5)
        for (T=0;T<500;T+=1) //I know, I know... shouldn't hard code this.
    {
//        beta=1/((float)T/300.0); // recalculate beta (used internally) based
//        on T-dep forloop

//        for (i=0;i<MCMegaSteps;i++)
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
*/
//            beta=1/((float)T/300.0);  

            // Do some MC moves!

//            initialise_lattice(); // RESET LATTICE!

            //#pragma omp parallel for //SEGFAULTS :) - non threadsafe code everywhere
            tic=time(NULL);
            for (k=0;k<MCMinorSteps;k++) //let's hope the compiler inlines this to avoid stack abuse. Alternatively move core loop to MC_move fn?
                MC_move();
            toc=time(NULL);
 
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
            outputlattice_dumb_terminal(); //Party like it's 1980
            outputlattice_dumb_terminal(); //2nd time so it will calibrate D parameter

//            radial_order_parameter(); // outputs directly to Terminal

            fprintf(stderr,"Efield: x %f y %f z %f | Dipole %f CageStrain %f K %f\n",Efield.x,Efield.y,Efield.z,Dipole,CageStrain,K);
            fprintf(stderr,"dipole_fraction: %f T: %d Landau: %f\n",dipole_fraction,T,landau_order());
//            fprintf(stdout,"Moves: %d CageStrain: %f T: %d Landau: %f\n",i*(MCMinorSteps/(X*Y*Z)),CageStrain,T,landau_order());
            fflush(stdout); // flush the output buffer, so we can live-graph / it's saved if we interupt
            fprintf(stderr,"MC Moves: %f MHz\n",1e-6*(double)(MCMinorSteps*X*Y*Z)/(double)(toc-tic));
                
            sprintf(name,"T_%04d_Strain_%f.log",T,CageStrain);
            lattice_potential_XYZ(name); 

            sprintf(name,"T_%04d_Strain_%f_MC-PNG_final.png",T,CageStrain);
            outputlattice_ppm_hsv(name);

            sprintf(name,"T_%04d_Strain_%f_MC-SVG_final.svg",T,CageStrain);
            outputlattice_svg(name);

            sprintf(name,"T_%04d_Strain_%f.png",T,CageStrain);
            outputpotential_png(name); //"final_pot.png");

            // Manipulate the run conditions depending on simulation time
            //        if (i==100) { DIM=3;}  // ESCAPE FROM FLATLAND
            //        if (i==200) { Efield.z=1.0;}      // relax back to nothing
            //        if (i==300) {Efield.z=0.0; Efield.x=1.0;}

       }
        /*
        // now data collection on equilibriated structure...

        P=0.0;

        for (i=0;i<10;i++)
        {
        P+=polarisation();
        for (k=0;k<MCMinorSteps;k++) //let's hope the compiler inlines this to avoid stack abuse. Alternatively move core loop to MC_move fn?
        MC_move();
        fprintf(stderr,","); 
        }
        // hard coded for loops for Hysterisis exploration
        P/=10;

        double maxfield=Efield.x;
        //    for (maxfield=10.0;maxfield<10.001;maxfield=maxfield+1.0)
        for (i=0;i<0;i++) // hysterisis loop counter
        { 
        for (Efield.x=maxfield;Efield.x>-maxfield;Efield.x-=0.0005)
        {
        fprintf(stderr,"-");
        for (k=0;k<MCMinorSteps;k++)
        MC_move();
        printf("T: %d Efield.x: %f Polar: %f\n",T,Efield.x,polarisation());
        }
        for (Efield.x=-maxfield;Efield.x<maxfield;Efield.x+=0.0005)
        {
        fprintf(stderr,"+");
        for (k=0;k<MCMinorSteps;k++)
        MC_move();
        printf("T: %d Efield.x: %f Polar: %f\n",T,Efield.x,polarisation());
        }
        }

        // P/=(float)MCMegaSteps; //average over our points
        P/=(float)X*Y;          // per lattice site
        // P/=-(float)Efield.x;     // by Electric Field
        // P*=Dipole;
        // See 6.5 (p 167) in Zangwill Modern Electrodynamics

        fprintf(stderr,"NORK! T: %d E: %f P: %f polarisation(per_site): %f\n",T,Efield.x,P,polarisation()/((float)X*Y));
        printf("T: %d Dipole: %f E: %f P: %f polarisation(per_site): %f\n",T,Dipole,Efield.x,P,polarisation()/((float)X*Y));
        */
    } 
    // OK; we're finished...

    fprintf(stderr,"\n");

    // Final data output / summaries.
//    outputlattice_ppm_hsv("MC-PNG_final.png");
//    outputlattice_svg("MC-SVG_final.svg");

    //    lattice_potential_log(log);
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

static double site_energy(int x, int y, int z, struct dipole *newdipole, struct dipole *olddipole)
{
    int dx,dy,dz=0;
    float d;
    double dE=0.0;
    struct dipole *testdipole, n;

    // Sum over near neighbours for dipole-dipole interaction
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

                testdipole=& lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z];

                n.x=(float)dx/d; n.y=(float)dy/d; n.z=(float)dz/d; //normalised diff. vector

                //True dipole like
                dE+=  Dipole * ( dot(newdipole,testdipole) - 3*dot(&n,newdipole)*dot(&n,testdipole) ) / (d*d*d)
                    - Dipole * ( dot(olddipole,testdipole) - 3*dot(&n,olddipole)*dot(&n,testdipole) ) / (d*d*d); 

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

    // interaction with strain of cage modelled as cos^2 function (low energy
    // is diagonal with MA ion along hypotenuse)
    //    dE += + K*cos(2*newangle)*cos(2*newangle)
    //          - K*cos(2*oldangle)*cos(2*oldangle);

    // This is to replicate nice cos^2 (angle) effect in dot products.
    // There must be a more sensible way - if only I could remember my AS
    // double-angle formulae!

    // along .x projection, squared
    n.x=1.0; n.y=0.0; n.z=0.0;
    dE +=   - K*fabs(dot(newdipole,&n))
        + K*fabs(dot(olddipole,&n));
    // along .y projection, squared
    n.x=0.0; n.y=1.0; n.z=0.0;
    dE +=   - K*fabs(dot(newdipole,&n))
        + K*fabs(dot(olddipole,&n));

    return(dE); 
}

static void MC_move()
{
    int x, y, z;
    float d;
    float dE=0.0;
    struct dipole newdipole, *olddipole;

    // Choose random dipole / lattice location

    x=rand_int(X);
    y=rand_int(Y);
    z=rand_int(Z);

    if (lattice[x][y][z].x==0.0 && lattice[x][y][z].y==0.0 && lattice[x][y][z].z==0.0) return; //dipole zero length .'. not present

    // random new orientation. 
    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well
    random_sphere_point(& newdipole);    
    //random_X_point(& newdipole);

    olddipole=& lattice[x][y][z];

    //calc site energy
    dE=site_energy(x,y,z, & newdipole,olddipole);

    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
    {
        lattice[x][y][z].x=newdipole.x;
        lattice[x][y][z].y=newdipole.y;
        lattice[x][y][z].z=newdipole.z;

        ACCEPT++;
    }
    else
        REJECT++;
}

