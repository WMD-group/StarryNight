/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

// Prototypes...
//  These are the backend functions; initialise_lattice is a pointer which
//  points to the chosen one.
void initialise_lattice_random();         // randomly sampled lattice
void initialise_lattice_ferroelectric();  // fully aligned lattice in --X-->
void initialise_lattice_buckled();        // 3D checkboard pattern in +-X,Y,Z
void initialise_lattice_antiferro_wall(); // anti-ferroelectric bi-partition domain; for domain wall creep
void initialise_lattice_ferro_wall();     // ferroelectric bi-partition domain; for domain wall creep
void initialise_lattice_antiferro_slip(); // anti-ferroelectric domains shifted by one unit to one another 
void initialise_lattice_spectrum();       // continuously varying angle dipoles; to test HSV colour output routines
void initialise_lattice_buckled();
void initialise_lattice_slab_delete();    // currently hardcoded to delete a slab of dipoles --> vacuum / surface calculations

void solid_solution();  // apply mix of dipoles / gaps; for Relaxor ferroelectrics...

void initialise_lattice_random()
{
    int x,y,z;
    float angle;

    //Random initial lattice
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                random_sphere_point(& lattice[x][y][z]);
            }
    //Print lattice
    /*
       for (i=0;i<X;i++)
       for (k=0;k<Y;k++)
       printf("\n %f %f %f %f",lattice[i][k].x,lattice[i][k].y,lattice[i][k].z,
       dot(&lattice[i][k],&lattice[i][k]));
       */  
}

void initialise_lattice_ferroelectric()
{
    int x,y,z;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            { lattice[x][y][z].x=1.0; lattice[x][y][z].y=0.0; lattice[x][y][z].z=0.0; }
}

void initialise_lattice_buckled()
{
    int x,y,z;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            { lattice[x][y][z].x=x%2; lattice[x][y][z].y=y%2; lattice[x][y][z].z=z%2; }
}

void initialise_lattice_antiferro_wall()
{
    int x,y,z;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            { 
                if (y<Y/2 ^ x>X/2) // bitwise XOR - to make checkerboard
                { lattice[x][y][z].x=(2.*((z+y)%2))-1.0; lattice[x][y][z].y=0.0; } // modulo arithmathic burns my brain
                else
                { lattice[x][y][z].x=0.0; lattice[x][y][z].y=(2.*((x+z)%2))-1.0; } 
                lattice[x][y][z].z=0.0; 
                //                fprintf(stderr,"Dipole: %d %d %d %f %f %f\n",x,y,z,lattice[x][y][z].x,lattice[x][y][z].y,lattice[x][y][z].z);
            }
}

void initialise_lattice_ferro_wall()
{
    int x,y,z;
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                if (x<X/2)
                    { lattice[x][y][z].x=0.0; lattice[x][y][z].y=-1.0; lattice[x][y][z].z=0.0; }
                else
                    { lattice[x][y][z].x=0.0; lattice[x][y][z].y= 1.0; lattice[x][y][z].z=0.0; }
            }
}

void initialise_lattice_antiferro_slip()
{
    int x,y,z;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            { 
                if (x<X/2)
                { lattice[x][y][z].x=(2.*((z+y)%2))-1.0; lattice[x][y][z].y=0.0; } // modulo arithmathic burns my brain
                else
                { lattice[x][y][z].x=(2.*((z+y+1)%2))-1.0; lattice[x][y][z].y=0.0; } // modulo arithmathic burns my brain
                lattice[x][y][z].z=0.0; 
                //                fprintf(stderr,"Dipole: %d %d %d %f %f %f\n",x,y,z,lattice[x][y][z].x,lattice[x][y][z].y,lattice[x][y][z].z);
            }
}

void initialise_lattice_spectrum()
{
    int x,y,z;
    float angle;

    // initial lattice on spectrum as test
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
                // continous set of dipole orientations to test colour output (should
                // appear as spectrum)
            {
                angle=2*M_PI*(x*X+y)/((float)X*Y); 
                lattice[x][y][z].x = sin(angle);
                lattice[x][y][z].y = cos(angle);
                lattice[x][y][z].z = 0.0;
            }
}

void initialise_lattice_slab_delete()
{
    int x,y,z;
    for (x=0;x<6;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                lattice[x][y][z].x=0.0;
                lattice[x][y][z].y=0.0;
                lattice[x][y][z].z=0.0;
            }
}

void solid_solution()
{
    int x,y,z;
    int i;
    float sample;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                // sample is on {0..1}
                sample=genrand_real1();
                // go through dipoles list; if sample < prevalence, chose this
                // dipole
                // Otherwise step through to next on list, taking away
                // prevalence from this random number
                for (i=0; sample>dipoles[i].prevalence; sample-=dipoles[i].prevalence, i++);
                fprintf(stderr,"SolidSoln: %d %d %d Chosing %f\n",x,y,z,dipoles[i].length);
                // set dipole length to sampled value
                lattice[x][y][z].length=dipoles[i].length;
            }
}

