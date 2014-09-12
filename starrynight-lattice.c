/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

// Prototypes...
void initialise_lattice();
void initialise_lattice_wall();
void initialise_lattice_slip();
void initialise_lattice_spectrum();
void initialise_lattice_buckled();
void initialise_lattice_slab_delete();

void initialise_lattice()
{
    int x,y,z;
    float angle;

    //Random initial lattice
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                if (genrand_real1()<dipole_fraction) //occupy fraction of sites...
                    random_sphere_point(& lattice[x][y][z]);
                else
                {lattice[x][y][z].x=0.0; lattice[x][y][z].y=0.0; lattice[x][y][z].z=0.0; }
            }
    //Print lattice
    /*
       for (i=0;i<X;i++)
       for (k=0;k<Y;k++)
       printf("\n %f %f %f %f",lattice[i][k].x,lattice[i][k].y,lattice[i][k].z,
       dot(&lattice[i][k],&lattice[i][k]));
       */  
}

void initialise_lattice_buckled()
{
    int x,y,z;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            { lattice[x][y][z].x=x%2; lattice[x][y][z].y=y%2; lattice[x][y][z].z=z%2; }
}

void initialise_lattice_wall()
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

void initialise_lattice_slip()
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


