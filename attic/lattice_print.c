/* Gorgophone
 * A Continuous-Time-Random-Walk Time of Flight simulator for thin film 
 * organic photovoltaic cells, using slithering-snake morphologies generated 
 * by Amphisbaena.
 * 
 * This file started 30th June 2005, by Jarvist Frost.
 * Based on previous code and algorithms in Jenny Nelson & Amanda Chatten's 
 * previous ToF simulator, adapted for snake considerations.
 * 
 * Molecular Electronic Materials and Devices
 * Experimental Solid State Physics
 * Blackett Laboratory
 * Imperial College, London
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

//double occtime[X][Y][Z]; //save occupation time on lattice
//#include "df3.c" //df3 povray density generation
#include "mt19937ar-cok.c" //Mersenne-Twister Random Number generator
#include "amphisbaena.h" //reuse lattice setup from Amphisbaena
#include "lattice_util.c" //reuse lattice utilities / print functions from Amphisbaena

#define MAX_HOPPERS 200 //number of hoppers simultaneously ToF'ing
#define R_MAX 2 //radius of sphere in which potential hops are considered
double BIAS=0.01; //0.01 //potential per site i.e. V/nm

#define RATE 1000000

#define LAMBDA 0.5 //0.5 //lattice reorganisation energy
#define K_BOLTZMANN 8.6173857e-5 //in units of eV
#define T 300.0 //in Kelvin

#define TOF_TIME 8000.0 //time at which to end ToF
#define TOFS 10 //number of ToF's to run to build up stable statistics

#define ntmax 150 //number of time bins
#define ntsample 12 //number of points to sample for moving window gradient calculation / kink detection
//must be an even number :)

#define MAX_RATES 36

int currentevents[ntmax];

double ratetables[X][Y][Z][MAX_RATES];

static double IR_0=(1.0/0.15); //0.15; //inverse of natural length for hopping
//equiv. to 1/r0B in JN code

double THERMAL=(1.0/(4*LAMBDA*K_BOLTZMANN*T));

double simtime=0.0; //current time of simulation
double timefac;

int current=0; //count of current flux 
int escapes=0; //number hoppers escaped electrode

int sphere_bits=0; //number of segments considered with calculation of hopping times
struct coord sphere[(2*R_MAX+1)*(2*R_MAX+1)*(2*R_MAX+1)]; //look-up table for offsets to map sphere

struct hopper 
{
    struct coord loc; //current location on lattice
    double esc; //time of escape
    struct coord dst; //destination of escape
    int hops; //count number of hops
    int id;
} hoppers [MAX_HOPPERS];

int compare_waits(a,b)
    struct hopper *a,*b;
{
    if ( a->esc < b->esc)
        return(-1);
    if(a->esc > b->esc)
        return (1);
    return (0);        
}

print_waits()
{
    int i;
    for (i=0;i<MAX_HOPPERS;i++)
        fprintf(stderr,"Hopper: %d EscTime: %f Loc: %d %d %d Dst: %d %d %d\n",hoppers[i].id,hoppers[i].esc,
                hoppers[i].loc.x,hoppers[i].loc.y,hoppers[i].loc.z,
                hoppers[i].dst.x,hoppers[i].dst.y,hoppers[i].dst.z);   
}

main()
{
    int timestart;
    char name[30];

    timestart=time(NULL);
    init_genrand(timestart); //seed with current time in seconds since 1970 Unix Epoch

    empty_lattice();

    fprintf(stderr,"Enter lattice data to load & ToF...\n");
    scanf("%s",&name);

    fprintf(stderr,"Loading Lattice Data : %s...\n",name);
    load_lattice_file(name);
    //   save_lattice_file("lattice_out.dat");
    fprintf(stderr,"Lattice Loaded.\nReprinting...\n");

    print_povray();
    //   gorgophone();
    //   fprintf(stderr,"Flipping latticce...\n");
    //   flip_lattice();
    //   fprintf(stderr,"Now doing Transient for other material...\n");
    //   gorgophone();  

    //   fprintf(stderr,"Time Taken: %d seconds\n",time(NULL)-timestart);

}


gorgophone()
{
    struct hopper temp;
    int i,tof_count=0,nt,sims,fastesthopper,unsorted_hoppers=0;
    double oldgrad, grad,quickest;
    double Sx,Sy,Sxx,Sxy;
    int NS=0,kinkpassed=0;
    double timebin;
    int tmptime;

    fprintf(stderr,"\n\tEntering Gorgophone, Time of Flight simulator...\n");
    printf("#Hoppers: %d\n#R_Max: %d\n#Bias: %f\n#Lambda: %f\n#T: %f\n#TOF_TIME: %f\n#TOFS: %d\n#IR_0: %f\n#Bins: %d\n",
            MAX_HOPPERS,R_MAX,BIAS, LAMBDA, T, 
            TOF_TIME, (1.0/IR_0), TOFS,ntmax);

    timefac=pow((TOF_TIME),1.0/ntmax);

    construct_shifts(); //fill table of dx,dy,dz for looking around sphere to avoid for loops   

    for (BIAS=0.01;BIAS<=0.1;BIAS+=0.01)
    {

        for (i=0;i<ntmax;i++)
            currentevents[i]=0.0;

        construct_ratetables();

        //   empty_occtime();

        for (sims=0;sims<TOFS;sims++)
        {

            simtime=0.0; current=0; escapes=0;
            fprintf(stderr,"Doing ToF simulation %d of %d\n",sims+1,TOFS);
            empty_perc();
            fprintf(stderr,"Perc lattice emptied...\n");
            expose_lattice(); //expose lattice to light, generating charge carriers
            fprintf(stderr,"Lattice Exposed to light, hoppers generated...\n");
            print_waits();
            qsort(hoppers,MAX_HOPPERS,sizeof(struct hopper),compare_waits);
            fprintf(stderr,"Sorted:\n");
            print_waits();
            //   print_lattice();      

            while (simtime<TOF_TIME)
            {
                /* The hopper queue is rather strange for reasons of super-speedy sortage :)
                 * Its a semi-sorted queue.
                 * That is to say, it starts fully sorted & then has the fastest hoppers [first of the queue] popped
                 * off and replaced with a random assortment.
                 * This random assortment + the first [fastest] member of the sorted queue is then looked over to find
                 * the fastest hopper.
                 * Periodically [when more than n% of the queue is unsorted], the queue is qsorted into its original state.
                 */

                if (unsorted_hoppers>180) //resort the whole thing...
                {
                    qsort(hoppers,MAX_HOPPERS,sizeof(struct hopper),compare_waits);
                    unsorted_hoppers=0;
                }

                quickest=10e20; fastesthopper=0;

                for (i=0;i<=unsorted_hoppers;i++) 
                    if (hoppers[i].esc<quickest)
                    {
                        quickest=hoppers[i].esc;
                        fastesthopper=i;
                    }

                if (hoppers[fastesthopper].esc>10e9) //all hoppers escaped
                    break;

                hop(fastesthopper); //hop the first hopper
                look_around_you(fastesthopper); //allow it to reset its wait time

                unsorted_hoppers++;
                if (unsorted_hoppers>=MAX_HOPPERS) //in case we're not sorting at all for this number of hoppers
                    unsorted_hoppers=MAX_HOPPERS-1;

                if (hoppers[fastesthopper].loc.z==0) //if we've reached the exit electrode
                {     
                    perc[hoppers[fastesthopper].loc.x][hoppers[fastesthopper].loc.y][hoppers[fastesthopper].loc.z]--; //remove hopper from old loc
                    //	     occtime[hoppers[0].loc.x][hoppers[0].loc.y][hoppers[0].loc.z]=-1.0; //set occtime for escaped carrier

                    hoppers[fastesthopper].loc.z=-1; //place off board
                    hoppers[fastesthopper].esc=10e10; //really big (equiv. infinite) escape time
                    escapes++;
                    fprintf(stderr,"Carrier Escape at time %f Current: %d randomwalks: %d\n",simtime,current,tof_count);

                    //	     nt=(int)((float)ntmax*(simtime/TOF_TIME));
                    //	     if (nt<0) nt=0;
                    //	     currentevents[nt]++;
                    //	     printf("Transient: logJ %f logt %f\n",log((float)escapes),log(simtime));
                }

                /*	
                    i=0;
                    while (hoppers[i].esc>hoppers[++i].esc && i<MAX_HOPPERS) //this is a rather dirty bubblesort to 
                //insert the hopper at correct loc in queue
                {
                //	     fprintf(stderr,"Bubble: %d ",i);
                temp= hoppers[i-1]; //this should be done with pointers... very inefficient currently.
                hoppers[i-1]= hoppers[i];
                hoppers[i]=temp;
                }
                */
                //	printf("Stay at %d\n",i);

                //	print_waits();
                //        print_lattice();


                tof_count++; //count number of CTRW's
                if (tof_count%1000000==0)
                    //	  printf(".\n");
                {
                    //	        generate_df3(tof_count/10000);
                    //	        empty_occtime();
                    fprintf(stderr,"tof: %d current: %d simtime: %f hops: %d\t",sims+1,current,simtime,tof_count);
                    fprintf(stderr,"Time Taken: %d s\n",time(NULL)-tmptime);
                    tmptime=time(NULL);
                    //	  print_lattice();

                }

            }
            fprintf(stderr,"\n");   

            printf("#Tof Sim %d Complete %d tof_count\n",sims+1,tof_count);
            /*     for (i=0;i<ntmax;i++)
                   {
                   if (nt>0) timebin=pow(timefac,i)-pow(timefac,i-1);
                   else timebin=pow(timefac,i);

                   printf("#Bin: %f Current: %f\n",
                   pow(timefac,i+0.5)-1.0
                   ,(float)currentevents[i]/(timebin*(sims+1)*MAX_HOPPERS));
                   }
                   */	
        }

        //   print_lattice();
        //   print_waits();

        for (i=0;i<ntmax;i++)
        {
            if (nt>0) timebin=pow(timefac,i)-pow(timefac,i-1);
            else timebin=pow(timefac,i);

            printf("%f %f %f\n",BIAS,
                    pow(timefac,i+0.5)-1.0
                    ,(float)currentevents[i]/(timebin*TOFS*MAX_HOPPERS));
        }


        //kink detection code copied+pasted from JN
        /*  for (nt=ntsample/2;nt<=ntmax-ntsample/2;nt++)
            {	
            Sx=0;Sy=0;Sxx=0;Sxy=0;NS=0;
            for (i=nt-ntsample/2;i<nt+ntsample/2;i++)
            {	     
            if (currentevents[i]>0.0)
            {	  
            NS++;
            Sx+=log(i*(TOF_TIME/ntmax));
            Sy+=log(currentevents[i]);
            Sxx+=log(i*(TOF_TIME/ntmax))*log(i*(TOF_TIME/ntmax));
            Sxy+=log(i*(TOF_TIME/ntmax))*log(currentevents[i]);
            }

            }

            oldgrad=grad;
            grad=(NS*Sxy-Sx*Sy)/(NS*Sxx-Sx*Sx);
        //	if (printallfiles) mobfile <<t[nt]<<'\t'<<grad<<endl;
        if (oldgrad>-1&&grad<=-1&&NS==ntsample&&kinkpassed==0) 
        {

        //	     tkink = (nt-1)*(TOF_TIME/ntmax)+(i*(TOF_TIME/ntmax))*(oldgrad+1)/(oldgrad-grad);
        printf("Kink at: %f with grad: %f\n",(nt-1)*(TOF_TIME/ntmax)+(i*(TOF_TIME/ntmax))*(oldgrad+1)/(oldgrad-grad),grad);
        // cout <<"kink at:"<<'\n'<<tkink<<"next pt:"<<t[nt]<<'\t'<<grad<<endl;
        kinkpassed=1;
        }
        }
        */

    }


    //   print_occtime_pnm();
    //   fprintf(stderr,"\nGenerating snakes.df3 povray density file...");

    //   generate_df3(tof_count);

    fprintf(stderr,"\n\tExit Gorgophone (Time of Flight simulator). %d random walk moves.\n",tof_count);   
}

int hop (int hopper)
{
    int nt;
    current+=(hoppers[hopper].loc.z-hoppers[hopper].dst.z); //change current according to hop

    perc[hoppers[hopper].loc.x][hoppers[hopper].loc.y][hoppers[hopper].loc.z]--; //move hopper from old loc
    //   perc[hoppers[hopper].dst.x][hoppers[hopper].dst.y][hoppers[hopper].dst.z]++; //into destination

    //fill in suitable timebin with current fluctuation
    //nt=(int)((float)ntmax*(simtime/TOF_TIME)); //arithmetic bins

    nt=(int)(log(simtime)/log(timefac));
    //   fprintf(stderr,"Simtime: %f\tnt: %d\n",simtime,nt);
    if (nt<0) nt=0;

    currentevents[nt]+=hoppers[hopper].loc.z-hoppers[hopper].dst.z;

    //   occtime[hoppers[hopper].loc.x][hoppers[hopper].loc.y][hoppers[hopper].loc.z]+=hoppers[hopper].esc-simtime; //count time occupied per lattice site

    hoppers[hopper].loc.x=hoppers[hopper].dst.x; //location becomes destination
    hoppers[hopper].loc.y=hoppers[hopper].dst.y;
    hoppers[hopper].loc.z=hoppers[hopper].dst.z;

    hoppers[hopper].hops++; //count number of hops of hopper

    simtime=hoppers[hopper].esc; //update global simulation time to that of current hop
    //   printf("h");
}


empty_perc()
{
    int x, y, z;
    for (x = 0; x < X; x++) //reset percolation / electrification lattice
        for (y = 0; y < Y; y++)
            for (z = 0; z < Z; z++)
                perc[x][y][z] = 0;
}

/*
   empty_occtime()
   {
   int x, y, z;
   for (x = 0; x < X; x++) //reset percolation / electrification lattice
   for (y = 0; y < Y; y++)
   for (z = 0; z < Z; z++)
   occtime[x][y][z] = 0.0;
   }
   */

expose_lattice()
{
    int i,x,y,z;
    for (i=0;i<MAX_HOPPERS;i++)
    {
        do 
        {
            x=rand_int(X);
            y=rand_int(Y);
            z=Z-1-rand_int(Z/10); //even distribution of charge carriers in first 1/10th of material
        }
        while (lattice[x][y][z]<0 && perc[x][y][z]==0); //while no carrier here already, and snake material present

        fprintf(stderr,"Hooper located at: (x,y,z) %d %d %d\n",x,y,z);
        perc[x][y][z]++; //put carrier on perc lattice
        hoppers[i].loc.x=x; hoppers[i].loc.y=y; hoppers[i].loc.z=z; //let hopper know its location
        hoppers[i].id=i; //identifiy hopper for tracking purposes
        hoppers[i].hops=0; //set hop count to zero

        look_around_you(i); //have the hopper choose its own destiny
    }
}

construct_shifts() //fill table of dx,dy,dz for looking around sphere to avoid for loops
{
    int dx,dy,dz;
    for (dx=-R_MAX;dx<=R_MAX;dx++) //this will be faster with a precomputed table of variables?
        for (dy=-R_MAX;dy<=R_MAX;dy++)
            for (dz=-R_MAX;dz<=R_MAX;dz++)
            {
                if ( ((dx*dx)+(dy*dy)+(dz*dz)) > R_MAX*R_MAX)
                    continue; //if outside our bounding sphere, skip on...
                if (dx==0 && dy==0 && dz==0)
                    continue; //no self hopping
                sphere[sphere_bits].x=dx; sphere[sphere_bits].y=dy; sphere[sphere_bits].z=dz;

                /*	    fprintf(stderr,"Sphere: %d dx: %d dy: %d dz: %d \t%d %d %d\n",
                        sphere_bits,
                        sphere[sphere_bits].x,
                        sphere[sphere_bits].y,
                        sphere[sphere_bits].z,dx,dy,dz);*/
                sphere_bits++;	    
            }   
}

construct_ratetables()
{
    int j,x,y,z,dx,dy,dz;
    double dE;
    double rate;

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                //	    fprintf(stderr,"\nTable %d %d %d\t",x,y,z);
                ratetables[x][y][z][0]=0.0; //used to store 'total rate'    	    
                ratetables[x][y][z][sphere_bits]=0.0; //used to store 'total rate'

                for (j=0;j<sphere_bits;j++) // ~32 loops for 5x5x5 sphere
                    // this is far faster than having 3 nested for loops over the 'sphere' of space we look in
                {
                    dx=sphere[j].x; dy=sphere[j].y; dz=sphere[j].z;

                    if ( (x+dx)>=X || (x+dx)<0 ||
                            (y+dy)>=Y || (y+dy)<0 ||
                            (z+dz)>=Z || (z+dz)<0 ||
                            lattice[x+dx][y+dy][z+dz]<0)
                        //attempting to step outside lattice or not snake material...
                    {
                        ratetables[x][y][z][j+1]=ratetables[x][y][z][j];

                        continue;
                    }

                    //	    fprintf(stderr,"Calculating hop to: (x,y,z) %d %d %d (dx,dy,dz) %d %d %d\n",x+dx,y+dy,z+dz,dx,dy,dz);

                    //SUPERIMPORTANTBIT
                    dE=-LAMBDA-dz*BIAS; //energy required to move
                    rate=RATE 
                        * exp( - ((double)dE*dE*THERMAL) )
                        * exp( - ( sqrt((double)((dx*dx)+(dy*dy)+(dz*dz)))*(double)IR_0) );

                    //	    fprintf(stderr,"[%d %d %d] \tdE: %f\t rate: %f\n",dx,dy,dz,dE,rate);

                    ratetables[x][y][z][j]+=rate;  //was 1/rate 20:30 Aug 18th 2005
                    //	    fprintf(stderr,"Wait: %f\n",1/rate);
                    ratetables[x][y][z][j+1]=ratetables[x][y][z][j];

                    //	    ratetables[x][y][z][sphere_bits+1]+=rate; //actual rate

                    //	    fprintf(stderr,"[%d %d %d]:",dx,dy,dz);	    
                    //	    fprintf(stderr,"%d:%f\t",j,rate);	    
                }
                //	    fprintf(stderr,"Tot:%d:%f",j,ratetables[x][y][z][j]);
            }   
}

look_around_you(int hopper)
{
    int x,y,z,i,trapcount=0;
    double chosen_rate;

    x=hoppers[hopper].loc.x; y=hoppers[hopper].loc.y; z=hoppers[hopper].loc.z;

    //   fprintf(stderr,"Look around you: (x,y,z) (%d,%d,%d) Totalrate: %lf\n",x,y,z,ratetables[x][y][z][sphere_bits]);

    do 
    {

        chosen_rate=ratetables[x][y][z][sphere_bits]*rand_float(); //choose which step to take
        //   fprintf(stderr,"Chosen rate: %f\t",chosen_rate);
        for (i=0;ratetables[x][y][z][i]<chosen_rate;i++);
        if (trapcount++>1000) //trapped by other charges
        {
            //set up false hop to self in order to give chance to recalculate hop once other charges have moved away
            hoppers[hopper].dst.x=x; hoppers[hopper].dst.y=y; hoppers[hopper].dst.z;

            break;


        }

    }
    while (perc[x+sphere[i].x][y+sphere[i].y][z+sphere[i].z]>0);



    hoppers[hopper].dst.x=x+sphere[i].x;
    hoppers[hopper].dst.y=y+sphere[i].y;
    hoppers[hopper].dst.z=z+sphere[i].z;

    perc[hoppers[hopper].dst.x]
        [hoppers[hopper].dst.y]
        [hoppers[hopper].dst.z]++; //ghost in the machine placed on lattice; intended dest

    //   fprintf(stderr,"Hopper: %d Hopper_id: %d Choosen option %d, dst: %d %d %d\n",hopper,hoppers[hopper].id,i,deltas[i].x,deltas[i].y,deltas[i].z);

    hoppers[hopper].esc=simtime-log(rand_float())/ratetables[x][y][z][sphere_bits]; //let hopper know when to escape...

    //   fprintf(stderr,"Simtime: %f Escape Time: %f\n",simtime,hoppers[hopper].esc);

}

/*
   print_occtime_pnm ()
   {
   double max=0;
   int x=0, y=0, z = 0;
   int maxpix=255*255;

   for (z = 0; z < Z; z++)
   for (y = 0; y < Y; y++)
   if (log(1.0+occtime[x][y][z])>max) max=log(1.0+occtime[x][y][z]);

   printf ("P3\n%d %d\n%d\n", Z, Y, maxpix);

   for (z = 0; z < Z; z++)
   {
   for (y = 0; y < Y; y++)
   {
   if (lattice[x][y][z] == -1)
   printf ("%d %d %d\t",0,0,maxpix);
   else
//printf("o");
printf ("%d %d %d\t", (int)(maxpix*(log(1.0+occtime[x][y][z])/max)),0,0);
}
printf ("\n");
}
printf ("\n\n");
}
*/
