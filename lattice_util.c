

void
empty_lattice ()
{
  int x, y, z;
  for (x = 0; x < X; x++)
    for (y = 0; y < Y; y++)
      for (z = 0; z < Z; z++)
	lattice[x][y][z] = -1; //0;
   //-1 is unoccupied
   //0 and above indicate occupied, with the number being the ID of the snake
}

void
print_lattice ()
{
  int x, y=0, z = 0;
      for (z = 0; z <Z; z++)

    {
      printf ("\n");
   for (x = 0; x < X; x++)
	//for (z = 0; z < Z; z++)
	if (lattice[x][y][z] == -1)
	  printf ("%c[37m.%c[0m",27,27);
	else
	 {	    
//	  printf("o");
	  printf ("%c[%d",27,31+(lattice[x][y][z]%7));
	  if (perc[x][y][z]>0)		      
	      printf(";7");
	  printf("m%c%c[0m",(lattice[x][y][z]%26)+'A',27);
	    
	 }
       
    }
  printf ("\n\n");
}

void
  print_xmakemol ()
{
   int s,seg;
   char mole;
   printf("%d\n\n",num_segments);//total number of snake atoms
   for (s=0;s<num_snakes;s++)
     {
//	printf("molecule\n"); //each snake is new molecule
	if (perc[snakes[s].oroborus[0].x][snakes[s].oroborus[0].y][snakes[s].oroborus[0].z]==1) //electric snake?
	  mole='C';
	else
	  mole='O';
	
	for (seg=0;seg<snakes[s].segs;seg++)
	  printf("%c %.1f %.1f %.1f\n",
		 mole,
		 1.5*(float)snakes[s].oroborus[seg].x,
		 1.5*(float)snakes[s].oroborus[seg].y,
		 1.5*(float)snakes[s].oroborus[seg].z);
//	printf("\n");
	
     }   
}

void
    print_povray (char * name)
{
   
      int s,seg;
   float r,g,b;
      char mole;
      FILE *fo;
   fo=fopen(name,"w");

      fprintf(stderr,"%d\n\n",num_segments);//total number of snake atoms
      for (s=0;s<num_snakes;s++)
     {
	r=g=b=0;
	//      printf("molecule\n"); //each snake is new molecule
         if (perc[snakes[s].oroborus[0].x][snakes[s].oroborus[0].y][snakes[s].oroborus[0].z]==1) //electric snake?
           g=1;
         else
           r=1;

	fprintf(fo,"sphere{< %f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
	                        (float)snakes[s].oroborus[snakes[s].head].x/Z,
	                        (float)snakes[s].oroborus[snakes[s].head].y/Z,
	                        (float)snakes[s].oroborus[snakes[s].head].z/Z,
	                        0.25/Z,
	                        2.0,2.0,2.0);
        fprintf(fo,"sphere{< %f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
	                        (float)snakes[s].oroborus[(snakes[s].head-1)%snakes[s].segs].x/Z,
	                        (float)snakes[s].oroborus[(snakes[s].head-1)%snakes[s].segs].y/Z,
	                        (float)snakes[s].oroborus[(snakes[s].head-1)%snakes[s].segs].z/Z,
	                        0.20/Z,
	                        0.0,0.0,2.0);

         for (seg=snakes[s].head;seg<(snakes[s].segs+snakes[s].head-1);seg++)
	  {

           fprintf(fo,"cylinder{< %f, %f, %f>, <%f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
                  (float)snakes[s].oroborus[seg%snakes[s].segs].x/Z,
                  (float)snakes[s].oroborus[seg%snakes[s].segs].y/Z,
                  (float)snakes[s].oroborus[seg%snakes[s].segs].z/Z,
		  (float)snakes[s].oroborus[(seg+1)%snakes[s].segs].x/Z,
		  (float)snakes[s].oroborus[(seg+1)%snakes[s].segs].y/Z,
		  (float)snakes[s].oroborus[(seg+1)%snakes[s].segs].z/Z,
		  0.1/Z,
		  r,g,b);
	     
	  fprintf(fo,"sphere{< %f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
		 (float)snakes[s].oroborus[seg%snakes[s].segs].x/Z,
		 (float)snakes[s].oroborus[seg%snakes[s].segs].y/Z,
		 (float)snakes[s].oroborus[seg%snakes[s].segs].z/Z,
		 0.15/Z,
		 r,g,b);
		 
 //      printf("\n");
	  } 
	
      }
   fclose(fo);
 } 


void
print_lattice_pnm ()
{
  printf ("P2\n%d %d\n%d\n", X, Y, num_snakes);

  int x, y, z = 0;
  for (x = 0; x < X; x++)
    {
      for (y = 0; y < Y; y++)
	{
	  if (lattice[x][y][z] == -1)
	    printf ("%d\t", num_snakes);
	  else
	    //printf("o");
	    printf ("%d\t", lattice[x][y][z]);
	}

      printf ("\n");
    }
  printf ("\n\n");
}

void save_lattice_file(char * name)
{
   int x,y,z,s,seg;
   FILE *fo;
   fo=fopen(name,"w");
   fprintf(fo,"%d %d %d %d\n",X,Y,Z,num_snakes);
   
// OK; need to access the snake data structure to paint in the dandelions

     for (x=0;x<X;x++)
       for (y=0;y<Y;y++)
	 {
	    fprintf(fo,"\n");
	    for (z=0;z<Z;z++)
	        fprintf (fo,"%d\t", lattice[x][y][z]);
	 }
    for (s=0;s<num_snakes;s++){
    fprintf(fo,"\n");
      for (seg=snakes[s].head;seg<(snakes[s].segs+snakes[s].head-1);seg++)
        fprintf(fo,"%d %d %d\t ",
                    snakes[s].oroborus[seg%snakes[s].segs].x, 
                       snakes[s].oroborus[seg%snakes[s].segs].y, 
                          snakes[s].oroborus[seg%snakes[s].segs].z);
    }

    fprintf(fo,"\n\n%d %f %f\n",TOTAL_SLITHERS,DENSITY,iE[1][1]);
    fclose(fo);
}

void load_lattice_file(char * name)
{
   int x,y,z,tmp;
   FILE *fo;
   fo=fopen(name,"r");
   fscanf(fo,"%d %d %d %d\n",&x,&y,&z,&num_snakes);
   if (x!=X || y!=Y || z!=Z )
     {
	
     fprintf(stderr,"CRASH! INPUT LATTICE PARAMS NOT COMPATIBLE. I READ %d %d %d, I EXPECTED %d %d %d!",x,y,z,X,Y,Z);
     exit(-1);
     }
   
     for (x=0;x<X;x++)
       for (y=0;y<Y;y++)
	 for (z=0;z<Z;z++)
	   fscanf(fo,"%d",&lattice[x][y][z]);
   fclose(fo);
}   

void print_lattice_pnm_file(int i)
{
 FILE *fo;
   char name[30]; //static buffers, security flaws, yada yada
     int x, y, z = 0;
   
   sprintf(name,"hydra_%d.pgm",i);
   fo=fopen(name,"w");   
   fprintf (fo,"P2\n%d %d\n%d\n", X, Y, num_snakes);
   

     for (y = 0; y < Y; y++)
     {
	      for (x = 0; x < X; x++)
	  {
	     
	               if (lattice[x][y][z] == -1)
	                   fprintf (fo,"%d\t", num_snakes);
	               else
	                   //printf("o");
	                   fprintf (fo,"%d\t", lattice[x][y][z]);
	  }
	
	
	      fprintf (fo,"\n");
     }
   
     fprintf (fo,"\n\n");
   fclose(fo);
   
}


void
print_snakes ()
{
  int i;
  for (i = 0; i < num_snakes; i++)
    {
      printf ("Id: %d \tSegments: %d \tHead: %d\n", snakes[i].id,
	      snakes[i].segs, snakes[i].head);

      printf ("Head Loc (x,y,z): %d %d %d\n",
	      snakes[i].oroborus[snakes[i].head].x,
	      snakes[i].oroborus[snakes[i].head].y,
	      snakes[i].oroborus[snakes[i].head].z);

    }
}

void
print_a_snake (int id)
{
  int i;
  printf ("Id: %d \tSegments: %d \tHead: %d\n", snakes[id].id,
	  snakes[id].segs, snakes[id].head);
  for (i = 0; i < snakes[id].segs; i++)
    printf ("Seg: %d (x,y,z) %d %d %d\n", i, snakes[id].oroborus[i].x,
	    snakes[id].oroborus[i].y, snakes[id].oroborus[i].z);
}



