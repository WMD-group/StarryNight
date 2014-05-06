starrynight: starrynight.c
	gcc -O4 -lm -lconfig -o starrynight starrynight.c

starrynight-mac-openmp: starrynight.c
	/usr/local/bin/gcc-4.8 -O4 -lm -lconfig -fopenmp -lgomp -o starrynight starrynight.c

parallel: starrynight
	seq 0 100 1000 | parallel  ./starrynight {}  | sort -k2 -g > T-dep.dat

superparallel: starrynight
	awk 'BEGIN{for (i=0;i<1000;i=i+20) { for (j=0.1;j<3;j=j+0.5) printf ("%f %f\n",i,j); }}' \
		| parallel --colsep ' ' ./starrynight {1} {2}  > aggregate.dat


all: starrynight

clean:
	rm starrynight *.pnm *.jpg *.gif *.avi *.svg *.png
