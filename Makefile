starrynight: starrynight.c
	gcc -O4 -lm -lconfig -o starrynight starrynight.c

starrynight-mac-openmp: starrynight
	/usr/local/bin/gcc-4.8 -O4 -lm -lconfig -fopenmp -lgomp -o starrynight starrynight.c

all: starrynight

clean:
	rm starrynight *.pnm *.jpg *.gif *.avi *.svg *.png
