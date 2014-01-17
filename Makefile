starrynight: starrynight.c
	gcc -O4 -lm -o starrynight starrynight.c

all: starrynight

clean:
	rm starrynight *.pnm *.jpg *.gif *.avi *.svg *.png
