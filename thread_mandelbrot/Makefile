CC=gcc
% :: %.c
	$(CC) -g -Wall -Wvla -Werror $< -o $@

all: thread_mandelbrot

thread_mandelbrot: thread_mandelbrot.c mandelbrot_functions.c read_ppm.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable thread_mandelbrot.c mandelbrot_functions.c read_ppm.c -o $@ -lm -lpthread

clean:
	rm thread_mandelbrot
