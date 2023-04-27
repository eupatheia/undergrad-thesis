CC=gcc
% :: %.c
	$(CC) -g -Wall -Wvla -Werror $< -o $@

all: mandelbulb raymarcher test

mandelbulb: mandelbulb.c read_write.c vec.c mat.c shaders.c ray_functions.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable mandelbulb.c read_write.c vec.c mat.c shaders.c ray_functions.c -o $@ -lm -lpthread

raymarcher: raymarcher.c read_write.c vec.c mat.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable raymarcher.c read_write.c vec.c mat.c -o $@ -lm -lpthread

test: test.c vec.c mat.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable test.c vec.c mat.c -o $@ -lm

clean:
	rm mandelbulb raymarcher test