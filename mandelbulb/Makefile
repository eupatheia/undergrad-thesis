CC=gcc
% :: %.c
	$(CC) -g -Wall -Wvla -Werror $< -o $@

all: mandelbulb raymarcher test

mandelbulb: mandelbulb.c read_write.c vec.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable mandelbulb.c read_write.c vec.c -o $@ -lm -lpthread

raymarcher: raymarcher.c read_write.c vec.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable raymarcher.c read_write.c vec.c -o $@ -lm -lpthread

test: test.c vec.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable test.c vec.c -o $@ -lm

clean:
	rm mandelbulb raymarcher test