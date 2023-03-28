CC=gcc
% :: %.c
	$(CC) -g -Wall -Wvla -Werror $< -o $@

all: raytracer

raytracer: raytracer.c read_write.c vec.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable raytracer.c read_write.c vec.c -o $@ -lm -lpthread

clean:
	rm raytracer