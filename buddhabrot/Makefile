CC=gcc
% :: %.c
	$(CC) -g -Wall -Wvla -Werror $< -o $@

all: buddhabrot color_buddhabrot

buddhabrot: buddhabrot.c read_ppm.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable buddhabrot.c read_ppm.c -o $@ -lm -lpthread

color_buddhabrot: color_buddhabrot.c read_ppm.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable color_buddhabrot.c read_ppm.c -o $@ -lm -lpthread

clean:
	rm buddhabrot color_buddhabrot
