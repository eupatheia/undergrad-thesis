CC=gcc
% :: %.c
	$(CC) -g -Wall -Wvla -Werror $< -o $@

all: julia

julia: julia.c read_ppm.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable julia.c read_ppm.c -o $@ -lm

clean:
	rm julia
