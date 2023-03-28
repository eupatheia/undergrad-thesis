CC=gcc
% :: %.c
	$(CC) -g -Wall -Wvla -Werror $< -o $@

all: raymarcher

raymarcher: raymarcher.c read_write.c vec.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable raymarcher.c read_write.c vec.c -o $@ -lm -lpthread

clean:
	rm raymarcher