CC=gcc
% :: %.c
	$(CC) -g -Wall -Wvla -Werror $< -o $@

all: flame transforms

flame: flame.c transforms.c read_write.c quicksort.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable flame.c transforms.c read_write.c quicksort.c -o $@ -lm -lpthread

transforms: transforms.c read_write.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable transforms.c read_write.c -o $@ -lm -lpthread

clean:
	rm flame transforms
