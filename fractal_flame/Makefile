CC=gcc
% :: %.c
	$(CC) -g -Wall -Wvla -Werror $< -o $@

all: flame functions

flame: flame.c functions.c read_write.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable flame.c functions.c read_write.c -o $@ -lm -lpthread

functions: functions.c read_write.c
	$(CC) -g -Wall -Wvla -Werror -Wno-unused-variable functions.c read_write.c -o $@ -lm -lpthread

clean:
	rm flame functions
