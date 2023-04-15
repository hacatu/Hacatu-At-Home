ALL: p833d

p833d: p833d.c
	gcc -Wextra -std=c2x -fno-strict-aliasing -flto -Ofast -march=native -mtune=native p833d.c -o p833d -lprimesieve -lnut -lcrater -lm

