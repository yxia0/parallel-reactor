SRC = src/simulation_configuration.c src/main.c src/simulation_support.c src/simulation_parallel.c
LFLAGS=-lm -fsanitize=address
CFLAGS=-Og -g -fsanitize=address -Wall -Wextra -Wpedantic -Wshadow

.PHONY: archer2 local build

archer2: CC=cc
archer2: build

cirrus: CC=gcc
cirrus: build

local: CC=mpicc
local: build

build: 
	$(CC) -o reactor $(SRC) $(CFLAGS) $(LFLAGS)

