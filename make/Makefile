CFLAGS := -O3 -ffast-math -Wall -Wextra -Werror -Wfatal-errors -Wconversion -std=c99 -fanalyzer

VPATH := ..

objs := cpl_bvh.o

.PHONY: all
all: ${objs}

cpl_bvh.o: cpl_bvh.c cpl_bvh.h

.PHONY: clean
clean:
	${RM} ${objs}
