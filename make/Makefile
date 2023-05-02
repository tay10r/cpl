CFLAGS := -O3 -ffast-math \
           -Wall -Wextra -Werror -Wfatal-errors -Wconversion -Wstrict-aliasing \
           -std=c89 -ansi -pedantic \
           -fanalyzer -ftree-vectorizer-verbose=2

VPATH := ..

#objs := cpl_bvh.o cpl_math.o
objs := cpl_math.o

.PHONY: all
all: ${objs}

cpl_bvh.o: cpl_bvh.c cpl_bvh.h

cpl_math.o: cpl_math.c cpl_math.h

.PHONY: clean
clean:
	${RM} ${objs}
