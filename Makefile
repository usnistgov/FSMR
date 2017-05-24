# the compiler: 
CC = gcc

# compiler flags:
CFLAGS  = -O3 -fopenmp
# linker flags:
LIBS = -lm -lgsl -lgslcblas
# the build target executable:
  
all: fit_jpd fit_rpd gen_jpd jpd2rpd

fit_jpd: fit_jpd.c
	$(CC) $(CFLAGS) -o ./bin/fit_jpd fit_jpd.c $(LIBS)

fit_rpd: fit_rpd.c
	$(CC) $(CFLAGS) -o ./bin/fit_rpd fit_rpd.c $(LIBS)

gen_jpd: gen_jpd.c
	$(CC) $(CFLAGS) -o ./bin/gen_jpd gen_jpd.c $(LIBS)

jpd2rpd: jpd2rpd.c
	$(CC) $(CFLAGS) -o ./bin/jpd2rpd jpd2rpd.c $(LIBS)

