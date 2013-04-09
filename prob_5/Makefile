# Makefile: CSE391 HW3 Problem 1
# Luby's maximal independent set algorithm

nbody: nbody.o particle.o
	g++ -Wall -fopenmp -o nbody nbody.o particle.o

nbody.o: nbody.C
	g++ -c nbody.C -fopenmp 

particle.o: particle.C
	g++ -c particle.C -fopenmp 

run:
	@ export OMP_NUM_THREADS=1
	@./nbody 

run2:
	@ export OMP_NUM_THREADS=2
	@./nbody

