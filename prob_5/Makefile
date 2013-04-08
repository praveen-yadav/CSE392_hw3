# Makefile: CSE391 HW3 Problem 1
# Luby's maximal independent set algorithm

nbody: nbody.C
			g++ nbody.C -fopenmp -o nbody

run:
	@ export OMP_NUM_THREADS=1
	@./nbody 

run2:
	@ export OMP_NUM_THREADS=2
	@./nbody

