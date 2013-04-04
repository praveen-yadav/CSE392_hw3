# Makefile: CSE391 HW3 Problem 1
# Luby's maximal independent set algorithm

mis_shared: mis_shared.C
			g++ mis_shared.C -fopenmp -o mis_shared

run:
	@ export OMP_NUM_THREADS=1
	@./mis_shared 

run2:
	@ export OMP_NUM_THREADS=2
	@./mis_shared

