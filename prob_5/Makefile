# Makefile: CSE391 HW3 Problem 1
# Luby's maximal independent set algorithm

nbody: nbody.o particle.o qtree.o point.o
	g++ -Wall -fopenmp -o nbody nbody.o particle.o qtree.o point.o

nbody.o: nbody.C qtree.h point.h particle.h
	g++ -c nbody.C -fopenmp 

particle.o: particle.C particle.h bits_numbers.h
	g++ -c particle.C -fopenmp 

qtree.o: qtree.C qtree.h bits_numbers.h
	g++ -c qtree.C -fopenmp 

point.o: point.C point.h
	g++ -c point.C -fopenmp 


run:
	@ export OMP_NUM_THREADS=1
	@./nbody 

run2:
	@ export OMP_NUM_THREADS=2
	@./nbody

clean:
	@rm *.o

euclid:
	@scp *.C *.h Makefile  euclid:/home/ysakamoto/work/prob_5
