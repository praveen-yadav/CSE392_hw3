#include <omp.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <bitset>

#include "../msort/msort.h"
#include "point.h"
#include "particle.h"
#include "qtree.h"


using namespace std;

extern int n_th;
extern int* threads_pointer;

int compar_vector( const particle left, const particle right)
{
	return left.mt_id<right.mt_id;
}

void build_tree( const double xmin,
				 const double ymin,
				 const int max_pts_per_node,
				 const int max_level,
				 point width,
				 vector<particle>* pts,
				 int np,
				 qtree* qt,
				 const int st_pt,
				 const int np_proc,
				 const int nt )
{
	point anch(xmin, ymin);
	
	qt->initialize_root( NULL, 0, anch, max_level,
						 max_pts_per_node, width, 0,
						 pts, np, 0, st_pt, np_proc );

	qt->insert_points(nt);
	
	// for(int i=0; i<4;i++){
	// 	for(int j=0; j<4; j++){
	// 		for(int k=0; k<4; k++){
	// 			cout<<qt.kids[i].kids[j].kids[k].idx.size()<<"+";
	// 		}
	// 	}
	// 	cout<<endl<<qt.kids[i].idx.size()<<endl;
	// }
	
	return;
}
	

// read from the file
int read_points( vector<particle>& pts,
				 const int np )
{
	ifstream file_in ("points.dat");
	if (!file_in.is_open()) return 1;

	pts.resize(np);
	
	for(int i=0; i<np; i++){
		file_in>>pts[i].x;
		file_in>>pts[i].y;
		file_in>>pts[i].mt_id;
	}
	
	file_in.close();
	
	return 0;
}


int write_points( const vector<particle>& pts,
				  const int np )
{
	ofstream ofile;
	ofile.open ("points.dat");

	for (int i=0; i<np; i++){
		ofile<<pts[i].x<<" "<<pts[i].y<<" "<<pts[i].mt_id<<" "<<pts[i].m<<endl;
	}
	ofile.close();

	return 0;
}

// not necessary?
void average_trees( qtree* qts, const int nt )
{
		
	return;
}

// evaluate the forces using quadtree
void evaluate_trees()
{
	
	
	
	return;
}


// main function
int main()
{
	// thread nesting enabled
	omp_set_nested(1);
	
	// number of threads
	const int nt = 1;

	// domain
	const double xmin = 0;
	const double xmax = 1;
	const double ymin = 0;
	const double ymax = 1;

	// number of poitns
	const int np = 2000000;
	// const int np = 2000;

	// mass
	const double mmin = 1.0;
	const double mrange = 10;
	
	// information necessary for qtree construction
	const int max_level = 10;
	const int max_pts_per_node = 1;

	// get domain range and grid size
	const double xrange = xmax-xmin;
	const double yrange = ymax-ymin;
	const double x_grid_size = (xmax-xmin)/pow(2.,max_level);
	
	const double y_grid_size = (ymax-ymin)/pow(2.,max_level);
	const point width(xrange, yrange);

	// generate points
	vector<particle> pts;
	pts.resize(np);
	#pragma omp parallel for default(none), shared(pts)
	for(int i=0; i<np; i++){
		if(i<np/30){
			pts[i].gen_coords_cluster(xmin, xrange, ymin, yrange,
				mmin, mrange, xmin+xrange/4,
				ymin+yrange/3, min(xrange,yrange)/5);
		}
		else if(i<np/10) {
			pts[i].gen_coords_cluster(xmin, xrange, ymin, yrange,
				mmin, mrange, xmin+3*xrange/4,
				ymin+2*yrange/3, min(xrange,yrange)/6);
		}
		else{
			pts[i].gen_coords_cluster(xmin, xrange, ymin, yrange,
				mmin, mrange,
				xmin+3*xrange/4, ymin,
				min(xrange,yrange)/6);
		}

		// pts[i].gen_coords(xmin, xrange, ymin, yrange);
	}
	
	// compute morton ids
	#pragma omp parallel for default(none), shared(pts)
	for(int i=0; i<np; i++)
		pts[i].get_morton_id(xmin, ymin, x_grid_size, y_grid_size, max_level);
	
	// parallel merge sort
	vector<particle> tmp;
	tmp.resize(np);
	mergesort<particle>(&pts[0], nt, np, &tmp[0], compar_vector);

	// read points data
	// read_points(pts, np);

	n_th=1;
	threads_pointer = new int*;
	threads_pointer[0] = 1;
	qtree* qt1 = new qtree;
	// parallel tree construction
	double start=omp_get_wtime();
	// #pragma omp parallel for shared(pts, np_local)
	for(int i=0; i<nt; i++){
		// build tree
		build_tree( xmin, ymin, max_pts_per_node, max_level, width,
		&pts, np, qt1, 0, np, 1);
	}
	double end=omp_get_wtime();
	cout<<"wall clock time = "<<end-start<<endl;

	n_th=2;
	qtree* qt2 = new qtree;
	// parallel tree construction
	start=omp_get_wtime();
	// #pragma omp parallel for shared(pts, np_local)
	for(int i=0; i<nt; i++){
		// build tree
		build_tree( xmin, ymin, max_pts_per_node, max_level, width,
		&pts, np, qt2, 0, np, 2);
	}
	end=omp_get_wtime();
	cout<<"wall clock time = "<<end-start<<endl;

	n_th=4;
	qtree* qt4 = new qtree;
	// parallel tree construction
	start=omp_get_wtime();
	// #pragma omp parallel for shared(pts, np_local)
	for(int i=0; i<nt; i++){
		// build tree
		build_tree( xmin, ymin, max_pts_per_node, max_level, width,
		&pts, np, qt4, 0, np, 4);
	}
	end=omp_get_wtime();
	cout<<"wall clock time = "<<end-start<<endl;

	n_th=8;
	qtree* qt8 = new qtree;
	// parallel tree construction
	start=omp_get_wtime();
	// #pragma omp parallel for shared(pts, np_local)
	for(int i=0; i<nt; i++){
		// build tree
		build_tree( xmin, ymin, max_pts_per_node, max_level, width,
		&pts, np, qt8, 0, np, 8);
	}
	end=omp_get_wtime();
	cout<<"wall clock time = "<<end-start<<endl;


	// evaluate_trees();
	
	// write_points(pts, np);

	delete qt1, qt2, qt4, qt8;

	// delete qt1;
	
	return 0;
}

