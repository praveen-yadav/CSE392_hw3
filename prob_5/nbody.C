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
#include "particle.h"
#include "qtree.h"

using namespace std;

int compar_vector( const particle left, const particle right)
{
	return left.mt_id<right.mt_id;
}

void build_tree( const double xmin,
				 const double ymin,
				 const int max_pts_per_node,
				 const int max_level,
				 point width,
				 particle* pts,
				 int np,
				 qtree* qt,
				 const int st_pt )
{
	point anch(xmin, ymin);
	
	qt->initialize_root(NULL, 0, anch, max_level,
						max_pts_per_node, width, 0, pts, np, st_pt );

	qt->create_kids();
	
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
	

int write_points( const vector<particle>& pts,
				  const int np )
{
	ofstream ofile;
	ofile.open ("points.dat");

	for (int i=0; i<np; i++)
		ofile<<pts[i].x<<" "<<pts[i].y<<" "<<pts[i].mt_id<<endl;

	ofile.close();

	return 0;
}

void average_trees( qtree* qts, const int nt )
{
	// for(int i=0; i<nt; i++)
	// 	qts[i].show_tree();

	
	for(int i=0; i<4; i++){
		for(int j=0; j<nt; j++){
			if(qts[j].kids[i].idx.size())
				cout<<j<<" "<<i<<" "<<qts[j].kids[i].gid<<endl;
		}
	}	
	
	
	return;
}

int main()
{
	// number of threads
	const int nt = 2;

	// domain
	const double xmin = 0;
	const double xmax = 1;
	const double ymin = 0;
	const double ymax = 1;

	// number of poitns
	const int np = 1000;

	// information necessary for qtree construction
	const int max_level = 4;
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
	for(int i=0; i<np; i++)
		pts[i].gen_coords(xmin, xrange, ymin, yrange);

	// compute morton ids
	#pragma omp parallel for default(none), shared(pts)
	for(int i=0; i<np; i++)
		pts[i].get_morton_id(xmin, ymin, x_grid_size, y_grid_size, max_level);
	
	// parallel merge sort
	vector<particle> tmp;
	tmp.resize(np);
	mergesort<particle>(&pts[0], nt, np, &tmp[0], compar_vector);
	
	// partition ids for parallel tree build
	int* np_local = new int[nt]; // length of partition
	int* st_pt = new int[nt]; // starting point of partition
	for(int i=0; i<nt-1; i++){
		np_local[i] = np/nt;
	}
	np_local[nt-1] = np-np/nt*(nt-1);
	st_pt[0] = 0;
	for(int i=1; i<nt; i++){
		st_pt[i] = st_pt[i-1]+np_local[i-1];
	}

	// pointer for quadtrees
	qtree* qts = new qtree[nt];

	// parallel tree construction
	double start=omp_get_wtime();
	// #pragma omp parallel for shared(pts, np_local)
	for(int i=0; i<nt; i++){
		// build tree
		build_tree( xmin, ymin, max_pts_per_node, max_level, width,
					&pts[st_pt[i]], np_local[i], &qts[i], st_pt[i]);
	}
	double end=omp_get_wtime();
	
	cout<<"wall clock time = "<<end-start<<endl;

	average_trees(qts, nt);
	
	// write_points(pts, np);

	delete[] qts;
	
	return 0;
}

