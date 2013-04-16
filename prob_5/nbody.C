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

// keep track of the number of threads at work
int n_th;

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

int write_boxes( const qtree* qt, ofstream& ofile, int parent )
{
	// only make ofstream instance for parent node
	if(parent){
		ofile.open ("boxes.dat");
	}

	// qt is not leaf
	if (!qt->isleaf){
		for(int i=0; i<4; i++){
			write_boxes(&(qt->kids[i]), ofile, 0);
		}
	}
	// qt is leaf
	else{
		ofile<<qt->level<<" "
			 <<qt->anchor.x<<" "<<qt->anchor.y<<" "
			 <<qt->width.x<<" "<<qt->width.y<<endl;
	}

	// only close ofstream for parent node
	if(parent){
		ofile.close();
	}

	return 0;
}

// not necessary?
void average_trees( qtree* qts, const int nt )
{
		
	return;
}

// get which kid has the node with id
qtree* get_local_node(const int id, qtree* qt, const int level )
{
	if ((level==0) || (qt->isleaf))
		return qt;

	
	for(int i=0; i<4; i++){
		for(int j=0; j<qt->kids[i].idx.size(); j++){
			if (id==qt->kids[i].idx[j]){
				return get_local_node(id, &(qt->kids[i]), level-1);
			}
		}
	}
	
	return NULL;
}


// check if the target and source are well-separated
int well_separated( qtree* target, qtree* source )
{
	double xdiff = abs((source->anchor.x) - (target->anchor.x));
	double ydiff = abs((source->anchor.y) - (target->anchor.y));

	if ((xdiff>target->width.x) || (ydiff>target->width.y) )
		return 1;
	
	return 0;
}

// evaluate the forces using quadtree
void evaluate_trees( const int id, const int level,
					 qtree* target, qtree* source, vector<particle>& pts)
{
	// get qtree pointer for the kid that has the id
	target = get_local_node(id, target, 1);

	// cout<<"x "<<pts[id].x<<" y "<<pts[id].y<<endl;
	// cout<<"gid "<<target->gid<<" level "<<target->level<<endl;

	// now compare the target with the kids of the source
	for(int i=0; i<4; i++){
		// not separated enough
		if(!well_separated(target, &(source->kids[i]))){
			// but the kid is a leaf
			if(source->kids[i].isleaf){
				cout<<target->gid<<" and "<<source->kids[i].gid
					<<" need direct evaluation on level "
					<<level+1<<endl;	// cout<<"leaf"<<endl;
			}
			// not leaf
			else{
				// cout<<"not separated enough"<<endl;
				// go one level down
				evaluate_trees(id, level+1, target, &(source->kids[i]), pts);
			}
		}
		else{
			cout<<target->gid<<" and "<<source->kids[i].gid
				<<" are sufficiently SEPARATED on level "
				<<level+1<<" "<<target->level<<endl;
		}
	}
		
	
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
	// const int np = 2000000;
	const int np = 5;

	// mass
	const double mmin = 1.0;
	const double mrange = 10;
	
	// information necessary for qtree construction
	const int max_level = 2;
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
				xmin+1*xrange/6, ymin,
				min(xrange,yrange)/3);
		}

		// pts[i].gen_coords(xmin, xrange, ymin, yrange);
	}
	
	// compute morton ids
	#pragma omp parallel for default(none), shared(pts)
	for(int i=0; i<np; i++)
		pts[i].get_morton_id(xmin, ymin, x_grid_size, y_grid_size, max_level);
	
	// parallel merge sort
	// vector<particle> tmp;
	// tmp.resize(np);
	// mergesort<particle>(&pts[0], nt, np, &tmp[0], compar_vector);

	// read points data
	// read_points(pts, np);

	n_th=1;
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

	
	evaluate_trees(0, 0, qt1, qt1, pts);
	cout<<"done"<<endl;
	
	ofstream ofs;
	write_boxes( qt1, ofs, 0 );

	delete qt1;

	
	// n_th=2;
	// qtree* qt2 = new qtree;
	// // parallel tree construction
	// start=omp_get_wtime();
	// // #pragma omp parallel for shared(pts, np_local)
	// for(int i=0; i<nt; i++){
	// 	// build tree
	// 	build_tree( xmin, ymin, max_pts_per_node, max_level, width,
	// 	&pts, np, qt2, 0, np, 2);
	// }
	// end=omp_get_wtime();
	// cout<<"wall clock time = "<<end-start<<endl;

	// n_th=4;
	// qtree* qt4 = new qtree;
	// // parallel tree construction
	// start=omp_get_wtime();
	// // #pragma omp parallel for shared(pts, np_local)
	// for(int i=0; i<nt; i++){
	// 	// build tree
	// 	build_tree( xmin, ymin, max_pts_per_node, max_level, width,
	// 	&pts, np, qt4, 0, np, 4);
	// }
	// end=omp_get_wtime();
	// cout<<"wall clock time = "<<end-start<<endl;

	// n_th=8;
	// qtree* qt8 = new qtree;
	// // parallel tree construction
	// start=omp_get_wtime();
	// // #pragma omp parallel for shared(pts, np_local)
	// for(int i=0; i<nt; i++){
	// 	// build tree
	// 	build_tree( xmin, ymin, max_pts_per_node, max_level, width,
	// 	&pts, np, qt8, 0, np, 8);
	// }
	// end=omp_get_wtime();
	// cout<<"wall clock time = "<<end-start<<endl;

	
	// write_points(pts, np);

	// delete qt1, qt2, qt4, qt8;

	
	return 0;
}

