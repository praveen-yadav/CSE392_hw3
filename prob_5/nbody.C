#include <omp.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <bitset>

// #include "../msort/msort.h"
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

// 2-norm
double norm2( const double x1, const double y1,
			  const double x2, const double y2)
{
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}


// kernel
double g_kernel( const double x1, const double y1,
				 const double x2, const double y2) 
{
	if ((x1!=x2) && (y1!=y2))
		return -1/(2*pi) * log( norm2(x1, y1, x2, y2) );
	
	return 0;
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
		ofile<<pts[i].x<<" "<<pts[i].y<<" "<<pts[i].mt_id<<" "
			 <<pts[i].m<<" "<<pts[i].u<<endl;
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
				return &(qt->kids[i]);
				return get_local_node(id, &(qt->kids[i]), level-1);
			}
		}
	}
	
	return NULL;
}

// get which kid has the node with id
qtree* get_local_node(const int id, qtree* qt, const vector<particle>& pts )
{
	if (qt->isleaf)
		return qt;
	
	for(int i=0; i<4; i++){
		if( pts[id].x > qt->kids[i].anchor.x &&
			pts[id].x < (qt->kids[i].anchor.x+qt->kids[i].width.x) &&
			pts[id].y > qt->kids[i].anchor.y &&
			pts[id].y < (qt->kids[i].anchor.y+qt->kids[i].width.y) )

			return &(qt->kids[i]);
			
	}
	
	return NULL;
}


// check if the target and source are well-separated
int well_separated( const double x_target,
					const double y_target,
					qtree* source )
{
	double xdiff = abs((source->anchor.x) - (x_target));
	double ydiff = abs((source->anchor.y) - (y_target));

	if ((xdiff>source->width.x) || (ydiff>source->width.y) )
		return 1;
	
	return 0;
}

// evaluate the potential at target using individual particles
double direct_evaluation( const int id, qtree* source,
						const vector<particle>& pts)
{
	double u=0.0;
	
	for(int i=0; i<source->idx.size(); i++){
		u += g_kernel( pts[id].x, pts[id].y,
					   pts[source->idx[i]].x, pts[source->idx[i]].y )
			* pts[source->idx[i]].m;
	}
		

	return u;
}

// evaluate the potential at a target using a box
double approximate_evaluation( const int id, qtree* source,
							   const vector<particle>& pts )
{
	double u = g_kernel( pts[id].x , pts[id].y,
						 source->centroid.x, source->centroid.y )
		* source->total_mass;

	return u;
}

// evaluate the forces using quadtree
void evaluate_trees( const int id,
					 qtree* source, vector<particle>& pts)
{
	// get qtree pointer of the kid that has the pts[id]
	// target = get_local_node(id, target, 1);
	// target = get_local_node(id, target, pts);

	
	// cout<<"x "<<pts[id].x<<" y "<<pts[id].y<<endl;
	// cout<<"gid "<<target->gid<<" level "<<target->level<<endl;

	// now compare the target with the kids of the source
	for(int i=0; i<4; i++){
		// not separated enough
		if(!well_separated(pts[id].x, pts[id].y, &(source->kids[i]))){
			// but the kid is a leaf
			if(source->kids[i].isleaf){
				// cout<<target->gid<<" and "<<source->kids[i].gid
					// <<" need direct evaluation on level "
					// <<level+1<<endl;	// cout<<"leaf"<<endl;
				pts[id].u += direct_evaluation(id, &(source->kids[i]), pts);
			}
			// not leaf
			else{
				// cout<<"not separated enough"<<endl;
				// go one level down
				evaluate_trees(id, &(source->kids[i]), pts);
			}
		}
		else{
			// cout<<target->gid<<" and "<<source->kids[i].gid
				// <<" are sufficiently SEPARATED on level "
				// <<level+1<<" "<<target->level<<endl;
			pts[id].u += approximate_evaluation(id, &(source->kids[i]), pts );
		}
	}
		
	
	return;
}


// main function
int main( int argc, char** argv )
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
	int np = 2000;
	if(argc>1) np = pow(2,atoi(argv[1]));
	// const int np = pow(2,22);//2000000;

	// const int np = 2000;

	// mass
	const double mmin = 1.0;
	const double mrange = 10;
	
	// information necessary for qtree construction
	const int max_level = 100;//2*log(np)/log(16);
	const int max_pts_per_node = max(1.0,np/(pow(2,max_level)));

	cout<<"max level: "<<max_level<<" "<<endl
		<<"max pts per node: "<<max_pts_per_node<<endl;
	
	// get domain range and grid size
	const double xrange = xmax-xmin;
	const double yrange = ymax-ymin;
	const double x_grid_size = (xmax-xmin)/pow(2.,max_level);
	
	const double y_grid_size = (ymax-ymin)/pow(2.,max_level);
	const point width(xrange, yrange);

	// for measuring wall clock time
	double start, end;
	
	// generate points
	vector<particle> pts;
	pts.resize(np);
	#pragma omp parallel for shared(pts)
	for(int i=0; i<np; i++){
		if(i<np/10){
			pts[i].gen_coords_cluster(xmin, xrange, ymin, yrange,
				mmin, mrange, xmin+xrange/4,
				ymin+yrange/3, min(xrange,yrange)/5);
		}
		else if(i<4*np/10) {
			pts[i].gen_coords_cluster(xmin, xrange, ymin, yrange,
				mmin, mrange, xmin+3*xrange/4,
				ymin+2*yrange/3, min(xrange,yrange)/6);
		}
		else{
			pts[i].gen_coords_cluster(xmin, xrange, ymin, yrange,
				mmin, mrange,
				xmin+xrange/4, ymin+3*yrange/4,
				min(xrange,yrange)/5);
		}

		// pts[i].gen_coords(xmin, xrange, ymin, yrange);
	}
	
	// compute morton ids
	start=omp_get_wtime();	
	#pragma omp parallel for default(none), shared(pts, np)
	for(int i=0; i<np; i++)
		pts[i].get_morton_id(xmin, ymin, x_grid_size, y_grid_size, max_level);
	end=omp_get_wtime();	
	cout<<"get_morton_id: "<<end-start<<endl;
	
	// parallel merge sort
	// vector<particle> tmp;
	// tmp.resize(np);
	// mergesort<particle>(&pts[0], nt, np, &tmp[0], compar_vector);

	// read points data
	// read_points(pts, np);

	cout<<endl;
	for(int i=0; i<4; i++){
		n_th=pow(2.0,i);
		cout<<"proc: "<<n_th<<endl;
			
		qtree* qt1 = new qtree;
		// parallel tree construction
		start=omp_get_wtime();
		for(int i=0; i<nt; i++){
			// build tree
			build_tree( xmin, ymin, max_pts_per_node, max_level, width,
						&pts, np, qt1, 0, np, 1);
		}
		end=omp_get_wtime();
		cout<<"build_tree: "<<end-start<<endl;

		start=omp_get_wtime();	
		#pragma omp parallel for shared(qt1, pts) num_threads(n_th)
		for(int i=0; i<np; i++)
			evaluate_trees(i, qt1, pts);
		end=omp_get_wtime();
		cout<<"evaluate_trees: "<<end-start<<endl;
	
		// ofstream ofs;
		// write_boxes( qt1, ofs, 1 );
		// write_points(pts, np);
	
		// delete qt1;
		
		cout<<endl;
	}
	
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

	
	// delete qt1, qt2, qt4, qt8;

	
	return 0;
}

