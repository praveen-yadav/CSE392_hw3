#include <omp.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdlib>
// #include <cmath>

using namespace std;

// Algorithm 36.1
// 1. Create a set S of candidates for I as follows. For each vertex v in
// parallel, include v \in S with probability 1/2d(v)
// 2. For each edge in E , if both its endpoints are in S , discard the one of
// lower degree; ties are resolved arbitrary (say by vertex number).
// The resulting set is I.

// random generator function:
int myrandom (int i) { return std::rand()%i;}

// E: sparse matrix to describe graph connectivity
// V: vertices
// I: independent set
// n: number of vertices
double mis_shared( int* val, int* col_ind, int* row_ptr,
				   int* V, vector<int>& I, int n)
{
	int* C = new int[n];
	int* r = new int[n];

	// copy vertices
	for (int i=0;i<n;i++){
		C[i]=V[i];
		r[i] = C[i];
	}
	srand(time(0));
	
	// assign random numbers to vertices
	random_shuffle(r,r+n,myrandom);

	for(int i=0; i<n; i++){
		cout<<r[i]<<endl;
	}

	// int r[6]={5,1,3,4,2,0};
	
	// number of elements searched
	int n_done=0;

	while(n_done!=n){

#pragma omp parallel for shared(I) num_threads(1)
		for(int i=0; i<n; i++){

			if(C[i]==-1)
				continue;
			
			cout<<"working on node "<<i<<endl;
			for(int j=0; j<n; j++)
				cout<<C[j]<<" ";
			cout<<endl;
		

			// get neibhors
			int n_nb = row_ptr[i+1]-row_ptr[i];
			int* nb = new int[n_nb];
			int k=0;
			for(int j=0; j<n_nb; j++){
				int i_nb =  col_ind[j+row_ptr[i]];
				if(C[i_nb]!=-1){
					nb[k] = i_nb;
					k++;
				}
			}
			n_nb=k;
			
			int fl=1;
			for(int j=0; j<n_nb; j++){
				if(r[i] > r[nb[j]]){
					fl = 0;
					break;
				}
			}

			// if the value is smaller than its neibhors
			if(fl){
#pragma omp critical
				{
					I.push_back(i);
				}
				C[i] = -1;
				for( int j=0; j<n_nb; j++)
					C[nb[j]]=-1;

				n_done += 1+n_nb;
			}
		}
	}
}



int main(int argc, char **argv)
{

	int val[18] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	int col_ind[18] = {1,2,3,0,3,0,3,4,5,0,1,2,5,2,5,2,3,4};
	int row_ptr[7] = {0,3,5,9,13,15,18};
	int n=6;
	int V[6] = {0,1,2,3,4,5};
	vector<int> I;
	
	mis_shared( val, col_ind, row_ptr,
				V, I, n);

	cout<<endl<<"result"<<endl;
	for(int i=0; i<I.size(); i++)
		cout<<I[i]<<endl;
	
	return 0;
  
}
