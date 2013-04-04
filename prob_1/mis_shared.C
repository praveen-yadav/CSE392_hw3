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
int mis_shared( const int* val, const int* col_ind, const int* row_ptr,
				const int* V, int* I, const int n)
{
	int* C = new int[n];
	int* r = new int[n];

	srand(time(0));

	// copy vertices
	for (int i=0;i<n;i++){
		C[i]=V[i];
		r[i] = rand()%(2*n);//C[i];
	}
	
	// assign random numbers to vertices
	// random_shuffle(r,r+n,myrandom);

	for(int i=0; i<n; i++){
		cout<<r[i]<<endl;
	}

	// int r[6]={5,1,3,4,2,0};
	
	// number of elements searched
	int n_done=0;

	// number of independent node
	int n_is=0;
	
	while(n_done<n){
		vector<int> nb;
		
#pragma omp parallel for private(nb) shared(I,C,r,n_done,row_ptr,col_ind,n_is) num_threads(1)
		for(int i=0; i<n; i++){
			

			// if( omp_get_thread_num()==0)
			cout<<"working on node "<<i<<" n_done "<<n_done<<endl;		
			if(C[i]==-1)
				continue;
			
			// get neibhors
			int n_nb = row_ptr[i+1]-row_ptr[i];
			// int* nb = new int[n_nb];
			int k=0;
			for(int j=0; j<n_nb; j++){
				int i_nb =  col_ind[j+row_ptr[i]];
				if(C[i_nb]!=-1){
					nb.push_back(i_nb);
					k++;
				}
			}
			// number of neighbors
			n_nb=k;
			int nb_size = nb.size();
			
			int fl=1;
			for(int j=0; j<n_nb; j++){
				if(r[i] > r[nb[nb_size-k+j+1]]){
					fl = 0;
					break;
				}
				else if(r[i] == r[nb[nb_size-k+j+1]] && i>nb[nb_size-k+j+1]){
					fl = 0;
					break;
				}
			}
			
			// if the value is smaller than its neibhors
			// or the values are the same but node number is smaller, 
			if(fl){
				// #pragma omp critical
				
				I[n_is]=i;
				n_is++;
				
				// add the node number too
				nb.push_back(i);
				
				//	#pragma omp atomic
				n_done += 1+n_nb;
			}
			
			
			for(int j=0; j<n; j++)
				cout<<C[j]<<" ";
			cout<<endl;
			
		}

		// sweep all the "used" nodes
		for( int j=0; j<nb.size(); j++){
			cout<<nb[j]<<endl;
			C[nb[j]]=-1;
		}
		
		
	}
	

return n_is;

}



int main(int argc, char **argv)
{

	int val[18] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	int col_ind[18] = {1,2,3,0,3,0,3,4,5,0,1,2,5,2,5,2,3,4};
	int row_ptr[7] = {0,3,5,9,13,15,18};
	int n=6;
	int V[6] = {0,1,2,3,4,5};
	int I[3];
	
	int n_is = mis_shared( val, col_ind, row_ptr,
						   V, I, n);

	cout<<endl<<"result"<<endl;
	for(int i=0; i<n_is; i++)
		cout<<I[i]<<endl;
	
	return 0;
  
}
