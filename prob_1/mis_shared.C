#include <omp.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <stack>
#include <cstdlib>
#include <cmath>

using namespace std;

// Algorithm 36.1
// 1. Create a set S of candidates for I as follows. For each vertex v in
// parallel, include v \in S with probability 1/2d(v)
// 2. For each edge in E , if both its endpoints are in S , discard the one of
// lower degree; ties are resolved arbitrary (say by vertex number).
// The resulting set is I.

// comparison function for sorting pairs
bool comp_pairs( const pair<int, int>& i, const pair<int, int>& j ) {
    if( i.first < j.first ) return true;
	else if( i.first == j.first) return i.second<j.second;
	else return false;
}

// create sparse symmetric matrix with 0 diagonal
int symm_matrix(const int n, const int ne, 	vector<pair <int, int> >& idx,
				vector<int>& col_ind, vector<int>& row_ptr)
{
	srand(time(0));

	int inc_rdm=n*n/ne-1;
	int n_curr = 0;
	int rn=0;

	cout<<"inc_rdm "<<inc_rdm<<endl;
	
  	while(rn < n){

		int inc = rand()%inc_rdm+1;
		n_curr += inc;
 
		if (n_curr>=n){
			rn++;
			n_curr=n_curr-n;
		}
		int cn = n_curr;
		//  n_curr=cn;

		if(rn>=n)
			break;

		//	  cout<<rn<<" "<<cn<<" "<<n_curr<<endl;
	  
		// only upper-triangular region
		int flag=1;
		if(rn<=cn){
			flag=0;
		}

		if(flag){
			pair<int, int> el_1(rn, cn);
			pair<int, int> el_2(cn, rn);
			idx.push_back(el_1);
			idx.push_back(el_2);
		}
	}


	cout<<"sorting..."<<endl;
	sort(idx.begin(), idx.end(), comp_pairs);
	cout<<"sorting done"<<endl;

	// for(int i=0; i<idx.size(); i++)
	// 	cout<<idx[i].first<<" "<<idx[i].second<<endl;

	col_ind.resize(idx.size(), 0);
	row_ptr.resize(idx[0].first+1,0);

	for(int i=0; i<idx.size()-1; i++){
		col_ind[i] = idx[i].second;
		if(idx[i].first<idx[i+1].first){
			int n_gap = idx[i+1].first-idx[i].first;
			for(int j=0; j<n_gap; j++)
				row_ptr.push_back(i+1);
		}
	}
	col_ind[idx.size()-1] = idx[idx.size()-1].second;

	// for the case that there are empty rows at the end of matrix
	for(int i=row_ptr.size(); i<n+1; i++)
		row_ptr.push_back(idx.size());

	// cout<<"col_ind"<<endl;
	// for(int i=0; i<col_ind.size(); i++)
	// 	cout<<col_ind[i]<<" ";
	// cout<<endl;

	// cout<<"row_ptr: "<<row_ptr.size()<<endl;
	// for(int i=0; i<row_ptr.size(); i++)
	// 	cout<<row_ptr[i]<<" ";
	// cout<<endl;

	return 0;
}


// write out the sparse matrix
int write_matrix(const int n,
				 const vector<int>& col_ind,
				 const vector<int>& row_ptr)
{
	ofstream file_out;
	file_out.open ("matrix.dat");

	if(!file_out.is_open()){
		return 1;
	}

	// matrix dimension
	file_out<<row_ptr.size()-1<<endl;
	// number of non-zero elements
	file_out<<col_ind.size()<<endl;
	
	for(int i=0; i<col_ind.size(); i++)
		file_out << col_ind[i] <<" ";
	file_out<<endl;
	
	for(int i=0; i<row_ptr.size(); i++)
		file_out << row_ptr[i] <<" ";
	file_out<<endl;
	
	file_out.close();
	
	return 0;
}

// read from the file
int read_matrix( vector<int>& col_ind,
				 vector<int>& row_ptr)
{
	string line;
	ifstream file_in ("matrix.dat");
	if (!file_in.is_open()) return 1;

	int ne;
	int n;

	// matrix dimension
	file_in>>n;
	// number of non-zero elements
	file_in>>ne;

	col_ind.resize(ne,0);
	row_ptr.resize(n+1,0);
	
	for(int i=0; i<ne; i++){
		file_in>>col_ind[i];
	}
	for(int i=0; i<n+1; i++){
		file_in>>row_ptr[i];
	}
	
	file_in.close();
	
	return 0;
}

// E: sparse matrix to describe graph connectivity
// V: vertices
// I: independent set
// n: number of vertices
int mis_shared_2( const vector<int>& col_ind,
				  const vector<int>& row_ptr,
				  const int* V, int* I, const int n,
				  const int nt )
{
	// copy V into C
	int* C = new int[n];
	//#pragma omp parallel for shared(C, V) num_threads(nt)
	for(int i=0; i<n; i++)
		C[i] = V[i];

	// generate random integers
	int* r = new int[n];
	srand(time(0));
	//#pragma omp parallel for shared(r) num_threads(nt)
	for (int i=0;i<n;i++){
		r[i] = rand()%(int(pow(n,4)));
		// r[i]=C[i];
	}


	
	// number of node removed
	int n_done = 0;
	// size of independent set
	int n_is = 0;

	int flag=0;
#pragma omp parallel num_threads(nt) shared(n_done, n_is, C, r, I)
	{
		int u=n/(nt)*omp_get_thread_num();
				
		while(n_done<n){

			
		// #pragma omp barrier
			if(omp_get_thread_num()==0){
			n_done = 0;
			// #pragma omp parallel for reduction(+:n_done)
			// flag=1;
			for(int j=0; j<n; j++){
				if(C[j]!=-1)
					flag=0;
				else
					n_done++;
			}
				// cout<<C[j]<<" ";
				// cout<<endl;
				// n_done += C[j];

			// cout<<n<<" "<<n_done<<endl;

			if(flag==0)
				cout<<n_done<<endl;
			
			if(n_done==n)
			flag=1;
		}

			if(flag==1){
				cout<<"break "<<omp_get_thread_num()<<endl;
				
				break;
			}

			
			// get random vertex
			// int u = (rand())%n;
			u = u%n;
			if (C[u]==-1){
				//#pragma omp atomic
				u++;
				continue;
			}

			// if (omp_get_thread_num()==1)
			// cout<<u<<endl;
			
			// get neibhors of u
			int n_nb = row_ptr[u+1]-row_ptr[u];
			int* nb = new int[n_nb];
			int k=0;
			//#pragma omp critical(about_C)
			{
				for(int j=0; j<n_nb; j++){
					int i_nb =  col_ind[j+row_ptr[u]];
					if(C[i_nb]!=-1){
						nb[k] = i_nb;
						k++;
					}
				}
			}
			// number of neighbors
			n_nb=k;
			
			// check if the node value is the locally smallest
			int fl=1;
			for(int j=0; j<n_nb; j++){
				if(r[u] > r[nb[j]]){
					fl = 0;
					break;
				}
				else if(r[u] == r[nb[j]] && u>nb[j]){
					fl = 0;
					break;
				}
			}

			// if the value is smaller than its neibhors
			// or the values are the same but node number is smallest, 
			if(fl){
#pragma omp critical(add_data)
				{
					// cout<<"u "<<u<<endl;
					I[n_is]=u;
				}

#pragma omp atomic
				n_is++;
// #pragma omp atomic
				// n_done += (1+n_nb);
									
				//#pragma omp critical(about_C)
				{
					C[u] = -1;
					for( int j=0; j<n_nb; j++){
						C[nb[j]]=-1;
					}
				}

			}

			
			delete[] nb;
			
			//#pragma omp atomic
			u++;

// #pragma omp barrier

			
		} // end while
		

			
	} // end parallel region

			
	delete[] C, r;
			
	return n_is;
	
}

// main
int main(int argc, char **argv)
{
	// number of vertices
	const int n =  1000;
	// number of edges
	const int ne = 4000;
	// number of threads
	// const int nt=4;
	
	vector<pair <int, int> > idx;
	vector<int> col_ind;
	vector<int> row_ptr;
	
	cout<<"creating symmetric sparse matrix"<<endl;
	symm_matrix(n, ne, idx, col_ind, row_ptr);
	cout<<"done!"<<endl;
	
	cout<<"writing matrix to file."<<endl;
	if(write_matrix(n, col_ind, row_ptr)){
		cout<<"file output failed."<<endl;
		return 1;
	}
	cout<<"done!"<<endl;

	cout<<"reading matrix from file."<<endl;
	if(read_matrix(col_ind, row_ptr)){
		cout<<"file input failed."<<endl;
		return 1;
	}
	cout<<"done!"<<endl;	
	
	// initialize vertices
	int V[n]; 
	for(int i=0; i<n; i++)
		V[i] = i;

	// initialize independent set 
	int I[n];

	// number of indepdnent nodes
	int n_is;
	
	cout<<"n_omp=1"<<endl;
	for(int p=0; p<1; p++){
	
		const double start=omp_get_wtime();
		n_is = mis_shared_2( col_ind, row_ptr,
								 V, I, n, 4);
		const double end=omp_get_wtime();
	
		cout<<"wall clock time = " <<end-start<<endl;
	}
	
	cout<<"n_omp=2"<<endl;
	for(int p=0; p<5; p++){
	
		const double start=omp_get_wtime();
		n_is = mis_shared_2( col_ind, row_ptr,
								 V, I, n, 2);
		const double end=omp_get_wtime();
	
		cout<<"wall clock time = " <<end-start<<endl;
	}

	cout<<"n_omp=4"<<endl;
	for(int p=0; p<5; p++){
	
		const double start=omp_get_wtime();
		n_is = mis_shared_2( col_ind, row_ptr,
								 V, I, n, 4);
		const double end=omp_get_wtime();
	
		cout<<"wall clock time = " <<end-start<<endl;
	}

	cout<<"n_omp=8"<<endl;
	for(int p=0; p<5; p++){
	
		const double start=omp_get_wtime();
		n_is = mis_shared_2( col_ind, row_ptr,
								 V, I, n, 8);
		const double end=omp_get_wtime();
	
		cout<<"wall clock time = " <<end-start<<endl;
	}

	cout<<"n_omp=16"<<endl;
	for(int p=0; p<5; p++){
	
		const double start=omp_get_wtime();
		n_is = mis_shared_2( col_ind, row_ptr,
								 V, I, n, 16);
		const double end=omp_get_wtime();
	
		cout<<"wall clock time = " <<end-start<<endl;
	}

	cout<<"n_omp=32"<<endl;
	for(int p=0; p<5; p++){
		//cout<<"generating maximum independent set."<<endl;
		const double start=omp_get_wtime();
		n_is = mis_shared_2( col_ind, row_ptr,
								 V, I, n, 32);
		const double end=omp_get_wtime();
		//	cout<<"done!"<<endl;	   
	
		cout<<"wall clock time = " <<end-start<<endl;

		//	cout<<endl<<"result "<<n_is<<endl;
		// for(int i=0; i<n_is; i++)
		// 	cout<<I[i]<<endl;
	
	}

	// cout<<endl<<"result "<<n_is<<endl;
	// for(int i=0; i<n_is; i++)
	// 	cout<<I[i]<<endl;

	return 0;
  
}
