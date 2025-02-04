#ifndef TE_H
#define TE_H
# include <iostream>
# include <vector>
# include <algorithm>
# include <sstream>
#include "cut.h"
#include <omp.h>
using namespace std;


//quick sort
void quicksort_index(double* a, int first, int last, int* b);

//read the file
void read_data_input(ifstream& infile_input, int&  last_snapshot, 
int& initial_snapshot);

void read_data_dump(ifstream& infile_input, int& num_atom, double* length, 
double* lengthh, double** pos_o, double** pos_h);

//void intialize the voronoi cell
void init_voronoi(vector<vert>& vert_te, vector<edge>& edge_te, 
vector<face>& face_te, vector<flist>& flist_te, const double* origin, int vlast=4, int elast=6,
int flast=4, int flistlast=12, double ranglim=31.0);

//synthesis
void voronoi_cell
(double* length, double* lengthh, double** pos_o, int index_origin, const int* order, int num_order, 
 double& vol, double& tot_surf);



#endif

