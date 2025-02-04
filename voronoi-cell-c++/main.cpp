#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#include "te.h"
#include "cut.h"
#include <omp.h>
using namespace std;
//#define OMP_NUM_THREADS 4

int main(){
   //allocate the space for the position of oxygen and hydrogen atoms
 double ** pos_o=new double *[1024];  for(int i0=0; i0<1024; i0++) pos_o[i0]=new double[3];
 double ** pos_h=new double *[2048]; for(int i0=0; i0<2048; i0++) pos_h[i0]=new double[3];
 int num_atom; double length[3], lengthh[3];
 double ** pos_near=new double* [1024]; for(int i=0; i<1024; i++) pos_near[i]=new double[3]; 
 int * index_near=new int [1024]; double* sorted_distance=new double [1024];
 double temp_r, vol, tot_surf; 
 int last_snapshot, initial_snapshot, m_1; 

//omp_set_num_threads(OMP_NUM_THREADS);
//omp_set_dynamic(0);
stringstream ss1;
for( int i=6; i<=6; i++){
  
    //read input file
    ss1<<"input"<<i; string s1=ss1.str(); 
    ifstream infile_input(s1.c_str(), ios_base::in);
    read_data_input( infile_input,  last_snapshot, initial_snapshot);
    cout<<last_snapshot<<" "<<initial_snapshot<<endl;
    infile_input.close(); ss1.str(""); ss1.clear();
    
   
     
    //read water-dump
    ss1<<"../../save"<<i<<"/water.dump"; s1=ss1.str(); cout<<s1<<endl;
    ifstream infile_dump(s1.c_str(),ios_base::in);
    int start_snapshot;
    //initialzation of start snapshot
    if(initial_snapshot < 100000) {start_snapshot=100000;
    for(int i=initial_snapshot; i <start_snapshot; i++){
    for(int j=1; j<=3081;j++)  infile_dump.ignore(1000, '\n');
                                                       }}
    else{start_snapshot=initial_snapshot;}
    //read data
    for(int i0=start_snapshot; i0<=start_snapshot; i0++){ 
    read_data_dump(infile_dump, num_atom, length, lengthh, pos_o, pos_h); 

  
#pragma omp parallel for default(none) private( m_1,vol, tot_surf) firstprivate(length, lengthh) shared(pos_o)  
 
 for(int i1=0; i1<1024; i1++){  m_1=0;
double ** pos_near=new double* [1024]; for(int i=0; i<1024; i++) pos_near[i]=new double[3];
 int * index_near=new int [1024]; double* sorted_distance=new double [1024];
//calculate the distance of neighbors
    for(int j1=0; j1<1023; j1++){
        if(j1 !=i1) { 
            for(int k1=0;k1<3;k1++) pos_near[m_1][k1]=pos_o[j1][k1]-pos_o[i1][k1]; }
            calc_distan(length, lengthh, pos_near[m_1], sorted_distance[m_1]);
            index_near[m_1]=j1; m_1+=1;
                                }
   quicksort_index(sorted_distance , 0, 1022, index_near);
  
//int p1=omp_get_num_threads(), p2=omp_get_thread_num();
    voronoi_cell( length, lengthh, pos_o, i1, index_near, 1023,  vol,  tot_surf);
    printf("%d  %lf  %lf  \n", i1, vol, tot_surf); 
    for(int i=0; i<1024;i++) delete[] pos_near[i]; delete[] pos_near;
    delete[] index_near; delete[] sorted_distance;   
    //cout<<i1<<" "<<vol<<" "<<tot_surf<<" "<<omp_get_thread_num()<<" "<<omp_get_num_threads()<<endl
    } 

    }





    
infile_dump.close();
}


    for(int i=0; i<1024;i++) { delete[] pos_o[i]; } delete[] pos_o;
    for(int i=0; i<2048;i++) { delete[] pos_h[i]; } delete[] pos_h;
   //for(int i=0; i<1024;i++) delete[] pos_near[i]; delete[] pos_near; 
    //delete[] index_near; delete[] sorted_distance;
    return 0;
}
