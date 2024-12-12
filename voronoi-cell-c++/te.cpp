# include "te.h"
# include "cut.h"
# include <cmath>
# include <vector>
# include <fstream>
# include <iostream>
using namespace std;
//.............quick sort
void quicksort_index(double* a, int first, int last, int* b) {
    if (first >= last) return; // Base case: no sorting needed
    double x = a[(first + last) / 2];
    int i = first;
    int j = last;
    while (true) {
        while (a[i] < x) i++;
        while (x < a[j]) j--;
        if (i >= j) break;
        swap(a[i], a[j]);
        // Swap corresponding elements in array b
        swap(b[i], b[j]);
        i++;
        j--;
    }
    // Recursive calls
    if (first < i - 1) quicksort_index(a, first, i - 1, b);
    if (j + 1 < last) quicksort_index(a, j + 1, last, b);
}

//.....read input data
void read_data_input(ifstream& infile_input, int&  last_snapshot, int& initial_snapshot)
{ infile_input.ignore(100,'\n');  infile_input>>last_snapshot; infile_input.ignore(100, '\n'); 
infile_input.ignore(100,'\n');  infile_input>>initial_snapshot;
}

//.....read dump data
void read_data_dump(ifstream& infile_dump, int& num_atom, double* length, 
double* lengthh, double** pos_o, double** pos_h){
for(int i1=1; i1<=3;i1++)  infile_dump.ignore(1000, '\n'); 
    infile_dump>>num_atom;
    infile_dump.ignore(100, '\n');
    infile_dump.ignore(100, '\n'); double vertex[3];
    for(int i1=0; i1<3; i1++) { infile_dump>>vertex[0]>>vertex[1]>>vertex[2];
    length[i1]=vertex[1]-vertex[0]; lengthh[i1]=length[i1]/(2.0);}
    infile_dump.ignore(100, '\n');
    infile_dump.ignore(100, '\n');
    int m_o=0, m_h=0, index_atom, index_type; double pos[3], vec[3];
    for(int i1=1; i1<=num_atom; i1++){
        infile_dump>>index_atom>>index_type>>
        pos[0]>>pos[1]>>pos[2]>>vec[0]>>vec[1]>>vec[2];
        switch (index_type){
            case 1: copy(pos, pos+3, pos_o[m_o]);  
                    //cout<<"o "<<" "<<pos_o[m_o][0]<<" "<<pos_o[m_o][1]<<" "<<pos_o[m_o][2]<<endl;
                    m_o+=1; break;
            case 2: copy(pos, pos+3, pos_h[m_h]); 
                    //cout<<"h "<<pos_h[m_h][0]<<" "<<pos_h[m_h][1]<<" "<<pos_h[m_h][2]<<endl; 
                    m_h+=1; break; 
        } } 
infile_dump.ignore(100, '\n');
};



//void intialize the voronoi cell
void init_voronoi(vector<vert>& vert_te, vector<edge>& edge_te, 
vector<face>& face_te, vector<flist>& flist_te, double* origin, int vlast, int elast,
int flast, int flistlast, double rangelim){
//default setting
int list_v_e[100]={1, 3, 6, 1, 2, 5, 2, 3, 4, 4, 5, 6};
int list_e_v[100]={1,2, 2, 3, 1, 3, 3, 4, 2, 4, 1,4};
int list_e_f[100]={1, 4, 1, 2, 1, 3, 2, 3, 2, 4, 3, 4};
int list_flist_e[100]={1, 2, 3, 2, 4, 5, 3, 4, 6, 6, 5, 1};
int list_flist_v[100]={1, 2, 3, 2, 3, 4, 1, 3, 4, 1, 4, 2 };
int list_f_vfar[100]={1, 2, 1, 1};
double temp_list[12]={
(sqrt(6.0)*rangelim)/(3.0), (-sqrt(2.0)*rangelim/(3.0)), -(rangelim/(3.0)),
-(sqrt(6.0)*rangelim)/(3.0), (-sqrt(2.0)*rangelim/(3.0)), -(rangelim/(3.0)), 
0.0*rangelim, (2*sqrt(2.0)*rangelim)/(3.0), -(rangelim/(3.0)),
0.0*rangelim, 0.0*rangelim, rangelim};
//vertex1-4
for(int i=0; i<vlast; i++){
    vert_te[i].posv[0]=temp_list[i*3+0]+origin[0]; 
    vert_te[i].posv[1]=temp_list[i*3+1]+origin[1]; 
    vert_te[i].posv[2]=temp_list[i*3+2]+origin[2];
    vert_te[i].posv[3]=rangelim; vert_te[i].stat=2;
    for(int i1=0; i1<3;i1++) vert_te[i].e[i1]=list_v_e[3*i+i1]-1;
    //cout<<"vertex: "<<i<<" "<<vert_te[i].posv[0]<<" "<<vert_te[i].posv[1]
    //<<" "<<vert_te[i].posv[2]<<endl;
}

//edge1-6
for(int i=0; i<elast; i++){
for(int i1=0; i1<2; i1++) {edge_te[i].v[i1]=list_e_v[i*2+i1]-1; 
edge_te[i].f[i1]=list_e_f[i*2+i1]-1;}
edge_te[i].stat=3;
//cout<<"edge: "<<i<<" "<<edge_te[i].v[0]<<" "<<edge_te[i].v[1]<<endl;
}
//f_te, flist_te, vert_te, flast, vol, and flist
for(int i=0;i<flistlast; i++) flist_te[i].link=i+1;
int m0=0;
for(int i=0; i<flast;i++){
    face_te[i].stat=3; face_te[i].vfar=list_f_vfar[i]-1; 
    face_te[i].fptr=m0;  face_te[i].dist=rangelim/(3.0);
    for(int i1=0; i1<3; i1++) { flist_te[m0].e=list_flist_e[m0]-1; 
    flist_te[m0].v=list_flist_v[m0]-1; m0=m0+1;
}
flist_te[m0-1].link=face_te[i].fptr;
//for(int j=0; j<3; j++) cout<<"link "<<flist_te[m0-3].link<<flist_te[m0-2].link
//<<flist_te[m0-1].link<<endl;
}
}




void voronoi_cell(double* length, double* lengthh, double** pos_o, 
int index_origin, int* order, int num_order, double& vol, double& tot_surf){
double rangelim=31.0, second[3];
int elast=6, vlast=4,  flast=4, flistlast=12;
double origin[3]; for(int i=0; i<3;i++) origin[i]=pos_o[index_origin][i];
vector<vert> vert_te(4); vector<edge> edge_te(6); vector<face> face_te(4); 
vector<flist> flist_te(12);
init_voronoi( vert_te, edge_te, face_te, flist_te, origin);
volume(face_te, flist_te, vert_te, flast, vol, tot_surf, length, lengthh);
//cout<<"initializeation "<<vol<<" "<<tot_surf<<endl;
bool near;
//for(int index_order=0; index_order<num_order; index_order++)
for(int index_order=0; index_order<num_order; index_order++){
    for(int i1=0; i1<3; i1++){second[i1]=pos_o[order[index_order]][i1]; }
    near=false; //double temp0=0.0;
    //calc_distan(length, lengthh, second, temp0);
    bisectplane( origin, second, length, lengthh, vert_te, edge_te, face_te , 
    flist_te, elast, vlast, flast, flistlast, near);
    int l1=face_te[flast-1].fptr; //cout<<"new added face"<<endl;
    /*
    while(flist_te[l1].link != face_te[flast-1].fptr)
    {cout<<"e "<<flist_te[l1].e<<" v "<<flist_te[l1].v<<" l "<<flist_te[l1].link<<endl; 
    l1=flist_te[l1].link;}
    cout<<"e "<<flist_te[l1].e<<" v "<<flist_te[l1].v<<" l "<<flist_te[l1].link<<endl;
    */
    
    /*
    cout<<" vol "<<vol<<tot_surf<<endl;
    for(int i1=0; i1<flast; i1++ ){if (face_te[i1].stat ==3) cout<<" "<<face_te[i1].dist;}
    cout<<endl;
    */
                              
                                                             }
 volume(face_te, flist_te, vert_te, flast, vol, tot_surf, length, lengthh);                                                       

 }
