#include "cut.h"
#include <iostream>
#include <cmath>
#include <omp.h>
//the function of struct 
//vert
using namespace std;

void calc_distan(const double * length, const double * lengthh, double * distan, double& tp){
    tp=0.0;
    for(int i=0; i<3;i++){
        if(std::abs(distan[i])> lengthh[i]) {
            if(distan[i]>0.0) distan[i]=distan[i]-length[i];
            else distan[i]=distan[i]+length[i]; 
            } 
            tp+=distan[i]*distan[i];
            }    
}

//cross and dot product
void crossproduct(const double* a, const double* b, double* c){
c[0]=a[1]*b[2]-a[2]*b[1];
c[1]=a[2]*b[0]-a[0]*b[2];
c[2]=a[0]*b[1]-a[1]*b[0];
}

void dotproduct(const double* a, const double* b, double& c, int length){
c=0.0; for(int i=0; i<length; i++) c+=a[i]*b[i];
}

//volume and surface
void volume(const vector<face>& face_te, const vector<flist>& flist_te, const vector<vert>& vert_te, 
int flast,  double& vol, double& tot_surf, const double* length, const double* lengthh){
vol=0.0; tot_surf=0.0; double surf;
for(int i0=0; i0<flast; i0++){
if (face_te[i0].stat == 3) { surface(face_te[i0], flist_te, vert_te, surf, length, lengthh); tot_surf+=surf;
vol+=(surf*face_te[i0].dist)/3.0;}
                             }
}

void surface(const face& face0, const vector<flist>& flist_te, const vector<vert>& vert_te, double& surf, 
const double* length, const double* lengthh){
surf=0.0; int link_index=flist_te[face0.fptr].link, m=0;
double pos0[4], pos1[4], tempr[4];
while( flist_te[link_index].link != face0.fptr) {
    for(int i=0; i<3; i++) 
    pos0[i]=vert_te[flist_te[link_index].v].posv[i]-vert_te[flist_te[face0.fptr].v].posv[i];
    link_index=flist_te[link_index].link;
    for(int i=0; i<3; i++) 
    pos1[i]=vert_te[flist_te[link_index].v].posv[i]-vert_te[flist_te[face0.fptr].v].posv[i];
    crossproduct(pos0, pos1, tempr);
    dotproduct(tempr, tempr, tempr[4]);
    surf+=(sqrt(tempr[4]))/(2.0);
}
}


//synthesis the process
void bisectplane( const double* pos0, const double* pos1, const double* length, const double* lengthh, 
vector<vert>& vert_te, vector<edge>& edge_te, vector<face>& face_te , 
vector<flist>& flist_te, int& elast, int& vlast, int& flast,int& flistlast, bool& near)
{ int vlast0=vlast, changed_face=0, m; near=false;
 judge_stat(pos0, pos1, length, lengthh, vert_te, edge_te, face_te, 
 flist_te, vlast, flast);
 //==========finding the new vert 
//!=====whether change the face of current Vc after the consideration of extra neighbor
find_new_vert(pos0, pos1, length, lengthh, vert_te, edge_te, elast, vlast, changed_face);
if(changed_face > 0 ){near=true; int elast0=elast, first_vertex; 
int first_cut_edge, second_cut_flist, second_vertex, second_cut_edge;
for(int i0=0; i0<flast; i0++){if (face_te[i0].stat ==2) {
//cout<<"face changed "<<i0<<endl;
cutting_face_accessory(face_te[i0], i0, edge_te, vert_te, flist_te, flistlast, elast, flast, 
first_vertex, first_cut_edge, second_cut_flist, second_vertex, second_cut_edge);} 
                              }
//!======the new flist to new surface flist_te new added  new face : stat fptr 
face f1; face_te.push_back(f1);
new_face_flist(flist_te, vert_te, edge_te, face_te[flast], vlast0, 
vlast, elast0, flistlast);
double temp0[4]; for(int i0=0; i0<3; i0++) temp0[i0]=pos1[i0]-pos0[i0]; 
calc_distan(length, lengthh, temp0, temp0[3]);
//new face:: dist
face_te[flast].dist=sqrt(temp0[3])/(2.0);
flistlast=flistlast+vlast-vlast0; flast=flast+1;
//change the stat of cutting face and cutting edge
for(int i0=0; i0<elast; i0++){if (edge_te[i0].stat == 2) edge_te[i0].stat=3; }
for(int i0=0; i0<flast; i0++){if (face_te[i0].stat == 2) face_te[i0].stat=3; }
                      }
}

//judge_vertex and status
void judge_stat(const double* pos0, const double* pos1, const double* length, const double* lengthh, 
vector<vert>& vert_te, vector<edge>& edge_te, vector<face>& face_te, 
vector<flist>& flist_te, int vlast, int flast){
bool vert_outside; int fs=0, l1;
//change the status of edge and vertex
for(int i0=0; i0<vlast; i0++){if(vert_te[i0].stat >1) 
{vert_outside=false;
judge_vertex(pos0, pos1, vert_te[i0].posv, length, lengthh, vert_outside );
if(vert_outside == true) {
    vert_te[i0].stat=1; //cout<<"vertex outside "<<i0<<endl;
    for(int i1=0; i1<3;i1++) edge_te[vert_te[i0].e[i1]].stat-=1;
                          }}}

//change the staus of face
for(int i0=0; i0<flast; i0++) { fs=0;
if(face_te[i0].stat > 1) { l1=face_te[i0].fptr;
if(flist_te[l1].link != face_te[i0].fptr)
{while(true){if(vert_te[flist_te[l1].v].stat <2) face_te[i0].stat=2; else fs+=1;
l1=flist_te[l1].link;  if(l1 ==face_te[i0].fptr) break; 
            }
}
if(fs<1) face_te[i0].stat=1;
                          }
   //cout<<"face "<<i0<<" "<<face_te[i0].stat<<endl;     
                              }
}


void judge_vertex( const double* pos0, const double* pos1, const double* pos2, 
const double* length, const double* lengthh, bool& vert_outside){
double distan[3], tempr, plu[3], temp0, temp1;
vert_outside=false; for(int i=0; i<3; i++)
{distan[i]=pos1[i]-pos0[i]; }
calc_distan(length, lengthh, distan, tempr);
for(int i=0; i<3;i++){plu[i]=distan[i]+pos0[i];}
dotproduct(pos0, pos0, temp0); dotproduct(plu, plu, temp1);
dotproduct(distan, pos2, tempr);
if ( tempr > ((temp1-temp0)/(2.0)) ) vert_outside=true;
//cout<<"sucessfully judge_vertex"<<endl;
}


//determing the new vertex pos and status
void find_new_vert(const double* pos0, const double* pos1, const double* length, const double* lengthh, 
vector<vert>& vert_te, vector<edge>& edge_te, int elast, int& vlast, int& changed_face){
double new_vertex[3], a, pos[4]; int i2, change; vert vert1; 
for(int i0=0; i0<elast; i0++){
    if (edge_te[i0].stat ==2){ changed_face=1; change=0;
    for(int i1=0; i1<2; i1++){ if(vert_te[edge_te[i0].v[i1]].stat==2 && change==0){
        if (i1==0) i2=1; else i2=0;
        new_vertex_coeff(pos0, pos1, vert_te[edge_te[i0].v[i1]].posv, 
        vert_te[edge_te[i0].v[i2]].posv, length, lengthh, a);
        //nerw_vertex
        for(int j1=0; j1<3; j1++) new_vertex[j1]=vert_te[edge_te[i0].v[i1]].posv[j1]+
        a*(vert_te[edge_te[i0].v[i2]].posv[j1]-vert_te[edge_te[i0].v[i1]].posv[j1]);
        vert1.stat=2; for(int j1=0; j1<3; j1++) vert1.posv[j1]=new_vertex[j1];
        vert1.e[0]=i0; //vert1.e[1]=vert1.e[2]=0; 
        //revised 
        vert1.e[1]=40000; vert1.e[2]=40000;
        //
        for(int j1=0; j1<3;j1++) pos[j1]=vert1.posv[j1]-pos0[j1];
        calc_distan(length, lengthh, pos, pos[3]); vert1.posv[3]=sqrt(pos[3]);
        change=change+1; vlast=vlast+1; 
        vert_te.push_back(vert1); edge_te[i0].v[i2]=vlast-1;
        //cout<<"edge "<<i0<<" changed "<<edge_te[i0].v[i2]<<" "
        //<<"remaining "<<edge_te[i0].v[i1]<<endl;
                                                                                }
                             }
                             }}
                                                                                }



void new_vertex_coeff(const double* pos0, const double* pos1, double* pos2, double* pos3, 
const double* length, const double* lengthh, double& a){
double distan[3], tempr, shift1[3], shift2[3];
for(int i=0; i<3;i++) { distan[i]=pos1[i]-pos0[i]; shift1[i]=pos3[i]-pos2[i];} 
calc_distan(length, lengthh, distan, tempr);
for(int i=0; i<3; i++) shift2[i]=pos0[i]+distan[i];
double temp1, temp0, b;
dotproduct(shift2, shift2, temp1);
dotproduct(pos0, pos0, temp0);
b=(temp1-temp0)/(2.0);
double temp2, temp3;
dotproduct(distan, shift1, temp2); dotproduct(distan, pos2, temp3);
a=(b-temp3)/(temp2+0.0);
//cout<<"new_vert coeff "<<a<<endl;
//cout<<"sucessfully new_vertex_coeff"<<endl;
}


// determine the new edge cut face and new face
// the revised face of VC due to the exisitance of new face
 void cutting_face_accessory( face& face_te0, int index_face, vector<edge>& edge_te, 
 vector<vert>& vert_te, vector<flist>& flist_te, int& flistlast, int& elast, int flast, 
 int& first_vertex, int& first_cut_edge, int& second_cut_flist, int& second_vertex, 
 int& second_cut_edge)
 { int firstlink=face_te0.fptr, fv=flist_te[face_te0.fptr].v;
if( vert_te[fv].stat==1) {while(vert_te[fv].stat ==1) {firstlink=flist_te[firstlink].link; 
fv=flist_te[firstlink].v;                             }}
face_te0.fptr=firstlink;
//finding the breaking point first_vertex
int link_index=firstlink, fe=flist_te[link_index].e; 
while(edge_te[fe].stat == 3){ link_index=flist_te[link_index].link; 
fe=flist_te[link_index].e;}

for(int i1=0; i1<2; i1++){if(edge_te[fe].v[i1] != flist_te[link_index].v) 
{first_vertex=edge_te[fe].v[i1]; first_cut_edge=fe; }      }
//!===== second_vertex
link_index=flist_te[link_index].link;
if(edge_te[flist_te[link_index].e].stat ==1){
    while(edge_te[flist_te[link_index].e].stat ==1){
        link_index=flist_te[link_index].link;
                                                   }}
second_cut_flist=flist_te[link_index].link;
//second vertex and second cuttin edge
for(int i1=0; i1<2; i1++) {
    if(edge_te[flist_te[link_index].e].v[i1] != flist_te[flist_te[link_index].link].v)
    { second_vertex=edge_te[flist_te[link_index].e].v[i1]; 
     second_cut_edge=flist_te[link_index].e;
    }}                                      
//adding new edge
edge edge1; edge1.stat=3; edge1.v[0]=first_vertex; edge1.v[1]=second_vertex; 
edge1.f[0]=index_face; edge1.f[1]=flast; edge_te.push_back(edge1); 
//!====add new edge of new vertex
for(int i1=0; i1<3;i1++){ 
if(vert_te[first_vertex].e[i1] ==40000) {vert_te[first_vertex].e[i1]=elast; break;}}
for(int i1=0; i1<3;i1++){ 
if(vert_te[second_vertex].e[i1] ==40000) {vert_te[second_vertex].e[i1]=elast; break;}}
//!==== create two new element of flist for face i0
flist flist1;  flist1.v=first_vertex; flist1.e=elast; flist1.link=flistlast+1; 
flist_te.push_back(flist1); 
flist flist2; flist2.v=second_vertex; flist2.e=second_cut_edge, 
flist2.link=second_cut_flist; flist_te.push_back(flist2);
//!========create the new loop by define the flist of first cut edge
link_index=face_te0.fptr;
while(edge_te[flist_te[link_index].e].stat  ==3 ) link_index=flist_te[link_index].link;
flist_te[link_index].link=flistlast;

link_index=face_te0.fptr;
/*
while( flist_te[link_index].link != face_te0.fptr ){cout<<"flist of face "
<<flist_te[link_index].v<<" "<<flist_te[link_index].e<<" "<<flist_te[link_index].link
<<endl; link_index=flist_te[link_index].link;
}

cout<<"flist of face "<<flist_te[link_index].v<<" "<<flist_te[link_index].e<<" "
<<flist_te[link_index].link<<endl;
*/
flistlast+=2; elast+=1;
}

// new face
void new_face_flist(vector<flist>& flist_te, vector<vert>& vert_te, vector<edge>& edge_te, 
face& face_te_new, int vlast0, int vlast, int elast0, int flistlast){
//the begining new added point of flist
int i0, m, tempi, sw; 
//cout<<"last "<<vlast0<<" "<<elast0<<endl;
i0=vlast0; m=1; 
face_te_new.fptr=flistlast; 
face_te_new.stat=3;
flist flist1;
for(int i1=0; i1<3; i1++){if(vert_te[i0].e[i1] >= elast0) {tempi=vert_te[i0].e[i1]; 
//cout<<"edge(begin) v "<<i0<<" e "<<tempi<<endl;
                                                           } }
while (m <= (vlast-vlast0)){
    if(m>1){ for(int i1=0; i1<3; i1++){ 
        if( vert_te[i0].e[i1]>= elast0 && vert_te[i0].e[i1] !=flist_te[flistlast+m-2].e)
        {tempi=vert_te[i0].e[i1]; //cout<<"next v "<<i0<<" e "<<tempi<<endl; 
        }      }
           }
    flist1.v=i0; flist1.e=tempi;
    //flist_te[flistlast+m-1].v=i0; flist_te[flistlast+m-1].e=tempi;
    if( m == (vlast-vlast0) )  flist1.link=flistlast; //flist_te[flistlast+m].link=flistlast+1; 
    else flist1.link=flistlast+m;
    flist_te.push_back(flist1);
    //flist_te[flistlast+m].link=flistlast+m+1;
    for(int i1=0; i1<2; i1++){  if(edge_te[tempi].v[i1] !=i0) sw=edge_te[tempi].v[i1];}
    i0=sw; m=m+1;
}   
}


