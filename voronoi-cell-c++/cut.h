#ifndef CUT_H
#define CUT_H
#include <vector>
using namespace std;

//................
//define struct of edge vert flist and face;
struct edge{ int f[2], v[2], stat; };
struct vert{ double posv[4]; int e[3], stat; };
struct flist{int e, v, link;};
struct face{double dist; int fptr, stat, vfar;};


//calculate the distance
void calc_distan(double* length, double* lengthh, double* distan, double& tp);

 
 //cross product
void crossproduct(double* a, double* b, double* c);
void dotproduct(double* a, double* b, double& c, int length=3);

//volume and surface
void volume(vector<face>& face_te, vector<flist>& flist_te, vector<vert>& vert_te, 
int flast, double& vol, double& tot_surf, double* length, double* lengthh);

void surface(face& face0, vector<flist>& flist_te, vector<vert>& vert_te, double& surf, 
double* length, double* lengthh);

//synthesis the process
void bisectplane( double* pos0, double* pos1, double* length, double* lengthh, 
vector<vert>& vert_te, vector<edge>& edge_te, vector<face>& face_te , 
vector<flist>& flist_te, int& elast, int& vlast, int& flast,int& flistlast, bool& near);

//judge_vertex and status
void judge_vertex( double* pos0, double* pos1,double* pos2,  
double* length, double* lengthh, bool& vert_outside);

void judge_stat(double* pos0, double* pos1, double* length, double* lengthh, 
vector<vert>& vert_te, vector<edge>& edge_te, vector<face>& face_te, 
vector<flist>& flist_te, int vlast, int flast);

//determing the new vertex pos and status
void new_vertex_coeff(double* pos0, double* pos1, double* pos2, double* pos3, 
double* length, double* lengthh, double& a); 

void find_new_vert(double* pos0, double* pos1, double* length, double* lengthh, 
vector<vert>& vert_te, vector<edge>& edge_te, int elast, int& vlast, int& changed_face);


// determine the new edge cut face and new face
// the revised face of VC due to the exisitance of new face
 void cutting_face_accessory( face& face_te0, int index_face, vector<edge>& edge_te, 
 vector<vert>& vert_te, vector<flist>& flist_te, int& flistlast, int& elast, int flast, 
 int& first_vertex, int& first_cut_edge, int& second_cut_flist, int& second_vertex, int& second_cut_edge);
 // new face
void new_face_flist(vector<flist>& flist_te, vector<vert>& vert_te, vector<edge>& edge_te, 
face& face_te_new, int vlast0, int vlast, int elast0, int flistlast);




#endif