#define pi 3.1415926535898

/*
    declarition of array,matrix and 3-d tensor allocate function
*/
double *dvec(int n);
float *fvec(int n);
int *ivec(int n);
double **dmat(int row,int col);
float **fmat(int row,int col);
int **imat(int row,int col);
double ***d3tensor(int row,int col,int deep);

/*
    Deallocate function
*/
void free_dvec(double *a);
void free_fvec(float *a);
void free_ivec(int *a);
void free_dmat(double **P,int row_p);
void free_fmat(float **a,int row);
void free_imat(int **a,int row);
void free_d3tensor(double ***P,int row,int col);

/*
    Matrix handling function
*/
void copyvec(double *a,double *b,int n);
void matmul(double **a,double **b,int rowa,int cola,int colb,double **c);
void copymat(double **a,double **b,int m,int n);
double dvecmin(double *a,int n);
void diff(double *a,double *h,int n);