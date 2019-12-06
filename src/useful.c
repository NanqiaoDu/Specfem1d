/*
    This file contains some very useful array allocating and deallocating function 
    which can be used in scientific computation.

    Written by Nanqiao Du, IGGCAS, 2018.12.20
*/

#include<stdlib.h>

double *dvec(int n){
    double *a=(double*)calloc(sizeof(double),n);
    return a;
}

float *fvec(int n){
    float *a=(float *)calloc(sizeof(float),n);
    return a;
}

int *ivec(int n){
    int *a=(int *)calloc(sizeof(int),n);
    return a;
}


float **fmat(int row,int col){
    int i;
    float **a=(float **)calloc(sizeof(float*),row);

    for(i=0;i<row;i++) a[i]=fvec(col);
    return a;
}

int **imat(int row,int col){
    int i;
    int **a=(int **)calloc(sizeof(int*),row);

    for(i=0;i<row;i++) a[i]=ivec(col);

    return a;
}

double **dmat(int row,int col){
    int i;
    double **a=(double **)calloc(sizeof(double*),row);

    for(i=0;i<row;i++) a[i]=dvec(col);
    return a;
}

double ***d3tensor(int row,int col,int deep){
    int i;
    double ***a=(double ***)calloc(sizeof(double**),row);

    for(i=0;i<row;i++) a[i]=dmat(col,deep);
    return a;
}

void free_dvec(double *a){
    free(a);
    a=NULL;
}
void free_fvec(float *a){
    free(a);
    a=NULL;
}

void free_ivec(int *a){
    free(a);
    a=NULL;
}

void free_dmat(double **P,int row_p){
    int i;
    for(i=0;i<row_p;i++) free_dvec(P[i]);
    free(P);
    P=NULL;
}

void free_fmat(float **a,int row){
    int i;
    for(i=0;i<row;i++) free_fvec(a[i]);
    free(a);
    a=NULL;
}

void free_imat(int **a,int row){
    int i;
    for(i=0;i<row;i++) free_ivec(a[i]);
    free(a);
    a=NULL;  
}

void free_d3tensor(double ***P,int row,int col){
    int i;
    for(i=0;i<row;i++) free_dmat(P[i],col);
    free(P);
    P=NULL;
}

void matmul(double **a,double **b,int rowa,int cola,int colb,double **c)
/*
    Matrix multiply c=a*b
    a(rowa,cola),b(cola,colb),c(rowa,colb)
*/
{
    int i,j,k;
    double s;

    for(i=0;i<rowa;i++)
    for(j=0;j<colb;j++){
        s=0;
        for(k=0;k<cola;k++) 
            s+=a[i][k]*b[k][j];
        c[i][j]=s;
    }
}

void copyvec(double *a,double *b,int n){
    int i;

    for(i=0;i<n;i++)
        b[i]=a[i];
}

void copymat(double **a,double **b,int m,int n)
/*
    copy a to b
*/
{
    int i,j;
    for(i=0;i<m;i++)
    for(j=0;j<n;j++) b[i][j]=a[i][j];
}

double dvecmin(double *a,int n){
    int i;
    double min=a[0];

    for(i=1;i<n;i++)
        if(min>a[i])
            min=a[i];
    
    return min;
}

void diff(double *a,double *h,int n){
    int i;

    for(i=0;i<n-1;i++)
        h[i]=a[i+1]-a[i];
}