#include"useful.h"
#include"lobgauss.h"
#include<stdio.h>

/*
    compute mass matrix for one element
    Parameters:
        h   : thickness of one element
        x,w : GLL points and weights
        rho : density at GLL points
        n   : number of points
    
    Returns:
        M   : mass matrix, a diagonal matrix, and we save it
              in a 1d array
*/
void mass_mat_one_element(double h,double *x,double *w,double *rho,int n,double *M)
{
    int i;
    for(i=0;i<n;i++)
        M[i]=w[i]*rho[i]*h/2;
}

/*
    compute stiff matrix for one element
    K_{ij}=\frac{2}{h}\sum_{k=0}^{n}\omega_{k}\mu_k \frac{dL_i}{dx_k}\frac{dL_j}{dx_k}

    Parameters:
        h   : thickness of one element
        x,w : GLL points and weights
        rho : density at GLL points
        n   : number of GLL points
    
    Returns:
        K   : stiff matrix 
*/
void stiff_mat_one_element(double h,double *x,double *w,double *mu,int n,double **K)
{
    int i,j,k;
    double s;
    double **L=dmat(n,n);

    // compute lagrange derivative
    lagrange_deriv(x,n,L);

    for(i=0;i<n;i++)
        for(j=0;j<n;j++){
            s=0.0;
            for(k=0;k<n;k++)
                s=s+w[k]*mu[k]*L[i][k]*L[j][k]*2/h;
            K[i][j]=s;
        }
    
    free_dmat(L,n);
}

/*
    compute global mass matrix

    Parameters:
        xe      ： grid points which delineate every elements shape(ne)
        h       : thickness of every elements
        ne      : number of grid points
        x0,w    : GLL points and weights    shape(n)
        rho     : density at GLL points     shape(ne)
        n       : number of GLL points
        M       : global mass matrix        shape(ng)
*/
void mass(double *xe,double *h,double *x0,double *w,double *rho,int ne,int n,double *M)
{
    double m[n];
    int i,j;
    int ng;

    ng=(ne-1)*(n-1)+1;   
    
    for(i=0;i<ng;i++)
        M[i]=0.0;
    
    for(i=0;i<ne-1;i++){
        mass_mat_one_element(h[i],x0,w,rho+i*(n-1),n,m);
        for(j=0;j<n;j++)
            M[i*(n-1)+j]+=m[j];
    }
}

/*
    compute global stiffness matrix

    Parameters:
        xe      ： grid points which delineate every elements
        h       : thickness of every elements
        ne      : number of grid points
        x0,w    : GLL points and weights
        n       : number of GLL points
        K       : global stiffness matrix
*/
void stiff(double *xe,double *h,double *x0,double *w,double *mu,int ne,int n,double **K)
{
    double **Ke;
    int i,j,k;
    int ng;

    ng=(ne-1)*(n-1)+1;
    Ke=dmat(n,n);

    for(i=0;i<ng;i++)
        for(j=0;j<ng;j++)
            K[i][j]=0.0;

    for(i=0;i<ne-1;i++){
        stiff_mat_one_element(h[i],x0,w,mu+i*(n-1),n,Ke);

        for(j=0;j<n;j++) 
            for(k=0;k<n;k++)    
                K[i*(n-1)+j][i*(n-1)+k]+=Ke[j][k];  
    }

    free_dmat(Ke,n);   
}