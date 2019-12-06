#include<stdio.h>
#include<stdlib.h>
#include"lobgauss.h"
#include"useful.h"
#include"specmat.h"
#include"finitediff.h"

int main(){

    double xmax,v,ep,rho0,mu0,f0;
    double dt,dx;
    int i,j,k;
    int n,ne,ng,nt;
    double *x0,*w,*xg,*xe,*h,*hg;
    double *rho,*mu;
    double *u;
    FILE *fp;

    // media parameters
    xmax=10000;
    rho0=2000;
    v=2500;
    mu0=v*v*rho0;
    
    // init elements parameters
    ne=250; // number of elements
    n=4;    // number of GLL points in each element
    ng=ne*(n-1)+1;  // total points
    ep=0.2;         // cfl constants
    nt=10000;       // time step evolved

    xe=dvec(ne+1); // grid points of each element
    xg=dvec(ng);   // grid points of every points in every element
    h=dvec(ne);    // thickness of every element
    rho=dvec(ng);   // density
    mu=dvec(ng);    // rigidity

    // Generate GLL points and weights
    x0=dvec(n);     
    w=dvec(n);      // GLL points and weights
    lobgauss(x0,w,n);

    // init element's metadata
    for(i=0;i<ne+1;i++)
        xe[i]=xmax/ne*i;

    diff(xe,h,ne+1);   
    
    for(i=0;i<ne;i++)
        for(j=0;j<n;j++)
            xg[i*(n-1)+j]=h[i]*(x0[j]+1)/2+xe[i];
    
    for(i=0;i<ng;i++){
        rho[i]=rho0;
        mu[i]=mu0;
    }

    // define dx and dt
    hg=dvec(ng-1);
    diff(xg,hg,ng);
    dx=dvecmin(hg,ng-1);
    dt=dx*ep/v;
    free_dvec(hg);

    // source parameters
    f0=5;

    // init wavefield
    u=dvec(ng);

    // init mass and stiff matrix
    double *M=dvec(ng);
    double **K=dmat(ng,ng);

    // compute global mass and stiffness matrix
    mass(xe,h,x0,w,rho,ne+1,n,M);
    stiff(xe,h,x0,w,mu,ne+1,n,K);

    // finite difference to extrapolate wavefield
    finitediff(u,M,K,ng,nt,dt,f0);

    // save results
    fp=fopen("out.txt","w");
    for(i=0;i<ng;i++)
        fprintf(fp,"%g %g \n",xg[i],u[i]);
    
    fclose(fp);

    // deallocate space
    free_dvec(x0);
    free_dvec(w);
    free_dvec(xg);
    free_dvec(xe);
    free_dvec(h);
    free_dvec(M);
    free_dvec(u);
    free_dmat(K,ng);

    system("python show.py");

    return 0;
}