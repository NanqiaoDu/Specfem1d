#include"source.h"
#include"useful.h"

void finitediff(double *u,double *M,double **K,int ng,int nt,
                double dt,double f0)
{

    double unew[ng],uold[ng],**tmp,**tmp1;
    double f[ng];
    int i,j;

    tmp=dmat(ng,1);
    tmp1=dmat(ng,1);

    for(i=0;i<nt;i++){
        f[ng/2]=rick(i*dt,f0);

        for(j=0;j<ng;j++)
            tmp[j][0]=u[j];
        matmul(K,tmp,ng,ng,1,tmp1);

        for(j=0;j<ng;j++)
            unew[j]=dt*dt*1/M[j]*(f[j]-tmp1[j][0])+2*u[j]-uold[j];

        copyvec(u,uold,ng);
        copyvec(unew,u,ng);
    }

    free_dmat(tmp,ng);
    free_dmat(tmp1,ng);
}