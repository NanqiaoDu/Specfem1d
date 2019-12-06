/*
    This program is used to generate Gauss-Lobatto-Legendre points and weights
    where GLL points are -1.0,1.0,plus the n-1 roots of pn^{\prime}(x),
    and weights are computed by using the relation:

        \omega_{i}=\frac{2}{(n-1)n p_{n-1}(x)^2}
*/

#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_interp.h>
#include"useful.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/*
    compute n-th order Legendre polynomial or its derivative

    Parameters:

        n       : order of Legendre polynomial
        x       : variable
        deri    : find derivative, could be 0 or 1(1-th order derivative)
*/
double pn(int n,double x,int deri)
{
    double fun;
    if(deri==0)
        fun=gsl_sf_legendre_Pl(n,x);
    else if(deri==1){

        double result[n+1],result_deri[n+1];

        gsl_sf_legendre_Pl_deriv_array(n,x,result,result_deri);
        fun=result_deri[n];
    }
    else{
        printf("deri should be 0 or 1, higher order derivative cannot be computed!\n");
        exit(0);
    }

    return fun;

}

/*
    use spline interpolate to find value of f(x0)
    with known values at 3 points x[0-2]

    Parameters:
        x,y : interpolate grids
        n   : length of x and y
        x0  : point to interpolate
*/
double neville(double *x,double *y,int n,double x0)
{   
    double fun;
    gsl_interp_accel *acc= gsl_interp_accel_alloc ();
    const gsl_interp_type *T=gsl_interp_polynomial;
    gsl_interp *interp=gsl_interp_alloc(T,n);
    gsl_interp_init(interp,x,y,n);

    fun=gsl_interp_eval(interp,x,y,x0,acc);

    gsl_interp_free(interp);
    gsl_interp_accel_free (acc);

    return fun;

}

/*
    find root of p_n^{\prime} in a interval [x0,x1] by using quadratic interpolation,
    where f(x0)*f(x1)<0
*/
double quadroot(double x0,double x1,int n)

{
    double f0,f1,fmid,xmid;
    double zero;
    double xarrc=pow(10.0,-10);
    double x[3],y[3];
    int count;

    xmid=(x0+x1)/2;
    f0=pn(n,x0,1);
    f1=pn(n,x1,1);
    fmid=pn(n,xmid,1);
    zero=x0;
    count=0;
    
    while(fabs(xmid-zero)>xarrc && (x1-x0>xarrc)){
        xmid=(x0+x1)/2;
        fmid=pn(n,xmid,1);
        if(f0<f1){
            x[0]=f0;
            x[1]=fmid;
            x[2]=f1;
            y[0]=x0;
            y[1]=xmid;
            y[2]=x1;
        }
        else{
            x[2]=f0;
            x[1]=fmid;
            x[0]=f1;
            y[2]=x0;
            y[1]=xmid;
            y[0]=x1;   
        }
        zero=neville(x,y,3,0.0);
        fmid=pn(n,zero,1);

        if(fabs(fmid)+fabs(f0)==fabs(fmid+f0)){
            x0=zero;
            f0=fmid;
        }
        else{
            x1=zero;
            f1=fmid;
        }

        count++;
        if(count==100){
            break;
        }
    } 
    return zero;
}
 // refine root in c1 and c2, search begin at x0
double refineroot(double x0,int n,double *c1,double *c2)
{
    double dc=0.001;
    double f1,f2;
    *c1=x0;
    *c2=x0+dc;

    f1=pn(n,*c1,1);
    f2=pn(n,*c2,1);

    while(f1*f2>0){
        *c1=*c2;
        *c2+=dc;
        f1=pn(n,*c1,1);
        f2=pn(n,*c2,1);
    }
}

/*
    Generate Gauss-Lobatto-Legendre points and weights
*/
void lobgauss(double *x,double *w,int n)
{
    double c1,c2;
    double dx=0.001;
    int i;

    x[0]=-1.0;
    x[n-1]=1.0;

    for(i=1;i<n-1;i++){

        refineroot(x[i-1]+dx,n-1,&c1,&c2);
        x[i]=quadroot(c1,c2,n-1);
    }

    for(i=0;i<n;i++)
        w[i]=2/(n*(n-1)*pow(pn(n-1,x[i],0),2));
}

// compute derivative of lagrange polynomial
void lagrange_deriv(double *x,int m,double **Lki)
{
    int n=m-1;
    int i,j,k;

    for(k=0;k<m;k++)
        for(i=0;i<m;i++)
            if(i==k && i==0)
                Lki[k][i]=-0.25*n*(n+1);
            else if(i==k && i==n)
                Lki[k][i]=0.25*n*(n+1);
            else if(i==k)
                Lki[k][i]=0.0;
            else
                Lki[k][i]=pn(n,x[i],0)/pn(n,x[k],0)/(x[i]-x[k]);
}