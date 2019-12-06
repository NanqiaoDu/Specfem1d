#include<math.h>
double rick(double t,double f0){
    double a,t0,c;

    a=4*f0;
    t0=1/f0;
    
    c=-2*a*(t-t0)*exp(-pow(a*(t-t0),2));
    
    return c;
}