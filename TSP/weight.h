//
//  weight.h
//  TSP
//
//  Created by Sora Sugiyama on 9/26/24.
//

#ifndef weight_h
#define weight_h

#include "TSP.h"
#include <cmath>
#include <algorithm>

namespace TSP{

double distEUC(std::vector<double>&u,std::vector<double>&v,int dim){
    double ret=0;
    for(int i=0;i<dim;i++)ret+=(u[i]-v[i])*(u[i]-v[i]);
    return sqrt(ret);
}

double distMAN(std::vector<double>&u,std::vector<double>&v,int dim){
    double ret=0;
    for(int i=0;i<dim;i++)ret+=std::abs(u[i]-v[i]);
    return ret;
}

double distMAX(std::vector<double>&u,std::vector<double>&v,int dim){
    double ret=0;
    for(int i=0;i<dim;i++)ret=std::max(ret,std::abs(u[i]-v[i]));
    return ret;
}

// TSPLIB GEO uses this exact (truncated) value of PI, not full precision.
const double PI=3.141592;
const double RRR=6378.388;

// Convert a TSPLIB "DDD.MM" coordinate (integer part = degrees, fraction =
// minutes) to radians, per the TSPLIB spec: PI*(deg + 5*min/3)/180.
static double geoRad(double x){
    const double deg=int(x);
    const double min=x-deg;
    return PI*(deg+5.0*min/3.0)/180.0;
}

double distGEO(std::vector<double>&u,std::vector<double>&v){
    const double lat1=geoRad(u[0]),lng1=geoRad(u[1]);
    const double lat2=geoRad(v[0]),lng2=geoRad(v[1]);

    const double q1=std::cos(lng1-lng2),q2=std::cos(lat1-lat2),q3=std::cos(lat1+lat2);
    return int(RRR*std::acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3))+1.0);
}

double distATT(std::vector<double>&u,std::vector<double>&v){
    const double xd=u[0]-v[0],yd=u[1]-v[1];
    double r=sqrt((xd*xd+yd*yd)/10.0);
    double t=int(r);
    if(t<r)r=r+1.0;
    return r;
}

void makeWeightMatrix(problem &X){
    const int n=X.dimension;
    X.W=std::vector<std::vector<double> >(n,std::vector<double>(n,0.0));

    // Resolve the weight type once instead of doing up to nine string
    // comparisons per cell inside the n^2 loop.
    enum Code{EUC2,EUC3,MAX2,MAX3,MAN2,MAN3,CEIL2,GEOc,ATTc,NONE};
    const std::string &t=X.WeightType;
    Code code=
        t=="EUC_2D"?EUC2 : t=="EUC_3D"?EUC3 :
        t=="MAX_2D"?MAX2 : t=="MAX_3D"?MAX3 :
        t=="MAN_2D"?MAN2 : t=="MAN_3D"?MAN3 :
        t=="CEIL_2D"?CEIL2 : t=="GEO"?GEOc :
        t=="ATT"?ATTc : NONE;
    if(code==NONE)return;

    auto &C=X.coord;
    // Full row-major fill (sequential writes stay cache-friendly). Symmetry is
    // deliberately not exploited here: mirroring into W[j][i] is a column write
    // strided by a whole row and thrashes the cache for large n, costing more
    // than the halved (cheap) distance evaluations it would save.
    for(int i=0;i<n;i++){
        std::vector<double>&Wi=X.W[i];
        std::vector<double>&Ci=C[i];
        for(int j=0;j<n;j++){
            switch(code){
                case EUC2: Wi[j]=std::round(distEUC(Ci,C[j],2)); break;
                case EUC3: Wi[j]=std::round(distEUC(Ci,C[j],3)); break;
                case MAX2: Wi[j]=std::round(distMAX(Ci,C[j],2)); break;
                case MAX3: Wi[j]=std::round(distMAX(Ci,C[j],3)); break;
                case MAN2: Wi[j]=std::round(distMAN(Ci,C[j],2)); break;
                case MAN3: Wi[j]=std::round(distMAN(Ci,C[j],3)); break;
                case CEIL2:Wi[j]=std::ceil (distEUC(Ci,C[j],2)); break;
                case GEOc: Wi[j]=distGEO(Ci,C[j]); break;
                case ATTc: Wi[j]=distATT(Ci,C[j]); break;
                default:   Wi[j]=0; break;
            }
        }
    }
}

double loss(int rout[],problem &X){
    double ret=0;
    const int n=X.dimension;
    for(int i=0;i<n;i++)ret+=X.W[rout[i]][rout[(i+1)%n]];
    return ret;
}
double loss(std::vector<int>&rout,problem &X){
    double ret=0;
    const int n=X.dimension;
    for(int i=0;i<n;i++)ret+=X.W[rout[i]][rout[(i+1)%n]];
    return ret;
}

}

#endif /* weight_h */
