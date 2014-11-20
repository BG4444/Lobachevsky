#include "point.h"
#include <string.h>

#include <gsl/gsl_linalg.h>

#include <limits>
#include <vector>
#include <math.h>
# define M_PI		3.14159265358979323846	/* pi */


#define ASSERT(a) (a);

bool Point::intervalIntersection(double a, double b, double c, double d,double& e,double& f)
{
    const double p=a<c?a:c;
    const double q=a<c?b:d;
    const double r=a<c?c:a;
    const double s=a<c?d:b;

    const double B=std::max(p,q);
    const double C=std::min(r,s);
    const double D=std::max(r,s);

    if(C<=B)
    {
        e=C;
        if(D>B)
        {
             f=B;
             return(true);
        }
        f=D;
        return(true);
    }
    return(false);
}

Point* Point::make1dPoint(double x)
{
    return(new Point(&x,1));
}

void Point::sub(Point *a, Point *b, double *res)
{
    ASSERT(a);
    ASSERT(b);
    ASSERT(res);
    ASSERT(a->ndims==b->ndims);
    for(int i=0;i<a->ndims;i++)
    {
        res[i]=a->coord[i]-b->coord[i];
    }
}


Point* Point::make2dPoint(double x, double y)
{
    double temp[]={x,y};
    return(new Point(temp,2));
}
Point* Point::make3dPoint(double x, double y,double z)
{
    double temp[]={x,y,z};
    Point* ret=new Point(3);
    for(int i=0;i<3;i++)
    {
        ret->coord[i]=temp[i];
    }
    return(ret);
}


void Point::mul(double m)
{
    for(int i=0;i<ndims;i++)
    {
        coord[i]*=m;
    }
}

Point::Point(double *crd,int ndims):ndims(ndims),red(0),green(0),blue(0)
{
    coord=new double[ndims];
    memcpy(coord,crd,ndims*sizeof(double));
}

Point::~Point()
{
    delete[] coord;
}


void Point::fill(double *dst, double* src,int D)
{
    memcpy(dst,src,D*sizeof(double));
    for(int i=D;i<3;i++)
    {
        dst[i]=0;
    }
}

double Point::colDet2(Point *a, Point *b)
{
    ASSERT(a->ndims==b->ndims);
    ASSERT(a->ndims==2);
    return(a->coord[0]*b->coord[1]-a->coord[1]*b->coord[0]);
}

Point::Point(int ndims, double v):ndims(ndims)
{
    coord=new double[ndims];
    for(int i=0;i<ndims;i++)
    {
        coord[i]=v;
    }
}

Point::Point(Point *a, int D):ndims(std::max(a->ndims,D))
{
    coord=new double[this->ndims];
    for(int i=0;i<a->ndims;i++)
    {
        coord[i]=a->coord[i];
    }
    for(int i=a->ndims;i<D;i++)
    {
        coord[i]=0;
    }
}

Point* Point::dropLast()
{
    ASSERT(ndims>1);
    Point* ret=new Point(ndims-1);
    for(int i=0;i<ndims-1;i++)
    {
        ret->coord[i]=coord[i];
    }
    return(ret);
}

Point::Point(Point *a, Point *b):ndims(a->ndims)
{
    coord=new double[ndims];
    sub(b,a,coord);
}

double Point::norm2b()
{
    return(sqrt(norm2b2()));
}

double Point::norm2b2()
{
    double ret=0;
    for(int i=0;i<ndims;i++)
    {
        ret+=pow(coord[i],2);
    }
    return(ret);
}


void Point::add(Point *pnt)
{
    ASSERT(pnt);
    ASSERT(ndims==pnt->ndims);
    for(int i=0;i<ndims;i++)
    {
        coord[i]+=pnt->coord[i];
    }
}

void Point::sub(Point *pnt)
{
    ASSERT(pnt);
    ASSERT(ndims==pnt->ndims);
    for(int i=0;i<ndims;i++)
    {
        coord[i]-=pnt->coord[i];
    }
}

double Point::scal(Point *a, Point *b)
{
    ASSERT(a->ndims==b->ndims);
    double ret=0;
    for(int i=0;i<a->ndims;i++)
    {
        ret+=a->coord[i]*b->coord[i];
    }
    return(ret);
}

void Point::ortho()
{
    std::swap(coord[0],coord[1]);
    coord[1]*=-1;
}

double Point::distanceToHyperplane(Point *target, Point *normal,Point* onPlane)
{
    const double D=-scal(normal,onPlane);
    const double c=fabs(scal(target,normal)+D);
    const double subdiv=normal->norm2b();
    ASSERT(subdiv>0);
    return(c/subdiv);
}

Point* Point::normal(Point* a)
{
    ASSERT(a->ndims==2);
    Point* ret=new Point(a);
    ret->ortho();
    return(ret);
}


Point* Point::mid(Point *a, Point *b)
{
    ASSERT(a->ndims==b->ndims);
    Point* ret=new Point(a->ndims);
    for(int i=0;i<a->ndims;i++)
    {
        ret->coord[i]=(a->coord[i])/2.0+(b->coord[i])/2.0;
    }
    return(ret);
}

double Point::dist2(Point *a, Point *b)
{
    double ret=0;
    ASSERT(a->ndims==b->ndims);
    for(int i=0;i<a->ndims;i++)
    {
        ret+=pow(a->coord[i]-b->coord[i],2);
    }
    ASSERT(ret<std::numeric_limits<double>::infinity());
    return(ret);
}

double Point::dist(Point *a, Point *b)
{
    return(sqrt(dist2(a,b)));
}

void Point::print()
{

    for(int i=0;i<ndims;i++)
    {
        printf("%12.10lf,",coord[i]);
    }
//    printf(") %llx",int64_t(bc));

}

void Point::rotate(double cost, double sint)
{
    ASSERT(ndims==2);
    const double x= ( coord[0]*cost+coord[1]*sint);
    const double y= (-coord[0]*sint+coord[1]*cost);
    coord[0]=x;
    coord[1]=y;
}

double Point::distanceToLine(Point *A, Point *B)
{
    ASSERT(A->ndims==B->ndims);
    ASSERT(A->ndims==2);
    Point napr(A,B);
    Point link(A,this);
    Point normal(A,B);
    normal.ortho();

    const double mainDet=colDet2(&napr,&normal);
    const double pDet   =colDet2(&link,&normal);

    const double t=pDet/mainDet;
    if((t>=0)&&(t<=1))
    {
        const double det1=colDet2(&napr,this);
        const double det2=colDet2(&napr,A);

        return(fabs(det1-det2)/napr.norm2b());
    }
    else
    {
        return(std::numeric_limits<double>::infinity());
    }
}
