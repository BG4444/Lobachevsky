#ifndef POINT_H
#define POINT_H

//#include <commondb.h>
//#include "signs.h"
//#include <containerelement.h>
//#include <GL/gl.h>
//#include <ibpp.h>

class Point
{
public:
    Point(double * crd,int ndims);
    Point(Point* a,Point* b);
    Point(int ndims,double v=0);
    Point* dropLast();
    Point(Point* a,int D=0);
    static Point* make3dPoint(double x,double y,double z);
    static Point* make2dPoint(double x,double y);
    static Point* make1dPoint(double x);
    ~Point();
    double norm2b();
    double norm2b2();
    static Point* mid(Point* a,Point* b);
    static double dist(Point* a,Point* b);
    static double dist2(Point* a,Point* b);
    void print();
    static Point* normal(Point* a);
    void ortho();
    void rotate(double cost,double sint);
    void add(Point* pnt);
    void sub(Point* pnt);
    static double scal(Point* a,Point* b);
    static double colDet2(Point*a,Point* b);
    void draw(double h);
    void mul(double m);
    double distanceToLine(Point* A,Point* B);
    static void setPoint(double *src,double h,int ndims);
    static void fill(double *dst, double* src,int D);
    static void sub(Point* a,Point* b,double* res);
    static Point* vecMul(Point** arg);
    static double distanceToHyperplane(Point* target, Point* normal, Point *onPlane);
    static bool intervalIntersection(double a, double b, double c, double d,double& e,double& f);


    int red;
    int green;
    int blue;
    double* coord;
    const int ndims;


private:
};

#endif // POINT_H
