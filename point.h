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
    Point(int ndims,double v=0);
    ~Point();
    void print();

    double* coord;
    const int ndims;


private:
};

#endif // POINT_H
