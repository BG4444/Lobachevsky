#include "point.h"
#include <string.h>

#include <gsl/gsl_linalg.h>

#include <limits>
#include <vector>
#include <math.h>
# define M_PI		3.14159265358979323846	/* pi */


#define ASSERT(a) (a);

Point::~Point()
{
    delete[] coord;
}


Point::Point(int ndims, double v):ndims(ndims)
{
    coord=new double[ndims];
    for(int i=0;i<ndims;i++)
    {
        coord[i]=v;
    }
}

void Point::print()
{

    for(int i=0;i<ndims;i++)
    {
        printf("%12.10lf,",coord[i]);
    }
//    printf(") %llx",int64_t(bc));

}
