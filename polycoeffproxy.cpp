#include "polycoeffproxy.h"

PolyCoeffProxy::PolyCoeffProxy()
{
}

void PolyCoeffProxy::print()
{
    for(int i=0;i<=power();i++)
    {
        printf("%lf\n",getCoeff(i).get_d());
    }
}
