#ifndef SAMPLEPOLY_H
#define SAMPLEPOLY_H

#include <polycoeffproxy.h>

class SamplePoly :public PolyCoeffProxy
{
public:
    SamplePoly();
    mpf_class getCoeff(int idx);
    int power();
};

#endif // SAMPLEPOLY_H
