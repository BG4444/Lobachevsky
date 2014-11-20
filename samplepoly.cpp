#include "samplepoly.h"

SamplePoly::SamplePoly()
{
}

mpf_class SamplePoly::getCoeff(int idx)
{
    static const mpf_class coeff[]={12,-8,1};
    return(coeff[idx]);
}

int SamplePoly::power()
{
    return(2);
}
