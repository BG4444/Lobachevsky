#ifndef POLYCOEFFPROXY_H
#define POLYCOEFFPROXY_H

#include <gmpxx.h>

class PolyCoeffProxy
{
public:
    PolyCoeffProxy();
    virtual mpf_class getCoeff(int idx)=0;
    virtual int power()=0;
    void print();
};

#endif // POLYCOEFFPROXY_H
