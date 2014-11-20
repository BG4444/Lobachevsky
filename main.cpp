#include <iostream>
#include <lobachevskysolver.h>
#include <samplepoly.h>
#include <point.h>
using namespace std;

int main()
{
    SamplePoly sample;
    Point* ret=LobachevskySolver::solve(&sample,1e-10,100);
    ret->print();
    delete ret;
    return 0;
}

