#include "lobachevskysolver.h"
#include <polycoeffproxy.h>
#include <point.h>

typedef mpf_class rational;
typedef int idx;

LobachevskySolver::LobachevskySolver()
{

}

Point* LobachevskySolver::solve(PolyCoeffProxy *prox,double meps,int nmaxiter)
{
    Point* ret=new Point(prox->power());
    rational* coeff[2];
    rational roots[prox->power()];
    for(idx i=0;i<2;i++)
    {
        coeff[i]=new rational[prox->power()+1];
    }
    for(idx i=0;i<prox->power();i++)
    {
        roots[i]=2;
    }
    for(idx i=0;i<=prox->power();i++)
    {
        coeff[0][i]=prox->getCoeff(i);
    }
    idx l=1;
    idx m=0;
    for(;m<nmaxiter;m++)
    {
        coeff[l][0]=coeff[1-l][0]*coeff[1-l][0];
        coeff[l][prox->power()]=coeff[1-l][prox->power()]*coeff[1-l][prox->power()];
        for(idx i=1;i<prox->power();i++)
        {
            const idx d=std::min(prox->power()-i,i);
            coeff[l][i]=coeff[1-l][i]*coeff[1-l][i];
            for(idx j=1;j<=d;j++)
            {

                idx k=(((1+j)%2)*2-1)<<1;
                coeff[l][i]+=k*coeff[1-l][i-j]*coeff[1-l][i+j];
            }
        }
        mpf_class eps=0;
        for(idx i=0;i<prox->power();i++)
        {
            const mpf_class znam=coeff[l][prox->power()-i];
            mpf_class tmp=coeff[l][prox->power()-1-i]/znam;
            tmp=abs(tmp);

            for(idx j=0;j<=m;j++)
            {
                tmp=sqrt(tmp);
            }
            const mpf_class ceps=abs(tmp-roots[i]);
            eps=std::max(eps,ceps);
            roots[i]=tmp;
        }
//        printf("%10.8lf\n",eps.get_d());
        if(eps<meps)
        {
            break;
        }
        l=1-l;
    }
    for(int i=0;i<prox->power();i++)
    {
        ret->coord[i]=roots[i].get_d();
    }
    for(int i=0;i<2;i++)
    {
        delete[] coeff[i];
    }

    return(ret);
}
