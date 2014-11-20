#ifndef LOBACHEVSKYSOLVER_H
#define LOBACHEVSKYSOLVER_H

class PolyCoeffProxy;
class Point;


class LobachevskySolver
{
public:
    LobachevskySolver();
    static Point* solve(PolyCoeffProxy* prox, double meps, int nmaxiter);
};

#endif // LOBACHEVSKYSOLVER_H
