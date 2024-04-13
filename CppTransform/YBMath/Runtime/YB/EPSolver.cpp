#include "EPSolver.h"

RTAns::RTAns(MyMat _r, Vec _t, double _resi)
{
    R = _r;
    T = _t;
    E = MyMat::SetByTCross(T.m[0], T.m[1], T.m[2]) * R;
    Resi = _resi;
}


