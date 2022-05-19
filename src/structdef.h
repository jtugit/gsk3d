#ifndef _STRUCTDEF_H_
#define _STRUCTDEF_H_
#include <iostream>
#include "parameters.h"


struct structg 
    {
        int face;
        int ig;  int jg;  int kg;
        int iw;  int jw;  int kw;
        double vx; double vy; double vz;
//        Vector3 vp;
//        double mass;
    };

struct structPar
{
    int face;
    int ig; int jg; int kg;
    double w000; double w100; double w010; double w110;
    double w001; double w101; double w011; double w111;
};
// Getface     
inline uint_64 Getface( double px_in, double py_in, double pz_in)
{
    return abs(px_in) > abs(py_in) ?
                abs(px_in) > abs(pz_in) ? 
                    px_in > 0 ? 0 : 3 :
                    pz_in > 0 ? 2 : 5 :
                abs(py_in) > abs(pz_in) ?
                    py_in > 0 ? 1 : 4 :
                    pz_in > 0 ? 2 : 5 ;
}
// inline weight calculation
inline void StructPar(struct structg &tempStr, struct structPar &tempStrPar)
{
    //
    tempStrPar.face = tempStr.face;
    tempStrPar.ig = tempStr.ig;
    tempStrPar.jg = tempStr.jg;
    tempStrPar.kg = tempStr.kg;
    //
    tempStrPar.w000 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (cellSize1 - tempStr.jw - 0.5) * (cellSize1 - tempStr.kw - 0.5);
    tempStrPar.w100 = 1.0 * (0.5 + tempStr.iw) * (cellSize1 - tempStr.jw - 0.5) * (cellSize1 - tempStr.kw - 0.5);
    tempStrPar.w010 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (1 + tempStr.jw - 0.5) * (cellSize1 - tempStr.kw - 0.5);
    tempStrPar.w110 = 1.0 * (0.5 + tempStr.iw) * (0.5 + tempStr.jw) * (cellSize1 - tempStr.kw - 0.5);
    tempStrPar.w001 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (cellSize1 - tempStr.jw - 0.5) * (0.5 + tempStr.kw);
    tempStrPar.w101 = 1.0 * (0.5 + tempStr.iw) * (cellSize1 - tempStr.jw - 0.5) * (0.5 + tempStr.kw);
    tempStrPar.w011 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (0.5 + tempStr.jw) * (0.5 + tempStr.kw);
    tempStrPar.w111 = 1.0 * (0.5 + tempStr.iw) * (0.5 + tempStr.jw) * (0.5 + tempStr.kw);
    //
    double totalWeight = tempStrPar.w000 + tempStrPar.w100 + tempStrPar.w010 + tempStrPar.w110 + tempStrPar.w001 + tempStrPar.w101 + tempStrPar.w011 + tempStrPar.w111;
    tempStrPar.w000 = tempStrPar.w000 / totalWeight;
    tempStrPar.w100 = tempStrPar.w100 / totalWeight;
    tempStrPar.w010 = tempStrPar.w010 / totalWeight;
    tempStrPar.w110 = tempStrPar.w110 / totalWeight;
    tempStrPar.w001 = tempStrPar.w001 / totalWeight;
    tempStrPar.w101 = tempStrPar.w101 / totalWeight;
    tempStrPar.w011 = tempStrPar.w011 / totalWeight;
    tempStrPar.w111 = tempStrPar.w111 / totalWeight;
}
#endif