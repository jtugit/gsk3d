/* JTU created on 3/18/2022 */
#include "vector3.h"
#include "module.h"
#include "module_0.h"
#include "module_1.h"
#include <cmath>

foot3::foot3()
{
    xx=0.0; yy=0.0; zz=0.0;
}

inline void mfd_line_trace(double xmag, double ymag, double zmag, 
    double &xf, double &yf,double &zf)
{
    double r, theta, phi, Lshell, sintheta, theta0;
    double R0=1.0+AltitudeMin/radius;

    SPHCAR_08(r, theta, phi, xmag, ymag, zmag, -1);

    sintheta=sin(theta);
    Lshell=r/(sintheta*sintheta);

    theta0=acos(sqrt(R0/Lshell));

    if (zmag < 0.0) theta0 += pi/2.0;

    SPHCAR_08(R0, theta0, phi, xf, yf, zf, 1);
}

foot3 ****calcFootprints(GridsPoints *****ptrArray)
{
    double xmag,ymag,zmag,xf,yf,zf;
    int face, i, j, k;

    foot3 ****footArray;
    footArray = new foot3 ***[totalFace];

//    RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,params[0],params[1],params[2]);

    for (face = 0; face < totalFace; face++) {
        footArray[face] = new foot3 **[fieldsGridsSize + 3];

        for (i = 0; i <= fieldsGridsSize + 2; i++) {
            footArray[face][i] = new foot3 *[fieldsGridsSize + 3];

            for (j = 0; j <= fieldsGridsSize + 2; j++) {
                footArray[face][i][j] = new foot3 [fieldsGridsSize + 3];

                for (k = 0; k <= fieldsGridsSize * grid_domain; k++) {
                    if (i >= 1 && j >= 1) {
                        footArray[face][i][j][k] = foot3();

                        //convert lengths to unit Re
                        xmag=ptrArray[face][i][j][k]->pos3.x()/radius;
                        ymag=ptrArray[face][i][j][k]->pos3.y()/radius;
                        zmag=ptrArray[face][i][j][k]->pos3.z()/radius;

                        mfd_line_trace(xmag,ymag,zmag,xf,yf,zf);

                        // field line footprints at 90 km. length in meter
                        footArray[face][i][j][k].Setx(xf*radius);
                        footArray[face][i][j][k].Sety(yf*radius);
                        footArray[face][i][j][k].Setz(zf*radius);
                    }
                }
            }
        }
    }

    return footArray;
}