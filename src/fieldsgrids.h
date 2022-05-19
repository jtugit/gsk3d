#ifndef _FIELDSGRIDS_H_
#define _FIELDSGRIDS_H_
/// Field grids information
#include <iostream>
#include "parameters.h"
#include "vector3.h"
#include <cmath>
using namespace std;

// set up grids container to contain class GridsPoints, the position 
// can be calculated but not stored
class GridsPoints
{
    public:
//************************************************************************
//************************************************************************
// initialize bx, by, bz using dipole field
//
//************************************************************************
//************************************************************************
inline void XYZtoB( Vector3 const& v)
{
    double r = sqrt(pow(v.x(),2.0) + pow(v.y(),2.0) + pow(v.z(),2.0));
    b3.Setx( 3.0 * dMoment * v.x() * v.z() / pow(r,5.0));
    b3.Sety( 3.0 * dMoment * v.y() * v.z() / pow(r,5.0));
    b3.Setz( dMoment * (3.0 * pow(v.z(),2.0) - pow(r,2.0)) / pow(r,5.0));    
}
//************************************************************************
//************************************************************************
// initialize vx, vy, vz using corotation assumption
// with pos3
//************************************************************************
//************************************************************************
inline void XYZtoVel()
{
    // update_type: 0-no current, 1-with current
    // move_type: 0-test particle, 1-not test particles
    // coordinate_sys: 0-lab coordinate, 1-rotating
    if( update_type == 0 && move_type == 0 && coordinate_rotate == 0)
    {
        Vector3 tempPos = pos3;
        Vector3 tempOmega = Vector3( 0.0, 0.0, omega_earth);
        tempPos.Setz( 0.0);
        v3 = tempOmega.CrossProduct(tempPos);
    } else if( update_type == 1 && move_type == 1)
    {
        Vector3 tempPos = pos3;
        Vector3 tempOmega = Vector3( 0.0, 0.0, omega_earth);
        tempPos.Setz( 0.0);
        v3 = tempOmega.CrossProduct(tempPos);
       
    } else
    {
      //  std::cout << " test " << "\n";
        v3 = Vector3( 0.0, 0.0, 0.0);
    }
    
}
//************************************************************************
//************************************************************************
// initialize vx, vy, vz of top boundary
//************************************************************************
//************************************************************************
inline void SetVe_Boundary( const Vector3& ve_in)
{
//    v3 = ve_in;
    ve3= ve_in;
}

//************************************************************************
//************************************************************************
// initialize vx, vy, vz of top boundary
//************************************************************************
//************************************************************************
inline void SetVel_Boundary( const Vector3& vel_in)
{
    v3 = vel_in;
    ve3= vel_in;
}
//************************************************************************
//************************************************************************
// initialize Ex, Ey, Ez using electron momentum equation 
// - vel X B
//************************************************************************
//************************************************************************
inline void XYZtoE()
{
    Vector3 temp; // for calculating grad(Pe) term
    // coordinate_sys: 0-at rest, 1-rotating
    if( pos3.norm() > 0 && coordinate_rotate == 0)
    {
//    temp = pos3.NormalizedVector().ScaleProduct( mi0_O * gravity); 
//    e3 = temp.PlusProduct( b3.CrossProduct(v3));
    e3 = b3.CrossProduct(v3);
    } else
    {
        e3 = Vector3(0.0, 0.0, 0.0);
    }
    
}

//************************************************************************
//************************************************************************
// initialize density
//************************************************************************
//************************************************************************
inline void VeltoE_topBoundary()
{
    e3 = b3.CrossProduct(v3);
}
//************************************************************************
//************************************************************************
// initialize density
//************************************************************************
//************************************************************************
inline void XYZtoDensity( int type)
{
    double initialDensity = 1.0e+7F;
    density_H = initialDensity;
    density_He= initialDensity;
    density_O = initialDensity;
}
inline void XYZtoDensity( )
{
/*    double scaleHeight = ikT / mi0_H / gravity;
    if( pos3.norm() > 0)
    density_H = N0_H * exp(-1 * (pos3.norm() - radius) / scaleHeight);

    scaleHeight = ikT / mi0_He / gravity;
    if( pos3.norm() > 0)
    density_He= N0_He * exp(-1 * (pos3.norm() - radius) / scaleHeight);
    
    scaleHeight = ikT / mi0_O / gravity;
    if( pos3.norm() > 0)
    density_O= N0_O * exp(-1 * (pos3.norm() - radius) / scaleHeight);  

    double r = pos3.norm() / radius;
    if( r > 0){
    density_H = N0_H / r * ( 1.0 - tanh( r - 6.5));
    density_He= N0_He / r * ( 1.0 - tanh( r - 6.5));
    density_O = N0_O / r * ( 1.0 - tanh( r - 6.5));
    }
*/
    // Set density for sin function
    //double PI = 3.1415926535897;
    double x = pos3.x();
    double y = pos3.y();
    double z = pos3.z();
    double longtitude = 0.0;
    double latitude = 0.0;
    double A = 0.5 * ( ne_density_max - ne_density_min);
    double A_average = 0.5 * ( ne_density_max + ne_density_min);
    double rho;
    double ratioH=0.0;
    double ratioHe=0.0;
    double ratioO=0.0;

    if( x == 0 && y == 0)
    { rho = A_average;}
    else if( x == 0 && y > 0)
    { longtitude = PI / 2.0;}
    else if( x == 0 && y < 0)
    { longtitude = PI / 2.0 * 3.0;}
    else if( x!= 0)
    {
    
    longtitude = atan( y / x);
    if( x<0) { longtitude = longtitude + PI;}
    }

    latitude = acos( z / sqrt( x*x + y*y + z*z));   
    
    double r = pos3.norm() / radius; // r is a non-unit value

    //double parameter0 = 0.25 * (1-tanh(r -7.25)) /(r-0.75);

    double parameter1 = 0.5 * ( 1.0 - tanh( r - 6.5)) / r;

    //double tempHight = ratioHight/radius + 1.0;
    
    if( update_type == 0)
    {
        if( r < LMin_centraldomain + 0.0001)
        {
            ratioH = ratioH_bot + ( ratioH2000 - ratioH_bot) / ( ratioHight - AltitudeMin/1000.0) * ( (r-1) * radius /1000.0 - AltitudeMin/1000.0);
            ratioO = ratioO_bot + ( ratioO2000 - ratioO_bot) / ( ratioHight - AltitudeMin/1000.0) * ( (r-1) * radius /1000.0 - AltitudeMin/1000.0);

            ratioHe = 1.0 - ratioH - ratioO;

            rho = ne_density_min;
            
            density_H = rho * ratioH * parameter1 ;
            density_He= rho *ratioHe * parameter1 ;
            density_O = rho * ratioO * parameter1 ;

      //     std::cout << rho << " " << ratioH << " " << parameter1 << " \n";
      //     std::cout << density_H << " " << density_He << " " << density_O << " \n";
        }

    } else 
    {
        //if( r < tempHight-0.001)
        if( r < LMin_maindomain + 0.01)
        {
        //     rho = A * sin( latitude) * sin( longtitude + PI / 2.0) * (1.0 - (r - LMin) / (tempHight - LMin)) + A_average;
            rho = A * sin(latitude) * sin(longtitude + PI / 2.0)  + A_average;
    
        }else
        {
            rho = ne_density_min;
        }

        // (1.0 - (r - LMin) / (LMax - LMin))
        //  try to make densities the same for surface at large r
        // notice * 0.1 factor
        //if( r < tempHight + 0.01)
        if( r < LMin_maindomain + 0.01)
        {
            ratioH = ratioH_bot + ( ratioH2000 - ratioH_bot) / ( ratioHight - AltitudeMin/1000.0) * ( (r-1) * radius /1000.0 - AltitudeMin/1000.0);
            ratioO = ratioO_bot + ( ratioO2000 - ratioO_bot) / ( ratioHight - AltitudeMin/1000.0) * ( (r-1) * radius /1000.0 - AltitudeMin/1000.0);

            ratioHe = 1.0 - ratioH - ratioO;


            density_H = rho * ratioH * parameter1 ;
            density_He= rho *ratioHe * parameter1 ;
            density_O = rho * ratioO * parameter1 ;
        }
        else if( r > LMax_maindomain - 0.01)
        {   
            density_H = rho * ratioH_top * parameter1;
            density_He= rho *ratioHe_top * parameter1;
            density_O = rho * ratioO_top * parameter1;

        } else
        {
        //    density_H = rho * initialDensityRate;
        //    density_He= rho * initialDensityRate;
        //    density_O = rho * initialDensityRate;

            ratioH = ratioH2000 + ( ratioH_top - ratioH2000) / ( (LMax-1)*radius/1000.0 - ratioHight/1000.0) * ( (r-1) * radius /1000.0 - ratioHight/1000.0);
            ratioO = ratioO2000 + ( ratioO_top - ratioO2000) / ( (LMax-1)*radius/1000.0 - ratioHight/1000.0) * ( (r-1) * radius /1000.0 - ratioHight/1000.0);

            ratioHe = 1.0 - ratioH - ratioO;

            density_H = rho * ratioH * parameter1 ;
            density_He= rho *ratioHe * parameter1 ;
            density_O = rho * ratioO * parameter1 ;
        }
    }
}
//************************************************************************
//************************************************************************
// initialize GradBNorm
//************************************************************************
//************************************************************************
inline void XYZtoGradBNorm()
{
    double r_xyz, r_xy;
    double sintheta, costheta, sinphi, cosphi;
    double fact_1, fact_2;

    r_xyz = pos3.norm();
    r_xy = pow( pow( pos3.x(), 2.0) + pow( pos3.y(), 2.0), 0.5);
    sintheta = r_xy / r_xyz;
    costheta = pos3.z() / r_xyz;
    sinphi = pos3.y() / r_xy;
    cosphi = pos3.x() / r_xy;

    fact_1 = - 3.0 * pow( 1+ 3* pow(costheta,2.0), 0.5) / pow( r_xyz, 4.0);
    fact_2 = - 3.0 * sintheta * costheta / pow( r_xyz, 4.0) / pow( 1+ 3* pow(costheta,2.0), 0.5);
        
    gradB3.Setx( fact_1 * sintheta * cosphi + fact_2 * costheta * cosphi);
    gradB3.Sety( fact_1 * sintheta * sinphi + fact_2 * costheta * sinphi);
    gradB3.Setz( fact_1 * costheta - fact_2 * sintheta);

}

//************************************************************************
//************************************************************************
// Initialization the pos3 for gridspoints
// giving face, i, j, k, all are int not int64
//************************************************************************
//************************************************************************

inline void InttoPos3( int face, int i, int j, int k)
{
    double px, py, pz;
    double temp[2];
    // 2.1 radial
    double L=1.01;
    if( grid_domain == 1) L = LMin * pow(10, logRatio *  k * cellSize1 );
    else if( grid_domain == 2)
    {
        if(k <= fieldsGridsSize)
            L = LMin + (LMid - LMin) * sinh(k * cellSize1 / grid_N1) / const_sinh1;
        else //if(k > fieldsGridsSize)
            L = LMid + (LMax - LMid) * sinh((k-fieldsGridsSize) * cellSize1 / grid_N2) / const_sinh2;
    } 
    else
    {
        std::cout << " grid_domain error \n";
        std::cin.get();
    }
    // 2.2 IgJg to ST note 0<ST<1
    temp[0] = (1.0 / fieldsGridsSize) * (i - 1);
    temp[1] = (1.0 / fieldsGridsSize) * (j - 1);
    // 2.3 ST to UV note -1<UV<1
    // tan(theta) = temp[i] = ST, note 0<ST<1
    temp[0] = tan(( temp[0] - 0.5) * pi / 2.0);
    temp[1] = tan(( temp[1] - 0.5) * pi / 2.0);

    // // faster access version
    // for (i=0; i<=1; i++)
    // {
    //     if (temp[i] >= 0.5) 
    //     {   
    //         temp[i]= (1.0/3.0) * (4.0*temp[i]*temp[i] - 1.0);
    //     }
    //     else
    //     {   
    //         temp[i]= (1.0/3.0) * (1.0 - 4.0*(1.0-temp[i])*(1.0-temp[i]));
    //     }
    // }
    // 2.4 UV to xyz 
    double kk = L * radius / sqrt(1.0 + temp[0] * temp[0] + temp[1] * temp[1]);
    switch (face)
    {
        case 0: px=1.0;           py=temp[0];     pz=temp[1]; break;
        case 1: px=-1.0*temp[0];  py=1.0;         pz=temp[1]; break;
        case 2: px=-1.0*temp[1];  py=temp[0];     pz=1.0;     break;
        case 3: px=-1.0;          py=-1.0*temp[0];pz=temp[1]; break;
        case 4: px=temp[0];       py=-1.0;        pz=temp[1]; break;
        default:px=temp[1];       py=temp[0];     pz=-1.0;    break;
    }
    px *= kk; py *= kk; pz *= kk;
    pos3 = Vector3( px, py, pz);

//    std::cout << L << " " << px << " " << py << " " << pz << std::endl;
}

//************************************************************************
//************************************************************************
// Reset parameters
//
//************************************************************************
//************************************************************************
inline void ResetParameters()
{
    density_H_cumu = 0.0;
    density_He_cumu = 0.0;
    density_O_cumu = 0.0;

    vH3_cumu = Vector3( 0.0, 0.0, 0.0);
    vHe3_cumu = Vector3( 0.0, 0.0, 0.0);
    vO3_cumu = Vector3( 0.0, 0.0, 0.0);
}

//************************************************************************
//************************************************************************
// Calculate weighting of density on grids as well as velocity
// weight: iw, jw, kw
// double number_in: number weight of each simualtion particle
// Vector3 vp_in: velocity of each simulation particle
//************************************************************************
//************************************************************************
inline void UpdateDueToWgt( int iw, int jw, int kw, double number_in, Vector3 vp_in, double ion_mass)
{
    int ionType_in;
    if( ion_mass == mi0_H) ionType_in = 1;
    if( ion_mass == mi0_He) ionType_in =4;
    if( ion_mass == mi0_O) ionType_in =16;
    switch (ionType_in)
    {
    case 1: {
    density_H  += number_in * iw * jw * kw / cellSize1/ cellSize1/ cellSize1; // acutally is number not number density
    vH3 = vH3.PlusProduct( Vector3( number_in * vp_in.x() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1,
                                     number_in * vp_in.y() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1,
                                    number_in * vp_in.z() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1 ));
    break;}
    case 4:{
    density_He += number_in * iw * jw * kw / cellSize1/ cellSize1/ cellSize1; // acutally is mass not density
    vHe3 = vHe3.PlusProduct( Vector3( number_in * vp_in.x() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1,
                                     number_in * vp_in.y() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1,
                                    number_in * vp_in.z() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1 ));
    break;}
    default:{
    density_O  += number_in * iw * jw * kw / cellSize1/ cellSize1/ cellSize1; // acutally is mass not density
    vO3 = vO3.PlusProduct( Vector3( number_in * vp_in.x() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1,
                                     number_in * vp_in.y() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1,
                                    number_in * vp_in.z() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1 ));
    break;}
    }
    
    v3 = v3.PlusProduct( Vector3( number_in * vp_in.x() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1,
                                  number_in * vp_in.y() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1,
                                  number_in * vp_in.z() * iw * jw * kw / cellSize1/ cellSize1/ cellSize1
                                ));
}

inline void UpdateDueToWgt( double weight, double number_in, Vector3 vp_in, double ion_mass)
{
    int ionType_in;
    if( ion_mass == mi0_H) ionType_in = 1;
    if( ion_mass == mi0_He) ionType_in =4;
    if( ion_mass == mi0_O) ionType_in =16;
    switch (ionType_in)
    {
    case 1:{
    density_H_cumu  += number_in * weight; // acutally is number not number density
    vH3_cumu = vH3_cumu.PlusProduct( Vector3( number_in * vp_in.x() * weight,
                                     number_in * vp_in.y() * weight,
                                    number_in * vp_in.z() * weight ));
    break;}
    case 4:{
    density_He_cumu += number_in * weight; // acutally is mass not density
    vHe3_cumu = vHe3_cumu.PlusProduct( Vector3( number_in * vp_in.x() * weight,
                                     number_in * vp_in.y() * weight,
                                    number_in * vp_in.z() * weight ));
    break;}
    default:{
    density_O_cumu  += number_in * weight; // acutally is mass not density
    vO3_cumu = vO3_cumu.PlusProduct( Vector3( number_in * vp_in.x() * weight,
                                     number_in * vp_in.y() * weight,
                                    number_in * vp_in.z() * weight ));
    break;}
    }
    
//   v3 = v3.PlusProduct( Vector3( number_in * vp_in.x() * weight,
//                                 number_in * vp_in.y() * weight,
//                                 number_in * vp_in.z() * weight
//                               ));
}

inline void UpdateDueToWgt( double weight, double number_in, Vector3 vp_in, int a)
{
    if( a ==0)
    {
        density_H_cumu  += number_in * weight; // acutally is number not number density
        vH3_cumu = vH3_cumu.PlusProduct( Vector3( number_in * vp_in.x() * weight,
                                     number_in * vp_in.y() * weight,
                                    number_in * vp_in.z() * weight ));
    }else if( a == 1)
    {
        density_He_cumu += number_in * weight; // acutally is mass not density
        vHe3_cumu = vHe3_cumu.PlusProduct( Vector3( number_in * vp_in.x() * weight,
                                     number_in * vp_in.y() * weight,
                                    number_in * vp_in.z() * weight ));
    }else if( a==2)
    {
        density_O_cumu  += number_in * weight; // acutally is mass not density
        vO3_cumu = vO3_cumu.PlusProduct( Vector3( number_in * vp_in.x() * weight,
                                     number_in * vp_in.y() * weight,
                                    number_in * vp_in.z() * weight ));
    }else
    {
        std::cout << " UpdateDueToWgt error \n";
    }
}

inline void UpdateDueToWgt( GridsPoints _gridsPoints, int a)
{
    if (a == 0)
    {
        density_H_cumu += _gridsPoints.density_H_cumu;
        vH3_cumu = vH3_cumu.PlusProduct( _gridsPoints.vH3_cumu);
    }
    else if (a == 1)
    {
        density_He_cumu += _gridsPoints.density_He_cumu;
        vHe3_cumu = vHe3_cumu.PlusProduct( _gridsPoints.vHe3_cumu);
    }
    else if (a == 2)
    {
        density_O_cumu += _gridsPoints.density_O_cumu;
        vO3_cumu = vO3_cumu.PlusProduct( _gridsPoints.vO3_cumu);
    }else
    {
        std::cout << " UpdateDueToWgt error \n";
    }
    //
}


// After all simulation are calculated once, the density means the total 
// mass at each grid points, and the v3 means the total momentum. 
inline void UpdateDueToWgt( GridsPoints***** ptrArray_in, double volume_in, int updateInfoPeriod_in)
{
    if ( density_H_cumu > 0.0 || density_He_cumu > 0.0 || density_O_cumu > 0.0)
    {
    //    v3 = vH3.PlusProduct(vHe3).PlusProduct(vO3).ScaleProduct(1.0 / (density_H + density_He + density_O) / updateInfoPeriod_in);
        v3 = vH3_cumu.PlusProduct(vHe3_cumu).PlusProduct(vO3_cumu).ScaleProduct(1.0 / (density_H_cumu + density_He_cumu + density_O_cumu));

 //       std::cout << vH3.x() << " " << vHe3.x() << " " << vO3.x() << " " <<  density_H << " " << density_He << " " << density_O << " " << volume_in << std::endl;
        if( density_H_cumu > 0.0)
        {
            vH3= vH3_cumu.ScaleProduct(1.0 / density_H_cumu);
        } else
        {
            vH3 = Vector3( 0.0, 0.0, 0.0);
        }

        if( density_He_cumu > 0.0)
        {
            vHe3= vHe3_cumu.ScaleProduct(1.0 / density_He_cumu);
        } else
        {
            vHe3 = Vector3( 0.0, 0.0, 0.0);
        }

        if( density_O_cumu > 0.0)
        {
            vO3= vO3_cumu.ScaleProduct(1.0 / density_O_cumu);
        } else
        {
            vO3 = Vector3( 0.0, 0.0, 0.0);
        }
        
        density_H = density_H_cumu / volume_in / updateInfoPeriod_in * normalized_N;
        density_He = density_He_cumu / volume_in / updateInfoPeriod_in * normalized_N;
        density_O = density_O_cumu / volume_in / updateInfoPeriod_in * normalized_N;
    } else
    {
        v3 = Vector3( 0.0, 0.0, 0.0);
    }
}


// update E from electron's momentum equation
// E = - Ve X B - grad Pe / N / qi0
// input grad Pe at grids

inline void updateE( Vector3 GradPe_in)
{
    if( update_type == 0)
    {
        if( density_H >0.0 || density_He > 0.0 || density_O > 0.0){
            e3 = b3.CrossProduct(v3).MinusProduct(GradPe_in.ScaleProduct(1.0 / (density_H + density_He + density_O ) / qi0));
           // e3 = GradPe_in.ScaleProduct( -1.0 / (density_H + density_He + density_O ) / qi0);
            
        } else
        {
            e3 = Vector3( 0.0, 0.0, 0.0);
        }
        
    } else
    {
        if( density_H > 0.0 || density_He > 0.0 || density_O > 0.0){
            e3 = b3.CrossProduct(ve3).MinusProduct(GradPe_in.ScaleProduct(1.0 / (density_H + density_He + density_O ) / qi0));
        } else
        {
            e3 = Vector3( 0.0, 0.0, 0.0);
        }
        
    }
    
}


// Update ve from ampere's law
// ve = vi - curl B / mu0 / qi0 / total number density
// input curl B
inline void updateve3( Vector3 curlB_in)
{
    if( density_H > 0.0 || density_He > 0.0 || density_O > 0.0){
    ve3 = v3.MinusProduct( curlB_in.ScaleProduct( 1.0 / (density_H + density_He + density_O ) / qi0 / mu0));
    } else
    {
        ve3 = Vector3( 0.0, 0.0, 0.0);
    }
    
}

// update B from Faraday's Law
// B = B - tstep * (curl E)
// input (curl E)
void updatedB( Vector3 E_in)
{

    dB3 = dB3.PlusProduct( E_in.ScaleProduct(-1 * tstep * updateInfoPeriod));
}


// return Vector3 gradB3
inline Vector3 GradB3()
{
    return gradB3;
}

// return Vector3 gradPe
inline Vector3 GradPe()
{
    return gradPe;
}


// return Vector3 e3
inline Vector3 E3()
{
    return e3;
}

inline Vector3 DE3()
{
    return dE3;
}

inline void SetE3( const Vector3& E3_other)
{
    e3.SetVector3( E3_other);
}

inline void SetdE3( const Vector3& E3_other)
{
    dE3.SetVector3( E3_other);
}

// return Vector3 b3 only
inline Vector3 B3_base()
{
    return b3;
}
// return Vector3 B3
inline Vector3 B3()
{ 
    return b3.PlusProduct(dB3);
}

// return Vector3 dB3
inline Vector3 DB3()
{
    return dB3;
}

inline void SetdB3( const Vector3& dB3_other)
{
    dB3 = dB3_other;
}

// return Vector3 pos3
inline Vector3 Pos3()
{
    return pos3;
}

// return density
inline double Density()
{
    return density_H + density_He + density_O;
}

// return density of H
inline double Density_H()
{
    return density_H;
}
inline double Density_H_cumu()
{
    return density_H_cumu;
}

// return density of H
inline double Density_He()
{
    return density_He;
}
// return density of H
inline double Density_He_cumu()
{
    return density_He_cumu;
}

// return density of H
inline double Density_O()
{
    return density_O;
}
// return density of H
inline double Density_O_cumu()
{
    return density_O_cumu;
}

// set density 
inline void Density_H( double den_in)
{
    density_H = den_in;
}
inline void PlusMassH3( const double other_mH3)
{
    density_H += other_mH3;
}
inline void Density_He( double den_in)
{
    density_He = den_in;
}
inline void PlusMassHe3( const double other_mHe3)
{
    density_He += other_mHe3;
}
inline void Density_O( double den_in)
{
    density_O = den_in;
}
inline void PlusMassO3( const double other_mO3)
{
    density_O += other_mO3;
}
// return velocity
inline Vector3 Vel3()
{
    return v3;
}

inline void SetVel3( Vector3 other_v3)
{
    v3 = other_v3;
}


inline Vector3 VelH3()
{
    return vH3;
}

inline Vector3 VelH3_cumu()
{
    return vH3_cumu;
}

inline void PlusVelH3( const Vector3& other_vH3)
{
    vH3 = vH3.PlusProduct(other_vH3);
}

inline Vector3 VelHe3()
{
    return vHe3;
}

inline Vector3 VelHe3_cumu()
{
    return vHe3_cumu;
}

inline void PlusVelHe3( const Vector3& other_vHe3)
{
    vHe3 = vHe3.PlusProduct( other_vHe3);
}

inline Vector3 VelO3()
{
    return vO3;
}

inline Vector3 VelO3_cumu()
{
    return vO3_cumu;
}

inline void PlusVelO3( const Vector3& other_vO3)
{
    vO3 = vO3.PlusProduct( other_vO3);
}

// return velocity e3
inline Vector3 Vel_e3()
{
    return ve3;
}

inline void SetVel_e3( const Vector3& other_v3)
{
//    ve3.SetVector3( other_v3);
    ve3 = other_v3;
}

inline void SetGradNormB( const Vector3& other_GradNormB)
{
    gradB3.SetVector3( other_GradNormB);
}

inline void SetGradPe( const Vector3& other_GradPe)
{
    gradPe.SetVector3( other_GradPe);
}
// return stopSign
inline int StopSign()
{
    return stopSign;
}
// return PSD_H
inline double PSD_H()
{
    return psd_H;
}
inline void SetPSD_H( double psd)
{
    psd_H = psd;
}
inline void SetPotential( double potential_in)
{
    potential = potential_in;
}
inline double Potential()
{
    return potential;
}
// set temperature
inline void SetTemperature( double temperature_in)
{
    temperature = temperature_in;
}
inline void SetTemperature( double a, double b, double c)
{
    temperature_H = a;
    temperature_He = b;
    temperature_O = c;
}
inline void SetTemperatureIons( double temp)
{
    temperature_H = temp;
    temperature_He = temp;
    temperature_O = temp;
}
inline double Temperature( )
{
    return temperature;
}
inline double Temperature_H()
{
    return temperature_H;
}
inline double Temperature_He()
{
    return temperature_He;
}
inline double Temperature_O()
{
    return temperature_O;
}
// set stopSign
inline void SetStopSign( int stopSign_in)
{
    stopSign = stopSign_in;
}
// Set grids
inline void CopyGridsPoints(const GridsPoints& other)
{  
    pos3 = Vector3( other.pos3);
    e3 = Vector3( other.e3);
    dE3 = Vector3( other.dE3);
    b3 = Vector3( other.b3);
    dB3= Vector3( other.dB3);
    v3 = Vector3( other.v3);
    vH3 = Vector3( other.vH3);
    vHe3 = Vector3( other.vHe3);
    vO3 = Vector3( other.vO3);
    vH3_cumu = Vector3( other.vH3_cumu);
    vHe3_cumu = Vector3( other.vHe3_cumu);
    vO3_cumu = Vector3( other.vO3_cumu);
    ve3= Vector3( other.ve3);
    gradB3 = Vector3( other.gradB3);
    gradPe = Vector3( other.gradPe);

    density_H = other.density_H;
    density_He = other.density_He;
    density_O = other.density_O;
    
    density_H_cumu = other.density_H_cumu;
    density_He_cumu = other.density_He_cumu;
    density_O_cumu = other.density_O_cumu;
    
    temperature = other.temperature;
    temperature = other.temperature;
    temperature = other.temperature;
    temperature = other.temperature;
    potential = other.potential;
    psd_H = other.psd_H;

    stopSign = other.stopSign;
}
// set
inline void ReadGridsPoints(double px_in, double py_in, double pz_in,
                            double ex_in, double ey_in, double ez_in,
                            double dEx_in, double dEy_in, double dEz_in,
                            double bx_in, double by_in, double bz_in,
                            double dBx_in,double dBy_in,double dBz_in,
                            double vx_in, double vy_in, double vz_in,
                            double vHx_in, double vHy_in, double vHz_in,
                            double vHex_in, double vHey_in, double vHez_in,
                            double vOx_in, double vOy_in, double vOz_in,
                            double vHx_in_cumu, double vHy_in_cumu, double vHz_in_cumu,
                            double vHex_in_cumu, double vHey_in_cumu, double vHez_in_cumu,
                            double vOx_in_cumu, double vOy_in_cumu, double vOz_in_cumu,
                            double vex_in, double vey_in, double vez_in,
                            double gradBx_in, double gradBy_in, double gradBz_in,
                            double gradPex_in, double gradPey_in, double gradPez_in,
                            double density_H_in,
                            double density_He_in,
                            double density_O_in,
                            double density_H_in_cumu,
                            double density_He_in_cumu,
                            double density_O_in_cumu,
                            double temperature_in,
                            double temperature_H_in,
                            double temperature_He_in,
                            double temperature_O_in,
                            double potential_in,
                            double psd_H_in,
                            int stopSign_in
                            )
{
    stopSign = stopSign_in;
    temperature = temperature_in;
    temperature_H = temperature_H_in;
    temperature_He = temperature_He_in;
    temperature_O = temperature_O_in;

    pos3 = Vector3( px_in, py_in, pz_in);
    e3 =   Vector3( ex_in, ey_in, ez_in);
    dE3 =  Vector3( dEx_in, dEy_in, dEz_in);
    b3 =   Vector3( bx_in, by_in, bz_in);
    dB3=   Vector3( dBx_in, dBy_in, dBz_in);
    v3 =   Vector3( vx_in, vy_in, vz_in);
    vH3=   Vector3( vHx_in, vHy_in, vHz_in);
    vHe3=  Vector3( vHex_in, vHey_in, vHez_in);
    vO3=   Vector3( vOx_in, vOy_in, vOz_in);
    vH3=   Vector3( vHx_in_cumu, vHy_in_cumu, vHz_in_cumu);
    vHe3=  Vector3( vHex_in_cumu, vHey_in_cumu, vHez_in_cumu);
    vO3=   Vector3( vOx_in_cumu, vOy_in_cumu, vOz_in_cumu);
    ve3=   Vector3( vex_in, vey_in, vez_in);
    gradB3 = Vector3( gradBx_in, gradBy_in, gradBz_in);
    gradPe = Vector3( gradPex_in, gradPey_in, gradPez_in);

    density_H = density_H_in;
    density_He= density_He_in;
    density_O = density_O_in;

    density_H_cumu = density_H_in_cumu;
    density_He_cumu= density_He_in_cumu;
    density_O_cumu = density_O_in_cumu;

    potential = potential_in;
    psd_H = psd_H_in;    
}
                            
//////////////////////////    
    // Constructors
    GridsPoints( double px_in, double py_in, double pz_in,
                 double ex_in, double ey_in, double ez_in,
                 double dEx_in, double dEy_in, double dEz_in,
                 double bx_in, double by_in, double bz_in,
                 double dBx_in,double dBy_in,double dBz_in,
                 double vx_in, double vy_in, double vz_in,
                 double vHx_in, double vHy_in, double vHz_in,
                 double vHex_in, double vHey_in, double vHez_in,
                 double vOx_in, double vOy_in, double vOz_in,
                 double vHx_cumu, double vHy_cumu, double vHz_cumu,
                 double vHex_cumu, double vHey_cumu, double vHez_cumu,
                 double vOx_cumu, double vOy_cumu, double vOz_cumu,
                 double vex_in, double vey_in, double vez_in,
                 double gradBx_in, double gradBy_in, double gradBz_in,
                 double gradPex_in, double gradPey_in, double gradPez_in,
                 double density_H_in, double density_He_in, double density_O_in,
                 double density_H_cumu, double density_He_cumu, double density_O_cumu,
                 double temperature_in,
                 double temperature_H,
                 double temperature_He,
                 double temperature_O,
                 double potential,
                 double psd_H,
                 int stopSign_in);
                 
    GridsPoints( const GridsPoints& other);

    GridsPoints();

    
//    private:
    Vector3 pos3;            // pos3: position for vector 3
    
    Vector3 e3;             // used for moving particles, averaged by dE3
    Vector3 dE3;             // used for propagating waves

    Vector3 b3;             // b3
    Vector3 dB3;            // perturbation of b3
    
    Vector3 v3;             // average ions velocity3
    Vector3 vH3;
    Vector3 vHe3;
    Vector3 vO3;

    Vector3 vH3_cumu;
    Vector3 vHe3_cumu;
    Vector3 vO3_cumu;

    Vector3 ve3;            // e velocity

    Vector3 gradB3;          // gradB in which B is norm of vector3 B
    Vector3 gradPe;         // gradPe

    double density_H;    // number density
    double density_He;   // number density
    double density_O;  

    double density_H_cumu;   
    double density_He_cumu;  
    double density_O_cumu;   

    double temperature;     // electron temperature: Te
    double temperature_H;
    double temperature_He;
    double temperature_O;

    double potential;   // E
    double psd_H;       // PSD power spectrel density at gyro of H

    int stopSign;
//    int face; int gi; int gj; int gk; //face, i, j, in fieldsgrids, radial
};
#endif