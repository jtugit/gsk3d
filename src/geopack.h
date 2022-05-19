#ifndef _GEOPACK_H_
#define _GEOPACK_H_

#include "parameters.h"
#include <cmath>

using namespace std;
/******************************************************************************
C
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
C
C  LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
*/
inline int SUN_08(int IYEAR,int IDAY,int IHOUR,int MIN,int ISEC,double &GST,double &SLONG,
           double &SRASN,double &SDEC)
{
//      IMPLICIT REAL*8 (A-H,O-Z)
    double FDAY, DJ, T, VL, GG, OBLIQ, SOB, SLP, SIND, COSD, SC;
    const double RAD=57.295779513e0;

    if(IYEAR<1901 || IYEAR>2099) return 0;
    FDAY=double(IHOUR*3600+MIN*60+ISEC)/86400.e0;
    DJ=365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY-0.5e0+FDAY;
    T=DJ/36525.e0;
    VL=fmod(279.696678e0+0.9856473354e0*DJ,360.e0);
    GST=fmod(279.690983e0+.9856473354e0*DJ+360.e0*FDAY+180.e0,360.e0)/RAD;
    GG=fmod(358.475845e0+0.985600267e0*DJ,360.e0)/RAD;
    SLONG=(VL+(1.91946e0-0.004789e0*T)*sin(GG)+0.020094e0*sin(2.e0*GG))/RAD;
    if (SLONG > 6.2831853e0) SLONG=SLONG-6.283185307e0;
    if (SLONG < 0.e0) SLONG=SLONG+6.283185307e0;
    OBLIQ=(23.45229e0-0.0130125e0*T)/RAD;
    SOB=sin(OBLIQ);
    SLP=SLONG-9.924e-5;
/*
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION DUE TO
C   EARTH'S ORBITAL MOTION
*/
    SIND=SOB*sin(SLP);
    COSD=sqrt(1.e0-SIND*SIND);
    SC=SIND/COSD;
    SDEC=atan(SC);
    SRASN=3.141592654e0-atan2(cos(OBLIQ)/SOB*SC,-cos(SLP)/COSD);

    return 0;
}

/*  CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICE VERSA
    (THETA AND PHI IN RADIANS).
*/
inline void SPHCAR_08(double &ra, double &thet, double &phis, double &x, double &y,
               double &z, int dir)
{
/*                 dir > 0           dir < 0
-----INPUT:   dir,R,THETA,PHI     dir,X,Y,Z
-----OUTPUT:      X,Y,Z           R,THETA,PHI

  NOTE: AT THE POLES (X=0 AND Y=0) WE ASSUME PHI=0 WHEN CONVERTING
        FROM CARTESIAN TO SPHERICAL COORDS (I.E., FOR J<0)
*/
    double sq; //, pi=acos(-1.0);

    if (dir < 0) {
        sq=x*x+y*y;
        ra=sqrt(sq+z*z);

        if (sq == 0.0) {
            phis=0.0;
            if (z >= 0.0) thet=0.0;
            else thet=pi;
        }
        else {
            sq=sqrt(sq);
            phis=atan2(y, x);
            thet=atan2(sq, z);
            if (phis < 0.0) phis=phis+2.0*pi;
            
        }
    }
    else {
        sq=ra*sin(thet);
        x=sq*cos(phis);
        y=sq*sin(phis);
        z=ra*cos(thet);
    }

    return;
}

/* CALCULATES LOCAL SPHERICAL FIELD COMPONENTS FROM THOSE IN CARTESIAN SYSTEM

-----INPUT:   X,Y,Z  - CARTESIAN COMPONENTS OF THE POSITION VECTOR
              BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD VECTOR
-----OUTPUT:  BR,BTHETA,BPHI - LOCAL SPHERICAL COMPONENTS OF THE FIELD VECTOR

  NOTE: AT THE POLES (THETA=0 OR THETA=PI) WE ASSUME PHI=0,
        AND HENCE BTHETA=BX, BPHI=BY
*/
inline void BCARSP_08(double X,double Y,double Z,double BX,double BY,double BZ,
        double &BR,double &BTHETA,double &BPHI)
{
    double RHO2, R, RHO, CPHI, SPHI, CT, ST;

    RHO2=X*X+Y*Y;
    R=sqrt(RHO2+Z*Z);
    RHO=sqrt(RHO2);

    if (RHO != 0.0) {
        CPHI=X/RHO;
        SPHI=Y/RHO;
    }
    else {
        CPHI=1.0;
        SPHI=0.0;
    }

    CT=Z/R;
    ST=RHO/R;

    BR=(X*BX+Y*BY+Z*BZ)/R;
    BTHETA=(BX*CPHI+BY*SPHI)*CT-BZ*ST;
    BPHI=BY*CPHI-BX*SPHI;
}

/* CALCULATES LOCAL CARTESIAN FIELD COMPONENTS FROM THOSE INSPHERICAL SYSTEM
    (YH. 2021)
-----INPUT:   theta, phi, r
              BR,BTHETA,BPHI - LOCAL SPHERICAL COMPONENTS OF THE FIELD VECTOR
-----OUTPUT:   BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD VECTOR

  NOTE: AT THE POLES (THETA=0 OR THETA=PI) WE ASSUME PHI=0,
        AND HENCE BTHETA=BX, BPHI=BY
*/
inline void BCARSPRR_08(double theta, double phi, double BR, double BTHETA,
                        double BPHI, double &BX, double &BY, double &BZ)
{
    // r_hat = sin(theta)*cos(phi)*x_hat + sin(theta)sin(phi)y_hat + cos(theta)z_hat;
    // theta_hat = cos(theta)cos(phi)x_hat + cos(theta)sin(phi)y_hat - sin(theta)z_hat;
    // phi_hat = - sin(phi)x_hat + cos(phi)y_hat;
    BX = BR* sin(theta) * cos(phi) + BTHETA * cos(theta) * cos(phi) - BPHI * sin(phi);
    BY = BR* sin(theta)*sin(phi) + BTHETA* cos(theta)*sin(phi) + BPHI * cos(phi);
    BZ = BR * cos(theta) - BTHETA * sin(theta);
}
/*=====================================================================================
C  1. PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN
C     SEVERAL COORDINATE SYSTEMS, MOST FREQUENTLY USED IN SPACE PHYSICS.
C
C  2. PREPARES COEFFICIENTS USED IN THE CALCULATION OF THE MAIN GEOMAGNETIC FIELD
C      (IGRF MODEL)
C
C  THIS SUBROUTINE SHOULD BE INVOKED BEFORE USING THE FOLLOWING SUBROUTINES:
C  IGRF_GEO_08, IGRF_GSW_08, DIP_08, GEOMAG_08, GEOGSW_08, MAGSW_08, SMGSW_08, GSWGSE_08,
c  GEIGEO_08, TRACE_08, STEP_08, RHAND_08.
C
C  THERE IS NO NEED TO REPEATEDLY INVOKE RECALC_08, IF MULTIPLE CALCULATIONS ARE MADE
C    FOR THE SAME DATE/TIME AND SOLAR WIND FLOW DIRECTION.
C
C-----INPUT PARAMETERS:
C
C     IYEAR   -  YEAR NUMBER (FOUR DIGITS)
C     IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
C     IHOUR -  HOUR OF DAY (00 TO 23)
C     MIN   -  MINUTE OF HOUR (00 TO 59)
C     ISEC  -  SECONDS OF MINUTE (00 TO 59)
C     VGSEX,VGSEY,VGSEZ - GSE (GEOCENTRIC SOLAR-ECLIPTIC) COMPONENTS OF THE OBSERVED
C                              SOLAR WIND FLOW VELOCITY (IN KM/S)
C
C  IMPORTANT: IF ONLY QUESTIONABLE INFORMATION (OR NO INFORMATION AT ALL) IS AVAILABLE
C             ON THE SOLAR WIND SPEED, OR, IF THE STANDARD GSM AND/OR SM COORDINATES ARE
C             INTENDED TO BE USED, THEN SET VGSEX=-400.0 AND VGSEY=VGSEZ=0. IN THIS CASE,
C             THE GSW COORDINATE SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM.
C
C             IF ONLY SCALAR SPEED V OF THE SOLAR WIND IS KNOWN, THEN SETTING
C             VGSEX=-V, VGSEY=29.78, VGSEZ=0.0 WILL TAKE INTO ACCOUNT THE ~4 degs
C             ABERRATION OF THE MAGNETOSPHERE DUE TO EARTH'S ORBITAL MOTION
C
C             IF ALL THREE GSE COMPONENTS OF THE SOLAR WIND VELOCITY ARE AVAILABLE,
C             PLEASE NOTE THAT IN SOME SOLAR WIND DATABASES THE ABERRATION EFFECT
C             HAS ALREADY BEEN TAKEN INTO ACCOUNT BY SUBTRACTING 29.78 KM/S FROM VYGSE;
C             IN THAT CASE, THE UNABERRATED (OBSERVED) VYGSE VALUES SHOULD BE RESTORED
C             BY ADDING BACK THE 29.78 KM/S CORRECTION. WHETHER OR NOT TO DO THAT, MUST
C             BE EITHER VERIFIED WITH THE DATA ORIGINATOR OR DETERMINED BY AVERAGING
C             VGSEY OVER A SUFFICIENTLY LONG TIME INTERVAL.
C
C-----OUTPUT PARAMETERS:  NONE (ALL OUTPUT QUANTITIES ARE PLACED
C                         INTO THE COMMON BLOCKS /GEOPACK1/ AND /GEOPACK2/)
C
C    OTHER SUBROUTINES CALLED BY THIS ONE: SUN_08
C
C    AUTHOR:  N.A. TSYGANENKO
C    DATE:    DEC.1, 1991
C
C    REVISION OF JANUARY 31, 2015:
C
C     The table of IGRF coefficients was extended to include those for the epoch 2015 (igrf-12)
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
C
C    REVISION OF NOVEMBER 30, 2010:
C
C     The table of IGRF coefficients was extended to include those for the epoch 2010
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
C
C    REVISION OF NOVEMBER 15, 2007: ADDED THE POSSIBILITY TO TAKE INTO ACCOUNT THE OBSERVED
C     DEFLECTION OF THE SOLAR WIND FLOW FROM STRICTLY RADIAL DIRECTION. TO THAT END, THREE
C     GSE COMPONENTS OF THE SOLAR WIND VELOCITY WERE ADDED TO THE INPUT PARAMETERS.
C
c    CORRECTION OF MAY 9, 2006:  INTERPOLATION OF THE COEFFICIENTS (BETWEEN
C     LABELS 50 AND 105) IS NOW MADE THROUGH THE LAST ELEMENT OF THE ARRAYS
C     G(105)  AND H(105) (PREVIOUSLY MADE ONLY THROUGH N=66, WHICH IN SOME
C     CASES CAUSED RUNTIME ERRORS)
c
C    REVISION OF MAY 3, 2005:
C     The table of IGRF coefficients was extended to include those for the epoch 2005
c       the maximal order of spherical harmonics was also increased up to 13
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
c
C    REVISION OF APRIL 3, 2003:
c    The code now includes preparation of the model coefficients for the subroutines
c    IGRF_08 and GEOMAG_08. This eliminates the need for the SAVE statements, used
c    in the old versions, making the codes easier and more compiler-independent.
*/
inline void RECALC_08(int IYEAR,int IDAY,int IHOUR,int IMIN,int ISEC,double VGSEX,double VGSEY,
               double VGSEZ)
{
    int    IY, N, N22, M, MN, MNN;
    double F1, F2, DT;
    double S, P, AA1, G_10, G_11, H_11;
    double SQ, SQQ, SQR;
    double S1, S2, S3, GST, SLONG, SRASN, SDEC;
    double DJ, T, OBLIQ, DZ1, DZ2, DZ3;
    double DY1, DY2, DY3;
    double V, DX1, DX2, DX3;
    double X1, X2, X3;
    double DIP1, DIP2, DIP3;
    double Y1, Y2, Y3, Y;
    double Z1, Z2, Z3;
    double EXMAGX, EXMAGY, EXMAGZ, EYMAGX, EYMAGY;

    //COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,
     //* SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33,
     //* E11,E21,E31,E12,E22,E32,E13,E23,E33
/*
C  THE COMMON BLOCK /GEOPACK1/ CONTAINS ELEMENTS OF THE ROTATION MATRICES AND OTHER
C   PARAMETERS RELATED TO THE COORDINATE TRANSFORMATIONS PERFORMED BY THIS PACKAGE
C
      COMMON /GEOPACK2/ G(105),H(105),REC(105)
C
C  THE COMMON BLOCK /GEOPACK2/ CONTAINS COEFFICIENTS OF THE IGRF FIELD MODEL, CALCULATED
C    FOR A GIVEN YEAR AND DAY FROM THEIR STANDARD EPOCH VALUES. THE ARRAY REC CONTAINS
C    COEFFICIENTS USED IN THE RECURSION RELATIONS FOR LEGENDRE ASSOCIATE POLYNOMIALS.
*/
    /*DIMENSION G65(105),H65(105),G70(105),H70(105),G75(105),H75(105),
           G80(105),H80(105),G85(105),H85(105),G90(105),H90(105),G95(105),
           H95(105),G00(105),H00(105),G05(105),H05(105),G10(105),H10(105),
           G15(105),H15(105),DG15(45),DH15(45)*/

    double G65[106]={0.0,0.e0,-30334.e0,-2119.e0,-1662.e0,2997.e0,1594.e0,1297.e0,
      -2038.e0,1292.e0,856.e0,957.e0,804.e0,479.e0,-390.e0,252.e0,
      -219.e0,358.e0,254.e0,-31.e0,-157.e0,-62.e0,45.e0,61.e0,8.e0,
      -228.e0,4.e0,1.e0,-111.e0,75.e0,-57.e0,4.e0,13.e0,-26.e0,-6.e0,
      13.e0,1.e0,13.e0,5.e0,-4.e0,-14.e0,0.e0,8.e0,-1.e0,11.e0,4.e0,
      8.e0,10.e0,2.e0,-13.e0,10.e0,-1.e0,-1.e0,5.e0,1.e0,-2.e0,-2.e0,
      -3.e0,2.e0,-5.e0,-2.e0,4.e0,4.e0,0.e0,2.e0,2.e0,0.e0,39*0.e0};

    double H65[106]={0.0,0.e0,0.e0,5776.e0,0.e0,-2016.e0,114.e0,0.e0,-404.e0,
      240.e0,-165.e0,0.e0,148.e0,-269.e0,13.e0,-269.e0,0.e0,19.e0,
      128.e0,-126.e0,-97.e0,81.e0,0.e0,-11.e0,100.e0,68.e0,-32.e0,-8.e0,
      -7.e0,0.e0,-61.e0,-27.e0,-2.e0,6.e0,26.e0,-23.e0,-12.e0,0.e0,7.e0,
      -12.e0,9.e0,-16.e0,4.e0,24.e0,-3.e0,-17.e0,0.e0,-22.e0,15.e0,7.e0,
      -4.e0,-5.e0,10.e0,10.e0,-4.e0,1.e0,0.e0,2.e0,1.e0,2.e0,6.e0,-4.e0,
      0.e0,-2.e0,3.e0,0.e0,-6.e0,39*0.e0};

    double G70[106]={0.0,0.e0,-30220.e0,-2068.e0,-1781.e0,3000.e0,1611.e0,1287.e0,
     -2091.e0,1278.e0,838.e0,952.e0,800.e0,461.e0,-395.e0,234.e0,
     -216.e0,359.e0,262.e0,-42.e0,-160.e0,-56.e0,43.e0,64.e0,15.e0,
     -212.e0,2.e0,3.e0,-112.e0,72.e0,-57.e0,1.e0,14.e0,-22.e0,-2.e0,
     13.e0,-2.e0,14.e0,6.e0,-2.e0,-13.e0,-3.e0,5.e0,0.e0,11.e0,3.e0,
     8.e0,10.e0,2.e0,-12.e0,10.e0,-1.e0,0.e0,3.e0,1.e0,-1.e0,-3.e0,
     -3.e0,2.e0,-5.e0,-1.e0,6.e0,4.e0,1.e0,0.e0,3.e0,-1.e0,39*0.e0};

    double H70[106]={0.0,0.e0,0.e0,5737.e0,0.e0,-2047.e0,25.e0,0.e0,-366.e0,
     251.e0,-196.e0,0.e0,167.e0,-266.e0,26.e0,-279.e0,0.e0,26.e0,
     139.e0,-139.e0,-91.e0,83.e0,0.e0,-12.e0,100.e0,72.e0,-37.e0,-6.e0,
     1.e0,0.e0,-70.e0,-27.e0,-4.e0,8.e0,23.e0,-23.e0,-11.e0,0.e0,7.e0,
     -15.e0,6.e0,-17.e0,6.e0,21.e0,-6.e0,-16.e0,0.e0,-21.e0,16.e0,6.e0,
     -4.e0,-5.e0,10.e0,11.e0,-2.e0,1.e0,0.e0,1.e0,1.e0,3.e0,4.e0,-4.e0,
     0.e0,-1.e0,3.e0,1.e0,-4.e0,39*0.e0};

    double G75[106]={0.0,0.e0,-30100.e0,-2013.e0,-1902.e0,3010.e0,1632.e0,1276.e0,
     -2144.e0,1260.e0,830.e0,946.e0,791.e0,438.e0,-405.e0,216.e0,
     -218.e0,356.e0,264.e0,-59.e0,-159.e0,-49.e0,45.e0,66.e0,28.e0,
     -198.e0,1.e0,6.e0,-111.e0,71.e0,-56.e0,1.e0,16.e0,-14.e0,0.e0,
     12.e0,-5.e0,14.e0,6.e0,-1.e0,-12.e0,-8.e0,4.e0,0.e0,10.e0,1.e0,
     7.e0,10.e0,2.e0,-12.e0,10.e0,-1.e0,-1.e0,4.e0,1.e0,-2.e0,-3.e0,
     -3.e0,2.e0,-5.e0,-2.e0,5.e0,4.e0,1.e0,0.e0,3.e0,-1.e0,39*0.e0};

    double H75[106]={0.0,0.e0,0.e0,5675.e0,0.e0,-2067.e0,-68.e0,0.e0,-333.e0,
     262.e0,-223.e0,0.e0,191.e0,-265.e0,39.e0,-288.e0,0.e0,31.e0,
     148.e0,-152.e0,-83.e0,88.e0,0.e0,-13.e0,99.e0,75.e0,-41.e0,-4.e0,
     11.e0,0.e0,-77.e0,-26.e0,-5.e0,10.e0,22.e0,-23.e0,-12.e0,0.e0,
     6.e0,-16.e0,4.e0,-19.e0,6.e0,18.e0,-10.e0,-17.e0,0.e0,-21.e0,
     16.e0,7.e0,-4.e0,-5.e0,10.e0,11.e0,-3.e0,1.e0,0.e0,1.e0,1.e0,3.e0,
     4.e0,-4.e0,-1.e0,-1.e0,3.e0,1.e0,-5.e0,39*0.e0};

    double G80[106]={0.0,0.e0,-29992.e0,-1956.e0,-1997.e0,3027.e0,1663.e0,1281.e0,
     -2180.e0,1251.e0,833.e0,938.e0,782.e0,398.e0,-419.e0,199.e0,
     -218.e0,357.e0,261.e0,-74.e0,-162.e0,-48.e0,48.e0,66.e0,42.e0,
     -192.e0,4.e0,14.e0,-108.e0,72.e0,-59.e0,2.e0,21.e0,-12.e0,1.e0,
     11.e0,-2.e0,18.e0,6.e0,0.e0,-11.e0,-7.e0,4.e0,3.e0,6.e0,-1.e0,
     5.e0,10.e0,1.e0,-12.e0,9.e0,-3.e0,-1.e0,7.e0,2.e0,-5.e0,-4.e0,
     -4.e0,2.e0,-5.e0,-2.e0,5.e0,3.e0,1.e0,2.e0,3.e0,0.e0,39*0.e0};

    double H80[106]={0.0,0.e0,0.e0,5604.e0,0.e0,-2129.e0,-200.e0,0.e0,-336.e0,
     271.e0,-252.e0,0.e0,212.e0,-257.e0,53.e0,-297.e0,0.e0,46.e0,
     150.e0,-151.e0,-78.e0,92.e0,0.e0,-15.e0,93.e0,71.e0,-43.e0,-2.e0,
     17.e0,0.e0,-82.e0,-27.e0,-5.e0,16.e0,18.e0,-23.e0,-10.e0,0.e0,
     7.e0,-18.e0,4.e0,-22.e0,9.e0,16.e0,-13.e0,-15.e0,0.e0,-21.e0,
     16.e0,9.e0,-5.e0,-6.e0,9.e0,10.e0,-6.e0,2.e0,0.e0,1.e0,0.e0,3.e0,
     6.e0,-4.e0,0.e0,-1.e0,4.e0,0.e0,-6.e0,39*0.e0};

    double G85[106]={0.0,0.e0,-29873.e0,-1905.e0,-2072.e0,3044.e0,1687.e0,1296.e0,
     -2208.e0,1247.e0,829.e0,936.e0,780.e0,361.e0,-424.e0,170.e0,
     -214.e0,355.e0,253.e0,-93.e0,-164.e0,-46.e0,53.e0,65.e0,51.e0,
     -185.e0,4.e0,16.e0,-102.e0,74.e0,-62.e0,3.e0,24.e0,-6.e0,4.e0,
     10.e0,0.e0,21.e0,6.e0,0.e0,-11.e0,-9.e0,4.e0,4.e0,4.e0,-4.e0,5.e0,
     10.e0,1.e0,-12.e0,9.e0,-3.e0,-1.e0,7.e0,1.e0,-5.e0,-4.e0,-4.e0,
     3.e0,-5.e0,-2.e0,5.e0,3.e0,1.e0,2.e0,3.e0,0.e0,39*0.e0};

    double H85[106]={0.0,0.e0,0.e0,5500.e0,0.e0,-2197.e0,-306.e0,0.e0,-310.e0,
     284.e0,-297.e0,0.e0,232.e0,-249.e0,69.e0,-297.e0,0.e0,47.e0,
     150.e0,-154.e0,-75.e0,95.e0,0.e0,-16.e0,88.e0,69.e0,-48.e0,-1.e0,
     21.e0,0.e0,-83.e0,-27.e0,-2.e0,20.e0,17.e0,-23.e0,-7.e0,0.e0,8.e0,
     -19.e0,5.e0,-23.e0,11.e0,14.e0,-15.e0,-11.e0,0.e0,-21.e0,15.e0,
     9.e0,-6.e0,-6.e0,9.e0,9.e0,-7.e0,2.e0,0.e0,1.e0,0.e0,3.e0,6.e0,
     -4.e0,0.e0,-1.e0,4.e0,0.e0,-6.e0,39*0.e0};

    double G90[106]={0.0,0.e0,-29775.e0,-1848.e0,-2131.e0,3059.e0,1686.e0,1314.e0,
          -2239.e0,  1248.e0,  802.e0,  939.e0, 780.e0, 325.e0,-423.e0,
            141.e0,  -214.e0,  353.e0,  245.e0,-109.e0,-165.e0, -36.e0,
             61.e0,    65.e0,   59.e0, -178.e0,   3.e0,  18.e0, -96.e0,
             77.e0,   -64.e0,    2.e0,   26.e0,  -1.e0,   5.e0,   9.e0,
              0.e0,    23.e0,    5.e0,   -1.e0, -10.e0, -12.e0,   3.e0,
              4.e0,     2.e0,   -6.e0,    4.e0,   9.e0,   1.e0, -12.e0,
              9.e0,    -4.e0,   -2.e0,    7.e0,   1.e0,  -6.e0,  -3.e0,
             -4.e0,     2.e0,   -5.e0,   -2.e0,   4.e0,   3.e0,   1.e0,
              3.e0,     3.e0,    0.e0,39*0.e0};

    double H90[106]={0.0,0.e0,  0.e0,5406.e0,   0.e0,-2279.e0,-373.e0,  0.e0,
           -284.e0,293.e0,-352.e0,   0.e0,  247.e0,-240.e0, 84.e0,
           -299.e0,  0.e0,  46.e0, 154.e0, -153.e0, -69.e0, 97.e0,
              0.e0,-16.e0,  82.e0,  69.e0,  -52.e0,   1.e0, 24.e0,
              0.e0,-80.e0, -26.e0,   0.e0,   21.e0,  17.e0,-23.e0,
             -4.e0,  0.e0,  10.e0, -19.e0,    6.e0, -22.e0, 12.e0,
             12.e0,-16.e0, -10.e0,   0.e0,  -20.e0,  15.e0, 11.e0,
             -7.e0, -7.e0,   9.e0,   8.e0,   -7.e0,   2.e0,  0.e0,
              2.e0,  1.e0,   3.e0,   6.e0,   -4.e0,   0.e0, -2.e0,
              3.e0, -1.e0,  -6.e0,39*0.e0};

    double G95[106]={0.0,0.e0,-29692.e0,-1784.e0,-2200.e0,3070.e0,1681.e0,1335.e0,
          -2267.e0,  1249.e0,  759.e0,  940.e0, 780.e0, 290.e0,-418.e0,
            122.e0,  -214.e0,  352.e0,  235.e0,-118.e0,-166.e0, -17.e0,
             68.e0,    67.e0,   68.e0, -170.e0,  -1.e0,  19.e0, -93.e0,
             77.e0,   -72.e0,    1.e0,   28.e0,   5.e0,   4.e0,   8.e0,
             -2.e0,    25.e0,    6.e0,   -6.e0,  -9.e0, -14.e0,   9.e0,
              6.e0,    -5.e0,   -7.e0,    4.e0,   9.e0,   3.e0, -10.e0,
              8.e0,    -8.e0,   -1.e0,   10.e0,  -2.e0,  -8.e0,  -3.e0,
             -6.e0,     2.e0,   -4.e0,   -1.e0,   4.e0,   2.e0,   2.e0,
              5.e0,     1.e0,    0.e0,  39*0.e0};

    double H95[106]={0.0,0.e0,  0.e0,5306.e0,  0.e0,-2366.e0,-413.e0,  0.e0,
           -262.e0,302.e0,-427.e0,  0.e0,  262.e0,-236.e0, 97.e0,
           -306.e0,  0.e0,  46.e0,165.e0, -143.e0, -55.e0,107.e0,
              0.e0,-17.e0,  72.e0, 67.e0,  -58.e0,   1.e0, 36.e0,
              0.e0,-69.e0, -25.e0,  4.e0,   24.e0,  17.e0,-24.e0,
             -6.e0,  0.e0,  11.e0,-21.e0,    8.e0, -23.e0, 15.e0,
             11.e0,-16.e0,  -4.e0,  0.e0,  -20.e0,  15.e0, 12.e0,
             -6.e0, -8.e0,   8.e0,  5.e0,   -8.e0,   3.e0,  0.e0,
              1.e0,  0.e0,   4.e0,  5.e0,   -5.e0,  -1.e0, -2.e0,
              1.e0, -2.e0,  -7.e0,39*0.e0};

    double G00[106]={0.0,0.e0,-29619.4e0,-1728.2e0,-2267.7e0,3068.4e0,1670.9e0,
          1339.6e0,  -2288.e0, 1252.1e0,  714.5e0, 932.3e0, 786.8e0,
            250.e0,   -403.e0,  111.3e0, -218.8e0, 351.4e0, 222.3e0,
          -130.4e0,  -168.6e0,  -12.9e0,   72.3e0,  68.2e0,  74.2e0,
          -160.9e0,    -5.9e0,   16.9e0,  -90.4e0,  79.0e0, -74.0e0,
              0.e0,    33.3e0,    9.1e0,    6.9e0,   7.3e0,  -1.2e0,
            24.4e0,     6.6e0,   -9.2e0,   -7.9e0, -16.6e0,   9.1e0,
             7.0e0,    -7.9e0,    -7.e0,     5.e0,   9.4e0,    3.e0,
           - 8.4e0,     6.3e0,   -8.9e0,   -1.5e0,   9.3e0,  -4.3e0,
            -8.2e0,    -2.6e0,    -6.e0,    1.7e0,  -3.1e0,  -0.5e0,
             3.7e0,      1.e0,     2.e0,    4.2e0,   0.3e0,  -1.1e0,
             2.7e0,    -1.7e0,   -1.9e0,    1.5e0,  -0.1e0,   0.1e0,
            -0.7e0,     0.7e0,    1.7e0,    0.1e0,   1.2e0,   4.0e0,
            -2.2e0,    -0.3e0,    0.2e0,    0.9e0,  -0.2e0,   0.9e0,
            -0.5e0,     0.3e0,   -0.3e0,   -0.4e0,  -0.1e0,  -0.2e0,
            -0.4e0,    -0.2e0,   -0.9e0,    0.3e0,   0.1e0,  -0.4e0,
             1.3e0,    -0.4e0,    0.7e0,   -0.4e0,   0.3e0,  -0.1e0,
             0.4e0,      0.e0,    0.1e0};

    double H00[106]={0.0,0.e0,   0.e0,5186.1e0,   0.e0,-2481.6e0,-458.0e0,   0.e0,
          -227.6e0,293.4e0,-491.1e0,   0.e0,  272.6e0,-231.9e0,119.8e0,
          -303.8e0,   0.e0,  43.8e0,171.9e0, -133.1e0, -39.3e0,106.3e0,
              0.e0,-17.4e0,  63.7e0, 65.1e0,  -61.2e0,   0.7e0, 43.8e0,
              0.e0,-64.6e0, -24.2e0,  6.2e0,    24.e0,  14.8e0,-25.4e0,
            -5.8e0,  0.0e0,  11.9e0,-21.5e0,    8.5e0, -21.5e0, 15.5e0,
             8.9e0,-14.9e0,  -2.1e0,  0.0e0,  -19.7e0,  13.4e0, 12.5e0,
            -6.2e0, -8.4e0,   8.4e0,  3.8e0,   -8.2e0,   4.8e0,  0.0e0,
             1.7e0,  0.0e0,   4.0e0,  4.9e0,   -5.9e0,  -1.2e0, -2.9e0,
             0.2e0, -2.2e0,  -7.4e0,  0.0e0,    0.1e0,   1.3e0, -0.9e0,
            -2.6e0,  0.9e0,  -0.7e0, -2.8e0,   -0.9e0,  -1.2e0, -1.9e0,
            -0.9e0,  0.0e0,  -0.4e0,  0.3e0,    2.5e0,  -2.6e0,  0.7e0,
             0.3e0,  0.0e0,   0.0e0,  0.3e0,   -0.9e0,  -0.4e0,  0.8e0,
             0.0e0, -0.9e0,   0.2e0,  1.8e0,   -0.4e0,  -1.0e0, -0.1e0,
             0.7e0,  0.3e0,   0.6e0,  0.3e0,   -0.2e0,  -0.5e0, -0.9e0};

    double G05[106]={0.0,0.e0,  -29554.6e0,-1669.0e0,-2337.2e0,3047.7e0,1657.8e0,
          1336.3e0,   -2305.8e0, 1246.4e0,  672.5e0, 920.6e0, 798.0e0,
           210.7e0,    -379.9e0,  100.0e0, -227.0e0, 354.4e0, 208.9e0,
          -136.5e0,    -168.1e0,  -13.6e0,   73.6e0,  69.6e0,  76.7e0,
          -151.3e0,     -14.6e0,   14.6e0,  -86.4e0,  79.9e0, -74.5e0,
            -1.7e0,      38.7e0,   12.3e0,    9.4e0,   5.4e0,   1.9e0,
            24.8e0,       7.6e0,  -11.7e0,   -6.9e0, -18.1e0,  10.2e0,
             9.4e0,     -11.3e0,   -4.9e0,    5.6e0,   9.8e0,   3.6e0,
            -6.9e0,       5.0e0,  -10.8e0,   -1.3e0,   8.8e0,  -6.7e0,
            -9.2e0,      -2.2e0,   -6.1e0,    1.4e0,  -2.4e0,  -0.2e0,
             3.1e0,       0.3e0,    2.1e0,    3.8e0,  -0.2e0,  -2.1e0,
             2.9e0,      -1.6e0,   -1.9e0,    1.4e0,  -0.3e0,   0.3e0,
            -0.8e0,       0.5e0,    1.8e0,    0.2e0,   1.0e0,   4.0e0,
            -2.2e0,      -0.3e0,    0.2e0,    0.9e0,  -0.4e0,   1.0e0,
            -0.3e0,       0.5e0,   -0.4e0,   -0.4e0,   0.1e0,  -0.5e0,
            -0.1e0,      -0.2e0,   -0.9e0,    0.3e0,   0.3e0,  -0.4e0,
             1.2e0,      -0.4e0,    0.8e0,   -0.3e0,   0.4e0,  -0.1e0,
             0.4e0,      -0.1e0,   -0.2e0};

    double H05[106]={0.0,0.e0,  0.0e0,5078.0e0,  0.0e0,-2594.5e0,-515.4e0,  0.0e0,
          -198.9e0,269.7e0,-524.7e0,  0.0e0,  282.1e0,-225.2e0,145.2e0,
          -305.4e0,  0.0e0,  42.7e0,180.3e0, -123.5e0, -19.6e0,103.9e0,
             0.0e0,-20.3e0,  54.8e0, 63.6e0,  -63.5e0,   0.2e0, 50.9e0,
             0.0e0,-61.1e0, -22.6e0,  6.8e0,   25.4e0,  10.9e0,-26.3e0,
            -4.6e0,  0.0e0,  11.2e0,-20.9e0,    9.8e0, -19.7e0, 16.2e0,
             7.6e0,-12.8e0,  -0.1e0,  0.0e0,  -20.1e0,  12.7e0, 12.7e0,
            -6.7e0, -8.2e0,   8.1e0,  2.9e0,   -7.7e0,   6.0e0,  0.0e0,
             2.2e0,  0.1e0,   4.5e0,  4.8e0,   -6.7e0,  -1.0e0, -3.5e0,
            -0.9e0, -2.3e0,  -7.9e0,  0.0e0,    0.3e0,   1.4e0, -0.8e0,
            -2.3e0,  0.9e0,  -0.6e0, -2.7e0,   -1.1e0,  -1.6e0, -1.9e0,
            -1.4e0,  0.0e0,  -0.6e0,  0.2e0,    2.4e0,  -2.6e0,  0.6e0,
             0.4e0,  0.0e0,   0.0e0,  0.3e0,   -0.9e0,  -0.3e0,  0.9e0,
             0.0e0, -0.8e0,   0.3e0,  1.7e0,   -0.5e0,  -1.1e0,  0.0e0,
             0.6e0,  0.2e0,   0.5e0,  0.4e0,   -0.2e0,  -0.6e0, -0.9e0};

    double G10[106]={0.0,0.00e0,-29496.57e0,-1586.42e0,-2396.06e0,3026.34e0,
           1668.17e0,  1339.85e0,-2326.54e0, 1232.10e0, 633.73e0, 
            912.66e0,   808.97e0,  166.58e0, -356.83e0,  89.40e0,
           -230.87e0,   357.29e0,  200.26e0, -141.05e0,-163.17e0,
             -8.03e0,    72.78e0,   68.69e0,   75.92e0,-141.40e0, 
            -22.83e0,    13.10e0,  -78.09e0,   80.44e0, -75.00e0,
             -4.55e0,    45.24e0,   14.00e0,   10.46e0,   1.64e0,
              4.92e0,    24.41e0,    8.21e0,  -14.50e0,  -5.59e0,
            -19.34e0,    11.61e0,   10.85e0,  -14.05e0,  -3.54e0,
              5.50e0,     9.45e0,    3.45e0,   -5.27e0,   3.13e0,
            -12.38e0,    -0.76e0,    8.43e0,   -8.42e0, -10.08e0,
             -1.94e0,    -6.24e0,    0.89e0,   -1.07e0,  -0.16e0,
              2.45e0,    -0.33e0,    2.13e0,    3.09e0,  -1.03e0,
             -2.80e0,     3.05e0,   -1.48e0,   -2.03e0,   1.65e0,
             -0.51e0,     0.54e0,   -0.79e0,    0.37e0,   1.79e0,
              0.12e0,     0.75e0,    3.75e0,   -2.12e0,  -0.21e0,
              0.30e0,     1.04e0,   -0.63e0,    0.95e0,  -0.11e0,
              0.52e0,    -0.39e0,   -0.37e0,    0.21e0,  -0.77e0,
              0.04e0,    -0.09e0,   -0.89e0,    0.31e0,   0.42e0,
             -0.45e0,     1.08e0,   -0.31e0,    0.78e0,  -0.18e0,
              0.38e0,     0.02e0,    0.42e0,   -0.26e0,  -0.26e0};

    double H10[106]={0.0,0.00e0,  0.00e0,4944.26e0,   0.00e0,-2708.54e0,
          -575.73e0,   0.00e0,-160.40e0, 251.75e0, -537.03e0, 0.00e0,
           286.48e0,-211.03e0, 164.46e0,-309.72e0,    0.00e0,44.58e0,
           189.01e0,-118.06e0,  -0.01e0, 101.04e0,    0.00e0,-20.90e0,
            44.18e0,  61.54e0, -66.26e0,   3.02e0,   55.40e0,  0.00e0,
           -57.80e0, -21.20e0,   6.54e0,  24.96e0,    7.03e0,-27.61e0,
            -3.28e0,   0.00e0,  10.84e0, -20.03e0,   11.83e0,-17.41e0, 
            16.71e0,   6.96e0, -10.74e0,   1.64e0,    0.00e0,-20.54e0,
            11.51e0,  12.75e0,  -7.14e0,  -7.42e0,    7.97e0,  2.14e0,
            -6.08e0,   7.01e0,   0.00e0,   2.73e0,   -0.10e0,  4.71e0,
             4.44e0,  -7.22e0,  -0.96e0,  -3.95e0,   -1.99e0, -1.97e0,
            -8.31e0,   0.00e0,   0.13e0,   1.67e0,   -0.66e0, -1.76e0,
             0.85e0,  -0.39e0,  -2.51e0,  -1.27e0,   -2.11e0, -1.94e0,
            -1.86e0,   0.00e0,  -0.87e0,   0.27e0,    2.13e0, -2.49e0,
             0.49e0,   0.59e0,   0.00e0,   0.13e0,    0.27e0, -0.86e0,
            -0.23e0,   0.87e0,   0.00e0,  -0.87e0,    0.30e0,  1.66e0,
            -0.59e0,  -1.14e0,  -0.07e0,   0.54e0,    0.10e0,  0.49e0,
             0.44e0,  -0.25e0,  -0.53e0,  -0.79e0};

    double G15[106]={0.0, 0.00,-29441.46,-1501.77,-2445.88,3012.20,1676.35,1350.33,
          -2352.26,  1225.85,  581.69,  907.42, 813.68, 120.49,-334.85,
             70.38,  -232.91,  360.14,  192.35,-140.94,-157.40,   4.30,
             69.55,    67.57,   72.79, -129.85, -28.93,  13.14, -70.85,
             81.29,   -75.99,   -6.79,   51.82,  15.07,   9.32,  -2.88,
              6.61,    23.98,    8.89,  -16.78,  -3.16, -20.56,  13.33,
             11.76,   -15.98,   -2.02,    5.33,   8.83,   3.02,  -3.22,
              0.67,   -13.20,   -0.10,    8.68,  -9.06, -10.54,  -2.01,
             -6.26,     0.17,    0.55,   -0.55,   1.70,  -0.67,   2.13,
              2.33,    -1.80,   -3.59,    3.00,  -1.40,  -2.30,   2.08,
             -0.79,     0.58,   -0.70,    0.14,   1.70,  -0.22,   0.44,
              3.49,    -2.09,   -0.16,    0.46,   1.23,  -0.89,   0.85,
              0.10,     0.54,   -0.37,   -0.43,   0.22,  -0.94,  -0.03,
             -0.02,    -0.92,    0.42,    0.63,  -0.42,   0.96,  -0.19,
              0.81,    -0.13,    0.38,    0.08,   0.46,  -0.35,  -0.36};

    double H15[106]={0.0, 0.00,   0.00,  4795.99,    0.00,-2845.41,-642.17, 0.00,
           -115.29, 245.04,  -538.70,    0.00,  283.54,-188.43, 180.95,
           -329.23,   0.00,    46.98,  196.98, -119.14,  15.98, 100.12,
              0.00, -20.61,    33.30,   58.74,  -66.64,   7.35,  62.41,
              0.00, -54.27,   -19.53,    5.59,   24.45,   3.27, -27.50,
             -2.32,   0.00,    10.04,  -18.26,   13.18, -14.60,  16.16,
              5.69,  -9.10,     2.26,    0.00,  -21.77,  10.76,  11.74,
             -6.74,  -6.88,     7.79,    1.04,   -3.89,   8.44,   0.00,
              3.28,  -0.40,     4.55,    4.40,   -7.92,  -0.61,  -4.16,
             -2.85,  -1.12,    -8.72,    0.00,    0.00,   2.11,  -0.60,
             -1.05,   0.76,    -0.20,   -2.12,   -1.44,  -2.57,  -2.01,
             -2.34,   0.00,    -1.08,    0.37,    1.75,  -2.19,   0.27,
              0.72,  -0.09,     0.29,    0.23,   -0.89,  -0.16,   0.72,
              0.00,  -0.88,     0.49,    1.56,   -0.50,  -1.24,  -0.10,
              0.42,  -0.04,     0.48,    0.48,   -0.30,  -0.43,  -0.71};

    double G20[106]={0.0, 0.0,  -29404.8, -1450.9, -2499.6, 2982.0, 1677.0, 1363.2,
          -2381.2,    1236.2,   525.7,   903.0,  809.5,   86.3, -309.4,
             48.0,    -234.3,   363.2,   187.8, -140.7, -151.2,   13.5,
             66.0,      65.5,    72.9,  -121.5,  -36.2,   13.5,  -64.7,
             80.6,     -76.7,    -8.2,    56.5,   15.8,    6.4,   -7.2,
              9.8,      23.7,     9.7,   -17.6,   -0.5,  -21.1,   15.3,
             13.7,     -16.5,    -0.3,     5.0,    8.4,    2.9,   -1.5,
             -1.1,     -13.2,     1.1,     8.8,   -9.3,  -11.9,   -1.9,
             -6.2,      -0.1,     1.7,    -0.9,    0.7,   -0.9,    1.9,
              1.4,      -2.4,    -3.8,     3.0,   -1.4,   -2.5,    2.3,
             -0.9,       0.3,    -0.7,    -0.1,    1.4,   -0.6,    0.2,
              3.1,      -2.0,    -0.1,     0.5,    1.3,   -1.2,    0.7,
              0.3,       0.5,    -0.3,    -0.5,    0.1,   -1.1,   -0.3,
              0.1,      -0.9,     0.5,     0.7,   -0.3,    0.8,    0.0,
              0.8,       0.0,     0.4,     0.1,    0.5,   -0.5,   -0.4};

    double H20[106]={0.0, 0.0,    0.0,    4652.5,     0.0, -2991.6, -734.6,    0.0,
            -82.1,  241.9,    -543.4,     0.0,   281.9, -158.4,  199.7,
           -349.7,    0.0,      47.7,   208.3,  -121.2,   32.3,   98.9,
              0.0,  -19.1,      25.1,    52.8,   -64.5,    8.9,   68.1,
              0.0,  -51.5,     -16.9,     2.2,    23.5,   -2.2,  -27.2,
             -1.8,    0.0,       8.4,   -15.3,    12.8,  -11.7,   14.9,
              3.6,   -6.9,       2.8,     0.0,   -23.4,   11.0,    9.8,
             -5.1,   -6.3,       7.8,     0.4,    -1.4,    9.6,    0.0,
              3.4,   -0.2,       3.6,     4.8,    -8.6,   -0.1,   -4.3,
             -3.4,   -0.1,      -8.8,     0.0,     0.0,    2.5,   -0.6,
             -0.4,    0.6,      -0.2,    -1.7,    -1.6,   -3.0,   -2.0,
             -2.6,    0.0,      -1.2,     0.5,     1.4,   -1.8,    0.1,
              0.8,   -0.2,       0.6,     0.2,    -0.9,    0.0,    0.5,
              0.0,   -0.9,       0.6,     1.4,    -0.4,   -1.3,   -0.1,
              0.3,   -0.1,       0.5,     0.5,    -0.4,   -0.4,   -0.6};

    double DG20[46]={0.0, 0.0,    5.7,     7.4,   -11.0,    -7.0,   -2.1, 2.2,  
              -5.9,    3.1,   -12.0,    -1.2,    -1.6,   -5.9,     5.2,
              -5.1,   -0.3,     0.5,    -0.6,     0.2,    1.3,     0.9,
              -0.5,   -0.3,     0.4,     1.3,    -1.4,    0.0,     0.9,
              -0.1,   -0.2,     0.0,     0.7,     0.1,   -0.5,    -0.8,
               0.8,    0.0,     0.1,    -0.1,     0.4,   -0.1,     0.4,
               0.3,   -0.1,     0.4}; 

    double DH20[46]={0.0, 0.0,    0.0,   -25.9,     0.0,   -30.2,  -22.4, 0.0,
               6.0,   -1.1,     0.5,     0.0,    -0.1,    6.5,     3.6,
              -5.0,    0.0,     0.0,     2.5,    -0.6,    3.0,     0.3,
               0.0,    0.0,    -1.6,    -1.3,     0.8,    0.0,     1.0,
               0.0,    0.6,     0.6,    -0.8,    -0.2,   -1.1,     0.1,
               0.3,    0.0,    -0.2,     0.6,    -0.2,    0.5,    -0.3,
              -0.4,    0.5,     0.0};

      IY=IYEAR;
/*
C  WE ARE RESTRICTED BY THE INTERVAL 1965-2020, FOR WHICH EITHER THE IGRF/DGRF COEFFICIENTS OR SECULAR VELOCITIES
c    ARE KNOWN; IF IYEAR IS OUTSIDE THIS INTERVAL, THEN THE SUBROUTINE USES THE
C      NEAREST LIMITING VALUE AND PRINTS A WARNING:

      IF(IY.LT.1965) THEN
       IY=1965
       WRITE (*,10) IYEAR,IY
      ENDIF

      IF(IY.GT.2020) THEN
       IY=2020
       WRITE (*,10) IYEAR,IY
      ENDIF

C  CALCULATE THE ARRAY REC, CONTAINING COEFFICIENTS FOR THE RECURSION RELATIONS,
C  USED IN THE IGRF SUBROUTINE FOR CALCULATING THE ASSOCIATE LEGENDRE POLYNOMIALS
C  AND THEIR DERIVATIVES:
*/
    for (N=1; N<=14; N++) {
        N22=2*N-1;
        N22=N22*(N22-2);
        for (M=1; M<=N; M++) {
            MN=N*(N-1)/2+M;
        }
        REC[MN]=double((N-M)*(N+M-2))/double(N22);
    }

//  INTERPOLATE BETWEEEN 1965 - 1970:
    if (IY<1970) {
        if (IY<1965) IY=1965;
        F2=(double(IY)+double(IDAY-1)/365.25e0-1965)/5.e0;
        F1=1.e0-F2;
        for (N=1; N<=105; N++) {
            G[N]=G65[N]*F1+G70[N]*F2;
            H[N]=H65[N]*F1+H70[N]*F2;
        }
    }

//  INTERPOLATE BETWEEN 1970 - 1975:
    else if (IY>=1970 && IY<1975) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-1970)/5.e0;
        F1=1.e0-F2;
        for (N=1; N<=105; N++) {
            G[N]=G70[N]*F1+G75[N]*F2;
            H[N]=H70[N]*F1+H75[N]*F2;
        }
    }

//  INTERPOLATE BETWEEN 1975 - 1980:

    else if (IY>=1975 && IY<1980) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-1975)/5.e0;
        F1=1.e0-F2;
        for(N=1; N<=105; N++) {
           G[N]=G75[N]*F1+G80[N]*F2;
           H[N]=H75[N]*F1+H80[N]*F2;
        }
    }

//  INTERPOLATE BETWEEN 1980 - 1985:
    else if (IY>=1980 && IY<1985) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-1980)/5.e0;
        F1=1.e0-F2;
        for (N=1; N<=105; N++) {
            G[N]=G80[N]*F1+G85[N]*F2;
            H[N]=H80[N]*F1+H85[N]*F2;
        }
    }

    //INTERPOLATE BETWEEN 1985 - 1990:
    else if (IY>=1985 && IY<1990) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-1985)/5.e0;
        F1=1.e0-F2;
        for (N=1; N<=105; N++) {
            G[N]=G85[N]*F1+G90[N]*F2;
            H[N]=H85[N]*F1+H90[N]*F2;
        }
    }

//  INTERPOLATE BETWEEN 1990 - 1995:

    else if (IY>=1990 && IY<1995) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-1990)/5.e0;
        F1=1.e0-F2;
        for (N=1; N<=105; N++) {
            G[N]=G90[N]*F1+G95[N]*F2;
            H[N]=H90[N]*F1+H95[N]*F2;
        }
    }

//  INTERPOLATE BETWEEN 1995 - 2000:

    else if (IY>=1995 && IY<2000) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-1995)/5.e0;
        F1=1.e0-F2;
        for (N=1; N<=105; N++) {
            G[N]=G95[N]*F1+G00[N]*F2;
            H[N]=H95[N]*F1+H00[N]*F2;
        }
    }

    //INTERPOLATE BETWEEN 2000 - 2005:

    else if (IY>=2000 && IY<2005) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-2000)/5.e0;
        F1=1.e0-F2;
        for (N=1; N<=105; N++) {
            G[N]=G00[N]*F1+G05[N]*F2;
            H[N]=H00[N]*F1+H05[N]*F2;
        }
    }

    //INTERPOLATE BETWEEN 2005 - 2010:

    else if (IY>=2005 && IY<2010) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-2005)/5.e0;
        F1=1.-F2;
        for (N=1; N<=105; N++) {
            G[N]=G05[N]*F1+G10[N]*F2;
            H[N]=H05[N]*F1+H10[N]*F2;
        }
    }

    //INTERPOLATE BETWEEN 2010 - 2015:

    else if (IY>=2010 && IY<2015) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-2010)/5.e0;
        F1=1.-F2;
        for (N=1; N<=105; N++) {
            G[N]=G10[N]*F1+G15[N]*F2;
            H[N]=H10[N]*F1+H15[N]*F2;
        }
    }

    else if (IY>=2015 && IY <2020) {
        F2=(double(IY)+double(IDAY-1)/365.25e0-2015)/5.e0;
        F1=1.-F2;
        for (N=1; N<=105; N++) {
            G[N]=G15[N]*F1+G20[N]*F2;
            H[N]=H15[N]*F1+H20[N]*F2;
        }
    }

    //EXTRAPOLATE BEYOND 2020:
    else {
        if (IY>2020) IY=2020;
        DT=double(IY)+double(IDAY-1)/365.25e0-2020.e0;
        for (N=1; N<=105; N++) {
            G[N]=G20[N];
            H[N]=H20[N];
            if (N>45) continue;
            G[N]=G[N]+DG20[N]*DT;
            H[N]=H[N]+DH20[N]*DT;
        }
    }

/*   COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
C   THEM BY SCHMIDT NORMALIZATION FACTORS:*/

    S=1.e0;
    for (N=2; N<=14; N++) {
        MN=N*(N-1)/2+1;
        S=S*double(2*N-3)/double(N-1);
        G[MN]=G[MN]*S;
        H[MN]=H[MN]*S;
        P=S;
        for (M=2; M<=N; M++) {
            AA1=1.e0;
            if (M==2) AA1=2.e0;
            P=P*sqrt(AA1*double(N-M+1)/double(N+M-2));
            MNN=MN+M-1;
            G[MNN]=G[MNN]*P;
            H[MNN]=H[MNN]*P;
        }
    }

    G_10=-G[2];
    G_11= G[3];
    H_11= H[3];

/*  NOW CALCULATE GEO COMPONENTS OF THE UNIT VECTOR EzMAG, PARALLEL TO GEODIPOLE AXIS:
C   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
         ST0 * CL0                ST0 * SL0                CT0 */

    SQ=G_11*G_11+H_11*H_11;
    SQQ=sqrt(SQ);
    SQR=sqrt(G_10*G_10+SQ);
    AAP[2]=-H_11/SQQ;       //SL0
    AAP[3]=-G_11/SQQ;       //CL0
    AAP[0]=SQQ/SQR;         //ST0
    AAP[1]=G_10/SQR;        //CT0
    AAP[5]=AAP[0]*AAP[3];   //STCL
    AAP[7]=AAP[0]*AAP[2];   //STSL
    AAP[6]=AAP[1]*AAP[2];   //CTSL
    AAP[4]=AAP[1]*AAP[3];   //CTCL

/*  NOW CALCULATE GEI COMPONENTS (S1,S2,S3) OF THE UNIT VECTOR S = EX_GSE
    POINTING FROM THE EARTH'S CENTER TO SUN */

    SUN_08(IY,IDAY,IHOUR,IMIN,ISEC,GST,SLONG,SRASN,SDEC);

    S1=cos(SRASN)*cos(SDEC);
    S2=sin(SRASN)*cos(SDEC);
    S3=sin(SDEC);

/*  NOW CALCULATE GEI COMPONENTS (DZ1,DZ2,DZ3) OF THE UNIT VECTOR EZGSE
C  POINTING NORTHWARD AND ORTHOGONAL TO THE ECLIPTIC PLANE, AS
C  (0,-SIN(OBLIQ),COS(OBLIQ)). FOR THE EPOCH 1978, OBLIQ = 23.44214 DEGS.
C  HERE WE USE A MORE ACCURATE TIME-DEPENDENT VALUE, DETERMINED AS:
*/
    DJ=double(365*(IY-1900)+(IY-1901)/4 +IDAY)
       -0.5+double(IHOUR*3600+IMIN*60+ISEC)/86400.e0;
    T=DJ/36525.e0;
    OBLIQ=(23.45229e0-0.0130125e0*T)/57.2957795e0;
    DZ1=0.e0;
    DZ2=-sin(OBLIQ);
    DZ3=cos(OBLIQ);

/*
C  NOW WE OBTAIN GEI COMPONENTS OF THE UNIT VECTOR EYGSE=(DY1,DY2,DY3),
C  COMPLETING THE RIGHT-HANDED SYSTEM. THEY CAN BE FOUND FROM THE VECTOR
C  PRODUCT EZGSE x EXGSE = (DZ1,DZ2,DZ3) x (S1,S2,S3):
*/
    DY1=DZ2*S3-DZ3*S2;
    DY2=DZ3*S1-DZ1*S3;
    DY3=DZ1*S2-DZ2*S1;
/*
C  NOW LET'S CALCULATE GEI COMPONENTS OF THE UNIT VECTOR X = EXGSW, DIRECTED ANTIPARALLEL
C  TO THE OBSERVED SOLAR WIND FLOW. FIRST, CALCULATE ITS COMPONENTS IN GSE:
*/
    V=sqrt(VGSEX*VGSEX+VGSEY*VGSEY+VGSEZ*VGSEZ);
    DX1=-VGSEX/V;
    DX2=-VGSEY/V;
    DX3=-VGSEZ/V;
/*
C  THEN IN GEI:
*/
    X1=DX1*S1+DX2*DY1+DX3*DZ1;
    X2=DX1*S2+DX2*DY2+DX3*DZ2;
    X3=DX1*S3+DX2*DY3+DX3*DZ3;
/*
C  NOW CALCULATE GEI COMPONENTS (DIP1,DIP2,DIP3) OF THE UNIT VECTOR DIP = EZ_SM = EZ_MAG,
C   ALIGNED WITH THE GEODIPOLE AND POINTING NORTHWARD FROM ECLIPTIC PLANE:
*/
    AAP[13]=cos(GST);       //CGST
    AAP[14]=sin(GST);       //SGST

    DIP1=AAP[5]*AAP[13]-AAP[7]*AAP[14];
    DIP2=AAP[5]*AAP[14]+AAP[7]*AAP[13];
    DIP3=AAP[1];
/*
C  THIS ALLOWS US TO CALCULATE GEI COMPONENTS OF THE UNIT VECTOR Y = EYGSW
C   BY TAKING THE VECTOR PRODUCT DIP x X AND NORMALIZING IT TO UNIT LENGTH:
*/
    Y1=DIP2*X3-DIP3*X2;
    Y2=DIP3*X1-DIP1*X3;
    Y3=DIP1*X2-DIP2*X1;
    Y=sqrt(Y1*Y1+Y2*Y2+Y3*Y3);
    Y1=Y1/Y;
    Y2=Y2/Y;
    Y3=Y3/Y;
/*
C   AND GEI COMPONENTS OF THE UNIT VECTOR Z = EZGSW = EXGSW x EYGSW = X x Y:
*/
    Z1=X2*Y3-X3*Y2;
    Z2=X3*Y1-X1*Y3;
    Z3=X1*Y2-X2*Y1;
/*
C   ELEMENTS OF THE MATRIX GSE TO GSW ARE THE SCALAR PRODUCTS:
C
C  E11=(EXGSE,EXGSW)  E12=(EXGSE,EYGSW)  E13=(EXGSE,EZGSW)
C  E21=(EYGSE,EXGSW)  E22=(EYGSE,EYGSW)  E23=(EYGSE,EZGSW)
C  E31=(EZGSE,EXGSW)  E32=(EZGSE,EYGSW)  E33=(EZGSE,EZGSW)
*/
    AAP[25]= S1*X1 +S2*X2 +S3*X3;   //E11
    AAP[28]= S1*Y1 +S2*Y2 +S3*Y3;   //E12
    AAP[31]= S1*Z1 +S2*Z2 +S3*Z3;   //E13
    AAP[26]=DY1*X1+DY2*X2+DY3*X3;   //E21
    AAP[29]=DY1*Y1+DY2*Y2+DY3*Y3;   //E22
    AAP[32]=DY1*Z1+DY2*Z2+DY3*Z3;   //E23
    AAP[27]=DZ1*X1+DZ2*X2+DZ3*X3;   //E31
    AAP[30]=DZ1*Y1+DZ2*Y2+DZ3*Y3;   //E32
    AAP[33]=DZ1*Z1+DZ2*Z2+DZ3*Z3;   //E33

/*
C   GEODIPOLE TILT ANGLE IN THE GSW SYSTEM: PSI=ARCSIN(DIP,EXGSW)
*/
    AAP[10]=DIP1*X1+DIP2*X2+DIP3*X3;     //SPS
    AAP[11]=sqrt(1.e0-AAP[10]*AAP[10]);  //CPS
    AAP[15]=asin(AAP[10]);               //PSI
/*
C   ELEMENTS OF THE MATRIX GEO TO GSW ARE THE SCALAR PRODUCTS:

C   A11=(EXGEO,EXGSW), A12=(EYGEO,EXGSW), A13=(EZGEO,EXGSW),
C   A21=(EXGEO,EYGSW), A22=(EYGEO,EYGSW), A23=(EZGEO,EYGSW),
C   A31=(EXGEO,EZGSW), A32=(EYGEO,EZGSW), A33=(EZGEO,EZGSW),

C   ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:
C
C  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
C  EXGSW=(X1,X2,X3),  EYGSW=(Y1,Y2,Y3),   EZGSW=(Z1,Z2,Z3)
C                                                    AND  THEREFORE:
*/
    AAP[16]=X1*AAP[13]+X2*AAP[14];     //A11
    AAP[19]=-X1*AAP[14]+X2*AAP[13];    //A12
    AAP[22]=X3;                        //A13
    AAP[17]=Y1*AAP[13]+Y2*AAP[14];     //A21
    AAP[20]=-Y1*AAP[14]+Y2*AAP[13];    //A22
    AAP[23]=Y3;                        //A23
    AAP[18]=Z1*AAP[13]+Z2*AAP[14];     //A31
    AAP[21]=-Z1*AAP[14]+Z2*AAP[13];    //A32
    AAP[24]=Z3;                        //A33
/*
C  NOW CALCULATE ELEMENTS OF THE MATRIX MAG TO SM (ONE ROTATION ABOUT THE GEODIPOLE AXIS);
C   THEY ARE FOUND AS THE SCALAR PRODUCTS: CFI=GM22=(EYSM,EYMAG)=(EYGSW,EYMAG),
C                                          SFI=GM23=(EYSM,EXMAG)=(EYGSW,EXMAG),
C    DERIVED AS FOLLOWS:
C
C IN GEO, THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
C  AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THEIR COMPONENTS ARE:
C  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
C            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
C            -ST0
C  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
C            -SL0*SIN(GST)+CL0*COS(GST)
C             0
C  NOW, NOTE THAT GEI COMPONENTS OF EYSM=EYGSW WERE FOUND ABOVE AS Y1, Y2, AND Y3,
C  AND WE ONLY HAVE TO COMBINE THESE QUANTITIES INTO SCALAR PRODUCTS:
*/
    EXMAGX=AAP[1]*(AAP[3]*AAP[13]-AAP[2]*AAP[14]);
    EXMAGY=AAP[1]*(AAP[3]*AAP[14]+AAP[2]*AAP[13]);
    EXMAGZ=-AAP[0];
    EYMAGX=-(AAP[2]*AAP[13]+AAP[3]*AAP[14]);
    EYMAGY=-(AAP[2]*AAP[14]-AAP[3]*AAP[13]);
    AAP[9]=Y1*EYMAGX+Y2*EYMAGY;                  //CFI
    AAP[8]=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ;        //SFI

    return;
}

/*========================================================================================
C
C    CONVERTS GEOGRAPHIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICE VERSA.
C
C                    J>0                       J<0
C-----INPUT:  J,XGEO,YGEO,ZGEO           J,XMAG,YMAG,ZMAG
C-----OUTPUT:    XMAG,YMAG,ZMAG           XGEO,YGEO,ZGEO
C
C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE GEOMAG_08 IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR AND/OR IDAY HAVE BEEN CHANGED
C
C  NO INFORMATION IS REQUIRED HERE ON THE SOLAR WIND VELOCITY, SO ONE
C  CAN SET VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0 IN RECALC_08.
C
C   LAST MOFIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C   AUTHOR:  N. A. TSYGANENKO
*/
inline void GEOMAG_08 (double &XGEO,double &YGEO,double &ZGEO,double &XMAG,double &YMAG,
                double &ZMAG,int J)
{
    if(J > 0) {
       XMAG=XGEO*AAP[4]+YGEO*AAP[6]-ZGEO*AAP[0];
       YMAG=YGEO*AAP[3]-XGEO*AAP[2];
       ZMAG=XGEO*AAP[5]+YGEO*AAP[7]+ZGEO*AAP[1];
    }
    else {
       XGEO=XMAG*AAP[4]-YMAG*AAP[2]+ZMAG*AAP[5];
       YGEO=XMAG*AAP[6]+YMAG*AAP[3]+ZMAG*AAP[7];
       ZGEO=ZMAG*AAP[1]-XMAG*AAP[0];
    }

    return;
}

/***********************************************************************
 function dayno
   get number of days in a year from Jan 1

   Author Jiannan Tu
   Date March 15, 2007, coverted to C, October 22, 2012
***********************************************************************/
inline int dayno(int yr, int mn, int dd)
{
    int iday=1, i;

    const int id[11]={31,59,90,120,151,181,212,243,273,304,334};

    if(mn ==1) {
        iday=dd;
        return iday;
    }
    else {
        for (i = 1; i < 11; i++) if (mn == i+2) iday=id[i]+dd;
    }

    if (mn > 2 && (yr % 4) == 0) iday=iday+1;

    return iday;
}

/**********************************************************************
 * find related point on the earth shell of polar (dipole)
 input: x, y, z 
    output: x_earth, y_earth, z_earth
    require r_earth < r_top, otherwise it has unknown numerical error
    set only north hemisphere
***********************************************************************/
inline void DipoleRelatedOnEarthPolar( double x, double y, double z, double &x_earth,
                                        double &y_earth, double &z_earth, double height)
{
    //double PI = acos(-1.0);
    double r_earth = radius + height;
    //
    double r_top = sqrt(x * x + y * y + z * z);
    double theta_top = acos(z / r_top);
    double phi_top;

    phi_top=atan2(y, x);
    /*if (x != 0.0)
    {
        phi_top = atan(y / x);
        if (x < 0)
        {
            phi_top = phi_top + PI;
        }
    }
    else if (x == 0.0 && y > 0.0)
    {
        phi_top = PI / 2.0;
    }
    else if (x == 0.0 && y < 0.0)
    {
        phi_top = PI / 2.0 * (-1.0);
    }
    else if (x == 0.0 && y == 0.0)
    {
        phi_top = 0.0; // special points
    }*/
    // find related point on the earth shell 
    double theta_earth = asin(sin(theta_top) * sqrt(r_earth / r_top));
    // 
    if( z < 0.0) theta_earth = PI - theta_earth;  // 
    double phi_earth = phi_top;
    x_earth = r_earth * sin(theta_earth) * cos(phi_earth);
    y_earth = r_earth * sin(theta_earth) * sin(phi_earth);
    z_earth = r_earth * cos(theta_earth);
}
#endif