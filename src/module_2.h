#ifndef _MODULE_H_2_
#define _MODULE_H_2_
#include<iostream>
#include "parameters.h"
#include <memory>
#include "particles.h"
#include <vector>
#include <cmath>
#include <limits>
#include <bitset>
#include "module_base.h"
#include "vector3.h"
#include <algorithm> 

using std::max;


//************************************************************************
// Calculate postcollision velocity ( Nanbu and Yonemura, JCP, 145, 639, 1998 and 1997)
// test model
inline void test_coll_vel_chng(int a, int b,
                                double dt, double nb,
                                double LnA,
                                GridsPoints***** ptrArry,
                                Particles& pa, Particles& pb)
{
    const double pi2 = 6.283185307; // 2pi
    // 0 - H; 1 - He; 2 - O
    const double mc[3][3] ={{0.50,0.20,0.058823529},
                            {0.80,0.50,0.20},
                            {0.94117647,0.80,0.50}};
    /*double mi0_a, mi0_b;
    if( a == 0)
    {   
        mi0_a = mi0_H;
    }else if( a == 1)
    {
        mi0_a = mi0_He;
    } else if( a == 2)
    {
        mi0_a = mi0_O;
    }
    if( b == 0)
    {   
        mi0_b = mi0_H;
    }else if( b == 1)
    {
        mi0_b = mi0_He;
    } else if( b == 2)
    {
        mi0_b = mi0_O;
    }*/
    //
    Vector3 va = pa.VelParticles();
    Vector3 vb = pb.VelParticles();
    //
    Vector3 gg = va.MinusProduct(vb); // gg
    Vector3 gx = Vector3(gg.x(), 0.0, 0.0); // g_parallel
    Vector3 gp = gg.MinusProduct(gx); // g_perpendicular
    double gab = gg.norm();
    double acosx = Cosx( a, b, dt, nb, gab, LnA);
    double sinx = sqrt( 1.0 - acosx*acosx);
    acosx = 1.0 - acosx;
    //
    double epsilon = dRand() * pi2;
    double cosep = cos(epsilon);
    double sinep = sin(epsilon);
    //
    double hhx = gp.norm() * cosep;
    double hhy = -(gg.y() * gg.x() * cosep + gab * gg.z() * sinep)/gp.norm();
    double hhz = -(gg.z() * gg.x() * cosep - gab * gg.y() * sinep)/gp.norm();
    Vector3 hh = Vector3(hhx, hhy, hhz);
    //
    Vector3 gh = gg.ScaleProduct(acosx).PlusProduct(hh.ScaleProduct(sinx));
    Vector3 va_new = Vector3( va.MinusProduct(gh.ScaleProduct( mc[b][a])));
    Vector3 vb_new = Vector3( vb.PlusProduct(gh.ScaleProduct(mc[a][b])));
    //
    pa.SetVelocity(va_new);
    pb.SetVelocity(vb_new);
    double energy = pa.VelParticles().norm2() + pb.VelParticles().norm2();
    std::cout << energy * mi0_H * 0.5 << " " << pa.VelParticles().norm2()* mi0_H * 0.5 << " " << pb.VelParticles().norm2()* mi0_H * 0.5;
    std::cout << " v " << pa.VelParticles().x() + pb.VelParticles().x() << " " << pa.VelParticles().y() + pb.VelParticles().y() << " " << pa.VelParticles().z() + pb.VelParticles().z();
    

}
//************************************************************************
// Calculate postcollision velocity ( Nanbu and Yonemura, JCP, 145, 639, 1998 and 1997)
// assume an average B field for each cell
inline void coll_vel_chng(int a, int b,
                          double dt, double nb,
                          double LnA,
                          Vector3& bVector3,
                          Particles& pa, Particles& pb)
{
    const double pi2 = 6.283185307; // 2pi
    // 0 - H; 1 - He; 2 - O
    const double mc[3][3] ={{0.50,0.20,0.058823529},
                            {0.80,0.50,0.20},
                            {0.94117647,0.80,0.50}};
    Vector3 b_unit = bVector3.NormalizedVector();
    // 
    double mi0_a, mi0_b;
    if( a == 0)
    {   
        mi0_a = mi0_H;
    }else if( a == 1)
    {
        mi0_a = mi0_He;
    } else if( a == 2)
    {
        mi0_a = mi0_O;
    }
    Vector3 va = pa.VelCollParticles( bVector3, mi0_a);
    if( b == 0)
    {   
        mi0_b = mi0_H;
    }else if( b == 1)
    {
        mi0_b = mi0_He;
    } else if( b == 2)
    {
        mi0_b = mi0_O;
    }
    Vector3 vb = pb.VelCollParticles( bVector3, mi0_b);
    //
    double proba = pb.WeightNi() / max( pa.WeightNi(), pb.WeightNi());
    double probb = pa.WeightNi() / max( pa.WeightNi(), pb.WeightNi());
    //
    Vector3 velGuiding_a = pa.VelParticles();
    Vector3 vel_para_a = b_unit.ScaleProduct( velGuiding_a.DotProduct(b_unit));
    Vector3 vel_drift_a = velGuiding_a.MinusProduct(vel_para_a);
    Vector3 velGuiding_b = pb.VelParticles();
    Vector3 vel_para_b = b_unit.ScaleProduct( velGuiding_b.DotProduct(b_unit));
    Vector3 vel_drift_b = velGuiding_b.MinusProduct(vel_para_b);
    //
    std::cout << " driftab " << vel_drift_a.norm() << " " << vel_drift_b.norm() << " ";
    //
    Vector3 gg = va.MinusProduct(vb); // gg
    Vector3 gx = Vector3(gg.x(), 0.0, 0.0); // g_parallel
    Vector3 gp = gg.MinusProduct(gx); // g_perpendicular
    double gab = gg.norm();
    double acosx = Cosx( a, b, dt, nb, gab, LnA);
    double sinx = sqrt( 1.0 - acosx*acosx);
    acosx = 1.0 - acosx;
    //
    double epsilon = dRand() * pi2;
    double cosep = cos(epsilon);
    double sinep = sin(epsilon);
    //
    double hhx = gp.norm() * cosep;
    double hhy = -(gg.y() * gg.x() * cosep + gab * gg.z() * sinep)/gp.norm();
    double hhz = -(gg.z() * gg.x() * cosep - gab * gg.y() * sinep)/gp.norm();
    Vector3 hh = Vector3(hhx, hhy, hhz);
    //
    Vector3 gh = gg.ScaleProduct(acosx).PlusProduct(hh.ScaleProduct(sinx));
    Vector3 va_new = Vector3( va.MinusProduct(gh.ScaleProduct( mc[b][a])));
    Vector3 vb_new = Vector3( vb.PlusProduct(gh.ScaleProduct(mc[a][b])));
    //
    std::cout << "+ " << (pa.VelCollParticles(bVector3, mi0_H).norm2() + pb.VelCollParticles(bVector3, mi0_H).norm2()) * mi0_H * 0.5 << " ";
    if( dRand() < proba)
    {
        pa.SetMuVelocity(va_new, velGuiding_a, bVector3, mi0_a);
    }
    if( dRand() < probb)
    {
        pb.SetMuVelocity(vb_new, velGuiding_b, bVector3, mi0_b);
    }
    std::cout << "- " << (pa.VelCollParticles(bVector3, mi0_H).norm2() + pb.VelCollParticles(bVector3, mi0_H).norm2()) * mi0_H * 0.5 << " ";

}
//************************************************************************
// Calculate postcollision velocity ( Nanbu and Yonemura, JCP, 145, 639, 1998 and 1997)
// not high efficiency
inline void coll_vel_chng(int a, int b,
                          double dt, double nb,
                          double LnA,
                          GridsPoints***** ptrArray_in,
                          Particles& pa, Particles& pb)
{
    const double pi2 = 6.283185307F; // 2pi
    // 0 - H; 1 - He; 2 - O
    //const double mc[3][3] ={{0.50F, 0.20F, 0.058823529F},
    //                        {0.80F, 0.50F, 0.20F},
    //                        {0.94117647F, 0.80F, 0.50F}};
    // calculate b3 at particles' position
    struct structg strg_in1 = {0,0,0,0,0,0,0, 0.0F, 0.0F, 0.0F};
    struct structg strg_in2 = {0,0,0,0,0,0,0, 0.0F, 0.0F, 0.0F};
    Vector3 tempb;
    Vector3 b3_pa, b3_pb;
    uint_64 posUint1, posUint2;
    // same with
    // inline Vector3 B3atParticles( GridsPoints***** ptrArray_in)
    posUint1 = pa.PosUint();
    strg_in1.face = posUint1 >> 61;
    for( int i = 0; i < fieldsGridsLevel; i++) 
    {
        strg_in1.ig = (strg_in1.ig << 1) + ((posUint1 >> (60   - i*3)) & 1);
        strg_in1.jg = (strg_in1.jg << 1) + ((posUint1 >> (60-1 - i*3)) & 1);
        strg_in1.kg = (strg_in1.kg << 1) + ((posUint1 >> (60-2 - i*3)) & 1);
    }
    for( int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg_in1.iw = (strg_in1.iw << 1) + ((posUint1 >> (60   - i*3)) & 1);
        strg_in1.jw = (strg_in1.jw << 1) + ((posUint1 >> (60-1 - i*3)) & 1);
        strg_in1.kw = (strg_in1.kw << 1) + ((posUint1 >> (60-2 - i*3)) & 1);
    } 
    Vector3 tempb1=ptrArray_in[strg_in1.face][strg_in1.ig+1][strg_in1.jg+1][strg_in1.kg]->B3();
    Vector3 tempb2=ptrArray_in[strg_in1.face][strg_in1.ig+2][strg_in1.jg+1][strg_in1.kg]->B3();
    Vector3 tempb3=ptrArray_in[strg_in1.face][strg_in1.ig+2][strg_in1.jg+2][strg_in1.kg]->B3();
    Vector3 tempb4=ptrArray_in[strg_in1.face][strg_in1.ig+1][strg_in1.jg+2][strg_in1.kg]->B3();
    Vector3 tempb5=ptrArray_in[strg_in1.face][strg_in1.ig+1][strg_in1.jg+1][strg_in1.kg+1]->B3();
    Vector3 tempb6=ptrArray_in[strg_in1.face][strg_in1.ig+2][strg_in1.jg+1][strg_in1.kg+1]->B3();
    Vector3 tempb7=ptrArray_in[strg_in1.face][strg_in1.ig+2][strg_in1.jg+2][strg_in1.kg+1]->B3();
    Vector3 tempb8=ptrArray_in[strg_in1.face][strg_in1.ig+1][strg_in1.jg+2][strg_in1.kg+1]->B3();
    double w1 = 1.0F- 1.0F*(strg_in1.iw +0.5F) * (strg_in1.jw +0.5F) * (strg_in1.kw +0.5F) / cellSize1/ cellSize1/ cellSize1;
    double w2 = 1.0F- 1.0F*(cellSize1- strg_in1.iw -0.5F)* (strg_in1.jw +0.5F) * (strg_in1.kw +0.5F) / cellSize1/ cellSize1/ cellSize1;
    double w3 = 1.0F- 1.0F*(cellSize1- strg_in1.iw -0.5F) * (cellSize1- strg_in1.jw -0.5F) * (strg_in1.kw +0.5F) / cellSize1/ cellSize1/ cellSize1;
    double w4 = 1.0F- 1.0F*(strg_in1.iw +0.5F) * (cellSize1- strg_in1.jw -0.5F)* (strg_in1.kw +0.5F)/ cellSize1/ cellSize1/ cellSize1;
    double w5 = 1.0F- 1.0F*(strg_in1.iw +0.5F) * (strg_in1.jw +0.5F) * (cellSize1- strg_in1.kw -0.5F) / cellSize1/ cellSize1/ cellSize1;
    double w6 = 1.0F- 1.0F*(cellSize1- strg_in1.iw -0.5F)* (strg_in1.jw +0.5F) * (cellSize1- strg_in1.kw -0.5F) / cellSize1/ cellSize1/ cellSize1;
    double w7 = 1.0F- 1.0F*(cellSize1- strg_in1.iw -0.5F) * (cellSize1- strg_in1.jw -0.5F) * (cellSize1- strg_in1.kw -0.5F) / cellSize1/ cellSize1/ cellSize1;
    double w8 = 1.0F- 1.0F*(strg_in1.iw +0.5F) * (cellSize1- strg_in1.jw -0.5F)* (cellSize1- strg_in1.kw -0.5F) / cellSize1/ cellSize1/ cellSize1; 
    tempb.Setx(tempb1.x()*w1 + tempb2.x()*w2 + tempb3.x()*w3 + tempb4.x()*w4 
                + tempb5.x()*w5 + tempb6.x()*w6 + tempb7.x()*w7 + tempb8.x()*w8);            
    tempb.Sety(tempb1.y()*w1 + tempb2.y()*w2 + tempb3.y()*w3 + tempb4.y()*w4 
                + tempb5.y()*w5 + tempb6.y()*w6 + tempb7.y()*w7 + tempb8.y()*w8);
    tempb.Setz(tempb1.z()*w1 + tempb2.z()*w2 + tempb3.z()*w3 + tempb4.z()*w4 
                + tempb5.z()*w5 + tempb6.z()*w6 + tempb7.z()*w7 + tempb8.z()*w8);
    b3_pa = Vector3(tempb);
    //
    posUint2 = pb.PosUint();
    strg_in2.face = posUint2 >> 61;
    for( int i = 0; i < fieldsGridsLevel; i++) 
    {
        strg_in2.ig = (strg_in2.ig << 1) + ((posUint2 >> (60   - i*3)) & 1);
        strg_in2.jg = (strg_in2.jg << 1) + ((posUint2 >> (60-1 - i*3)) & 1);
        strg_in2.kg = (strg_in2.kg << 1) + ((posUint2 >> (60-2 - i*3)) & 1);
    }
    for( int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg_in2.iw = (strg_in2.iw << 1) + ((posUint2 >> (60   - i*3)) & 1);
        strg_in2.jw = (strg_in2.jw << 1) + ((posUint2 >> (60-1 - i*3)) & 1);
        strg_in2.kw = (strg_in2.kw << 1) + ((posUint2 >> (60-2 - i*3)) & 1);
    } 
    tempb1=ptrArray_in[strg_in2.face][strg_in2.ig+1][strg_in2.jg+1][strg_in2.kg]->B3();
    tempb2=ptrArray_in[strg_in2.face][strg_in2.ig+2][strg_in2.jg+1][strg_in2.kg]->B3();
    tempb3=ptrArray_in[strg_in2.face][strg_in2.ig+2][strg_in2.jg+2][strg_in2.kg]->B3();
    tempb4=ptrArray_in[strg_in2.face][strg_in2.ig+1][strg_in2.jg+2][strg_in2.kg]->B3();
    tempb5=ptrArray_in[strg_in2.face][strg_in2.ig+1][strg_in2.jg+1][strg_in2.kg+1]->B3();
    tempb6=ptrArray_in[strg_in2.face][strg_in2.ig+2][strg_in2.jg+1][strg_in2.kg+1]->B3();
    tempb7=ptrArray_in[strg_in2.face][strg_in2.ig+2][strg_in2.jg+2][strg_in2.kg+1]->B3();
    tempb8=ptrArray_in[strg_in2.face][strg_in2.ig+1][strg_in2.jg+2][strg_in2.kg+1]->B3();
    w1 = 1.0F- 1.0F*(strg_in2.iw +0.5F) * (strg_in2.jw +0.5F) * (strg_in2.kw +0.5F) / cellSize1/ cellSize1/ cellSize1;
    w2 = 1.0F- 1.0F*(cellSize1- strg_in2.iw -0.5F)* (strg_in2.jw +0.5F) * (strg_in2.kw +0.5F) / cellSize1/ cellSize1/ cellSize1;
    w3 = 1.0F- 1.0F*(cellSize1- strg_in2.iw -0.5F) * (cellSize1- strg_in2.jw -0.5F) * (strg_in2.kw +0.5F) / cellSize1/ cellSize1/ cellSize1;
    w4 = 1.0F- 1.0F*(strg_in2.iw +0.5F) * (cellSize1- strg_in2.jw -0.5F)* (strg_in2.kw +0.5F)/ cellSize1/ cellSize1/ cellSize1;
    w5 = 1.0F- 1.0F*(strg_in2.iw +0.5F) * (strg_in2.jw +0.5F) * (cellSize1- strg_in2.kw -0.5F) / cellSize1/ cellSize1/ cellSize1;
    w6 = 1.0F- 1.0F*(cellSize1- strg_in2.iw -0.5F)* (strg_in2.jw +0.5F) * (cellSize1- strg_in2.kw -0.5F) / cellSize1/ cellSize1/ cellSize1;
    w7 = 1.0F- 1.0F*(cellSize1- strg_in2.iw -0.5F) * (cellSize1- strg_in2.jw -0.5F) * (cellSize1- strg_in2.kw -0.5F) / cellSize1/ cellSize1/ cellSize1;
    w8 = 1.0F- 1.0F*(strg_in2.iw +0.5F) * (cellSize1- strg_in2.jw -0.5F)* (cellSize1- strg_in2.kw -0.5F) / cellSize1/ cellSize1/ cellSize1; 
    tempb.Setx(tempb1.x()*w1 + tempb2.x()*w2 + tempb3.x()*w3 + tempb4.x()*w4 
                + tempb5.x()*w5 + tempb6.x()*w6 + tempb7.x()*w7 + tempb8.x()*w8);            
    tempb.Sety(tempb1.y()*w1 + tempb2.y()*w2 + tempb3.y()*w3 + tempb4.y()*w4 
                + tempb5.y()*w5 + tempb6.y()*w6 + tempb7.y()*w7 + tempb8.y()*w8);
    tempb.Setz(tempb1.z()*w1 + tempb2.z()*w2 + tempb3.z()*w3 + tempb4.z()*w4 
                + tempb5.z()*w5 + tempb6.z()*w6 + tempb7.z()*w7 + tempb8.z()*w8);
    b3_pb = Vector3(tempb);
     //
    double mi0_a, mi0_b;
    if( a == 0)
    {   
        mi0_a = mi0_H;
    }else if( a == 1)
    {
        mi0_a = mi0_He;
    } else if( a == 2)
    {
        mi0_a = mi0_O;
    }
    Vector3 va = pa.VelCollParticles( b3_pa, mi0_a);
    Vector3 va_drift = pa.VelDriftParticles(b3_pa);
    if( b == 0)
    {   
        mi0_b = mi0_H;
    }else if( b == 1)
    {
        mi0_b = mi0_He;
    } else if( b == 2)
    {
        mi0_b = mi0_O;
    }
    Vector3 vb = pb.VelCollParticles( b3_pb, mi0_b);
    Vector3 vb_drift = pb.VelDriftParticles(b3_pb);
    // 
    //double proba = pb.WeightNi() / max( pa.WeightNi(), pb.WeightNi());
    //double probb = pa.WeightNi() / max( pa.WeightNi(), pb.WeightNi());
    //
    Vector3 gg = va.MinusProduct(vb); // gg
    Vector3 gx = Vector3(gg.x(), 0.0, 0.0); // g_parallel
    Vector3 gp = gg.MinusProduct(gx); // g_perpendicular
    double gab = gg.norm();
    double acosx = Cosx( a, b, dt, nb, gab, LnA);
    if( acosx > 1.0 || acosx < - 1.0)
    {
        std::cout << " Cosx > 1.0 in coll_vel_chng ";
        exit(0);
    }
    //
    ////acosx = dRand();
    //
    double sinx = sqrt( 1.0 - acosx*acosx);
    acosx = 1.0 - acosx;
    //
    double epsilon = dRand() * pi2;
    double cosep = cos(epsilon);
    double sinep = sin(epsilon);
    //
    double hhx = gp.norm() * cosep;
    double hhy = -(gg.y() * gg.x() * cosep + gab * gg.z() * sinep)/gp.norm();
    double hhz = -(gg.z() * gg.x() * cosep - gab * gg.y() * sinep)/gp.norm();
    Vector3 hh = Vector3(hhx, hhy, hhz);
    //
    double m_a = pa.WeightNi() * mi0_a;
    double m_b = pb.WeightNi() * mi0_b;
    Vector3 gh = gg.ScaleProduct(acosx).PlusProduct(hh.ScaleProduct(sinx));
    //Vector3 va_new = Vector3( va.MinusProduct(gh.ScaleProduct( mc[b][a])));
    //Vector3 vb_new = Vector3( vb.PlusProduct(gh.ScaleProduct(mc[a][b])));
    Vector3 va_new = Vector3( va.MinusProduct(gh.ScaleProduct(m_b / (m_a + m_b))));
    Vector3 vb_new = Vector3( vb.PlusProduct(gh.ScaleProduct(m_a / (m_a + m_b))));
    //
    ////std::cout << "g " << gh.norm() << " acos " << acosx << " hh " << hh.norm() << "\n";
    ////std::cout << " cab " << 3.692628212d-15 << " " << 1.90d-1 << " " << 0.9064176D-01 << " " << 0.1658835D+01 << " " << 2.0D-1 << "\n";;
    ////Vector3 ba = pa.B3atParticles(ptrArray_in);
    ////Vector3 bb = pb.B3atParticles(ptrArray_in);
    ////std::cout<< "energy0 " << pa.VelCollParticles(ba, mi0_a).norm2()* mi0_a * pa.WeightNi()<< " " 
    ////                        << pb.VelCollParticles(bb, mi0_b).norm2() * mi0_b * pb.WeightNi()<< " " 
    ////                        << pa.VelCollParticles(ba, mi0_a).norm2()* mi0_a * pa.WeightNi()+ pb.VelCollParticles(bb, mi0_b).norm2() * mi0_b * pb.WeightNi()<< " \n";
    ////std::cout<< "energy1 " << va.norm2() * mi0_a * pa.WeightNi()<< " " 
    ////                        << vb.norm2() * mi0_b * pb.WeightNi()<< " " 
    ////                        << va.norm2() * mi0_a * pa.WeightNi() + vb.norm2() * mi0_b * pb.WeightNi()<< " \n";
    ////std::cout<< "energy2 " << va_new.norm2() * mi0_a * pa.WeightNi()<< " " 
    ////                        << vb_new.norm2() * mi0_b * pb.WeightNi()<< " " 
    ////                        << va_new.norm2() * mi0_a * pa.WeightNi()+ vb_new.norm2() * mi0_b * pb.WeightNi()<<" \n";
    ////std::cout << " momentum0 x " << va.x() * m_a + vb.x() * m_b << " y " << va.y() * m_a + vb.y() * m_b << " z " << va.z() * m_a + vb.z() * m_b << "\n";  
    ////std::cout << " momentum1 x " << va_new.x() * m_a + vb_new.x() * m_b << " y " << va_new.y() * m_a + vb_new.y() * m_b << " z " << va_new.z() * m_a + vb_new.z() * m_b << "\n";   
    //
    //if( dRand() < proba)
    //{
    //    pa.SetMuVelocity(va_new, va_drift, b3_pa, mi0_a);
    //}
    //if( dRand() < probb)
    //{
    //    pb.SetMuVelocity(vb_new, vb_drift, b3_pb, mi0_b);
    //}
    pa.SetMuVelocity(va_new, va_drift, b3_pa, mi0_a);
    pb.SetMuVelocity(vb_new, vb_drift, b3_pb, mi0_b);
    ////ba = pa.B3atParticles(ptrArray_in);
    ////bb = pb.B3atParticles(ptrArray_in);
    ////std::cout<< "energy3 " << pa.VelCollParticles(ba, mi0_a).norm2()* mi0_a * pa.WeightNi()<< " " 
    ////                        << pb.VelCollParticles(bb, mi0_b).norm2() * mi0_b * pb.WeightNi()<< " " 
    ////                        << pa.VelCollParticles(ba, mi0_a).norm2()* mi0_a * pa.WeightNi()+ pb.VelCollParticles(bb, mi0_b).norm2() * mi0_b * pb.WeightNi()<< " \n";
    ////
}

//************************************************************************
// Calculate postcollision velocity ( Nanbu and Yonemura, JCP, 145, 639, 1998 and 1997)
inline void coll_vel_chng(int a, int b,
                          double dt, double nb,
                          double LnA,
                          Particles& pa, Particles& pb)
{
    // const double pi2 = 6.283185307; // 2pi
    // 0 - H; 1 - He; 2 - O
    //const double mc[3][3] ={{0.50F, 0.20F, 0.058823529F},
    //                        {0.80F, 0.50F, 0.20F},
    //                        {0.94117647F, 0.80F, 0.50F}};
    // 
    double mi0_a, mi0_b;
    if( a == 0)
    {   
        mi0_a = mi0_H;
    }else if( a == 1)
    {
        mi0_a = mi0_He;
    } else //if( a == 2)
    {
        mi0_a = mi0_O;
    }
    if( b == 0)
    {   
        mi0_b = mi0_H;
    }else if( b == 1)
    {
        mi0_b = mi0_He;
    } else //if( b == 2)
    {
        mi0_b = mi0_O;
    }
    //
    //double proba = pb.WeightNi() / max( pa.WeightNi(), pb.WeightNi());
    //double probb = pa.WeightNi() / max( pa.WeightNi(), pb.WeightNi());
    //
    Vector3 va = pa.VelParticles();
    Vector3 vb = pb.VelParticles();
    Vector3 gg = va.MinusProduct(vb); // gg
    Vector3 gx = Vector3(gg.x(), 0.0, 0.0); // g_parallel @ x
    Vector3 gp = gg.MinusProduct(gx); // g_perpendicular @ x
    double gab = gg.norm();
    double acosx = Cosx( a, b, dt, nb, gab, LnA);
    if( acosx > 1.0 || acosx < - 1.0)
    {
        std::cout << " Cosx > 1.0 in coll_vel_chng ";
        exit(0);
    }
    double sinx = sqrt( 1.0 - acosx*acosx);
    acosx = 1.0 - acosx;
    //
    double epsilon = dRand() * 6.283185307; // 2pi
    double cosep = cos(epsilon);
    double sinep = sin(epsilon);
    //
    double hhx = gp.norm() * cosep;
    double hhy = -(gg.y() * gg.x() * cosep + gab * gg.z() * sinep)/gp.norm();
    double hhz = -(gg.z() * gg.x() * cosep - gab * gg.y() * sinep)/gp.norm();
    Vector3 hh = Vector3(hhx, hhy, hhz);
    //
    double m_a = pa.WeightNi() * mi0_a;
    double m_b = pb.WeightNi() * mi0_b;
    Vector3 gh = gg.ScaleProduct(acosx).PlusProduct(hh.ScaleProduct(sinx));
    Vector3 va_new = Vector3( va.MinusProduct(gh.ScaleProduct(m_b / (m_a + m_b))));
    Vector3 vb_new = Vector3( vb.PlusProduct(gh.ScaleProduct(m_a / (m_a + m_b))));    
    // 
    pa.SetVelocity(va_new);
    pb.SetVelocity(vb_new);
}

#endif