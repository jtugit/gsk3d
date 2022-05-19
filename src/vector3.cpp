#include <iostream>
#include "mathutil.h"
#include "parameters.h"
#include "vector3.h"
 
Vector3::Vector3()
{
    v_x=0.0;
    v_y=0.0;
    v_z=0.0;
}
Vector3::Vector3(double x, double y, double z) 
{
    v_x = x;
    v_y = y;
    v_z = z;
}
Vector3::Vector3( const Vector3& v3b)
{
    v_x = v3b.v_x;
    v_y = v3b.v_y;
    v_z = v3b.v_z;
}