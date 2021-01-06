/*   Zodiac - molecular modelling environment. www.zeden.org
    Copyright (C) 2008  Nicola Zonta (nicola.zonta(at)zeden.org)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "maths.h"
#include <iostream>



#define DELTA 0.000000001f

using namespace std;
quaternion::quaternion (double w, double x, double y, double z) {
	_x = x;
	_y = y;
	_z = z;
	_w = w;
}


quaternion::quaternion () {
	_x = 1.;
	_y = 0.;
	_z = 0.;
	_w = 0.;
}


vect::vect (double xi, double yi, double zi) : OpenBabel::vector3 (xi, yi, zi) {

}

vect::vect (const OpenBabel::vector3& v) {vect (v.x(), v.y(), v.z ()); }

void vect::null () {
    Set (0.0, 0.0, 0.0);
}

void vect::normalise () {
    double mod = module (); 
    if (mod < DELTA) {
   //     cerr << "normalising null vector"<<endl;
        null ();
    }
    else {
        x () /= mod;
        y () /= mod;
        z () /= mod;
    }
}

void vect::multiply (double f) {
  //  if (f < DELTA) {
   //     cerr << "multiplying vector by 0"<<endl;
   //     null ();
   // }
  //  else {
        Set (x () * f, y () * f, z () * f);
   // }
}
    
void vect::scale_to (double f) {
    normalise ();
    multiply (f);
}



double square_distance (vect  coord1, vect  coord2) {
    return (coord1.x()-coord2.x())*(coord1.x()-coord2.x()) + (coord1.y()-coord2.y())*(coord1.y()-coord2.y()) + (coord1.z()-coord2.z())*(coord1.z()-coord2.z());
}



vect circumcenter (vect v1, vect v2, vect v3) {
	vect A = v1;
	vect B = v2;
	vect C = v3;
	double a2 = square_distance(B, C);
	double b2 = square_distance(C, A);
	double c2 = square_distance(A, B);
	double cx = C.x() * c2*(a2+b2-c2);
	double cy = C.y() * c2*(a2+b2-c2);
	double cz = C.z() * c2*(a2+b2-c2);
	double ax = A.x() * a2*(c2+b2-a2);	
	double ay = A.y() * a2*(c2+b2-a2);
	double az = A.z() * a2*(c2+b2-a2);
	double bx = B.x() * b2*(a2+c2-b2);
	double by = B.y() * b2*(a2+c2-b2);
	double bz = B.z() * b2*(a2+c2-b2);
	

	double den = 2 * (a2*b2 + a2*c2 + b2*c2)-(a2*a2 + b2*b2 + c2*c2);
	if (den == 0.f) den = 0.0001f;
	vect out = vect ((ax+bx+cx)/den, (ay+by+cy)/den, (az+bz+cz)/den);
	return out;	
}

void interpolate (vect v1, vect v2, vect v3, vect& i1, vect &i2) {
	vect d1 = subtract (v1, v2);
	vect d2 = subtract (v3, v2);
	d1.normalise ();
	d2.normalise ();
	if (dot_product(d1, d2) < -0.999) {
		i1 = mean (v1, v2);
		i2 = mean (v2, v3);
	}
	else {
	vect c = circumcenter (v1, v2, v3);
	vect m1 = mean (v1, v2);
	vect m2 = mean (v2, v3);
	vect r1 = subtract (m1, c);
	vect r2 = subtract (m2, c);
	r1.normalise();
	r2.normalise();
	float r = dist (c, v1);

	r1.multiply(r);
	r2.multiply(r);
	i1 = sum (c, r1);
	i2 = sum (c, r2);
	}
}

double dist (vect v1, vect v2) {
    vect v = subtract (v1, v2); 
    return v.module ();
}


vect cross_product (vect v1, vect v2) {
    vect vv ;

    vv.x() =   v1.y()*v2.z() - v1.z()*v2.y() ;
    vv.y() = - v1.x()*v2.z() + v1.z()*v2.x() ;
    vv.z() =   v1.x()*v2.y() - v1.y()*v2.x() ;

    return ( vv ) ;
  
}

