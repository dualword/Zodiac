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

