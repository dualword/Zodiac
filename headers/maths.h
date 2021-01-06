
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

#ifndef MATHS_H
#define MATHS_H

#include <math.h>
#include "obabel_includes.h"


//using namespace OpenBabel;
//struct vector3;
//typedef OpenBabel::vector3 vect;

class vect : public OpenBabel::vector3 {
    public:
    vect (double x = 0., double y = 0., double z = 0.);
    vect (const OpenBabel::vector3& );
 //   double x;
 //   double y;
 //   double z;
   
    inline double square_module () {return length_2 ();}
    inline double module () {return length ();};
    inline void trunc_at (double f) {if (module () > f) scale_to (f);}
    void multiply (double f);
    void null ();
    void normalise (); 
    void scale_to (double f);

};




class quaternion {
	public:
	quaternion (double w, double x, double y, double z);
	quaternion ();
	inline double x () {return _x;}
	inline double y () {return _y;}
	inline double z () {return _z;}
	inline double w () {return _w;}
	private:
	double _x, _y, _z, _w;
	
};
double square_distance (vect coord1, vect coord2);
inline double dot_product (vect v1, vect v2) {return dot (v1, v2);};
vect cross_product (vect v1, vect v2);

inline vect sum (vect v1, vect v2) {return vect (v1.x()+v2.x(), v1.y()+v2.y(), v1.z()+v2.z());};
inline vect subtract (vect v1, vect v2) {return vect (v1.x()-v2.x(), v1.y()-v2.y(), v1.z()-v2.z());};
inline vect mean (vect v1, vect v2) {vect s = sum (v1, v2); s.multiply (0.5f); return s;}
vect circumcenter (vect v1, vect v2, vect v3);
double dist (vect v1, vect v2);
void interpolate (vect v1, vect v2, vect v3, vect& i1, vect &i2);





#endif
