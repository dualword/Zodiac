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


#include "constants.h"
#include <math.h>

#ifdef WIN32
#include <float.h>
#define isnan _isnan
#endif // WIN32


color average_3_colors (float score, color c1, color c2, color c3, float begin, float mid, float end) {
float mid_width = 0.1 * (mid - begin);
float mid_width2 = 0.1 * (end - mid);
if (mid_width2 < mid_width) mid_width = mid_width2;



float q = score;
if (q <= begin) {
	return c1;
}
else if (q >= end) {
	return c3;
}
else if (q>= mid - mid_width && q<=  mid + mid_width) {
	return c2;
	
}
else if (q<= mid) {
	float perc = (mid - q - mid_width) / (mid - begin - mid_width);
	float b_red = c1.redF ();
	float b_green = c1.greenF ();
	float b_blue = c1.blueF ();
	float b_alpha = c1.alphaF ();
	float m_red = c2.redF ();
	float m_green = c2.greenF ();
	float m_blue = c2.blueF ();
	float m_alpha = c2.alphaF ();
	
	float r = (b_red  * perc + m_red   * (1-perc)) ;
	float g = (b_green* perc + m_green * (1-perc)) ;
	float b = (b_blue * perc + m_blue  * (1-perc)) ;
	float a = (b_alpha * perc + m_alpha  * (1-perc)) ;
	return color (r, g, b, a);
	
}
else {
	float perc = (q - mid -mid_width) / (end -mid_width- mid);
	float e_red = c3.redF ();
	float e_green = c3.greenF ();
	float e_blue = c3.blueF ();
	float e_alpha = c3.alphaF ();
	float m_red = c2.redF ();
	float m_green = c2.greenF ();
	float m_blue = c2.blueF ();
	float m_alpha = c2.alphaF ();
	
	float r = (m_red  * perc + e_red   * (1-perc)) ;
	float g = (m_green* perc + e_green * (1-perc)) ;
	float b = (m_blue * perc + e_blue  * (1-perc)) ;
	float a = (m_alpha * perc + e_alpha  * (1-perc)) ;
	return color (r, g, b, a);
	
}

}





float angle (vect coord1, vect coord2, vect coord3) {
    vect v1 = subtract (coord1, coord2);
    vect v2 = subtract (coord3, coord2);
    float cosine = dot_product (v1, v2)/(v1.module () * v2.module ());
    if (cosine < -1) cosine = -1;
    else if (cosine > 1) cosine = 1;
    return (acos (cosine)) * 180 / PI;
	
    
}


float signed_angle (vect coord1, vect coord2, vect coord3) {
    vect v1 = subtract (coord1, coord2);
    vect v2 = subtract (coord3, coord2);
	vect c = cross_product (v1, v2);
    float angle = atan2(c.module (), dot_product(v1, v2))* 180 / PI;
	vect reference (0., 0., 1.);
    return dot_product (c, reference) < 0.f ? -angle : angle;
}
	

double dihedral (vect coord1, vect coord2, vect coord3, vect coord4) {
    vect v1 = subtract (coord1, coord2);
    vect v2 = subtract (coord3, coord2);
    vect n1 = cross_product (v1, v2);
    vect v3 = v2;
    v3.multiply (-1.0f);
    vect v4 = subtract (coord4, coord3);
    vect n2 = cross_product (v3, v4);
    vect orig;
    orig.null ();
	float mult = 1.f;
	if (dot_product (v2, cross_product(n1, n2)) < 0.) mult = -1.f;
    return mult * angle (n1, orig, n2);
}




double wilson (vect coord1, vect coord2, vect coord3, vect coord4) {
    vect v1 = subtract (coord1, coord2);
    vect v2 = subtract (coord3, coord2);
    vect n1 = cross_product (v1, v2);
    float a = n1.x();
    float b = n1.y();
    float c = n1.z();
    float d = -1.f * (coord2.x() * n1.x() + coord2.y() * n1.y() + coord2.z() * n1.z());
    double  dis = (a*coord4.x() + b*coord4.y() +c*coord4.z() +d) / n1.module ();
    double sin = dis / dist (coord4, coord2);
    if (sin < -1.f) sin = -1.f;
    if (sin > 1.f) sin = 1.f;
    return asin (sin);
}

double  wilson (float  *coord1, float  *coord2, float  *coord3, float  *coord4){
    float  v1x, v1y, v1z, v2x, v2y, v2z, n1[3], a, b, c, d;

    v1x = coord1[0]-coord2[0];
    v1y = coord1[1]-coord2[1];
    v1z = coord1[2]-coord2[2];
    

    v2x = coord3[0]-coord2[0];
    v2y = coord3[1]-coord2[1];
    v2z = coord3[2]-coord2[2];

    n1[0] = v1y*v2z - v1z*v2y;
    n1[1] = v1z*v2x - v1x*v2z;
    n1[2] = v1x*v2y - v1y*v2x;

    a = n1[0];
    b = n1[1];
    c = n1[2];
    d = -1* (coord2[0]*n1[0]+coord2[1]*n1[1]+coord2[2]*n1[2]);
    double  dist = (a*coord4[0] + b*coord4[1] +c*coord4[2] +d)/sqrt (a*a +b*b+c*c);
    double si = dist/distance (coord4, coord2);
    if (si > 1.) si = 1.;
    else if (si < -1.) si = -1.;
    return asin (si);
}


void rotate_around_vector (vect &point, vect torque, vect center, float angle) {
    torque.normalise ();
    vect dir = subtract (point, center);
    float sa =sin(angle);
    float ca =cos(angle);
    float u = torque.x();
    float v = torque.y();
    float w = torque.z();
    float x = dir.x();
    float y = dir.y();   
    float z = dir.z();

    float ux = u*x;
    float uy = u*y;
    float uz = u*z;
    float vx = v*x;
    float vy = v*y;
    float vz = v*z;
    float wx = w*x; 
    float wy = w*y;
    float wz = w*z;

    point.x()= u*(ux+vy+wz)+(x*(v*v+w*w)-u*(vy+wz))*ca+(-wy+vz)*sa+center.x();
    point.y()= v*(ux+vy+wz)+(y*(u*u+w*w)-v*(ux+wz))*ca+(wx-uz)*sa+center.y();
    point.z()= w*(ux+vy+wz)+(z*(u*u+v*v)-w*(ux+vy))*ca+(-vx+uy)*sa+center.z();
    
}


void rotate_to_plane (vect up_point, vect right_point, vect cent_point) {
    vect pointA, pointB, pointC;
    pointA = subtract (up_point, cent_point);
    pointB = subtract (right_point, cent_point);
    pointC = cent_point;
    float rangle;
    float A, B, C;
    glTranslatef (pointC.x(), pointC.y(), pointC.z());
    //equation of plane Ax + By + Cz = 0

    if (!pointA.z() && !pointB.z()) { //all three points on z=0 plane...  we're already there

    }
    else if (pointA.x() && pointB.x() && (pointA.y()*pointB.x()-pointA.x()*pointB.y())){
    //equation (A/C)*x + (B/C)*y +z = 0
        C = 1.f;
        B = (pointA.x()*pointB.z()-pointB.x()*pointA.z())/(pointA.y()*pointB.x()-pointA.x()*pointB.y());
        A = (-pointA.z()-B*pointA.y())/pointA.x();
        //rotation axis A*x + B*y
        //perpendicular aA + bB = 0


        //perpendicular -B*x + A*y = 0
        //point x=1 y=B/A
        //point on the plane z= -(A + B*B/A)
        float xp, yp, zp;
        xp = 1.f;
        yp = B/A;
        zp = -(A + B*B/A);

        float dist =  sqrt(xp*xp+ yp*yp+zp*zp);
        float s = zp / dist;
        if (s > 1.f) s= 1.f;
        else if (s < -1.f) s = -1.f;


        rangle = asin (s);
//        cout<<rangle*180/PI;
   //     cout<<(-A/B*xp+yp)*zp;
        if (rangle<0) {
            rangle = -rangle;

        }
    //        if (yp > -A/B*xp);
     //   rangle = PI+rangle;
        glRotatef (-rangle*180/PI, B, -A, 0); 



    }

    vect O;
    O.null ();
    vect rot_vec (B, -A, 0.f);
    vect rot_vec2 (0.f, 0.f, 1.f);

    vect a_point = pointA;
    vect b_point = pointB;
  
    rotate_around_vector (a_point, rot_vec, O, rangle);
//    cout << a_point [0]<<" "<< a_point [1]<<" "<< a_point [2]<<" "<<rangle<<endl;


    float dis = a_point.module ();
    float sang;
    float start_ang;
    float c = a_point.x()/dis;
    if (c > 1.f) c= 1.f;
    else if (c < -1.f) c = -1.f;
    sang = acos (c);


    if (a_point.y()<0.f) start_ang = sang;

 //   cout <<sang<<endl;
    else start_ang = -sang;
    glRotatef (-start_ang*180/PI, 0, 0, 1);


    rotate_around_vector (b_point, rot_vec, O, rangle);
    rotate_around_vector (b_point, rot_vec2, O, start_ang);
  //  cout << b_point [0]<<" "<< b_point [1]<<" "<< b_point [2]<<" "<<start_ang<<endl;
    if (b_point.y()>0)     glRotatef (180, 1, 0, 0);


}




void components (vect vec, vect ref, vect &parallel, vect &normal){
    assert (!isnan (vec.module ()));
    assert (!isnan (ref.module ()));
	ref.normalise ();
    double mod = dot_product (vec, ref);
//    cerr << vec.x() << " "<< vec.y() << " "<< vec.z() << " "<< ref.x() << " "<< ref.y() << " "<< ref.z() << " "<< endl;
    ref.multiply (mod);
    parallel = ref;
    normal = subtract (vec, parallel);
    assert (!isnan (parallel.module ()));
    assert (!isnan (normal.module ()));
}







vect torque (vect force, vect coords, vect center) {
    vect pos = subtract (coords, center);
    return cross_product (pos, force);
}


void axis_angle_to_quaternion (vect axis, float angle, float *quaternion) {
    axis.normalise ();
    float s = sin(angle/2);
    quaternion[0] = cos(angle/2);
    quaternion[1] = axis.x() * s;
    quaternion[2] = axis.y() * s;
    quaternion[3] = axis.z() * s;
}

quaternion axis_angle_to_quaternion (vect axis, float angle) { //angle in rads
    axis.normalise ();
    float s = sin(angle/2);
	quaternion q (cos(angle/2), axis.x() * s,  axis.y() * s, axis.z() * s);
	return q;
}

void axis_angle_to_quaternion (float *axis, float angle, float *quaternion) { //angle in rads

    float ax[3];
    ax[0] = axis[0];
    ax[1] = axis[1];
    ax[2] = axis[2];
    float ax_mod = sqrt (ax[0]*ax[0]+ax[1]*ax[1]+ax[2]*ax[2]);
    if (ax_mod && angle) {
        ax[0] /= ax_mod;
        ax[1] /= ax_mod;
        ax[2] /= ax_mod;

   //     cout << "axis "<<ax[0]<<" "<<ax[1]<<" "<<ax[2]<< endl;
        float s = sin(angle/2);
        quaternion[0] = cos(angle/2);
        quaternion[1] = ax[0] * s;
        quaternion[2] = ax[1] * s;
        quaternion[3] = ax[2] * s;

    }

    else {
        quaternion[0] = 1;
        quaternion[1] = 0;
        quaternion[1] = 0;
        quaternion[1] = 0;
    }
}


quaternion yaw_pitch_roll_to_quaternion (float yaw, float pitch, float roll) {
	quaternion x = axis_angle_to_quaternion (vect(1.,0.,0.), roll);
	quaternion y = axis_angle_to_quaternion (vect(0.,1.,0.), pitch);
	quaternion z = axis_angle_to_quaternion (vect(0.,0.,1.), yaw);



 x = multiply_quaternions ( x, y );
 y = multiply_quaternions (x, z);

 return x;
}


void multiply_quaternions (float *q1, float *q2, float *prod) {

    prod[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
    prod[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
    prod[2] = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1];
    prod[3] = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0];



}


quaternion multiply_quaternions (quaternion q1, quaternion q2) {
	double w = q1.w ()*q2.w () - q1.x ()*q2.x () - q1.y ()*q2.y () - q1.z ()*q2.z ();
	double x = q1.w ()*q2.x () + q1.x ()*q2.w () + q1.y ()*q2.z () - q1.z ()*q2.y ();
	double y = q1.w ()*q2.y () - q1.x ()*q2.z () + q1.y ()*q2.w () + q1.z ()*q2.x ()	;
	double z = q1.w ()*q2.z () + q1.x ()*q2.y () - q1.y ()*q2.x () + q1.z ()*q2.w ()	;
	return quaternion (w, x, y, z);
	
}



void invert_matrix_9 (double *m, double *i) {
    double m11, m12, m13, m21, m22, m23, m31, m32, m33;
    m11  = m[0];
    m12  = m[1];
    m13  = m[2];
    m21  = m[3];
    m22  = m[4];
    m23  = m[5];
    m31  = m[6];
    m32  = m[7];
    m33  = m[8];
	double det = m11*m22*m33 + m12*m23*m31 + m13*m31*m32 - m11*m23*m32 - m12*m21*m33 - m13*m22*m31;
	if (det == 0) det = 0.0000001f;
	i[0] = (m22*m33 - m23*m32) / det;
	i[1] = (m13*m32 - m12*m33) / det;
	i[2] = (m12*m23 - m13*m22) / det;
	i[3] = (m23*m31 - m21*m33) / det;
	i[4] = (m11*m33 - m13*m31) / det;
	i[5] = (m13*m21 - m11*m23) / det;
	i[6] = (m21*m32 - m22*m31) / det;
	i[7] = (m12*m31 - m11*m32) / det;
	i[8] = (m11*m22 - m12*m21) / det;

}


void map_vector_on_vector_quaternion (vect v1, vect v2, float *quaternion) {
	v1.normalise();
	v2.normalise();
	vect axis = cross_product(v1, v2);
	float angle = acos (dot_product(v1, v2));
	axis_angle_to_quaternion (axis, angle, quaternion);
	
}


quaternion map_vector_on_vector_quaternion (vect v1, vect v2) {
	quaternion q (0., 1., 0., 0.);
	v1.normalise();
	v2.normalise();
	vect axis = cross_product(v1, v2);
	double dotproduct = dot_product(v1, v2);
	if (dotproduct >1.) dotproduct = 1.;
	if (dotproduct <-1.) dotproduct = -1.;
	float angle = acos (dotproduct);
	q = axis_angle_to_quaternion (axis, angle);
	return q;
	
}



vect rotate_vector_using_matrix_9 (vect v, double *m) {
    double m11, m12, m13, m21, m22, m23, m31, m32, m33;
    m11  = m[0];
    m12  = m[1];
    m13  = m[2];
    m21  = m[3];
    m22  = m[4];
    m23  = m[5];
    m31  = m[6];
    m32  = m[7];
    m33  = m[8];
    vect out;
    out.x() = m11*v.x() + m12*v.y() + m13*v.z();
    out.y() = m21*v.x() + m22*v.y() + m23*v.z();
    out.z() = m31*v.x() + m32*v.y() + m33*v.z();
    return out;
}


vect rotate_vector_using_matrix_16 (vect v, double *m) {
    double m11, m12, m13, m21, m22, m23, m31, m32, m33; // w, x, y, z;
    m11  = m[0];
    m12  = m[1];
    m13  = m[2];
    m21  = m[4];
    m22  = m[5];
    m23  = m[6];
    m31  = m[7];
    m32  = m[8];
    m33  = m[9];
    vect out;
    out.x() = m11*v.x() + m12*v.y() + m13*v.z();
    out.y() = m21*v.x() + m22*v.y() + m23*v.z();
    out.z() = m31*v.x() + m32*v.y() + m33*v.z();
    return out;
}

vect rotate_vector_using_quaternion (vect v, quaternion q) {
	double m11, m12, m13, m21, m22, m23, m31, m32, m33, w, x, y, z;
    w = q.w ();
    x = q.x ();
    y = q.y ();
    z = q.z ();
    m11  = 1 - 2 * ( y*y + z*z );
    m12  =     2 * ( x*y - z*w );
    m13  =     2 * ( x*z + y*w );
    m21  =     2 * ( x*y + z*w );
    m22  = 1 - 2 * ( x*x + z*z );
    m23  =     2 * ( y*z - x*w );
    m31  =     2 * ( x*z - y*w );
    m32  =     2 * ( y*z + x*w );
    m33  = 1 - 2 * ( x*x + y*y );
    vect out;
    out.x() = m11*v.x() + m12*v.y() + m13*v.z();
    out.y() = m21*v.x() + m22*v.y() + m23*v.z();
    out.z() = m31*v.x() + m32*v.y() + m33*v.z();
    return out;
}



vect rotate_vector_using_quaternion (vect v, float *q) {
    double m11, m12, m13, m21, m22, m23, m31, m32, m33, w, x, y, z;
    w = q[0];
    x = q[1];
    y = q[2];
    z = q[3];
    m11  = 1 - 2 * ( y*y + z*z );
    m12  =     2 * ( x*y - z*w );
    m13  =     2 * ( x*z + y*w );
    m21  =     2 * ( x*y + z*w );
    m22  = 1 - 2 * ( x*x + z*z );
    m23  =     2 * ( y*z - x*w );
    m31  =     2 * ( x*z - y*w );
    m32  =     2 * ( y*z + x*w );
    m33  = 1 - 2 * ( x*x + y*y );
    vect out;
    out.x() = m11*v.x() + m12*v.y() + m13*v.z();
    out.y() = m21*v.x() + m22*v.y() + m23*v.z();
    out.z() = m31*v.x() + m32*v.y() + m33*v.z();
    return out;
}


void rotate_vector_using_quaternion (float *v, float *q, float *out) {
    float m11, m12, m13, m21, m22, m23, m31, m32, m33, w, x, y, z;
    w = q[0];
    x = q[1];
    y = q[2];
    z = q[3];
/*    m11 = 1-2*y*y-2*z*z;
    m12 = 2*x*y-2*w*z;
    m13 = 2*x*z+2*w*y;
    m21 = 2*x*y+2*w*z;
    m22 = 1-2*x*x-2*z*z;
    m23 = 2*y*z-2*w*x;
    m31 = 2*x*z-2*w*y;
    m32 = 2*y*z-2*w*x;
    m33 = 1-2*x*x-2*y*y;
*/


    m11  = 1 - 2 * ( y*y + z*z );
    m12  =     2 * ( x*y - z*w );
    m13  =     2 * ( x*z + y*w );
    m21  =     2 * ( x*y + z*w );
    m22  = 1 - 2 * ( x*x + z*z );
    m23  =     2 * ( y*z - x*w );
    m31  =     2 * ( x*z - y*w );
    m32  =     2 * ( y*z + x*w );
    m33  = 1 - 2 * ( x*x + y*y );








    out[0] = m11*v[0] + m12*v[1] + m13*v[2];
    out[1] = m21*v[0] + m22*v[1] + m23*v[2];
    out[2] = m31*v[0] + m32*v[1] + m33*v[2];
 
}

void normalize_quaternion (float *q) {
    float mod = sqrt (q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
//    cout <<q[0]<<"   "<<q[0]/mod<<"      "<<mod<<endl;
    q[0] /= mod;
    q[1] /= mod;
    q[2] /= mod;
    q[3] /= mod;
}



double string_to_double (string s) {
	istringstream ss (s);
	double d;
	ss >> d;
	return d;
};


int string_to_int (string s) {
	istringstream ss (s);
	int i;
	ss >> i;
	return i;
}
string int_to_string (int i) {
	stringstream ss;
	ss << i;
	return ss.str ();
}

string double_to_string (double d) {
	stringstream ss;
	ss << d;
	return ss.str ();
}


color mean (color c1, color c2) {
	return color ((float)(c1.redF()+c2.redF())*0.5f,(float) (c1.greenF()+c2.greenF())*0.5f,(float) (c1.blueF()+c2.blueF())*0.5f,(float) (c1.alphaF()+c2.alphaF())*0.5f);
}
color::color () : QColor ()
{
    setRedF (1.f);
    setGreenF (1.f);
    setBlueF (1.f);
    setAlphaF (1.f);
}

color::color (float r, float g, float b, float a) :  QColor ()
{
    setRedF (r);
    setGreenF (g);
    setBlueF (b);
    setAlphaF (a);
}


color::color (int r, int g, int b, int a) : QColor (r, g, b, a) {}




