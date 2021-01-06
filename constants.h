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


#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>
#include <sstream>
#include <iostream>
#include <qwidget.h>
#include <vector>
#include <QtOpenGL/qgl.h>
#include <qpainter.h>
#include "maths.h"
#include <assert.h>


#define WIREFRAME 0
#define STICK 1
#define CPK 2
#define BALLANSTICK 3
#define BALLANDLINE 4

#define AROMATIC_RING 0
#define KEKULE 1
#define AROMATIC_BONDS 2


using namespace std;
using namespace Qt;
using namespace OpenBabel;

class color;







typedef struct {
	int type, only_to, excluding;
    float intensity;
} color_mask;







class color : public QColor {
    public:
    color ();
    color (float r, float g, float b, float a = 1.f);
    color (int r, int g, int b, int a = 255);
    inline void add (float r, float g, float b) {float cr = redF (); setRedF (cr +r);      float cg = greenF (); setGreenF (cg +g);     float cb = blueF (); setBlueF (cb +b);}

};







const    GLfloat specular[] = { 0.8f, 0.8f, 0.8f, 1.0f };
const    GLfloat diffuse[] = { 0.8f, 0.8f, 0.8f, 0.0f };

const string TITLE = "ZODIAC     v - ";
const string VERSION = "0.4 beta release";
/*
#ifdef WIN32
const string PARAMETER_DIR = "C:\\Documents and Settings\\Nick Avis\\My Documents\\Visual Studio 2005\\Projects\\zodiac\\zodiac\\parameters\\";
#else // WIN32
const string PARAMETER_DIR = "./parameters/";
#endif // WIN32
*/
//const QColor ERR_COL = QColor (255,205,105);
const QPalette::ColorRole ERR_COL = QPalette::AlternateBase;

//const Qt::BackgroundMode base = Qt::PaletteBase;

const QPalette::ColorRole base = QPalette::Base;
const float PI = 3.1415926f;

const float FF_FAR_NONBONDED_CUTOFF = 8.0;
const float FF_NEAR_NONBONDED_CUTOFF = 2.0;
const float STICK_RAD = 0.1f;
const float DOUBLE_BOND_INTER_DISTANCE = 0.15f;
const float AROMATIC_BOND_INTER_DISTANCE = 0.15f;
const int STICK_PRECISION = 15;
const int SPHERE_PRECISION = 15;
const float SPHERE_RADIUS = 0.1f;
const float VDW_SCALE = 0.2f;
const int VDW_PRECISION = 15;
const float DOUBLE_BOND_STICK_RADIUS_SCALE = 0.7f;
const float SURFACE_RESOLUTION = 3.0;

const int AROMATIC_DISPLAY_STYLE = KEKULE;

const int FOG_BEGIN = 100;

const color BACKGROUND_COLOR = color ( 1.0f, 1.0f, 1.0f, 1.0f );
const color WATER_COLOR  ( 0.0f, 1.0f, 0.0f, 0.5f );
const color BINDINGSITE_COLOR  ( 0.3f, 0.0f, 0.6f, 0.2f );
const color CURRENT_COLOR  ( 0.3f, 0.3f, 0.3f, 1.0f );
const color SELECT_COLOR  ( 0.5f, 0.1f, 0.6f, 0.5f );

const float TIMER_INTERVAL = 0.5; //msecs


double square_distance (float *coord1, float *coord2);
//double distance (float *coord1, float *coord2);
float angle (vect coord1, vect coord2, vect coord3);
//double angle (float *coord1, float *coord2, float *coord3);
//double angle (float *coord1, vector <float> coord2, float *coord3);
double dihedral (vect coord1, vect coord2, vect coord3, vect coord4);
//double dihedral (float *coord1, float *coord2, float *coord3, float *coord4);

double wilson (vect coord1, vect coord2, vect coord3, vect coord4);
//double wilson (float *coord1, float *coord2, float *coord3, float *coord4);


vect torque (vect force, vect coords, vect center);
//vector<float> torque (vector<float> force, float *coords, vector<float> center);
//vector<float> torque (float *force, float *coords, vector <float> center);  
//vector<float> torque (vector<float> force, float *coords, float *center); 
//vector<float> torque (float *force, float *coords, float *center); 


void rotate_around_vector (vect &point, vect torque, vect center, float angle);
//void rotate_around_vector (float *point, vector <float> torque, vector <float> center, float angle);
//void rotate_around_vector (float *point, float *torque, float *center, float angle);
//void rotate_around_vector (float *point, float *torque, vector <float> center, float angle);

void rotate_to_plane (vect pointA, vect pointb, vect center);

//void cross_product (float *vec1, float *vec2, float *product);
//void components (vector<float> vec, vector<float> ref, vector<float> &parall, vector<float> &normal);
//void components (float *vec, float *ref, float *parall, float *normal);
void components (vect vec, vect ref, vect &parall, vect &normal);


void axis_angle_to_quaternion (vect axis, float angle, float *quaternion);
//void axis_angle_to_quaternion (float *axis, float angle, float *quaternion);
void multiply_quaternions (float *q1, float *q2, float *product);
vect rotate_vector_using_quaternion (vect v, float *q);
//void rotate_vector_using_quaternion (float *v, float *q, float *out);
void normalize_quaternion (float *q);


//minimization

const int MAX_ITERATIONS = 1000;
const float MIN_ENERGY = 0.005f;
const float DX = 0.0005f;
const float STEP_SIZE = 0.0005f;













#endif
