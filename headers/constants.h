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


#include <iomanip>
#include <queue>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <qwidget.h>
#include <QtOpenGL/qgl.h>
#include <qpainter.h>
#include "maths.h"
#include <math.h>
#include <assert.h>
#include <QTime>





#include "obabel_includes.h"


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

static QTime last_time;



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
const string VERSION = "0.6.5 beta";
/*
#ifdef WIN32
const string PARAMETER_DIR = "C:\\Documents and Settings\\Nick Avis\\My Documents\\Visual Studio 2005\\Projects\\zodiac\\zodiac\\parameters\\";

#else // WIN32
const string PARAMETER_DIR = "./parameters/";
#endif // WIN32
*/

const QPalette::ColorRole ERR_COL = QPalette::AlternateBase;


const QPalette::ColorRole base = QPalette::Base;
//const double M_PI = 3.1415;
const double PI = M_PI;

const float FF_FAR_NONBONDED_CUTOFF = 5.5;
const float FF_NEAR_NONBONDED_CUTOFF = 1.5;
const float DOUBLE_BOND_INTER_DISTANCE = 0.15f;
const float AROMATIC_BOND_INTER_DISTANCE = 0.15f;



const float DOUBLE_BOND_STICK_RADIUS_SCALE = 0.7f;
const float SURFACE_RESOLUTION = 3.0;

const int AROMATIC_DISPLAY_STYLE = KEKULE;

const int FOG_BEGIN = 0;

//const color BACKGROUND_COLOR = color ( 1.0f, 1.0f, 1.0f, 1.0f );
//const color WATER_COLOR  ( 0.0f, 1.0f, 0.0f, 0.5f );
//const color BINDINGSITE_COLOR  ( 0.3f, 0.0f, 0.6f, 0.2f );
//const color CURRENT_COLOR  ( 0.3f, 0.3f, 0.3f, 1.0f );
const color SELECT_COLOR  ( 19, 3, 82, 255 );

const float TIMER_INTERVAL = 0.5; //msecs



color average_3_colors (float score, color c1, color c2, color c3, float begin, float mid, float end) ;


//double distance (float *coord1, float *coord2);
float angle (vect coord1, vect coord2, vect coord3);

float signed_angle (vect coord1, vect coord2, vect coord3);

//double angle (float *coord1, float *coord2, float *coord3);
//double angle (float *coord1, vector <float> coord2, float *coord3);
double dihedral (vect coord1, vect coord2, vect coord3, vect coord4);
//double dihedral (float *coord1, float *coord2, float *coord3, float *coord4);

double wilson (vect coord1, vect coord2, vect coord3, vect coord4);
//double wilson (float *coord1, float *coord2, float *coord3, float *coord4);


vect torque (vect force, vect coords, vect center);
color mean (color c1, color c2);


void rotate_around_vector (vect &point, vect torque, vect center, float angle);


void rotate_to_plane (vect pointA, vect pointb, vect center);


void components (vect vec, vect ref, vect &parall, vect &normal);
quaternion map_vector_on_vector_quaternion (vect v1, vect v2);
void axis_angle_to_quaternion (vect axis, float angle, float *quaternion);
quaternion yaw_pitch_roll_to_quaternion (float yaw, float pitch, float roll);
quaternion axis_angle_to_quaternion(vect axis, float angle);
void multiply_quaternions (float *q1, float *q2, float *product);
quaternion multiply_quaternions (quaternion q1, quaternion q2);
//vect rotate_vector_using_quaternion (vect v, float *q);
vect rotate_vector_using_quaternion(vect v, quaternion q);
vect rotate_vector_using_matrix_16 (vect v, double *m);
vect rotate_vector_using_matrix_9 (vect v, double *m);
void invert_matrix_9 (double *mat, double *inverse); 
void normalize_quaternion (double *q);


double string_to_double (string s);
int string_to_int (string s);
string int_to_string (int i);
string double_to_string (double d);


//minimisation

const int MAX_ITERATIONS = 1000;
const double MIN_ENERGY = 0.005f;
const double DX = 0.0005f;
const double STEP_SIZE = 0.0005f;













#endif
