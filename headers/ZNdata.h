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

#ifndef DATA_H
#define DATA_H


//haptics stuff
#ifdef HAPI
#include <HAPI/AnyHapticsDevice.h>
#include <HAPI/HapticForceField.h>
#include <H3DUtil/Threads.h>
#endif //HAPI

#include "constants.h"
#include "MMFF.h"
#include "PLP.h"
#include "chemscore.h"
#include "ZNmolecule.h"
#include "ddwin.h"
#include <qapplication.h>
#include "minimize.h"
#include <QUndoStack>
#include "actions.h"
#include "iodevice.h"

#include <sstream>
#include "molsketch_helium/mainwindow.h"


class Myline;
class MylineF;

#ifdef HAPI
using namespace HAPI;
#endif //HAPI
class Actions;
class Data;

class DDWin;
class Minimize;



class DataVar {
	public:
	DataVar () {};
	string name;
	virtual string write_string () = 0 ;
	virtual void load (string buffer) = 0 ;
};

class StringDataVar : public DataVar {
	private:
	string _value, _default;
	public:
	StringDataVar (string nam, string def) {name = nam; _default = def; _value = def;};
	string write_string () {return name+" "+_value;};
	void load (string buffer);
};

class ColorDataVar : public DataVar {
	private:
	color _value, _default;
	public:
	color *value_ptr () {return &_value;};
	ColorDataVar (string nam, color def) {name = nam; _default = def; _value = def;};
	string write_string () {return name+" "+int_to_string (_value.red ())+" "+int_to_string (_value.green ())+" "+int_to_string (_value.blue ())+" "+int_to_string (_value.alpha ());      };
	void load (string buffer);

};

class DoubleDataVar : public DataVar {
	private:
	double _value, _default;
	public:
	double *value_ptr () {return &_value;};
	DoubleDataVar (string nam, double def) {name = nam; _default = def; _value = def;};
	string write_string () {return name+" "+double_to_string (_value);      };
	void load (string buffer);

};


class Data : public QObject {
Q_OBJECT
public:
#ifdef HAPI
  AnyHapticsDevice *haptic_device;
#endif //HAPI
	int ncounter;

 double current_force_x;
 double current_force_y;
 double current_force_z;
 double current_position_x;
 double current_position_y;
 double current_position_z;
 double current_yaw, current_pitch, current_roll, last_pitch, last_roll, last_yaw;

QReadWriteLock *haptic_position_lock;

    Data (QApplication *master);//, DDWin *ddw);
 
    void set_ddwin (DDWin *ddw);
 //   void set_FF (MMFF *mmff, TriposFF *triposff, PLP *plp);
    QApplication *qapp;
    DDWin *ddwin;
	Actions *actions;
    MMFF *mmff;
    Minimize *minimize ;
	MainWindow *twodwin;

	color *background_color;
	double *quality_scale, *stick_radius, *vdw_scale, *sphere_radius, *double_bond_stick_scale, *double_bond_separation, *line_width;
	double *backbone_tube_helix_a, *backbone_tube_helix_b, *backbone_tube_helix_c, *backbone_tube_sheet_a, *backbone_tube_sheet_b, *backbone_tube_sheet_c, *backbone_tube_random_a, *backbone_tube_random_b, *backbone_tube_random_c;   
 //   ZNMolecule *protein;
	
//    vector <ZNMolecule*> ligands;
//    vector <vector <float> > waters;

    QUndoStack * undo_stack;

//	vector <IODevice *> iodevices;

    color constant_color;
    color charge_begin_color, charge_end_color;
    color score_begin_color, score_mid_color, score_end_color;
	vector <DataVar *> vars;

    double charge_begin, charge_end;
    double score_begin, score_mid, score_end;
	
	float total_energy_haptic;


    ZNMolecule* check_mol2 (string filename);
	void load_defaults ();
	void load_preferences ();
	void write_preferences ();

	void load_haptic_device ();
	//void haptic_callback ();
	string plants_exe;



private:	
	void _load_var (string& buffer, string reference, string &target); 	

private slots:
	void from_2D_to_3D ();
};

void *haptic_callback (void *dat);
bool is_db_extended (ZNMolecule *mol);


#endif
