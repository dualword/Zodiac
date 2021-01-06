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
#ifndef THREAD_H
#define THREAD_H

#include <QMutex>
#include <QThread>
#include "graphical_object.h"
#include "minimize.h"
#include "ddwin.h"
#include "obabel_includes.h"
#include "wiimote.h"


class DDWin; class Minimize;
class Thread : public QThread   {
    Q_OBJECT
    public:
    Thread (QObject *parent);
    bool is_running;
	virtual void stop ();
	string _name;
	string _doing;
	int _min_step;
	int _max_step;
	int _current_step;
  //  protected:
//      void run();
};

class MapThread : public Thread {
    Q_OBJECT
    public:
    MapThread (QObject *parent, Map *mp, DDWin *ddw);
    Map * map;
    float alph, res;    
    DDWin *ddwin;


    protected:
    void run ();
    void compute_map_data ();
};

class SurfaceThread : public Thread {
    Q_OBJECT
    public:
    SurfaceThread (QObject *parent, Surface *surf, DDWin *ddw);
    Surface *surface;
    float alph, res;    
    DDWin *ddwin;


    protected:
    void run ();
    void compute_surface_data ();
};


class SystematicConformationThread : public Thread {
    Q_OBJECT
public:
    SystematicConformationThread (QObject *parent, ZNMolecule *mol, Database *dat, DDWin *ddw);
    ZNMolecule *molecule;
	Database *database;
    DDWin *ddwin;
	
	
protected:
    void run ();
	void generate_conformation (vector <double> &dofs, int level);
};

class StochasticConformationThread : public Thread {
    Q_OBJECT
public:
    StochasticConformationThread (QObject *parent, ZNMolecule *mol, Database *dat, DDWin *ddw);
    ZNMolecule *molecule;
	Database *database;
    DDWin *ddwin;
	void set_results_number (int n) {_results_number = n;};
	void set_time_limit (int n) {_time_limit = n;};
	
protected:
    void run ();
	int _results_number, _time_limit;
};


class DockingThread : public Thread {
    Q_OBJECT
public:
    DockingThread (QObject *parent, ZNMolecule *mol, Database *dat, DDWin *ddw);
    ZNMolecule *molecule;
	Database *database;
    DDWin *ddwin;
	void set_results_number (int n) {_results_number = n;};
	void set_time_limit (int n) {_time_limit = n;};	
	void set_bindingsite_radius (float r) {_bindingsite_radius = r;};
	void set_bindingsite_center (vect c) {_bindingsite_center = c;};
	
protected:
    void run ();
	int _results_number, _time_limit;
	vect _bindingsite_center;
	float _bindingsite_radius;
};


class MinimiseThread : public Thread {
    Q_OBJECT
    public:
    MinimiseThread (QObject *parent, DDWin *ddw);
    ZNMolecule *molecule;
    DDWin *ddwin;

    void inline set_molecule (ZNMolecule *mol) {molecule = mol;};


    signals:
    void ask_redraw (ZNMolecule *mol);

	private:

    protected:
    void run ();
};

class DatabaseMinimiseThread : public Thread {
	Q_OBJECT
	public:
	DatabaseMinimiseThread (QObject *parent, Database *dat, DDWin *ddw);
	Database * database;
	DDWin *ddwin;
	protected:
	void run ();	
	};

class HapticWorkerThread;
class HapticForceThread : public Thread {
	Q_OBJECT
	public:
	HapticForceThread (QObject *parent, DDWin *ddw);
	DDWin *ddwin;
	    protected:
    void run ();
};

class HapticThread : public Thread {
    Q_OBJECT
    public:
    HapticThread (QObject *parent, DDWin *ddw);
	void stop ();
	bool save_result, last_is_user_selected, user_save_result, is_worse, saving;
	float save_E_threshold, mult;
    ZNMolecule *molecule;
	Database *results;
    DDWin *ddwin;
	Minimize *minimise;
	vector <float> viscous_force_x, viscous_force_y, viscous_force_z; 
    void inline set_molecule (ZNMolecule *mol) {molecule = mol;};
    int haptic_dof_mode;
    bool automove, color_by_score;
    void initialise (Minimize * min);

	int next_atom_to_check, internal_interaction_counter, number_of_atoms;
	vector <Atom *> atoms;
	vect last_center;
	//vector <vect> forces;
	//vector <double> scores; 


	
	QMutex atoms_mutex, internal_interactions_mutex, nonbonded_interactions_mutex;

	
	
    signals:

    void ask_color_by_score (ZNMolecule *mol);


    protected:
    void run ();
	
	private:
	vect total_torque;
	double total_energy;
	int ncycles;
	int number_of_threads;
	vector <HapticWorkerThread *> worker_threads;

public:
		vector <ForceFieldInteraction *> internal_interactions;
		//vector <ForceFieldInteraction *> restrains;	
private:
	vector <HBond> Hbonds;
	queue <ForceFieldInteraction *> nonbonded_interactions;

	
	public:
		bool have_to_work_on_internal_interactions ();
		bool have_to_work_on_nonbonded_interactions ();
		bool are_atoms_finished ();
		bool are_internal_interactions_finished ();
		bool is_end_of_cycle ();
		void compute_one_internal_interaction ();	
		void compute_one_nonbonded_interaction ();
		void compute_nonbonded_interactions_for_atom (int a);
		void compute_internal_interaction (int i);
		void compute_restrain (int i);
		void load_internal_interactions ();
		void load_interactions_for_new_atom ();
		void end_of_calculation_for_atom (Atom *at);
		void end_of_cycle ();
};


class HapticWorkerThread : public Thread {
    public:
			
		enum HapticWorkerType { INTERNAL, INTERACTION, SINGLE};
		HapticWorkerThread (HapticThread *ht, int i, HapticWorkerType typ);

		HapticWorkerType thread_type;
    protected:
		void run ();
	private:
		HapticThread *master;
		 int number;

};



class GridThread : public Thread {
    Q_OBJECT
    public:
    GridThread (QObject *parent, Grid *gri, DDWin *ddw);
    Grid *grid;
    DDWin *ddwin;


//    signals:
//    void ask_redraw (ZNMolecule *mol);
//    void ask_color_by_score (ZNMolecule *mol);


    protected:
    void run ();

};


class HeadTrackingThread : public Thread {
    Q_OBJECT
    public:
    HeadTrackingThread (QObject *parent, DDWin *ddw);
    DDWin *ddwin;



//    signals:
//    void ask_redraw (ZNMolecule *mol);
//    void ask_color_by_score (ZNMolecule *mol);

    Wiimote *wiimote;

    protected:
    void run ();
//    int lastx, lasty;
    signals:
	void head_moved (int x, int y);	

};



class WiimoteTrackingThread : public Thread {
    Q_OBJECT
    public:
    WiimoteTrackingThread (QObject *parent, DDWin *ddw);
    DDWin *ddwin;

    Wiimote *wiimote;
    void switch_to_mode (int i);

    protected:
    int mode;
    bool is_moving, is_rotating;
    void run ();
    float lastx, lasty, lastd;
    vect last_orientation;
    signals:
	void move_camera (float x, float y, float z);	
	void move_target (float x, float y, float z);
	void rotate_world (double x1, double y1, double z1, double x2, double y2, double z2);
	void rotate_target (double x1, double y1, double z1, double x2, double y2, double z2);
};




#endif //THREAD_H
