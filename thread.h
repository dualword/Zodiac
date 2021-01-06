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
#include <openbabel/forcefield.h>




#include "wiimote.h"


class DDWin; class Minimize;
class Thread : public QThread   {
    Q_OBJECT
    public:
    Thread (QObject *parent);
    bool is_running;
	virtual void stop ();
  //  protected:
//      void run();
};



class SurfaceThread : public Thread {
    Q_OBJECT
    public:
    SurfaceThread (QObject *parent, Surface *surf, DDWin *ddw);
    Surface *surface;
    float a, res;    
    DDWin *ddwin;


    protected:
    void run ();
    void compute_surface_data ();
};





class MinimiseThread : public Thread {
    Q_OBJECT
    public:
    MinimiseThread (QObject *parent, DDWin *ddw);
    Molecule *molecule;
    DDWin *ddwin;

    void inline set_molecule (Molecule *mol) {molecule = mol;};


    signals:
    void ask_redraw (Molecule *mol);

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

class HapticThread : public Thread {
    Q_OBJECT
    public:
    HapticThread (QObject *parent, DDWin *ddw);
	void stop ();
    Molecule *molecule;
    DDWin *ddwin;
	Minimize *minimise;
    void inline set_molecule (Molecule *mol) {cerr <<" setting mol "<< endl; molecule = mol;};
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
    void ask_redraw (Molecule *mol);
    void ask_color_by_score (Molecule *mol);


    protected:
    void run ();
	
	private:
	int number_of_threads;
	vector <HapticWorkerThread *> worker_threads;
	vector <ForceFieldInteraction *> internal_interactions;
	queue <ForceFieldInteraction *> nonbonded_interactions;

	
	public:
		bool have_to_work_on_internal_interactions ();
		bool have_to_work_on_nonbonded_interactions ();
		bool are_atoms_finished ();
		bool are_internal_interactions_finished ();
		bool is_end_of_cycle ();
		void compute_one_internal_interaction ();	
		void compute_one_nonbonded_interaction ();
		void load_internal_interactions ();
		void load_interactions_for_new_atom ();
		void end_of_calculation_for_atom (Atom *at);
		void end_of_cycle ();
};


class HapticWorkerThread : public Thread {
    public:
		HapticWorkerThread (HapticThread *ht, int i);
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
//    void ask_redraw (Molecule *mol);
//    void ask_color_by_score (Molecule *mol);


    protected:
    void run ();

};


class HeadTrackingThread : public Thread {
    Q_OBJECT
    public:
    HeadTrackingThread (QObject *parent, DDWin *ddw);
    DDWin *ddwin;
  //  cwiid_wiimote_t *wiimote;	/* wiimote handle */


//    signals:
//    void ask_redraw (Molecule *mol);
//    void ask_color_by_score (Molecule *mol);

    cwiid_wiimote_t *wiimote;

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
  //  cwiid_wiimote_t *wiimote;	/* wiimote handle */


//    signals:
//    void ask_redraw (Molecule *mol);
//    void ask_color_by_score (Molecule *mol);

    cwiid_wiimote_t *wiimote;

    protected:
    bool is_moving;
    void run ();
    float lastx, lasty, lastd;
    signals:
	void move_camera (float x, float y, float z);	

};




#endif //THREAD_H
