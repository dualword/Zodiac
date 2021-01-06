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

#include <QThread>
#include "graphical_object.h"
#include "minimize.h"
#include "ddwin.h"
#include <openbabel/forcefield.h>


class DDWin; class Minimize;
class Thread : public QThread   {
    Q_OBJECT
    public:
    Thread (QObject *parent);

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

 

class HapticThread : public Thread {
    Q_OBJECT
    public:
    HapticThread (QObject *parent, DDWin *ddw);
    Molecule *molecule;
    DDWin *ddwin;

    void inline set_molecule (Molecule *mol) {molecule = mol;};
    bool is_running;
    int haptic_dof_mode;
    bool automove, color_by_score;
    void initialise (Minimize * min);
    void stop ();

    signals:
    void ask_redraw (Molecule *mol);
    void ask_color_by_score (Molecule *mol);


    protected:
    void run ();
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


#endif //THREAD_H
