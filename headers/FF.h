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

#ifndef FF_H
#define FF_H
#include "ZNmolecule.h"
#include "constants.h"
#include "cutoffGrid.h"


class ForceFieldInteraction {
    public:
    ForceFieldInteraction ();
    Atom *at1, *at2;
    virtual float value () = 0;
    virtual double derive (Atom *at, double *val = NULL);
	virtual double derive_y (Atom *at, double *val = NULL);
	virtual double derive_z (Atom *at, double *val = NULL);
    virtual void set_forces (bool score = FALSE,  double mult = 1.);
	inline virtual bool isHbond () {return false;};
	inline virtual bool isElectrostatic () {return false;};
    private:
};


class ElasticRestrain : public ForceFieldInteraction {
public:
	ElasticRestrain ();
    double k;
	double dist0;
    float value ();
};


class ForceField {
public:
    ForceField ();
    bool is_initialised;





    virtual ~ForceField();
    void clear ();
    void load_mol (ZNMolecule *mol);
    void load_environment (vector<ZNMolecule *> envir, ZNMolecule *mol, vect cent = vect (0., 0., 0.), double rad = 0);
	virtual void load_grids (vect cent, double rad) {};
    void update ();

    void initialize_internal (ZNMolecule *mol, vector<ZNMolecule *> envir);
    void initialize_interaction (ZNMolecule *mol, vector<ZNMolecule *> envir, vect cent = vect (0., 0., 0.), float rad = 0.);
    void initialize (ZNMolecule *mol, vector<ZNMolecule *> envir);
    virtual double compute_total_energy ()=0;
    float total_energy;

//    int getPLPtype(Atom* at);
  //  int getPLPinteractiontype (int I, int J);

    virtual void compute_forces ()=0;

    virtual void clear_nonbonded_interactions ()=0;
    virtual void load_nonbonded_interactions ()=0;
	virtual void load_nonbonded_interactions_for_atom (Atom *at, queue <ForceFieldInteraction*> *q) {};
    virtual void load_internal_interactions (vector <ForceFieldInteraction*> *v) {};

	vector <ForceFieldInteraction *> restrains;
    ZNMolecule *target_mol;

protected:

    vector <Atom *> environment;





    cutoffGrid<Atom*>* near_grid;
    cutoffGrid<Atom*>* far_grid;
};














#endif
