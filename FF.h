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
#include "molecule.h"
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
    virtual void set_forces (bool score = FALSE);
    private:
};


class ForceField {
public:
    ForceField ();
    bool is_initialised;





    virtual ~ForceField();
    void clear ();
    void load_mol (Molecule *mol);
    void load_environment (vector<Molecule *> envir, Molecule *mol);
    void update ();

    void initialize_internal (Molecule *mol, vector<Molecule *> envir);
    void initialize_interaction (Molecule *mol, vector<Molecule *> envir);
    void initialize (Molecule *mol, vector<Molecule *> envir);
    virtual double compute_total_energy ()=0;
    float total_energy;

//    int getPLPtype(Atom* at);
  //  int getPLPinteractiontype (int I, int J);

    virtual void compute_forces ()=0;

    virtual void clear_nonbonded_interactions ()=0;
    virtual void load_nonbonded_interactions ()=0;
	virtual void load_nonbonded_interactions_for_atom (Atom *at, queue <ForceFieldInteraction*> *q)=0;
    virtual void load_internal_interactions (vector <ForceFieldInteraction*> *v)=0;


    Molecule *target_mol;

protected:

    vector <Atom *> environment;





    cutoffGrid<Atom*>* near_grid;
    cutoffGrid<Atom*>* far_grid;
};














#endif
