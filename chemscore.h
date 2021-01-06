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


#ifndef CHEMSCORE_H
#define CHEMSCORE_H



#include "molecule.h"
#include <vector>
#include "cutoffGrid.h"
#include "constants.h"
#include "FF.h"
#include <openbabel/obiter.h>


float f (float x, float x1, float x2);

class ChemscoreHBInteraction : public ForceFieldInteraction
{
    public:
	Atom *at1, *at2;
    Atom * root;
    float value ();
};


class ChemscoreLiInteraction : public ForceFieldInteraction
{   
    public:
	Atom *at1, *at2;
    bool metal;
    float value ();
};


class ChemscoreClInteraction : public ForceFieldInteraction
{
    public:
	Atom *at1, *at2;
    int type;
    float value ();
};






class Chemscore : public ForceField {

public:

    Chemscore ();
    void clear_nonbonded_interactions ();
    void load_mol (Molecule *mol);
    void load_environment (vector<Molecule *> envir, Molecule *mol);
    void load_internal_interactions ();
    void load_nonbonded_interactions ();
    void update ();


    inline double compute_total_energy () {return 0.;};
    void initialize_mol (Molecule *mol);


    void compute_forces ();
//    void compute_force (ChemscoreHBInteraction *hbint);
//    void compute_force (ChemscoreLiInteraction *liint);
//    void compute_force (ChemscoreClInteraction *clint);
 //   float compute_HB_interaction (ChemscoreHBInteraction *hbint);
 //   float compute_Li_interaction (ChemscoreLiInteraction *hbint);
 //   float compute_Cl_interaction (ChemscoreClInteraction *clint);

private:
    vector<ChemscoreClInteraction *> ClInteractions;
    vector<ChemscoreHBInteraction *> HBInteractions;
    vector<ChemscoreLiInteraction *> LiInteractions;
  
    int getChemscoretype (Atom *at);


};

#endif
