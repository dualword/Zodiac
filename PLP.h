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
#ifndef PLP_H
#define PLP_H



#include "molecule.h"
#include <vector>
#include "cutoffGrid.h"
#include "constants.h"
#include "FF.h"


class PLPInteraction : public ForceFieldInteraction
{
    public:
	Atom *at1, *at2;
    vect root;
    int I, J, type;
    float value ();
    
}; 

class PLP : public ForceField {

public:
    PLP ();

 //   void type ();

    void clear_nonbonded_interactions ();
    void load_mol (Molecule *mol);
    void load_environment (vector<Molecule *> envir, Molecule *mol);
    void load_internal_interactions ();
    void load_nonbonded_interactions ();
    void update ();

    void initialize (Molecule *mol, vector<Molecule *> envir);
    inline double compute_total_energy () {return 0.;};
    void initialize_mol (Molecule *mol);
    int getPLPtype(Atom* at);
    int getPLPinteractiontype (int I, int J);

    void compute_forces ();
 //   void compute_force (PLPInteraction *plpint);
 //   float compute_interaction (PLPInteraction *plpint);


private:
    vector<PLPInteraction *> NBInteractions;

};

#endif
