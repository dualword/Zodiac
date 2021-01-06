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


#ifndef BUILDER_H
#define BUILDER_H

#include "molecule.h"
#include "maths.h"
#include "constants.h"
#include "ZNdata.h"
#include "ddwin.h"


class DDWin;
class vect; 
class Builder {
    public:
    Builder (DDWin *ddw);
    
    void set_magic_pencil_atomic_number (int atomn);
    void save_Hs (Atom *at, vector <Atom *> &prev_ats, vector <Bond *> &prev_bonds);
    void add_Hs (Atom *at, vector <Atom *> &atoms, vector <Bond *> &bonds);
    Atom* new_atom (vect coord, int atomn);
    Bond *new_bond (Atom *at1, Atom *at2, int order=1);
    void add_atom (int atnum);
    void add_bond (int order);
    Atom * last_magic_pencil_atom, *start_magic_pencil_atom;
    void set_aromatic (Ring *ring);
    void set_non_aromatic (Ring *ring);
    void set_bond (Bond *bo, int order);
    void mutate_atom_to (Atom *at, unsigned int atomnum = 6 );
    void redefine_mol (Molecule * mol);
    Atom * add_atom_bonded_to (vect coord, int atmnum, Atom *at);
    int magic_pencil_atomic_number;

    void delete_atom (Atom *at);
    void delete_bond (Bond *bo);

    void delete_Hs (Atom *at);
    void delete_Hs (Molecule *mol);


    DDWin * ddwin;

    bool find_target ();

  
    void add_H (Atom *at,  vector <Atom *> &add_atoms, vector <Bond *> &add_bonds);  
    void delete_Hs (Atom *at, vector <Atom *> &del_atoms, vector <Bond *> &del_bonds); 

    void one_H (Atom *center, vect &coord2);
    void two_tetrahedral_coordinates (Atom *root1, Atom *root2, Atom *center, vect &coord2, vect &coord3);
    void three_coordinates (Atom * root, Atom *center, vect &coord2, vect &coord3);
    void four_coordinates (Atom * root, Atom *center, vect &coord2, vect &coord3, vect &coord4);




    private:
    int mcounter;
    
};



#endif
