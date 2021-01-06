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
#ifndef MOLECULE_H
#define MOLECULE_H


#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "maths.h"
#include "constants.h"
#include <openbabel/mol.h>
#include <openbabel/generic.h>
using namespace std;


class Database; class Fragment; class vect; class Molecule;


typedef OpenBabel::OBAtom Atom;
typedef OpenBabel::OBBond Bond;

typedef OpenBabel::OBRing Ring;
typedef OpenBabel::OBResidue Resid;

typedef OpenBabel::OBPairTemplate< double > DoubleData;
typedef OpenBabel::OBPairTemplate< bool > BoolData;
typedef OpenBabel::OBPairTemplate< int > IntData;
typedef OpenBabel::OBPairTemplate< color > ColorData;
typedef OpenBabel::OBPairTemplate< vect > VectorData;

vect find_mass_center (vector<Atom*>& invec);
color get_color (Atom *at);
void set_color (Atom *at, color col);
int get_ds (Atom *at);
void set_ds (Atom *at, int);
int get_ds (Bond *b);
void set_ds (Bond *b, int);
bool get_visible (Atom *at);
void set_visible (Atom *at, bool);
bool get_visible (Bond *bo);
void set_force (Atom *at, vect);
vect get_force (Atom *at);
void set_score (Atom *at, double);
double get_score (Atom *at);
void set_old_score (Atom *at, double);
double get_old_score (Atom *at);


bool get_selected (Atom *at);
void set_selected (Atom *at, bool);
bool get_selected (Bond *b);
void set_selected (Bond *b, bool);
bool get_sad (Atom *at);
void set_sad (Atom *at, bool);


void mend_coordinates (Molecule *mol);
int get_MMFFtype (Atom *at);
bool IsInSameRing(Atom* a, Atom* b);



double get_vdw (Atom *at);

int CountBonds (Atom *);


class Molecule : public OpenBabel::OBMol {
    public:

    Molecule ();
    int line_list, backbone_list1, backbone_list2, stick_list;
	int atoms_default_ds, bonds_default_ds;

    bool selection, multi;
    bool needs_redraw, needs_recolor;
    vect center, min_corner, max_corner;
     
    int bonded_to (Atom *, int, int);

    void ZNinit ();
    void ZNSetConformers ();
	void ZNAddHydrogens ();
    void ZNinit_atom (Atom *at);
    void ZNinit_bond (Bond *bo);
    bool ZNAddHydrogens (Atom *at);
    Atom * ZNAddAtom (Atom *at);
    Bond * ZNAddBond (Bond *bo);

    void find_limits ();
    void find_center ();
    color get_color_mw (Atom *at);
    void  set_color_mw (Atom *at);
	int get_ds_from_neighbour (Atom *at);
	int get_ds_from_neighbour (Bond *bo);
    void add_atom_bonded_to (Atom *to_add, Bond *bond, Atom *partner);
    void add_atom_bonded_to (vect coords, int atomnum, Atom *at);

    bool RemoveBond (Bond *bo);
    bool RemoveAtom (Atom *at);
//    bool RemoveHydrogen (Atom *at);

};






/*
class Molecule {
public:
    Molecule ();
//    Molecule (Molecule const& mol);
    ~Molecule ();

    Molecule & operator=(const Molecule& mol);
    int readMultiMOL2(string filename, ifstream& file, bool& continueRead);
    void copy_from (Molecule *mol);
  //  unsigned int aromatic_display_style;
    void find_residues ();
    void find_rings ();
    void find_fragments ();
    void find_limits ();
    void find_center ();
    void find_bound ();
    void number_atoms ();
    void number_bonds ();
    Bond *find_bond (Atom *at1, Atom *at2);
    void add_atom_bonded_to (Atom *to_add, Atom *partner);
    void add_atom_bonded_to (Atom *to_add, Bond *bond, Atom *partner);
    Atom *add_atom_bonded_to (vect coords, int atomnum, Atom *at);

    void del (Bond *bo);
    void del (Atom *at);

    void add (Bond *bo);
    void add (Atom *at);

    void remove (Bond *bo);
    void remove (Atom *at);

    bool close_ring (Bond *bo);
    Atom* findOtherAtom (Bond* bond, Atom* atom);
    void findBonds (Atom *atom, vector<Bond*>& out);
    vector<string> res_strings;
    vect center;
    vect min_corner;
    vect max_corner;

    void find_kekule ();

    bool needs_redraw, needs_recolor;


//private:
    bool multi;
    bool selection;
    string name, description; 
 //   unsigned int atomCount, bondCount, substructCount;
    vector<Atom*> atoms;
    vector<Bond*> bonds;
    vector<Resid*> residues;
    vector<Ring*> rings;
    vector<Fragment*> fragments;
    bool valid;
    
    int line_list, stick_list, surface_list, backbone_list1, backbone_list2;

};
*/


class Selection : public Molecule {
    public:
    Selection ();
    Atom * ZNAddAtomToSelection (Atom *at);
    Bond * ZNAddBondToSelection (Bond *bo);
    void deselect ();
    void add_mol (Molecule *mol);
    vector <Molecule *> molecules;
};





/*
class Atom {
    public:
        Atom ();
        void find_mol2_type ();
        float charge, MMFFstartcharge, vdw, score, old_score;
        string MMFFstring;
        vect coordinates;// ref_coords;
        unsigned int ID, atomicNumber, mol2Type, MMFFtype,  displayStyle;        
        string  name, substructureName, atomType;
        bool visited;
        bool isProtein, isLigand, inBackBone, isWater, mol2BACKBONE, mol2DICT, mol2DIRECT, mol2ESSENTIAL, no2bond; //no2bond: for kekule structures this atom can't get a double bond (like pyrrole N)
        unsigned int number, substructureNumber;
        Resid *residue;
        Fragment *fragment;
        bool isSp2 ();
        unsigned int bonded_to (int bondtype=-2, int atnumb=-2);
        unsigned int bonded_to (string MMFFstring);
        bool bonded_to (Atom *at);
        vector<Ring*> in_ring;
        vector<Atom*> bound;
        vector<Bond*> bonds;
        color col;
        void getVdw ();
        void set_color_mw ();
        color get_color_mw ();
   //     void set_color_charge ();
   //     void set_color_score ();
        bool is_aromatic ();
        bool is_in_ring (Ring *);
        bool sphere_already_drawn, visible, selected;
        vect force;


};
*/
/*
class Resid {
public:
    Resid ();
    int backbone_style;
    bool backbone_visible;
    color backbone_color;
    void refine_backbone ();
    vector<Atom*> atoms;
    vector<Bond*> bonds;
    Atom *CA;
    Atom *NB;
    Atom *CB;
    bool hasCA, hasNB, hasCB;
    Molecule *molecule;

    vector <vect> backbone_list;
    
 //   Resid *prev;
 //   Resid *foll;
    unsigned int number;
    string name;   
};

*/
class Fragment {
public:
    float *rotation_quaternion;
    vect translation_vector;
    vect center;
    Fragment ();
    vector<Atom*> atoms;
    vector<float *> internal_coordinates;
    Molecule *molecule;

    unsigned int number;
    string name;   
};

/*
class Bond {
public:
    Bond ();
    unsigned int atomID [2], fr_dir;
    Atom* atomPTR [2];
    Bond *fr_parent;
    vector <Ring*> in_ring;
    unsigned int number,mol2Type, kekule, displayStyle;
//    bool hasCoplanar0, hasCoplanar1;
    bool fr_visited;
    Atom *coplanarAtom0, *coplanarAtom1;
    bool is_in_ring (Ring *);
    bool is_aromatic ();
    bool is_rotatable ();
    bool conjugated_to (unsigned int order, Ring *ring);
    void set_order (int ord);
};
*/
/*
class Ring {
public: 
    int kekule_type; //marks the order in which kekule structures are found
    Ring ();
    vector<Atom*> atoms;
    vector<Bond*> bonds;
    Atom *alpha_atom;
    bool aromatic, IM_CAT, N5ANION;
    vect center;
    unsigned int count (int an);
    void test_aromaticity ();
    void test_kekule_aromaticity ();
    void find_kekule ();
    void find_kekule_pseudoaromatic ();
    void complete_kekule ();
    bool all_aromatic ();

};
*/


    vect find_mass_center (vector<Atom*>& atv);
#endif

