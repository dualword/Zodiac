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


#include "constants.h"
#include "obabel_includes.h"
#include <QMutex>
#include <QReadWriteLock>



//molecule data
const unsigned int MINIMISATION_DATA_INDEX = OBGenericDataType::CustomData0;
const unsigned int MOLECULAR_FRAGMENT_DATA_INDEX = OBGenericDataType::CustomData1;
const unsigned int MOLECULAR_DISPLAY_DATA_INDEX = OBGenericDataType::CustomData2;
const unsigned int GEOMETRY_DATA_INDEX = OBGenericDataType::CustomData3;


//atom data
const unsigned int ATOMIC_FRAGMENT_DATA_INDEX = OBGenericDataType::CustomData0;
const unsigned int ATOMIC_DISPLAY_DATA_INDEX = OBGenericDataType::CustomData1;

//bond data
const unsigned int BOND_DISPLAY_DATA_INDEX = OBGenericDataType::CustomData0;

//residue data
const unsigned int BACKBONE_DATA_INDEX = OBGenericDataType::CustomData0;


typedef OpenBabel::OBAtom Atom;
typedef OpenBabel::OBBond ZNBond;
class Fragment;


using namespace std;

typedef struct {
	vect v1;
	vect v2;
	float perc;
} HBond;

class rotorConnection  {
public:
	rotorConnection () {bond = NULL; from_atom = NULL; to_atom = NULL; to_fragment = NULL; to_keep = false; dihedral = 0.;}
	~rotorConnection () {};
	ZNBond *bond;
	Atom *from_atom, *to_atom;
	Fragment *to_fragment;
	bool to_keep;
	double dihedral;
	void set_dihedral (double angle);
};




class SurfVertex
	{
    public:
		  SurfVertex () : col (color (255,255,255)) {};
		inline vect& GetVector () {return coordinates;};
		vect coordinates;
		color col;
		vect normal;
		int n;
	}; 



class Database;  class vect; class ZNMolecule;



typedef OpenBabel::OBRing Ring;
typedef OpenBabel::OBResidue Resid;


//typedef OpenBabel::OBPairTemplate <QMutex *>	MutexData;


class MinimisationData:public OBGenericData  {
public:
	MinimisationData ();
	MinimisationData (const MinimisationData &dat) 
//        : OBGenericData( *((OBGenericData::OBGenericData *)&dat) )  
        : OBGenericData( dat ) 
    {
        _force = dat._force; 
        _back_force = dat._back_force; 
        _score = dat._score; 
        _back_score = dat._back_score;
    };
	virtual OBGenericData* Clone(OBBase* parent) const	{return new MinimisationData(*this);};
	void set_force_value (vect force) {_force = force;};
	vect get_force_value () {return _force;};
	void set_back_force_value (vect force) 
{
//#pragma omp parallel 
//#pragma omp critical
{
_back_force = force;
}
};
	vect get_back_force_value () {return _back_force;};
	void set_score_value (double score) {_score = score;};
	double get_score_value () {return _score;};
	void set_back_score_value (double score) 
{
//#pragma omp parallel 
//#pragma omp critical
{
_back_score = score;
}
};
	double get_back_score_value () {return _back_score;};
private:
	vect _force;
	vect _back_force;
	double _score, _back_score;
};


class AtomicFragmentData:public OBGenericData  {
public:	
	AtomicFragmentData ();
	AtomicFragmentData (const AtomicFragmentData &dat)
//        : OBGenericData( *(OBGenericData::OBGenericData *)&dat )   
        : OBGenericData( dat )   
    {_fragment= NULL; _kinematic_chain_visited = false;};
	virtual OBGenericData* Clone(OBBase* parent) const	{return new AtomicFragmentData(*this);};
	Fragment *get_fragment () {return _fragment;};
	void set_fragment (Fragment *frag) {_fragment = frag;};
	bool get_kinematic_chain_visited () {return _kinematic_chain_visited;};
	void set_kinematic_chain_visited (bool kcv) {_kinematic_chain_visited = kcv;};
private:
	Fragment *_fragment;
	bool _kinematic_chain_visited;
};



class MolecularFragmentData:public OBGenericData {
public:
	MolecularFragmentData ();
	MolecularFragmentData (const MolecularFragmentData &dat)
//        : OBGenericData( *(OBGenericData::OBGenericData *)&dat )   
        : OBGenericData( dat )   
    {_fragments.clear (); _connections.clear ();};
	virtual OBGenericData* Clone(OBBase* parent) const	{return new MolecularFragmentData(*this);};
	Fragment *get_fragment (int i) {int n = _fragments.size (); if (i<n && i>=0) return _fragments[i]; else return NULL;}
	void set_fragment_list (vector <Fragment *> fragment_list) {_fragments = fragment_list;};
	vector <Fragment *> &get_fragment_list () {return _fragments;};
	void add_rotor (rotorConnection *connection) {_connections.push_back (connection);}
	void clear_rotors () {_connections.clear ();}
	vector <rotorConnection *> &get_rotors () {return _connections;}; 
private:
	vector <Fragment *> _fragments;
	vector <rotorConnection *> _connections;
};





class GeometryData:public OBGenericData {
public:
	GeometryData ();
	GeometryData(const GeometryData &dat) 
//        : OBGenericData( *(OBGenericData::OBGenericData *)&dat )  
        : OBGenericData( dat )  
    {_center = dat._center; _min_corner = dat._min_corner; _max_corner = dat._max_corner; lock = new QReadWriteLock; _backbone_perceived = false; _fragments_perceived = false;}
	virtual OBGenericData* Clone(OBBase* parent) const	{return new GeometryData(*this);};
	vect get_center () {return _center;};
	void set_center (vect v) {_center = v;};
	vect get_min_corner () {return _min_corner;};
	vect get_max_corner () {return _max_corner;};
	void set_min_corner (vect v) {_min_corner = v;};
	void set_max_corner (vect v) {_max_corner = v;};
	QReadWriteLock *lock;
		bool _backbone_perceived, _fragments_perceived;
	
private:
	vect _center, _min_corner, _max_corner;


};

class MolecularDisplayData:public OBGenericData {
public:
	MolecularDisplayData ();
	MolecularDisplayData (const MolecularDisplayData &dat)
//        : OBGenericData( *(OBGenericData::OBGenericData *)&dat )   
        : OBGenericData( dat )   
    {
		_atom_display_style = dat._atom_display_style;
		_bond_display_style = dat._bond_display_style;
		_backbone_display_style = dat._backbone_display_style;
		_line_list = dat._line_list; _backbone_list1 = dat._backbone_list1;
		_backbone_list2 = dat._backbone_list2; _stick_list = dat._stick_list; _needs_redraw = true; _needs_backbone_redraw = true;
		_clippable = dat._clippable;
	};
	virtual OBGenericData* Clone(OBBase* parent) const	{return new MolecularDisplayData(*this);};
	int get_atoms_display_style () {return _atom_display_style;};
	void set_atoms_display_style (int i) {_atom_display_style = i;};
	int get_bonds_display_style () {return _bond_display_style;};
	void set_bonds_display_style (int i) {_bond_display_style = i;};
	int get_backbone_display_style () {return _backbone_display_style;};
	void set_backbone_display_style (int i) {_backbone_display_style = i;};
	int get_line_list () {return _line_list;};
	int get_backbone_list1 () {return _backbone_list1;};
	int get_backbone_list2 () {return _backbone_list2;};
	int get_stick_list () {return _stick_list;};
	void set_line_list (int i) {_line_list = i;};
	void set_backbone_list1 (int i) {_backbone_list1 = i;};
	void set_backbone_list2 (int i) {_backbone_list2 = i;};
	void set_stick_list (int i) {_stick_list = i;};
	void set_needs_redraw (bool b) {_needs_redraw =  b;};
	bool get_needs_redraw () {return _needs_redraw;}; 
	void set_needs_backbone_redraw (bool b) {_needs_backbone_redraw =  b;};
	bool get_needs_backbone_redraw () {return _needs_backbone_redraw;}; 

	void set_clippable (bool b) {_clippable =  b;};
	bool get_clippable () {return _clippable;}; 
	

	
private:
	int _atom_display_style;
	int _bond_display_style;
	int _backbone_display_style;
	int _line_list; 
	int _backbone_list1;
	int _backbone_list2; 
	int _stick_list;
	bool _needs_redraw, _clippable, _needs_backbone_redraw;
	
};

class AtomicDisplayData:public OBGenericData {
public:
	AtomicDisplayData ();
	AtomicDisplayData(const AtomicDisplayData &dat) 
//        : OBGenericData( *(OBGenericData::OBGenericData *)&dat )  
        : OBGenericData( dat )  
    {
        _visible = dat._visible; _selected = dat._selected; _sphere_already_drawn = false; _display_style = dat._display_style;
		_color = dat._color;
	}
	virtual OBGenericData* Clone(OBBase* parent) const	{return new AtomicDisplayData(*this);};
	bool get_visible () {return _visible;};
	bool get_selected () {return _selected;};
	bool get_sad () {return _sphere_already_drawn;};
	int get_display_style () {return _display_style;};
	color get_color () {return _color;};
	
	void set_visible (bool b) {_visible = b;};
	void set_display_style (int i) {_display_style = i;};
	void set_selected(bool b) {_selected = b;};
	void set_sad (bool b) {_sphere_already_drawn = b;}; 
	void set_color (color c) {_color = c;};
private:
	bool _visible, _selected, _sphere_already_drawn;
	int _display_style;
	color _color;
};


class BondDisplayData:public OBGenericData {
public:
	BondDisplayData ();
	BondDisplayData(const BondDisplayData &dat) 
//        : OBGenericData( *(OBGenericData::OBGenericData *)&dat )   
        : OBGenericData( dat )   
    {
        _display_style = dat._display_style;
    }
	virtual OBGenericData* Clone(OBBase* parent) const	{return new BondDisplayData(*this);};
	int get_display_style () {return _display_style;};
	void set_display_style (int i) {_display_style = i;};

private:
	int _display_style;

};


class BackboneData:public OBGenericData {
public:
	BackboneData ();

	
	color col;
	vector <vect> backbone_points;
	vect backbone_dir;
	unsigned int ss;
	
};


vect get_coordinates (SurfVertex *v);

//access atomic data

vect get_coordinates (Atom *at);
void set_coordinates (Atom *at, vect v);
void sum_to_coordinates (Atom *at, vect v);

vect find_mass_center (vector<Atom*>& invec);
color get_color (Atom *at);
color get_color (Resid *res);
void set_color (Resid *res, color c);
void set_color (Atom *at, color col);
int get_ds (Atom *at);
void set_ds (Atom *at, int);
bool get_visible (Atom *at);
void set_visible (Atom *at, bool);
void set_force (Atom *at, vect);
vect get_force (Atom *at);
void lock_force_mutex (Atom *at);
void unlock_force_mutex (Atom *at);
void set_back_force (Atom *at, vect);
vect get_back_force (Atom *at);
void flush_forces (Atom *at);
void flush_scores (Atom *at);
void set_score (Atom *at, double);
double get_score (Atom *at);
void set_back_score (Atom *at, double);
double get_back_score (Atom *at);




bool get_selected (Atom *at);
void set_selected (Atom *at, bool);
bool get_sad (Atom *at);
void set_sad (Atom *at, bool);



bool get_kinematic_chain_visited (Atom *at);
void set_kinematic_chain_visited (Atom *at, bool b);
void set_fragment (Atom *at, Fragment *fr);
Fragment *get_fragment (Atom *at);


Resid *get_previous_residue (Resid *res, int gap = 1);
Resid *get_following_residue (Resid *res, int gap = 1);

Atom *get_O (Resid *res);
Atom *get_CA (Resid *res);
Atom *get_C (Resid *res);
Atom *get_N (Resid *res);
double get_phi (Resid *res);
double get_psi (Resid *res);
bool is_helix (Resid *res);
bool is_sheet (Resid *res);
bool is_random (Resid *res);

void find_secondary_structure(ZNMolecule *mol);
void find_secondary_structure (Resid *res);

void color_backbone_ss (ZNMolecule *mol, color hel, color sheet, color random);
void color_backbone_color (ZNMolecule *mol, color c);
void set_ss (Resid *res, int);
int get_ss (Resid *res);
void set_backbone_direction (Resid *res, vect v);
void set_backbone_points (Resid *res, vector <vect> points);
vect get_backbone_direction (Resid *res);
vect get_start_reference (Resid *res);
vect get_end_reference (Resid *res);

vector <vect> get_backbone_points (Resid *res);
void find_backbone_data (ZNMolecule *mol);
void set_backbone_style (Atom *at, int n);
void find_backbone_points (Resid *res) ;
void find_backbone_direction (Resid *res);
vector <vect> smooth_list (vector <vect> lis);
void add_guide_points_to_backbone (Resid *res, vector <vect> &lis);
int get_MMFFtype (Atom *at);
bool IsInSameRing(Atom* a, Atom* b);


double get_vdw (Atom *at);

int CountBonds (Atom *);
bool is_polar (Atom *);
Atom* root_at (Atom *);


//access bond data
int get_ds (ZNBond *b);
void set_ds (ZNBond *b, int);
bool get_visible (ZNBond *bo);
bool get_selected (ZNBond *b);
void set_selected (ZNBond *b, bool);

//access molecule data
void set_needs_redraw (ZNMolecule *, bool b);
bool get_needs_redraw (ZNMolecule *);
void set_needs_backbone_redraw (ZNMolecule *, bool b);
bool get_needs_backbone_redraw (ZNMolecule *);
void set_clippable (ZNMolecule *, bool b);
bool get_clippable (ZNMolecule *);

void set_fragment_list (ZNMolecule *mol, vector <Fragment *> fragment_list);
void add_rotor (ZNMolecule *mol, rotorConnection *connection);
void clear_rotors (ZNMolecule *mol);
void build_kinematic_chain (ZNMolecule *mol);
void mend_coordinates (ZNMolecule *mol);
vector <rotorConnection *> &get_rotors (ZNMolecule *mol);
void set_dihedrals (ZNMolecule *mol, vector <double> dihedrals);
void build_molecule_from_dofs (ZNMolecule *mol, vector <float> dihedrals, bool move = false);

vect get_root_fragment_center (ZNMolecule *mol) ;
vector <Fragment *> get_fragments(ZNMolecule *mol);
void initialise_dihedrals (ZNMolecule *mol, bool complete = false);
vector <double> get_dihedrals(ZNMolecule *mol);
void finalise_molecule(ZNMolecule *mol);
void set_center (ZNMolecule *mol, vect v);
void set_min_corner (ZNMolecule *mol, vect v);
void set_max_corner (ZNMolecule *mol, vect v);
vect get_min_corner (ZNMolecule *mol);
vect get_max_corner (ZNMolecule *mol);
vect get_center (ZNMolecule *mol);

void find_limits (ZNMolecule *mol);
void find_center (ZNMolecule *mol);

int get_atoms_display_style (ZNMolecule *mol);
int get_bonds_display_style (ZNMolecule *mol);
int get_backbone_display_style (ZNMolecule *mol);

void set_atoms_display_style (ZNMolecule *mol, int i);
void set_bonds_display_style (ZNMolecule *mol, int i);
void set_backbone_display_style (ZNMolecule *mol, int i);

void lock_geometry_for_read (ZNMolecule *mol);
void lock_geometry_for_write (ZNMolecule *mol);
void unlock_geometry (ZNMolecule *mol);


float very_fast_RMSD (ZNMolecule *a, ZNMolecule *b);
int get_line_list (ZNMolecule *mol);
int get_backbone_list1 (ZNMolecule *mol);
int get_backbone_list2 (ZNMolecule *mol);
int get_stick_list (ZNMolecule *mol);
void set_display_lists (ZNMolecule *mol, int ll, int bl1, int bl2, int sl);
void translate_molecule (ZNMolecule *mol, vect v);
void rotate_molecule (ZNMolecule *mol, quaternion q, vect v);
void rotate_atom (Atom *at, quaternion q, vect v);
bool is_bsheet (Resid *res1, Resid *res2) ;
Resid *find_sheet_partner (Resid *res, ZNMolecule *mol) ;

void set_backbone_perceived (ZNMolecule *mol, bool b);
bool get_backbone_perceived (ZNMolecule *mol);
void set_fragments_perceived (ZNMolecule *mol, bool b);
bool get_fragments_perceived (ZNMolecule *mol);

void molecule_has_changed (ZNMolecule *mol);

ZNMolecule *sum (ZNMolecule *mol1, ZNMolecule *mol2);

class ZNMolecule : public OpenBabel::OBMol {
    public:

    ZNMolecule ();
	ZNMolecule (const ZNMolecule &);
	ZNMolecule &operator =(const ZNMolecule &mol);
	ZNMolecule &operator +=(const ZNMolecule &mol);	

	
    bool selection, multi;
    bool needs_recolor;

     
    int bonded_to (Atom *, int, int);

    void ZNinit ();
    void ZNSetConformers ();
	void ZNAddHydrogens (double ph = 7.4);
    void ZNinit_atom (Atom *at);
    void ZNinit_bond (ZNBond *bo);
	void ZNinit_residue (Resid *re);
    bool ZNAddHydrogens (Atom *at);
    Atom * ZNAddAtom (Atom *at);
    ZNBond * ZNAddBond (ZNBond *bo);


    color get_color_mw (Atom *at);
    void  set_color_mw (Atom *at);
	int get_ds_from_neighbour (Atom *at);
	int get_ds_from_neighbour (ZNBond *bo);
    void add_atom_bonded_to (Atom *to_add, ZNBond *bond, Atom *partner);
    void add_atom_bonded_to (vect coords, int atomnum, Atom *at);

    bool RemoveBond (ZNBond *bo);
    bool RemoveAtom (Atom *at);
	
	


};

class Selection : public ZNMolecule {
    public:
    Selection ();
    Atom * ZNAddAtomToSelection (Atom *at);
    ZNBond * ZNAddBondToSelection (ZNBond *bo);
    void deselect ();
    void add_mol (ZNMolecule *mol);

	void extend_to_residue ();
	void extend_to_fragment ();
	void extend_to_neighbours ();
	void select_atom (Atom *at);
	void select_atoms (vector <Atom *> to_add);
	vector <ZNMolecule *>&get_molecules () {return molecules;};
	private:
	    vector <ZNMolecule *> molecules;
};


class Fragment {
public:
	Fragment ();
	void set_rotation (quaternion q) {_rotation_quaternion = q;};
	quaternion get_rotation () {return _rotation_quaternion;};
	void set_translation (vect v) {_translation_vector = v;};
	vect get_tranlation () {return _translation_vector;};
	vect get_center () {return _center;};
	void set_number (unsigned int n) {_number = n;};
	unsigned int get_number () {return _number;};
	unsigned int number_of_connections () {return _rotor_connections.size ();};
	rotorConnection *get_rotor_connection (unsigned int i) {return _rotor_connections[i];};
	void add_atom (Atom *at) {_atoms.push_back (at); _reference_coordinates.push_back (get_coordinates(at)); _last_coordinates.push_back (get_coordinates (at));};
	void add_rotor_connection (rotorConnection *rc) {_rotor_connections.push_back (rc);};
	unsigned int size () {return _atoms.size ();};
	void set_visited (bool b) {_visited = b;};
	bool get_visited () {return _visited;};
	vector <Atom *>& get_atoms_list () {return _atoms;}
	void clean_connections ();
	void quaternion_rotate (vect center, quaternion quat);
	void find_children () {_children.clear (); _children = all_children();}
	void set_center () {_center = find_mass_center (_atoms);}
	void reset_coordinates (bool complete = false) {
		for (unsigned int i = 0; i < _atoms.size (); i++) {
			_last_coordinates[i] = _reference_coordinates[i];
			if (complete) set_coordinates (_atoms[i], _reference_coordinates[i]);
			}
		};
	vector <Fragment *> get_children () {return _children;};
	vector <Fragment *> all_children ();	
	vector <Fragment *> get_first_order_children ();
	vector <rotorConnection *> get_connections () {return _rotor_connections;};
private:
	vector <rotorConnection *> _rotor_connections;
    quaternion _rotation_quaternion;
    vect _translation_vector;
    vect _center;
	bool _visited;
    vector<Atom*> _atoms;
	vector<vect> _reference_coordinates, _last_coordinates;
	vector<Fragment *> _children;
    unsigned int _number;   

};

/*
class Ring {
public: 
    int kekule_type; //marks the order in which kekule structures are found
    Ring ();
    vector<Atom*> atoms;
    vector<ZNBond*> bonds;
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

