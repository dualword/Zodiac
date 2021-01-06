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
#ifndef MMFF_H
#define MMFF_H

#include "molecule.h"
#include <vector>
#include "cutoffGrid.h"
#include "constants.h"
#include "FF.h"
#include <QFile>
#include <QTextStream>

typedef struct {
	int type, I, J;
    double increment;
} MMFFciParameter;

typedef struct {
    string parent, strin;
} MMFFHParameter;

typedef struct {
	string strin, description;
    int number;
} MMFFatypeParameter;


typedef struct {
	int type, atomicNumber, valence, crd, mltb;
    bool lonePairs, linear, aromatic, multipleBond;
} MMFFatomParameter;


typedef struct {
   string old_type, arom_type;
    int at_number, ring_size, L5;
    bool IMCAT, N5ANION;
} MMFFaromaticParameter;


typedef struct {
	int type, I, J, K;
    double ka, theta0;
} MMFFabParameter;


class MMFFbsInteraction : public ForceFieldInteraction {
    public:

    double r0, kb;
    int type;
    float value ();
};


class MMFFabInteraction : public ForceFieldInteraction {
    public:
 //   MMFFabInteraction ();
	Atom *at3;
    double ka, theta0;
    int type;
    bool linear;
    float value ();
    void set_forces ();
};


class MMFFsbInteraction : public ForceFieldInteraction {
    public:
	Atom *at3;
    double theta0, kijk, kkji, r0ij, r0kj;
    float value ();
    void set_forces ();
};

class MMFFopInteraction : public ForceFieldInteraction {
    public:
	Atom *at3, *at4;
    double koop;
    float value ();
    void set_forces ();
};


class MMFFtoInteraction : public ForceFieldInteraction {
    public:
	Atom *at3, *at4;
    double v1, v2, v3;
    int type;
    float value ();
    void set_forces ();
};

class MMFFvwInteraction : public ForceFieldInteraction {
    public:
    double e, r0;
    float value ();
};

class MMFFelInteraction : public ForceFieldInteraction {
    public:
    double scale;
    float value ();
};



typedef struct {
	int type, I, J;
    double kb, r0;
} MMFFbsParameter;



typedef struct {
	int type, I, J, K;
    double kijk, kkji;
} MMFFsbParameter;

typedef struct {
	int I, J, K, L;
    double koop;
} MMFFopParameter;




typedef struct {
	int type, I, J, K, L;
    double V1, V2, V3;
} MMFFtoParameter;





typedef struct {
	int type, I, J;
    double alpha, N, A, G;
    int DA;
} MMFFvwParameter;


class MMFF  : public ForceField{

public:
    MMFF ();

  //  void type ();


    void clear_nonbonded_interactions ();
    void clear_internal_interactions ();


    void load_internal_interactions ();
    void load_nonbonded_interactions ();



    double compute_total_energy ();
    void initialize_mol (Molecule *mol);
    void compute_partial_charges (Molecule *mol);
    void get_strings_mol (Molecule *mol);
    int getMMFFtype(Atom* at);
    string getMMFFcarbonstring (Atom* at);
    string getMMFFeteroatomstring (Atom* at);
    string getMMFFHstring (Atom* at);
    double getMMFFstartcharge (Atom* at);
    double compute_charge (Atom *at);


    void compute_forces ();

 //   void compute_total_van_der_waals_force (vector<float>& force, vector<float>& torque);
 //   void compute_total_electrostatic_force (vector<float>& force, vector<float>& torque);


    double compute_bond_stretchings ();
    double compute_angle_bendings ();
    double compute_stretch_bend_interactions ();
    double compute_out_of_plane_bendings ();
    double compute_electrostatic_interactions ();
    double compute_torsion_interactions ();
    double compute_van_der_waals_interactions ();

//    float compute_total_van_der_waals_energy ();
//    float compute_total_electrostatic_energy ();

  //  Molecule *target_mol;


private:
  //  vector <Atom *> environment;
  //  cutoffGrid<Atom*>* env_el_grid;
  //  cutoffGrid<Atom*>* env_vw_grid;


    vector<float> compute_electrostatic_force_vector (MMFFelInteraction *elint); 
    vector<float> compute_van_der_waals_force_vector (MMFFvwInteraction *vwint); 



    void compute_bond_stretching_forces ();
    void compute_van_der_waals_forces ();
    void compute_angle_bending_forces ();
    void compute_stretch_bend_interaction_forces ();
    void compute_out_of_plane_bending_forces ();
    void compute_electrostatic_forces ();
    void compute_torsion_forces ();
    void compute_nonbonded_van_der_waals_forces ();
    void compute_nonbonded_electrostatic_forces ();




/*
    double compute_bond_stretching (MMFFbsInteraction *bsint);
    double compute_van_der_waals_interaction (MMFFvwInteraction *vwint);
    double compute_angle_bending (MMFFabInteraction *abinter);
    double compute_stretch_bend_interaction (MMFFsbInteraction *sbint);
    double compute_out_of_plane_bending (MMFFopInteraction *opint);
    double compute_electrostatic_interaction (MMFFelInteraction *elint);
    double compute_torsion_interaction (MMFFtoInteraction *toint);
*/





/*
    void compute_bond_stretching_force (MMFFbsInteraction *bsint);
    void compute_van_der_waals_force (MMFFvwInteraction *vwint);
    void compute_angle_bending_force (MMFFabInteraction *abinter);
    void compute_stretch_bend_interaction_force (MMFFsbInteraction *sbint);
    void compute_out_of_plane_bending_force (MMFFopInteraction *opint);
    void compute_electrostatic_force (MMFFelInteraction *elint);
    void compute_torsion_force (MMFFtoInteraction *toint);
*/


/*
    float compute_bond_stretching_force_module (MMFFbsInteraction *bsint);
    float compute_van_der_waals_force_module (MMFFvwInteraction *vwint);
    float compute_angle_bending_force_module (MMFFabInteraction *abinter);
    float compute_stretch_bend_interaction_force_module (MMFFsbInteraction *sbint);
    float compute_out_of_plane_bending_force_module (MMFFopInteraction *opint);
    float compute_electrostatic_force_module (MMFFelInteraction *elint);
    float compute_torsion_force_module (MMFFtoInteraction *toint);
*/



    vector <MMFFciParameter*> ciParameters;
    vector <MMFFaromaticParameter*> aromaticParameters;
    vector <MMFFatypeParameter*> atypeParameters;
    vector <MMFFHParameter*> HParameters;
    vector <MMFFatomParameter*> atomParameters;
    vector <MMFFabParameter*> abParameters;
    vector <MMFFbsParameter*> bsParameters;
    vector <MMFFsbParameter*> sbParameters;
    vector <MMFFsbParameter*> sb2Parameters;
    vector <MMFFtoParameter*> toParameters;
    vector <MMFFvwParameter*> vwParameters;
    vector <MMFFopParameter*> opParameters;

    vector <MMFFbsInteraction*> bsInteractions;
    vector <MMFFabInteraction*> abInteractions;
    vector <MMFFtoInteraction*> toInteractions;
    vector <MMFFvwInteraction*> vwInteractions;
    vector <MMFFopInteraction*> opInteractions;
    vector <MMFFelInteraction*> elInteractions;
    vector <MMFFsbInteraction*> sbInteractions;

    vector <MMFFvwInteraction*> vwNBInteractions;
    vector <MMFFelInteraction*> elNBInteractions;

    int load_parameters ();
    int load_ci_parameters ();
    int load_H_parameters ();
    int load_atype_parameters ();
    int load_atom_parameters ();
    int load_ab_parameters ();
    int load_bs_parameters ();
    int load_sb_parameters ();
    int load_sb2_parameters ();
    int load_to_parameters ();
    int load_vw_parameters ();
    int load_op_parameters ();
    int load_aromatic_parameters ();

    void empiric_to_parameters (Atom *at1, Atom *at2, Atom *at3, Atom *at4, float& v1, float &v2, float &v3);

    int get_bond_type (int I, int J);


    double get_charge_increment (int I, int J);

    double get_vdw_e (Atom *at1, Atom *at2);
    double get_vdw_r0 (Atom *at1, Atom *at2); 
    double get_bs_kb (int I, int J);
    double get_bs_r0 (int I, int J);
    double get_ab_ka (int I, int J, int K);
    double get_ab_theta0 (int I, int J, int K);
    double get_charge_rho (int I);
    double get_charge_theta (int I);
    int get_pt_row (int an);

    string get_aromatic_string (Atom *at); 
    int get_atype_line (int I); 
    int get_ab_pl (int type, int I, int J, int K);    
    int get_bs_pl (int type, int I, int J);
    int get_sb_pl (int I, int J, int K);
    int get_sb2_pl (int RI, int RJ, int RK);
    int get_to_pl (int type, int I, int J, int K, int L);
    int get_vw_pl (int I);
    int get_op_pl (int I, int J, int K, int L);
    bool get_linear (int J);
	bool get_sbmb (int J);
	bool get_arom (int J);
	
	int getBondType (Atom *a, Atom *b);
	int getAngleType (Atom *a, Atom *b, Atom *c);



};















#endif
