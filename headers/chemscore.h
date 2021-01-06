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


#define ZN_CS_DONOR 0
#define ZN_CS_ACCEPTOR 1
#define ZN_CS_BOTH 2
#define ZN_CS_POLAR 3
#define ZN_CS_NONPOLAR 4
#define ZN_CS_METAL 5

#include "FF.h"
#include "obabel_includes.h"
#include "datagrid.h"

float f (float x, float x1, float x2);

class ChemscoreHBInteraction : public ForceFieldInteraction
{
    public:
    Atom * root;
    float value ();
		inline bool isElectrostatic () {return true;};
		inline bool isHbond () {return true;};	
};


class ChemscoreLiInteraction : public ForceFieldInteraction
{   
    public:
    bool metal;
    float value ();
};


class ChemscoreClInteraction : public ForceFieldInteraction
{
    public:
    int type;
    float value ();
};






class Chemscore : public ForceField {

public:

    Chemscore ();
    void clear_nonbonded_interactions ();
    void load_mol (ZNMolecule *mol);
	void load_grids (vect cent, double rad);
 //   void load_internal_interactions ();
	
	void load_nonbonded_interactions_for_atom (Atom *a, queue <ForceFieldInteraction *> *queue = 0) {

		ZNMolecule *mol = (ZNMolecule *) a -> GetParent ();
		objectList<Atom*>* nbAtoms = far_grid->getNeighborObjects(get_coordinates (&*a));
		if (nbAtoms) {
			vector <Atom *> neighbours = nbAtoms->objects;  
		//	cerr <<"load nonbonded inter"<<neighbours.size ()<<endl;
			for (unsigned int j=0; j<neighbours.size (); j++) {
				Atom * root = NULL;
				bool hb = false;
				bool li = false;
				bool me = false;
				int ltype = getChemscoretype (&*a);
				int ptype = getChemscoretype (neighbours[j]);
				if (ltype == ZN_CS_DONOR) {
					if (ptype == ZN_CS_ACCEPTOR) {
						hb = true;
						OBBondIterator i;
						root = a -> BeginNbrAtom (i);
					}
				}
				else if (ltype == ZN_CS_ACCEPTOR) {
					if (ptype == ZN_CS_DONOR) {
						hb = true;
						OBBondIterator i;
						root = neighbours[j] -> BeginNbrAtom (i);
					}
					else if (ptype == ZN_CS_METAL) {
						me = true;
					}
				}
				else if (ltype == ZN_CS_NONPOLAR) {
					
					if (ptype == ZN_CS_NONPOLAR) {
						li = true;
					}
				}

				if (hb) {
					ChemscoreHBInteraction *hbint = new ChemscoreHBInteraction;
					hbint->at1 = &*a;
					hbint->at2=neighbours[j];
					hbint->root = root;
					if (queue) queue -> push (hbint);
					else HBInteractions.push_back (hbint);  
				}    
				else if (li || me) {
					ChemscoreLiInteraction *liint = new ChemscoreLiInteraction;
					liint->at1 = &*a;
					liint->at2 = neighbours[j];
					liint->metal = me;
					if (queue) queue -> push (liint);
					else LiInteractions.push_back (liint);  
				}
				//clash
	//			if (neighbours[j] -> GetAtomicNum ()==1 || a -> GetAtomicNum () ==1) continue;
				ChemscoreClInteraction *clint = new ChemscoreClInteraction;
				clint->at1 = &*a;
				clint->at2=neighbours[j];
				if (ptype == ZN_CS_METAL && ltype == ZN_CS_ACCEPTOR) clint->type = 1;
				else if ((clint->at1->GetAtomicNum () ==7 || clint->at1->GetAtomicNum () ==8) && (clint->at2->GetAtomicNum () ==7 || clint->at2->GetAtomicNum () ==8)) {
					if (mol -> bonded_to (clint->at1,1, 1) || mol -> bonded_to (clint->at2,1,1)) clint->type = 0;
				}
				else clint->type = 2;
				if (queue) queue -> push (clint);
				else ClInteractions.push_back (clint);  
				
				
				
			}        
		}
		
		
	}
	
    void load_nonbonded_interactions () {
		FOR_ATOMS_OF_MOL (a, target_mol) {
		load_nonbonded_interactions_for_atom (&*a);
		}
	};
    void update ();

	
    double compute_total_energy () {
		double E = 0., Ehb = 0., Ecl = 0., Eli = 0., Ehb2 = 0.,  Ecl2 = 0., Eli2 = 0.;
		
		for (unsigned int i = 0; i < HBInteractions.size (); i++) {
			E += HBInteractions [i] ->value ();
	//		Ehb += HBInteractions [i] ->value ();
		}
		for (unsigned int i = 0; i < ClInteractions.size (); i++) {
			E += ClInteractions [i] ->value ();
	//		Ecl += ClInteractions [i] ->value ();
		}
		for (unsigned int i = 0; i < LiInteractions.size (); i++) {
			E += LiInteractions [i] ->value ();
	//		Eli += LiInteractions [i] ->value ();
		}
		return E;
		
	/*	
		FOR_ATOMS_OF_MOL (a, target_mol) {
			vect v = get_coordinates (&*a);
			int atype = getChemscoretype (&*a);
			if (atype == ZN_CS_NONPOLAR) { 
				
				Eli2 += LiGrid -> get_value (v.x(), v.y(), v.z());
			}
			else if (atype == ZN_CS_ACCEPTOR) {
				objectList<Atom*>* nbAtoms = far_grid ->getNeighborObjects(v);
				if (nbAtoms) {
					ChemscoreHBInteraction *hbint = new ChemscoreHBInteraction;
					hbint->at1 = &*a;
					vector <Atom *> neighbours = nbAtoms->objects;  
					
					for (unsigned int j=0; j<neighbours.size (); j++) {
						if (getChemscoretype (neighbours[j]) == ZN_CS_DONOR) {
							OBBondIterator i;
							hbint ->root = neighbours[j] -> BeginNbrAtom (i);
							hbint ->at2 = neighbours[j];
							Ehb2 += hbint -> value ();
						}
					}
					delete hbint;
				}
			}
			else if (atype == ZN_CS_DONOR) {
				objectList<Atom*>* nbAtoms = far_grid ->getNeighborObjects(v);
				if (nbAtoms) {
					ChemscoreHBInteraction *hbint = new ChemscoreHBInteraction;
					hbint->at1 = &*a;
					vector <Atom *> neighbours = nbAtoms->objects;  
					OBBondIterator i;
					hbint ->root = a -> BeginNbrAtom (i);
					for (unsigned int j=0; j<neighbours.size (); j++) {
						if (getChemscoretype (neighbours[j]) == ZN_CS_ACCEPTOR) {
							hbint ->at2 = neighbours[j];
							Ehb2 += hbint -> value ();
						}
					}
					delete hbint;
				}
			}
			
			else {}
			if (a ->GetAtomicNum () !=1) {
				Ecl2 += ClGrid ->get_value (v.x(), v.y(), v.z());
			}
		}	
	//	cerr <<"HB = "<<Ehb<<"             Cl = "<<Ecl<<"               Li = "<<Eli<<endl;
	//	cerr <<"HB = "<<Ehb2<<"             Cl = "<<Ecl2<<"               Li = "<<Eli2<<endl;
		
		return Ehb2* + Ecl2 + Eli2;	
		*/
	};
	void initialize_mol (ZNMolecule *mol);
	
			
    void compute_forces ();
	
	

    vector<ChemscoreClInteraction *> ClInteractions;
    vector<ChemscoreHBInteraction *> HBInteractions;
    vector<ChemscoreLiInteraction *> LiInteractions;
	
	float Electrostatic_potential (vect v) {
			float E = 0.; 
		objectList<Atom*>* nbAtoms = far_grid ->getNeighborObjects(v);
		if (nbAtoms) {
					vector <Atom *> neighbours = nbAtoms->objects;
			for (unsigned int j=0; j<neighbours.size (); j++) {
				float q = neighbours[j] -> GetPartialCharge ();
				
				if (q) {

					float r = dist (v, get_coordinates (neighbours[j]));
					float e = 332.0716;
					E += q*e/(r);
				}
			}
		}	
		return E;		 
	};
	
	
		float Donorvalue (vect v) {
			float E = 0.; 
		objectList<Atom*>* nbAtoms = far_grid ->getNeighborObjects(v);
		if (nbAtoms) {
					vector <Atom *> neighbours = nbAtoms->objects;
			for (unsigned int j=0; j<neighbours.size (); j++) {
				int ptype = getChemscoretype (neighbours[j]);

				if (ptype == ZN_CS_ACCEPTOR) {

					float r = dist (v, get_coordinates (neighbours[j]));
				    float dghbond = -3.34f;
					float r0 = 1.85f;
					float drab = r-r0;
					if (drab < 0) drab = -drab;
					float dr1 = 0.25f;
					float dr2 = 0.65f;

					E += dghbond * f(drab, dr1, dr2);
				}
			}
		}	
		return E;		 
	};

	
	
	
	
		float Acceptorvalue (vect v) {
			float E = 0.; 
		objectList<Atom*>* nbAtoms = far_grid ->getNeighborObjects(v);
		if (nbAtoms) {
					vector <Atom *> neighbours = nbAtoms->objects;
			for (unsigned int j=0; j<neighbours.size (); j++) {
				int ptype = getChemscoretype (neighbours[j]);

				if (ptype == ZN_CS_DONOR) {

					float r = dist (v, get_coordinates (neighbours[j]));
				    float dghbond = -3.34f;
					float r0 = 1.85f;
					float drab = r-r0;
					if (drab < 0) drab = -drab;
					float dr1 = 0.25f;
					float dr2 = 0.65f;

					E += dghbond * f(drab, dr1, dr2);
				}
			}
		}	
		return E;		 
	};
	
	
	
	
	float Livalue (vect v) {
		float E = 0.;   
		
		Atom *dummy = new Atom;
		dummy -> SetVector (v); 
		objectList<Atom*>* nbAtoms = far_grid ->getNeighborObjects(v);
		if (nbAtoms) {
			ChemscoreLiInteraction *liint = new ChemscoreLiInteraction;
			liint ->metal = false;
			liint->at1 = dummy;
			vector <Atom *> neighbours = nbAtoms->objects;  
			
			for (unsigned int j=0; j<neighbours.size (); j++) {
				int ptype = getChemscoretype (neighbours[j]);
				if (ptype == ZN_CS_NONPOLAR) {
					
					
					liint->at2 = neighbours[j];
					E += liint -> value ();
					
				}
			}
			delete liint;
			
		} 
		delete dummy;
	//	cerr << E << "energy" << endl;
		return E;
		
	};
	

		float Clvalue (vect v) {
		float E = 0.;   
		
		Atom *dummy = new Atom;
		dummy -> SetVector (v); 
		objectList<Atom*>* nbAtoms = far_grid ->getNeighborObjects(v);
		if (nbAtoms) {
			ChemscoreClInteraction *clint = new ChemscoreClInteraction;
			clint ->type = 2;
			clint->at1 = dummy;
			vector <Atom *> neighbours = nbAtoms->objects;  
			
			for (unsigned int j=0; j<neighbours.size (); j++) {
						
					if (neighbours[j] ->GetAtomicNum () !=1) {
						clint->at2 = neighbours[j];
						E += clint -> value ();
					}
					
			}
			delete clint;
			
		} 
		delete dummy;
	//	cerr << E << "energy" << endl;
		return E;
		
	};

	
		private:
    int getChemscoretype (Atom *at);
	
	DataGrid *LiGrid, *ClGrid;
	
	
};

#endif
