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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <queue>
#include "molecule.h"
#include "constants.h"
#include <math.h>



vect get_coordinates (Atom *at) {
	lock_force_mutex(at);
	vect v = (vect &) at -> GetVector ();
	unlock_force_mutex(at);
	return v;
}

vect get_coordinates (SurfVertex *ver) {
	vect v = (vect &) ver -> GetVector ();
	return v;
}



void set_coordinates (Atom *at, vect v) {
	lock_force_mutex(at);
	at -> SetVector (v);
	unlock_force_mutex(at);
}

void sum_to_coordinates (Atom *at, vect v) {
	lock_force_mutex(at);
	vect v1 = (vect &) at -> GetVector ();
	vect v2 = sum (v, v1);
	at -> SetVector (v2);
	unlock_force_mutex(at);	
}



void mend_coordinates (Molecule *mol) {
	double x = 0.;
	int x1 = 1;
	double x2 = 0.;
	FOR_ATOMS_OF_MOL (a, mol) {
		x += 0.1;
		x2 += 0.2;
		x1 *= -1;
		vect v (x+(x2*x1), x*x1, x);
		a -> SetVector (v);
	}
}



vect find_mass_center (vector<Atom*>& invec){

    unsigned int i;
    vect midCoo;

	
	unsigned int numMid = 0;
	for (i = 0; i < invec.size(); i++) {
	    vect a = sum (midCoo, get_coordinates (invec[i]));
        midCoo = a;
		numMid++;
		
	}
	
    midCoo.multiply (1.0f/(float) numMid);
	
/*
	float radius = 0.0;
	for (i = 0; i < invec.size(); i++) {

			float r2 = (invec[i]-> GetVector ()[0] - midCoo[0]) *(invec[i]-> GetVector ()[0] - midCoo[0]) +
						(invec[i]-> GetVector ()[1] - midCoo[1]) *(invec[i]-> GetVector ()[1] - midCoo[1]) +
						(invec[i]-> GetVector ()[2] - midCoo[2]) *(invec[i]-> GetVector ()[2] - midCoo[2]);
			if (r2 > radius) {
				radius = r2;
			}
	}
	
*/

    return midCoo;

}

double get_vdw (Atom *at) {
    int atn = at -> GetAtomicNum ();
    double rad = etab.GetVdwRad (atn);
    return rad;
}

bool get_visible (Atom *at) {
    BoolData *d = (BoolData *) at -> GetData ("visible");
    assert (d);
    return d -> GetGenericValue ();
}


bool get_visible (Bond *bo) {
    return (get_visible (bo -> GetBeginAtom ()) && (get_visible (bo -> GetEndAtom ())));
}

void set_visible (Atom *at, bool vis) {
    BoolData *d = (BoolData *) at -> GetData ("visible");
    assert (d);
    d -> SetValue (vis);
}


color get_color (Atom *at) {
    ColorData *d = (ColorData *) at -> GetData ("color");
    assert (d);
    return d -> GetGenericValue ();
}


void set_color (Atom *at, color col) {
    ColorData *d = (ColorData *) at -> GetData ("color");
    assert (d);
    d -> SetValue (col);
}

vect get_force (Atom *at) {
	VectorData *d = (VectorData *) at -> GetData ("f");
	assert (d);
	return d -> GetGenericValue ();
}

void set_force (Atom *at, vect v) {
	VectorData *d = (VectorData *) at -> GetData ("f");
	assert (d);
	d -> SetValue (v);
}

vect get_back_force (Atom *at) {
	VectorData *d = (VectorData *) at -> GetData ("back_force");
	assert (d);
	return d -> GetGenericValue ();
}

void set_back_force (Atom *at, vect v) {
	VectorData *d = (VectorData *) at -> GetData ("back_force");
	assert (d);
	d -> SetValue (v);
}

void flush_forces (Atom *at) {
	vect f = get_back_force (at);
	set_force (at, f);
	vect n (0., 0., 0.);
	set_back_force (at, n);
}

void flush_scores (Atom *at) {
	double f = get_back_score (at);
//	cerr << f << endl;
	set_score (at, f);
	set_back_score (at, 0.);
}


void lock_force_mutex (Atom *at) {
	MutexData *d = (MutexData *) at -> GetData ("force_mutex");
	assert (d);
	d -> GetGenericValue () ->lock ();
}

void unlock_force_mutex (Atom *at) {
	MutexData *d = (MutexData *) at -> GetData ("force_mutex");
	assert (d);
	d -> GetGenericValue () ->unlock ();
}



int get_MMFFtype (Atom *atom) {
	Molecule *mol = (Molecule *) atom -> GetParent ();
    OBBond *bond;
    int oxygenCount, nitrogenCount, sulphurCount, doubleBondTo;
    ////////////////////////////////
    // Aromatic Atoms
    ////////////////////////////////
    if (atom->IsAromatic()) {
      if (atom->IsInRingSize(5)) {
        bool IsAromatic = false;
        vector<OBAtom*> alphaPos, betaPos;
        vector<OBAtom*> alphaAtoms, betaAtoms;

	if (atom->IsSulfur()) {
	  return 44; // Aromatic 5-ring sulfur with pi lone pair (STHI)
	}
	if (atom->IsOxygen()) {
	  return 59; // Aromatic 5-ring oxygen with pi lone pair (OFUR)
	}
	if (atom->IsNitrogen()) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
	    if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
	      return 82; // N-oxide nitrogen in 5-ring alpha position, 
	                 // N-oxide nitrogen in 5-ring beta position, 
			 // N-oxide nitrogen in other 5-ring  position, 
			 // (N5AX, N5BX, N5OX) 
	    }
	  }
	}
        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (!((mol -> GetBond(atom, &*nbr))->IsAromatic()) || !nbr->IsInRingSize(5))
            continue;
   
	  if (IsInSameRing(atom, &*nbr)) {
	    alphaPos.push_back(&*nbr);
	  }
          
	  FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
            if (nbrNbr->GetIdx() == atom->GetIdx())
              continue;
            if (!((mol -> GetBond(&*nbr, &*nbrNbr))->IsAromatic()) || !nbrNbr->IsInRingSize(5))
              continue;
             
            IsAromatic = true;
	    
	    if (IsInSameRing(atom, &*nbrNbr)) {
              betaPos.push_back(&*nbrNbr);
	    }
          }
        }
	
	if (IsAromatic) {
          
	  
	  for (unsigned int i = 0; i < alphaPos.size(); i++) {
  	    if (alphaPos[i]->IsSulfur()) {
              alphaAtoms.push_back(alphaPos[i]);
	    } else if (alphaPos[i]->IsOxygen()) {
              alphaAtoms.push_back(alphaPos[i]);
	    } else if (alphaPos[i]->IsNitrogen() && (alphaPos[i]->GetValence() == 3)) {
	      bool IsNOxide = false;
	      FOR_NBORS_OF_ATOM (nbr, alphaPos[i]) {
	        if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
		  IsNOxide = true;
		}
	      }

	      if (!IsNOxide) {
                alphaAtoms.push_back(alphaPos[i]);
	      }
	    }
          }
          for (unsigned int i = 0; i < betaPos.size(); i++) {
            if (betaPos[i]->IsSulfur()) {
              betaAtoms.push_back(betaPos[i]);
            } else if (betaPos[i]->IsOxygen()) {
              betaAtoms.push_back(betaPos[i]);
            } else if (betaPos[i]->IsNitrogen() && (betaPos[i]->GetValence() == 3)) {
	      bool IsNOxide = false;
              FOR_NBORS_OF_ATOM (nbr, betaPos[i]) {
	        if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
		  IsNOxide = true;
		}
	      }

	      if (!IsNOxide) {
                betaAtoms.push_back(betaPos[i]);
	      }
            }
          }
	  if (!betaAtoms.size()) {
	    nitrogenCount = 0;
	    FOR_NBORS_OF_ATOM (nbr, atom) {
	 //     cout << "BOSum=" << nbr->BOSum() << endl;
	      if (nbr->IsNitrogen() && (nbr->GetValence() == 3)) {
	        if ((nbr->BOSum() == 4) && nbr->IsAromatic()) {
	          nitrogenCount++;
	        } else if ((nbr->BOSum() == 3) && !nbr->IsAromatic()) {
	          nitrogenCount++;
	        }
	      }
	    }
	    if (nitrogenCount >= 2) {
	      return 80; // Aromatic carbon between N's in imidazolium (CIM+)
	    }
          }
	  if (!alphaAtoms.size() && !betaAtoms.size()) {
            if (atom->IsCarbon()) {
	      // there is no S:, O:, or N:
              // this is the case for anions with only carbon and nitrogen in the ring
              return 78; // General carbon in 5-membered aromatic ring (C5)
	    } else if (atom->IsNitrogen()) {
	      if (atom->GetValence() == 3) {
	        // this is the N: atom
	        return 39; // Aromatic 5 ring nitrogen with pi lone pair (NPYL)
	      } else {
	        // again, no S:, O:, or N:
	        return 76; // Nitrogen in 5-ring aromatic anion (N5M)
	      }
	    }
          }
	  if (alphaAtoms.size() == 2) {
	    if (atom->IsCarbon() && IsInSameRing(alphaAtoms[0], alphaAtoms[1])) {
	      if (alphaAtoms[0]->IsNitrogen() && alphaAtoms[1]->IsNitrogen()) {
	        if ((alphaAtoms[0]->GetValence() == 3) && (alphaAtoms[1]->GetValence() == 3)) {
	          return 80; // Aromatic carbon between N's in imidazolium (CIM+)
		}
	      }
	    }
	  }
	  if (alphaAtoms.size() && !betaAtoms.size()) {
	    if (atom->IsCarbon()) {
	      return 63; // Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
	    } else if (atom->IsNitrogen()) {
	      if (atom->GetValence() == 3) {
	        return 81; // Posivite nitrogen in 5-ring alpha position (N5A+)
	      } else {
	        return 65; // Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
	      }
	    }
	  }
	  if (!alphaAtoms.size() && betaAtoms.size()) {
	    if (atom->IsCarbon()) {
	      return 64; // Aromatic 5-ring C, beta to N:, O:, or S: (C5B)
	    } else if (atom->IsNitrogen()) {
	      if (atom->GetValence() == 3) {
	        return 81; // Posivite nitrogen in 5-ring beta position (N5B+)
	      } else {
	        return 66; // Aromatic 5-ring N, beta to N:, O:, or S: (N5B)
	      }
	    }
	  }
	  if (alphaAtoms.size() && betaAtoms.size()) {
	    for (unsigned int i = 0; i < alphaAtoms.size(); i++) {
	      for (unsigned int j = 0; j < betaAtoms.size(); j++) {
	        if (!IsInSameRing(alphaAtoms[i], betaAtoms[j])) {
		  if (atom->IsCarbon()) {
                    return 78; // General carbon in 5-membered aromatic ring (C5)
	          } else if (atom->IsNitrogen()) {
	            return 79; // General nitrogen in 5-membered aromatic ring (N5)
		  }
		}
	      }
	    }
	    for (unsigned int i = 0; i < alphaAtoms.size(); i++) {
              if (alphaAtoms[i]->IsSulfur() || alphaAtoms[i]->IsOxygen()) {
	        if (atom->IsCarbon()) {
	          return 63; // Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
	        } else if (atom->IsNitrogen()) {
	          return 65; // Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
		}
              }
	    }
	    for (unsigned int i = 0; i < betaAtoms.size(); i++) {
              if (betaAtoms[i]->IsSulfur() || betaAtoms[i]->IsOxygen()) {
	        if (atom->IsCarbon()) {
	          return 64; // Aromatic 5-ring C, beta to N:, O:, or S: (C5B)
	        } else if (atom->IsNitrogen()) {
	          return 66; // Aromatic 5-ring N, beta to N:, O:, or S: (N5B)
		}
              }
	    }
	    
	    if (atom->IsCarbon()) {
              return 78; // General carbon in 5-membered aromatic ring (C5)
	    } else if (atom->IsNitrogen()) {
	      return 79; // General nitrogen in 5-membered aromatic ring (N5)
	    }
	  }
        }
      }
    
      if (atom->IsInRingSize(6)) {
	
	if (atom->IsCarbon()) {
          return 37; // Aromatic carbon, e.g., in benzene (CB)
	} else if (atom->IsNitrogen()) {
	  FOR_NBORS_OF_ATOM (nbr, atom) {
	    if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
	      return 69; // Pyridinium N-oxide nitrogen (NPOX)
	    }
	  }
	  
	  if (atom->GetValence() == 3) {
	    return 58; // Aromatic nitrogen in pyridinium (NPD+)
	  } else {
	    return 38; // Aromatic nitrogen with sigma lone pair (NPYD)
	  }
	}
      }
    }
    
    ////////////////////////////////
    // Hydrogen
    ////////////////////////////////
    if (atom->GetAtomicNum() == 1) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (nbr->IsCarbon()) {
	  return 5; // Hydrogen attatched to carbon (HC)
	}
        if (nbr->GetAtomicNum() == 14) {
	  return 5; // Hydrogen attatched to silicon (HSI)
	}
        if (nbr->IsOxygen()) {
	  if (nbr->BOSum() == 3) {
	    if (nbr->GetValence() == 3) {
	      return 50; // Hydrogen on oxonium oxygen (HO+)
	    } else {
	      return 52; // Hydrogen on oxenium oxygen (HO=+)
	    }
	  }
	  
	  int hydrogenCount = 0;
          FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	    if (nbrNbr->IsHydrogen()) {
	      hydrogenCount++;
	      continue;
	    }
	    if (nbrNbr->IsCarbon()) {
              if (nbrNbr->IsAromatic()) {
		return 29; // phenol
              }
                     
	      FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
	        if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
		  continue;

		bond = mol -> GetBond(&*nbrNbr, &*nbrNbrNbr);
	        if (bond->IsDouble()) {
		  if (nbrNbrNbr->IsOxygen()) {
		    return 24; // Hydroxyl hydrogen in carboxylic acids (HOCO)
		  }
		  if (nbrNbrNbr->IsCarbon() || nbrNbrNbr->IsNitrogen()) {
		    return 29; // Enolic or phenolic hydroxyl hydrogen,
		               // Hydroxyl hydrogen in HO-C=N moiety (HOCC, HOCN)
		  }
		}
	      }
	    }
	    if (nbrNbr->IsPhosphorus()) {
	      return 24; // Hydroxyl hydrogen in H-O-P moiety (HOP)
	    }
	    if (nbrNbr->IsSulfur()) {
	      return 33; // Hydrogen on oxygen attached to sulfur (HOS)
	    }
	  
	  }
	  if (hydrogenCount == 2) {
	    return 31; // Hydroxyl hydrogen in water (HOH)
	  }

	  return 21; // Hydroxyl hydrogen in alcohols, Generic hydroxyl hydrogen (HOR, HO)
	}
        if (nbr->IsNitrogen()) {
	  switch (get_MMFFtype(&*nbr)) {
	    case 81:
	      return 36; // Hydrogen on imidazolium nitrogen (HIM+)
	    case 68:
	      return 23; // Hydrogen on N in N-oxide (HNOX)
	    case 67:
	      return 23; // Hydrogen on N in N-oxide (HNOX)
	    case 62:
	      return 23; // Generic hydrogen on sp3 nitrogen, e.g., in amines (HNR)
	    case 56:
	      return 36; // Hydrogen on guanimdinium nitrogen (HGD+)
	    case 55:
	      return 36; // Hydrogen on amidinium nitrogen (HNN+)
	    case 43:
	      return 28; // Hydrogen on NSO, NSO2, or NSO3 nitrogen, Hydrogen on N triply bonded to C (HNSO, HNC%)
	    case 39:
	      return 23; // Hydrogen on nitrogen in pyrrole (HPYL)
	    case 8:
	      return 23; // Generic hydrogen on sp3 nitrogen, e.g., in amines, Hydrogen on nitrogen in ammonia (HNR, H3N)
	  }

	  if (nbr->BOSum() == 4) {
	    if (nbr->GetValence() == 2) {
	      return 28; // Hydrogen on N triply bonded to C (HNC%)
	    } else {
	      return 36; // Hydrogen on pyridinium nitrogen, Hydrogen on protonated imine nitrogen (HPD+, HNC+)
	    }
	  }

	  if (nbr->GetValence() == 2) {
	    FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	      if (nbrNbr->IsHydrogen())
	        continue;
              
	      bond = mol -> GetBond(&*nbr, &*nbrNbr);
	      if (bond->IsDouble()) {
	        if (nbrNbr->IsCarbon() || nbrNbr->IsNitrogen()) {
	          return 27; // Hydrogen on imine nitrogen, Hydrogen on azo nitrogen (HN=C, HN=N) 
	        }

		return 28; // Generic hydrogen on sp2 nitrogen (HSP2)
	      }
	    }
	  }
	  
	  FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	    if (nbrNbr->IsHydrogen())
	      continue;
	    
	    if (nbrNbr->IsCarbon()) {
	      if (nbrNbr->IsAromatic()) {
	        return 28; // deloc. lp pair
	      }
	      
	      FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
	        if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
	          continue;
              
	        bond = mol -> GetBond(&*nbrNbr, &*nbrNbrNbr);
	        if (bond->IsDouble()) {
	          if (nbrNbrNbr->IsCarbon() || nbrNbrNbr->IsNitrogen() || nbrNbrNbr->IsOxygen() || nbrNbrNbr->IsSulfur()) {
	            return 28; // Hydrogen on amide nitrogen, Hydrogen on thioamide nitrogen,
		               // Hydrogen on enamine nitrogen, Hydrogen in H-N-C=N moiety (HNCO, HNCS, HNCC, HNCN)
	          }
	        }
	      }
	    }
	    if (nbrNbr->IsNitrogen()) {
              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
	        if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
	          continue;
              
	        bond = mol -> GetBond(&*nbrNbr, &*nbrNbrNbr);
	        if (bond->IsDouble()) {
	          if (nbrNbrNbr->IsCarbon() || nbrNbrNbr->IsNitrogen()) {
	            return 28; // Hydrogen in H-N-N=C moiety, Hydrogen in H-N-N=N moiety (HNNC, HNNN)
	          }
	        }
	      }
	    }
	    if (nbrNbr->IsSulfur()) {
              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
	        if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
	          continue;
              
	        if (nbrNbrNbr->IsOxygen() || (nbrNbrNbr->GetValence() == 1)) {
	          return 28; // Hydrogen on NSO, NSO2 or NSO3 nitrogen (HNSO)
	        }
	      }
	    }
	  }
              
	  return 23; // Generic hydrogen on sp3 nitrogen e.g., in amines,
	             // Hydrogen on nitrogen in pyrrole, Hydrogen in ammonia,
		     // Hydrogen on N in N-oxide (HNR, HPYL, H3N, HNOX)
	}
        if (nbr->IsSulfur() || nbr->IsPhosphorus()) {
	  return 71; // Hydrogen attached to sulfur, Hydrogen attached to >S= sulfur doubly bonded to N,
	             // Hydrogen attached to phosphorus (HS, HS=N, HP)
	}
      }
    }

    ////////////////////////////////
    // Lithium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 3) {
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 92; // Lithium cation (LI+)
      }
    }
 
    ////////////////////////////////
    // Carbon
    ////////////////////////////////
    if (atom->GetAtomicNum() == 6) {
      // 4 neighbours
      if (atom->GetValence() == 4) {
        if (atom->IsInRingSize(3)) {
          return 22; // Aliphatic carbon in 3-membered ring (CR3R)
	} 
	
	if (atom->IsInRingSize(4)) {
          return 20; // Aliphatic carbon in 4-membered ring (CR4R)
	}
        
	return 1; // Alkyl carbon (CR)
      }
      // 3 neighbours
      if (atom->GetValence() == 3) {
        int N2count = 0;
	int N3count = 0;
        oxygenCount = sulphurCount = doubleBondTo = 0;

	FOR_NBORS_OF_ATOM (nbr, atom) {
	  bond = mol -> GetBond(&*nbr, atom);
	  if (bond->IsDouble()) {
	    doubleBondTo = nbr->GetAtomicNum();
          }

	  if (nbr->GetValence() == 1) {
	    if (nbr->IsOxygen()) {
	      oxygenCount++;
	    } else if (nbr->IsSulfur()) {
	      sulphurCount++;
	    }
	  } else if (nbr->GetValence() == 3) {
	    if (nbr->IsNitrogen()) {
	      N3count++;
	    }
	  } else if ((nbr->GetValence() == 2) && bond->IsDouble()) {
	    if (nbr->IsNitrogen()) {
	      N2count++;
	    }
	  }
	}
	if ((N3count >= 2) && (doubleBondTo == 7) && !N2count) {
	  // N3==C--N3
	  return 57; // Guanidinium carbon, Carbon in +N=C-N: resonance structures (CGD+, CNN+)
	}
	if ((oxygenCount == 2) || (sulphurCount == 2)) {
	  // O1-?-C-?-O1 or S1-?-C-?-S1
	  return 41; // Carbon in carboxylate anion, Carbon in thiocarboxylate anion (CO2M, CS2M)
        }
	if (atom->IsInRingSize(4) && (doubleBondTo == 6)) {
	        return 30; // Olefinic carbon in 4-membered ring (CR4E)
        }
        if ((doubleBondTo ==  7) || (doubleBondTo ==  8) || 
	    (doubleBondTo == 15) || (doubleBondTo == 16)) {
	  // C==N, C==O, C==P, C==S
	  return 3; // Generic carbonyl carbon, Imine-type carbon, Guanidine carbon,
	            // Ketone or aldehyde carbonyl carbon, Amide carbonyl carbon,
		    // Carboxylic acid or ester carbonyl carbon, Carbamate carbonyl carbon,
		    // Carbonic acid or ester carbonyl carbon, Thioester carbonyl (double
		    // bonded to O or S), Thioamide carbon (double bonded to S), Carbon
		    // in >C=SO2, Sulfinyl carbon in >C=S=O, Thiocarboxylic acid or ester 
		    // carbon, Carbon doubly bonded to P (C=O, C=N, CGD, C=OR, C=ON, COO,
		    // COON, COOO, C=OS, C=S, C=SN, CSO2, CS=O, CSS, C=P)
        }
	
	return 2; // Vinylic Carbon, Generic sp2 carbon (C=C, CSP2)
	
      }
      // 2 neighbours
      if (atom->GetValence() == 2) {
        return 4; // Acetylenic carbon, Allenic caron (CSP, =C=)
      }
      // 1 neighbours
      if (atom->GetValence() == 1) {
        return 60; // Isonitrile carbon (C%-)
      }
    }

    ////////////////////////////////
    // Nitrogen
    ////////////////////////////////
    if (atom->GetAtomicNum() == 7) {
      // 4 neighbours
      if (atom->GetValence() == 4) {
	FOR_NBORS_OF_ATOM (nbr, atom) {
	  if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
	    return 68; // sp3-hybridized N-oxide nitrogen (N3OX)
          }
	}

	return 34; // Quaternary nitrogen (NR+)
      }
      // 3 neighbours
      if (atom->GetValence() == 3) {
        if (atom->BOSum() == 4) {
	  oxygenCount = nitrogenCount = doubleBondTo = 0;
	  
	  FOR_NBORS_OF_ATOM (nbr, atom) {
	    if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
	      oxygenCount++;
	    }
	    if (nbr->IsNitrogen()) {
	      bond = mol -> GetBond(&*nbr, atom);
	      if (bond->IsDouble()) {
	        doubleBondTo = 7;
              }
	    }
	    if (nbr->IsCarbon()) {
	      bond = mol -> GetBond(&*nbr, atom);
	      if (bond->IsDouble()) {
	        FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	          if (nbrNbr->IsNitrogen() && (nbrNbr->GetValence() == 3)) {
		    nitrogenCount++;
		  }
	        }
              }
	    }
	  }

	  if (oxygenCount == 1) {
	    return 67; // sp2-hybridized N-oxide nitrogen (N2OX)
	  }
	  if (oxygenCount >= 2) {
	    return 45; // Nitrogen in nitro group, Nitrogen in nitrate group (NO2, NO3)
	  }

          if (nitrogenCount == 1) {
	    return 54; // Iminium nitrogen (N+=C)
	  }
          if (nitrogenCount == 2) {
	    return 55; // Either nitrogen in N+=C-N: (NCN+)
	  }
          if (nitrogenCount == 3) {
	    return 56; // Guanidinium nitrogen (NGD+)
	  }
	  
	  if (doubleBondTo == 7) {
	    return 54; // Positivly charged nitrogen doubly bonded to nitrogen (N+=N)
	  }
	}
	
	if (atom->BOSum() == 3) {
	  bool IsAmide = false;
	  bool IsSulfonAmide = false;
	  bool IsNNNorNNC = false;
	  int tripleBondTo = 0;
	  doubleBondTo = 0;
	  
	  FOR_NBORS_OF_ATOM (nbr, atom) {
	    if (nbr->IsSulfur() || nbr->IsPhosphorus()) {
	      oxygenCount = 0;
	      
	      FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
		if (nbrNbr->IsOxygen() && (nbrNbr->GetValence() == 1)) {
		  oxygenCount++;
		}
	      }
              if (oxygenCount >= 2) {
	        IsSulfonAmide = true;
		//return 43; // Sulfonamide nitrogen (NSO2, NSO3)
	      }
	    }
	  }
	
	  FOR_NBORS_OF_ATOM (nbr, atom) {
	    if (nbr->IsCarbon()) {
	      FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	        bond = mol -> GetBond(&*nbr, &*nbrNbr);
		if (bond->IsDouble() && (nbrNbr->IsOxygen() || nbrNbr->IsSulfur())) {
	          IsAmide = true;
		  //return 10; // Amide nitrogen, Thioamide nitrogen (NC=O, NC=S)
		}
	      }
	    }
	  }
	  
	  FOR_NBORS_OF_ATOM (nbr, atom) {
	    if (nbr->IsCarbon()) {
	      int N2count = 0;
	      int N3count = 0;
	      oxygenCount = sulphurCount = 0;

	      FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	        bond = mol -> GetBond(&*nbr, &*nbrNbr);
		if (bond->IsDouble()) {
		  doubleBondTo = nbrNbr->GetAtomicNum();
		}
		if (bond->IsAromatic()) {
		  if ((nbrNbr->GetAtomicNum() == 7) || (nbrNbr->GetAtomicNum() == 6)) {
		    doubleBondTo = nbrNbr->GetAtomicNum();
		  }
		}
		if (bond->IsTriple()) {
		  tripleBondTo = nbrNbr->GetAtomicNum();
		}
		if (nbrNbr->IsNitrogen() && (nbrNbr->GetValence() == 3)) {
	          int nbrOxygen = 0;
		  FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
		    if (nbrNbrNbr->IsOxygen()) {
		      nbrOxygen++;
		    }
		  }
		  if (nbrOxygen < 2) {
		    N3count++;
		  }
		}
		if (nbrNbr->IsNitrogen() && (nbrNbr->GetValence() == 2) && (bond->IsDouble() || bond->IsAromatic())) {
		  N2count++;
		}
		if (nbrNbr->IsAromatic()) {
		  if (nbrNbr->IsOxygen()) {
		    oxygenCount++;
		  }
		  if (nbrNbr->IsSulfur()) {
		    sulphurCount++;
		  }
                }
	      }
	      if (N3count == 3) {
	       return 56; // Guanidinium nitrogen (NGD+)
	      }
	
	      if (!IsAmide && !IsSulfonAmide && !oxygenCount && !sulphurCount && nbr->IsAromatic()) {
	        return 40;
	      }

	      if ((N3count == 2) && (doubleBondTo == 7) && !N2count) {
	        return 55; // Either nitrogen in N+=C-N: (NCN+)
	      }
            }

	    if (nbr->IsNitrogen()) {
	      nitrogenCount = 0;
	      FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	        bond = mol -> GetBond(&*nbr, &*nbrNbr);
		if (bond->IsDouble()) {
		  if (nbrNbr->IsCarbon()) {
		    oxygenCount = sulphurCount = 0;
		    FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
		      if (nbrNbrNbr->IsOxygen()) {
		        oxygenCount++;
		      }
		      if (nbrNbrNbr->IsSulfur()) {
		        sulphurCount++;
		      }
		      if (nbrNbrNbr->IsSulfur()) {
		        nitrogenCount++;
	              }
		    }
		    if (!oxygenCount && !sulphurCount && (nitrogenCount == 1)) {
		      bool bondToAromC = false;
	              FOR_NBORS_OF_ATOM (nbr2, atom) {
		        if (nbr2->IsAromatic() && nbr2->IsCarbon() && nbr2->IsInRingSize(6)) {
			  bondToAromC = true;
			}
		      }
		      if (!bondToAromC) {
		        IsNNNorNNC = true;
		      }
		    }
		  }
		  if (nbrNbr->IsNitrogen()) {
    	            bool bondToAromC = false;
	            FOR_NBORS_OF_ATOM (nbr2, atom) {
	              if (nbr2->IsAromatic() && nbr2->IsCarbon() && nbr2->IsInRingSize(6)) {
			bondToAromC = true;
		      }
		    }
		    if (!bondToAromC) {
	              IsNNNorNNC = true;
		    }
		  }
		}
	      }
	    }	    
          }
          
          if (IsSulfonAmide) {
	    return 43; // Sulfonamide nitrogen (NSO2, NSO3)
	  }
	  if (IsAmide) {
	    return 10; // Amide nitrogen, Thioamide nitrogen (NC=O, NC=S)
	  }
 
	  if ((doubleBondTo ==  6) || (doubleBondTo == 7) ||(doubleBondTo == 15) || (tripleBondTo == 6)) {
	    return 40; // Enamine or aniline nitrogen (deloc. lp), Nitrogen in N-C=N with deloc. lp,
	               // Nitrogen in N-C=N with deloc. lp, Nitrogen attached to C-C triple bond
	               // (NC=C, NC=N, NC=P, NC%C)
	  }
	  if (tripleBondTo == 7) {
	    return 43; // Nitrogen attached to cyano group (NC%N)
	  }
	  if (IsNNNorNNC) {
	    return 10; // Nitrogen in N-N=C moiety with deloc. lp
		       // Nitrogen in N-N=N moiety with deloc. lp (NN=C, NN=N)
	  }
	
	  return 8; // Amine nitrogen (NR)
	}
      }
      // 2 neighbours
      if (atom->GetValence() == 2) {
        if (atom->BOSum() == 4) {
	  FOR_NBORS_OF_ATOM (nbr, atom) {
	    bond = mol -> GetBond(&*nbr, atom);
	    if (bond->IsTriple()) {
	      return 61; // Isonitrile nitrogen (NR%)
	    }
	  }

	  return 53; // Central nitrogen in C=N=N or N=N=N (=N=)
	} 
	
	if (atom->BOSum() == 3) {
	  doubleBondTo = 0;

	  FOR_NBORS_OF_ATOM (nbr, atom) {
	    bond = mol -> GetBond(&*nbr, atom);
	    if (nbr->IsOxygen() && bond->IsDouble() && (nbr->GetValence() == 1)) {
	      return 46; // Nitrogen in nitroso group (N=O)
	    }
	    if ((nbr->IsCarbon() || nbr->IsNitrogen()) && bond->IsDouble()) {
	      return 9; // Iminie nitrogen, Azo-group nitrogen (N=C, N=N)
	    }
          }
	  FOR_NBORS_OF_ATOM (nbr, atom) {
	    if (nbr->IsSulfur()) {
	      oxygenCount = 0;
	      
	      FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
		if (nbrNbr->IsOxygen() && (nbrNbr->GetValence() == 1)) {
		  oxygenCount++;
		}
	      }
              if (oxygenCount >= 2) {
	        return 43; // Sulfonamide nitrogen (NSO2, NSO3)
	      }
            }
          }	
	} 
	
	if (atom->BOSum() == 2) {
	  oxygenCount = sulphurCount = 0;

	  FOR_NBORS_OF_ATOM (nbr, atom) {
	    if (nbr->IsSulfur()) {
	      FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	        if (nbrNbr->IsOxygen() && (nbrNbr->GetValence() == 1)) {
		  oxygenCount++;
		}
	      }
              if (oxygenCount == 1) {
	        return 48; // Divalent nitrogen replacing monovalent O in SO2 group (NSO)
	      }
	    }
	  }

	  return 62; // Anionic divalent nitrogen (NM)
	} 
      }
      // 1 neighbours
      if (atom->GetValence() == 1) {
       	FOR_NBORS_OF_ATOM (nbr, atom) {
	  bond = mol -> GetBond(&*nbr, atom);
	  if (bond->IsTriple()) {
	    return 42; // Triply bonded nitrogen (NSP)
	  }
	  if (nbr->IsNitrogen() && (nbr->GetValence() == 2)) {
	    return 47; // Terminal nitrogen in azido or diazo group (NAZT)
          }
	}
      }
    }

    ////////////////////////////////
    // Oxygen
    ////////////////////////////////
    if (atom->GetAtomicNum() == 8) {
      // 3 neighbours
      if (atom->GetValence() == 3) {
	return 49; // Oxonium oxygen (O+)
      }
      // 2 neighbours
      if (atom->GetValence() == 2) {
        int hydrogenCount = 0;
	FOR_NBORS_OF_ATOM (nbr, atom) {
	  if (nbr->IsHydrogen()) {
	    hydrogenCount++;
	  }
	}

	if (hydrogenCount == 2) {
	  // H--O--H
	  return 70; // Oxygen in water (OH2)
	}
	if (atom->BOSum() == 3) {
	  return 51; // Oxenium oxygen (O=+)
	}
        
	return 6; // Generic divalent oxygen, Ether oxygen, Carboxylic acid or ester oxygen,
	          // Enolic or phenolic oxygen, Oxygen in -O-C=N- moiety, Divalent oxygen in
		  // thioacid or ester, Divalent nitrate "ether" oxygen, Divalent oxygen in
		  // sulfate group, Divalent oxygen in sulfite group, One of two divalent
		  // oxygens attached to sulfur, Divalent oxygen in R(RO)S=O, Other divalent
		  // oxygen attached to sulfur, Divalent oxygen in phosphate group, Divalent
		  // oxygen in phosphite group, Divalent oxygen (one of two oxygens attached
		  // to P), Other divalent oxygen (-O-, OR, OC=O, OC=C, OC=N, OC=S, ONO2, 
		  // ON=O, OSO3, OSO2, OSO, OS=O, -OS, OPO3, OPO2, OPO, -OP)

	// 59 ar
      }
      // 1 neighbour
      if (atom->GetValence() == 1) {
        oxygenCount = sulphurCount = 0;
        
	FOR_NBORS_OF_ATOM (nbr, atom) {
	  bond = mol -> GetBond(&*nbr, atom);

	  if (nbr->IsCarbon() || nbr->IsNitrogen()) {
            FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	      if (nbrNbr->IsOxygen() && (nbrNbr->GetValence() == 1)) {
	        oxygenCount++;
              }
	      if (nbrNbr->IsSulfur() && (nbrNbr->GetValence() == 1)) {
	        sulphurCount++;
	      }
	    }
	  }
          // O---H
	  if (nbr->IsHydrogen()) {
	    return 35;
	  }
          // O-?-C
	  if (nbr->IsCarbon()) {
	    if (oxygenCount == 2) {
	      // O-?-C-?-O
	      return 32; // Oxygen in carboxylate group (O2CM)
            }
	    if (bond->IsSingle()) { 
	      // O--C
	      return 35; // Oxide oxygen on sp3 carbon, Oxide oxygen on sp2 carbon (OM, OM2)
	    } else { 
	      // O==C
	      return 7; // Generic carbonyl oxygen, Carbonyl oxygen in amides,
	                // Carbonyl oxygen in aldehydes and ketones, Carbonyl
			// oxygen in acids or esters (O=C, O=CN, O=CR, O=CO)
	    }
	  }
          // O-?-N
	  if (nbr->IsNitrogen()) {
	    if (oxygenCount >= 2) { 
	      // O-?-N-?-O
	      return 32; // Oxygen in nitro group, Nitro-group oxygen in nitrate,
	                 // Nitrate anion oxygen (O2N, O2NO, O3N)
            }
	    if (bond->IsSingle()) { 
	      // O--N
	      return 32; // Oxygen in N-oxides (ONX)
	    } else { 
	      // O==N
	      return 7; // Nitroso oxygen (O=N)
	    }
	  }
          // O-?-S
	  if (nbr->IsSulfur()) {
	    if (sulphurCount == 1) { 
	      // O1-?-S-?-S1
	      return 32; // Terminal oxygen in thiosulfinate anion (OSMS)
            }
	    if (bond->IsSingle()) { 
	      // O--S
	      return 32; // Single terminal oxygen on sulfur, One of 2 terminal O's on sulfur, 
	                 // One of 3 terminal O's on sulfur, Terminal O in sulfate anion, 
			 // (O-S, O2S, O3S, O4S)
	    } else { 
	      // O==S
	      return 7; // Doubly bonded sulfoxide oxygen, O=S on sulfur doubly bonded 
	                // to, e.g., C (O=S, O=S=)
	    }
	  }

	  return 32; // Oxygen in phosphine oxide, One of 2 terminal O's on sulfur, 
	             // One of 3 terminal O's on sulfur, One of 4 terminal O's on sulfur, 
                     // Oxygen in perchlorate anion (OP, O2P, O3P, O4P, O4Cl)
	}
      }
    }
    
    ////////////////////////////////
    // Flourine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 9) {
      // 1 neighbour
      if (atom->GetValence() == 1) {
        return 11; // Fluorine (F)
      }
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 89; // Fluoride anion (F-)
      }
    }
    
    ////////////////////////////////
    // Sodium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 11) {
      return 93; // Sodium cation (NA+)
    }
    
    ////////////////////////////////
    // Magnesium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 12) {
      return 99; // Dipositive magnesium cation (MG+2)
    }
 
    ////////////////////////////////
    // Silicon
    ////////////////////////////////
    if (atom->GetAtomicNum() == 14) {
      return 19; // Silicon (SI)
    }
 
    ////////////////////////////////
    // Phosphorus
    ////////////////////////////////
    if (atom->GetAtomicNum() == 15) {
      if (atom->GetValence() == 4) {
        return 25; // Phosphate group phosphorus, Phosphorus with 3 attached oxygens,
	           // Phosphorus with 2 attached oxygens, Phosphine oxide phosphorus,
		   // General tetracoordinate phosphorus (PO4, PO3, PO2, PO, PTET)
      }
      if (atom->GetValence() == 3) {
        return 26; // Phosphorus in phosphines (P)
      }
      if (atom->GetValence() == 2) {
        return 75; // Phosphorus doubly bonded to C (-P=C)
      }
    }
    
    ////////////////////////////////
    // Sulfur
    ////////////////////////////////
    if (atom->GetAtomicNum() == 16) {
      // 4 neighbours
      if (atom->GetValence() == 4) {
        return 18; // Sulfone sulfur, Sulfonamide sulfur, Sulfonate group sulfur,
	           // Sulfate group sulfur, Sulfur in nitrogen analog of sulfone 
		   // (SO2, SO2N, SO3, SO4, SNO)
      }
      // 3 neighbours
      if (atom->GetValence() == 3) {
        oxygenCount = sulphurCount = doubleBondTo = 0;

	FOR_NBORS_OF_ATOM (nbr, atom) {
	  bond = mol -> GetBond(&*nbr, atom);
	  if (bond->IsDouble()) {
	    doubleBondTo = nbr->GetAtomicNum();
          }

	  if (nbr->GetValence() == 1) {
	    if (nbr->IsOxygen()) {
	      oxygenCount++;
	    } else if (nbr->IsSulfur()) {
	      sulphurCount++;
	    }
	  } 
	}

	if (oxygenCount == 2) {
	  if (doubleBondTo == 6) {
	    return 18; // Sulfone sulfur, doubly bonded to carbon (=SO2)
	  }
	  return 73; // Sulfur in anionic sulfinate group (SO2M)
	}
	if (oxygenCount && sulphurCount)
	  return 73; // Tricoordinate sulfur in anionic thiosulfinate group (SSOM)
  
	//if ((doubleBondTo == 6) || (doubleBondTo == 8))
	  return 17; // Sulfur doubly bonded to carbon, Sulfoxide sulfur (S=C, S=O)
      }
      // 2 neighbours
      if (atom->GetValence() == 2) {
        doubleBondTo = 0;

	FOR_NBORS_OF_ATOM (nbr, atom) {
	  if (nbr->IsOxygen()) {
	    bond = mol -> GetBond(&*nbr, atom);
	    if (bond->IsDouble()) {
	      doubleBondTo = 8;
            }
	  }
	}

        if (doubleBondTo == 8)
	  return 74; // Sulfinyl sulfur, e.g., in C=S=O (=S=O)
	
	return 15; // Thiol, sulfide, or disulfide sulfor (S)
      }
      // 1 neighbour
      if (atom->GetValence() == 1) {
        sulphurCount = doubleBondTo = 0;

	FOR_NBORS_OF_ATOM (nbr, atom) {
	  FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
	    if (nbrNbr->IsSulfur() && (nbrNbr->GetValence() == 1)) {
	      sulphurCount++;
	    }
	  }
	  bond = mol -> GetBond(&*nbr, atom);
	  if (bond->IsDouble()) {
	    doubleBondTo = nbr->GetAtomicNum();
          }
	}

	if ((doubleBondTo == 6) && (sulphurCount != 2)) {
	  return 16; // Sulfur doubly bonded to carbon (S=C)
        }

        return 72; // Terminal sulfur bonded to P, Anionic terminal sulfur,
	           // Terminal sulfur in thiosulfinate group (S-P, SM, SSMO)
      }

      // 44 ar
    }
    
    ////////////////////////////////
    // Clorine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 17) {
      // 4 neighbour
      if (atom->GetValence() == 4) {
        oxygenCount = 0;
        
	FOR_NBORS_OF_ATOM (nbr, atom) {
	  if (nbr->IsOxygen()) {
	    oxygenCount++;
	  }
	}
	if (oxygenCount == 4)
          return 77; // Perchlorate anion chlorine (CLO4)
      }
      // 1 neighbour
      if (atom->GetValence() == 1) {
        return 12; // Chlorine (CL)
      }
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 90; // Chloride anion (CL-)
      }
    }
    
    ////////////////////////////////
    // Potasium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 19) {
      return 94; // Potasium cation (K+)
    }
    
    ////////////////////////////////
    // Calcium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 20) {
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 96; // Dipositive calcium cation (CA+2)
      }
    }
 
    ////////////////////////////////
    // Iron
    ////////////////////////////////
    if (atom->GetAtomicNum() == 26) {
      return 87; // Dipositive iron (FE+2)
      return 88; // Tripositive iron (FE+3)
    }
    
    ////////////////////////////////
    // Copper
    ////////////////////////////////
    if (atom->GetAtomicNum() == 29) {
      return 97; // Monopositive copper cation (CU+1)
      return 98; // Dipositive copper cation (CU+2)
    }
    
    ////////////////////////////////
    // Zinc
    ////////////////////////////////
    if (atom->GetAtomicNum() == 30) {
      return 95; // Dipositive zinc cation (ZN+2)
    }
 
    ////////////////////////////////
    // Bromine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 35) {
      // 1 neighbour
      if (atom->GetValence() == 1) {
        return 13; // Bromine (BR)
      }
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 91; // Bromide anion (BR-)
      }
    }
 
    ////////////////////////////////
    // Iodine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 53) {
      // 1 neighbour
      if (atom->GetValence() == 1) {
        return 14; // Iodine (I)
      }
    }
 


    return 0;
}	
	
	
    bool IsInSameRing(Atom* a, Atom* b) {
	Molecule *mol = (Molecule *) a -> GetParent ();
    bool a_in, b_in;
    vector<OBRing*> vr;
    vr = mol -> GetSSSR();
    
    vector<OBRing*>::iterator i;
    vector<int>::iterator j;
    
    for (i = vr.begin();i != vr.end();i++) {
      a_in = false;
      b_in = false;
      for(j = (*i)->_path.begin();j != (*i)->_path.end();j++) {
        if ((unsigned)(*j) == a->GetIdx())
          a_in = true;
        if ((unsigned)(*j) == b->GetIdx())
          b_in = true;
      }
      
      if (a_in && b_in)
        return true;
    }
    
    return false;
  };



int get_ds (Atom *at) {
    IntData *d = (IntData *) at -> GetData ("display style");
    assert (d);
    int out = d -> GetGenericValue ();
    if (out == -1) {
    	Molecule *mol = (Molecule *) at -> GetParent ();
    	return mol -> atoms_default_ds;
   	}
   	else return out;
}


void set_ds (Atom *at, int ds) {
    IntData *d = (IntData *) at -> GetData ("display style");
    assert (d);
    d -> SetValue (ds);
}


int get_ds (Bond *b) {
    IntData *d = (IntData *) b -> GetData ("display style");
    assert (d);
    int out = d -> GetGenericValue ();
    if (out == -1) {
    	Molecule *mol = (Molecule *) b -> GetParent ();
    	return mol -> bonds_default_ds;
   	}
   	else return out;
}


void set_ds (Bond *b, int ds) {
    IntData *d = (IntData *) b -> GetData ("display style");
    assert (d);
    d -> SetValue (ds);
}


bool get_selected (Atom *at) {
    BoolData *d = (BoolData *) at -> GetData ("selected");
    assert (d);
    return d -> GetGenericValue ();
}


void set_selected (Atom *at, bool s) {
    BoolData *d = (BoolData *) at -> GetData ("selected");
    assert (d);
    d -> SetValue (s);
}


bool get_selected (Bond *b) {
    bool s1 = get_selected (b -> GetBeginAtom ());
    bool s2 = get_selected (b -> GetEndAtom ());
    return (s1 && s2);
}


void set_selected (Bond *b, bool s) {
    set_selected (b -> GetBeginAtom (), s);
    set_selected (b -> GetEndAtom (), s);
}




double get_score (Atom *at) {
	DoubleData *d = (DoubleData *) at -> GetData ("score");
	assert (d);
	return d -> GetGenericValue ();
}

void set_score (Atom *at, double dd) {
	DoubleData *d = (DoubleData *) at -> GetData ("score");
	assert (d);
	d -> SetValue (dd);
}


double get_back_score (Atom *at) {
	DoubleData *d = (DoubleData *) at -> GetData ("back_score");
	assert (d);
	return d -> GetGenericValue ();
}

void set_back_score (Atom *at, double dd) {
	DoubleData *d = (DoubleData *) at -> GetData ("back_score");
	assert (d);
	d -> SetValue (dd);
}



int CountBonds (Atom *at) {
    int count = 0;
    OBBond *bond;
    OBBondIterator i;
    for (bond = at -> BeginBond(i);bond;bond = at -> NextBond(i))
        count++;

    return(count);
}





Molecule::Molecule () : OBMol () {
    multi = false;
    selection = false;
    atoms_default_ds = 0;
    bonds_default_ds = 1;

}

bool get_sad (Atom *at) {
    BoolData *d = (BoolData *) at -> GetData ("sphere already drawn");
    assert (d);
    return d -> GetGenericValue ();
}


void set_sad (Atom *at, bool s) {
  //  cerr << "set sad" ;
  //  cerr << at -> GetAtomicNum ();
    BoolData *d = (BoolData *) at -> GetData ("sphere already drawn");
    assert (d);
    d -> SetValue (s);
}





void Molecule::find_center () {
    int n = 0;
    double x, y, z;
    x = 0.;
    y = 0.;
    z = 0.;

    FOR_ATOMS_OF_MOL (a, this) {
        vect coor = get_coordinates (&*a);

        n++;
        x += coor. x ();
        y += coor. y ();
        z += coor. z ();
    }
    x /= n;
    y /= n;
    z /= n;
  
    this -> center = vect (x, y, z);
}   


void Molecule::find_limits () {
    double xm, ym, zm, xM, yM, zM;
    xm = ym = zm = 1000000.;
    xM = yM = zM = -1000000.;
	double xt=0., yt=0., zt=0.;
    FOR_ATOMS_OF_MOL (a, this) {

        vect coor = get_coordinates (&*a);
        xt = coor.x();
        yt = coor.y();
        zt = coor.z();
        if (xt < xm) xm = xt;
        if (xt > xM) xM = xt;
        if (yt < ym) ym = yt;
        if (yt > yM) yM = yt;
        if (zt < zm) zm = zt;
        if (zt > zM) zM = zt;
     }
	 if (xm > xM || ym > yM || zm > zM) cerr << "error in calculating limits" <<endl;
     min_corner = vect (xm, ym, zm);
     max_corner = vect (xM, yM, zM);

}

void Molecule::ZNinit_bond (Bond *b) {
    if (!b -> HasData ("display style")) {
        IntData *display_style = new IntData;
        display_style -> SetAttribute ("display style");
        display_style -> SetValue (-1);
        b -> SetData (display_style);
    }
}

void Molecule::ZNinit_atom (Atom *a) {

        if (!a -> HasData ("score")) {
            DoubleData *score = new DoubleData;
            score -> SetAttribute ("score");
            a -> SetData (score);
        }

        if (!a -> HasData ("back_score")) {
            DoubleData *back_score = new DoubleData;
            back_score -> SetAttribute ("back_score");
            a -> SetData (back_score);
        }

        if (!a -> HasData ("display style")) {
            IntData *display_style = new IntData;
            display_style -> SetAttribute ("display style");
            display_style -> SetValue (-1);
            a -> SetData (display_style);
        }

        if (!a -> HasData ("color")) {
            ColorData *col = new ColorData;
            col -> SetAttribute ("color");
            a -> SetData (col);
        }


        if (!a -> HasData ("f")) {
            VectorData *force = new VectorData;
            force -> SetAttribute ("f");
            a -> SetData (force);
        }
		
		if (!a -> HasData ("back_force")) {
            VectorData *force = new VectorData;
            force -> SetAttribute ("back_force");
            a -> SetData (force);
        }
	
	if (!a -> HasData ("force_mutex")) {
		MutexData *mut = new MutexData;
		QMutex *mutex = new QMutex;
		mut -> SetValue (mutex);
		mut -> SetAttribute ("force_mutex");
		a -> SetData (mut);
	}


        if (!a -> HasData ("selected")) {
            OpenBabel::OBPairTemplate< bool > *selected= new OpenBabel::OBPairTemplate< bool>;
            selected -> SetAttribute ("selected");
            selected -> SetValue (false);
            a -> SetData (selected);
        }


        if (!a -> HasData ("sphere already drawn")) {
            BoolData *sad= new BoolData;
            sad -> SetAttribute ("sphere already drawn");
            sad -> SetValue (false);
            a -> SetData (sad);
        }

        if (!a -> HasData ("visible")) {
            OpenBabel::OBPairTemplate< bool > *visible = new OpenBabel::OBPairTemplate< bool >;
            visible -> SetAttribute ("visible");
            visible -> SetValue (true);
            a -> SetData (visible);
        }


        set_color_mw (&*a) ;


   //     a -> displayStyle = 0;
    //    stringstream ss;
   //     ss <<(GetEnergy ());
   //     cerr <<ss.str ()<<endl;
 //       double d = 0.000003;
  //      cerr << "d "<<d;

    //        a -> vdw = etab.GetVdwRad (a -> GetAtomicNum ());
        //    cerr << a->vdw;
  //       a -> col = a -> get_color_mw () ;
}



void Molecule::ZNinit () {
    needs_redraw = false; needs_recolor = false;
	cerr << HasNonZeroCoords () << endl;
	if (!HasNonZeroCoords ()) mend_coordinates (this);

    FOR_ATOMS_OF_MOL (a, this) {
        ZNinit_atom (&*a);
     }
     FOR_BONDS_OF_MOL (b, this) { 

        ZNinit_bond (&*b);
    }
}

Bond * Molecule::ZNAddBond (Bond *bond) {
    int first = bond -> GetBeginAtomIdx();
    int second = bond -> GetEndAtomIdx();
    int order = bond -> GetBO();
    int flags = bond -> GetFlags();
    int insertpos = -1;




    if (first == second)
      return(NULL);

    //    BeginModify();

    if ((unsigned)first <= NumAtoms() && (unsigned)second <= NumAtoms()
        && !GetBond(first, second))
      //atoms exist and bond doesn't
      {
        if (!bond)
          {
            EndModify();
            return(NULL);
          }

        OBAtom *bgn,*end;
        bgn = GetAtom(first);
        end = GetAtom(second);
        if (!bgn || !end)
          {
            obErrorLog.ThrowError(__FUNCTION__, "Unable to add bond - invalid atom index", obDebug);
            return(NULL);
          }
        bond->Set(_nbonds,bgn,end,order,flags);
        bond->SetParent(this);

        //set aromatic flags if it has the appropriate order
        if (order == 5)
          {
            bond->SetAromatic();
            bgn->SetAromatic();
            end->SetAromatic();
          }

#define OBBondIncrement 100
        if (_vbond.empty() || _nbonds+1 >= _vbond.size())
          {
            _vbond.resize(_nbonds+OBBondIncrement);
            vector<OBBond*>::iterator i;
            for (i = _vbond.begin(),i+=(_nbonds+1);i != _vbond.end();++i)
              *i = (OBBond*)NULL;
          }
#undef  OBBondIncrement

        _vbond[_nbonds] = (OBBond*)bond;
        _nbonds++;

        if (insertpos == -1)
          {
            bgn->AddBond(bond);
            end->AddBond(bond);
          }
        else
          {
            if (insertpos >= static_cast<int>(bgn->GetValence()))
              bgn->AddBond(bond);
            else //need to insert the bond for the connectivity order to be preserved
              {    //otherwise stereochemistry gets screwed up
                vector<OBBond*>::iterator bi;
                bgn->BeginNbrAtom(bi);
                bi += insertpos;
                bgn->InsertBond(bi,bond);
              }
            end->AddBond(bond);
          }
      }
    else //at least one atom doesn't exist yet - add to bond_q
      SetData(new OBVirtualBond(first,second,order,flags));

    //    EndModify();
    return(bond);

}
void Molecule::ZNSetConformers () {
    int n = NumConformers ();
    if (n) {
    cerr <<  "ZNSetConformers" << endl;
    //clear out the multiconformer data
        vector<double*>::iterator k;
        for (k = _vconf.begin();k != _vconf.end();++k)
          delete [] *k;
    _vconf.clear();
    }

    

    AddConformer (_c);
}

  bool Molecule::ZNAddHydrogens(Atom *atom)
  {
   // cerr << "1" << endl;
    OBAtom *h;

    //count up number of hydrogens to add
    int hcount,count=0;
    vector<pair<OBAtom*,int> > vhadd;
 //   cerr << "2" << endl;
    hcount = atom->GetImplicitValence() - atom->GetValence();

    //Jan 05 Implicit valency now left alone; use spin multiplicity for implicit Hs
    int mult = atom->GetSpinMultiplicity();
  //  cerr << "3" << endl;
    if(mult==2) //radical
      hcount-=1;
    else if(mult==1 || mult==3) //carbene
      hcount-=2;

    if (hcount < 0)
      hcount = 0;
    if (hcount)
      {
        vhadd.push_back(pair<OBAtom*,int>(atom,hcount));
        count += hcount;
      }

    if (count == 0)
      return(true);
 //   cerr << "4" << endl;

 /*   //realloc memory in coordinate arrays for new hydroges
    double *tmpf;
    vector<double*>::iterator j;
    for (j = _vconf.begin();j != _vconf.end();++j)
      {
        cerr << "ok "<< endl;
        tmpf = new double [(NumAtoms()+count)*3+10];
        cerr << *j << endl;
        memcpy(tmpf,(*j),sizeof(double)*NumAtoms()*3);
        cerr << "ok "<< endl;
        delete []*j;
        cerr << "ok "<< endl;
        *j = tmpf;
        cerr << "ok "<< endl;
      }
*/
    IncrementMod();
 //   cerr << "5" << endl;
    int m,n=0;
    vector3 v;
    vector<pair<OBAtom*,int> >::iterator k;
    double hbrad = etab.CorrectedBondRad(1,0);

    for (k = vhadd.begin();k != vhadd.end();++k)
      {
   //     cerr << "ok" << endl;
        atom = k->first;
        double bondlen = hbrad+etab.CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
        for (m = 0;m < k->second;++m)

          {
     //   cerr << "ok1" << endl;
  //          for (n = 0;n < NumConformers();++n)


//              {
    //    cerr << "ok2" << endl;
                SetConformer(n);
      //          cerr <<"vx"<< endl;
                atom->GetNewBondVector(v,bondlen);

   //             _c[(NumAtoms())*3]   = v.x();
   //             _c[(NumAtoms())*3+1] = v.y();
  //              _c[(NumAtoms())*3+2] = v.z();
        //        cerr << v.x ()<< "vx";
  //            }
            h = NewAtom();
            h->SetType("H");
            h->SetAtomicNum(1);
            AddBond(atom->GetIdx(),h->GetIdx(),1);
            h->SetVector(v);
          }
      }

    DecrementMod();
    SetConformer(0);

    //reset atom type and partial charge flags
    //_flags &= (~(OB_PCHARGE_MOL|OB_ATOMTYPES_MOL));

    return(true);
  }


Atom * Molecule::ZNAddAtom (Atom *atom) {

    //    BeginModify();
    ZNinit_atom (atom);
    OBAtom *obatom = CreateAtom();
    obatom = atom;
  //  cerr << "color" << obatom -> HasData ("color") << endl;
    obatom->SetIdx(_natoms+1);
    obatom->SetParent(this);


#define OBAtomIncrement 100

    if (_vatom.empty() || _natoms+1 >= _vatom.size())
      {
        _vatom.resize(_natoms+OBAtomIncrement);
        vector<OBAtom*>::iterator j;
        for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();++j)
          *j = (OBAtom*)NULL;
      }
#undef OBAtomIncrement

    _vatom[_natoms] = (OBAtom*)obatom;
    _natoms++;

    if (HasData(OBGenericDataType::VirtualBondData))
      {
        /*add bonds that have been queued*/
        OBVirtualBond *vb;
        vector<OBGenericData*> verase;
        vector<OBGenericData*>::iterator i;
        for (i = BeginData();i != EndData();++i)
          if ((*i)->GetDataType() == OBGenericDataType::VirtualBondData)
            {
              vb = (OBVirtualBond*)*i;
              if (static_cast<unsigned int>(vb->GetBgn()) > _natoms ||
                  static_cast<unsigned int>(vb->GetEnd()) > _natoms)
                continue;
              if (obatom->GetIdx() == static_cast<unsigned int>(vb->GetBgn())
                  || obatom->GetIdx() == static_cast<unsigned int>(vb->GetEnd()))
                {
                  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
                  verase.push_back(*i);
                }
            }

        if (!verase.empty())
          DeleteData(verase);
      }

    //    EndModify();

    return(obatom);

}



int Molecule::bonded_to (Atom *at, int bondtype, int atnumb) {
    int out = 0;

    
    if (bondtype!=5) {
        FOR_NBORS_OF_ATOM (n, at) {
            int an = n -> GetAtomicNum ();
            if (an == atnumb || atnumb == -2) {
                Bond *bond = at -> GetBond (&*n);
                if (bond -> GetBO () == bondtype || bondtype == -2) {
                    out += 1;
                }
            }
            
        }
    }
    return out;
}



bool Molecule::RemoveAtom (Atom *atom) {
////    if (atom->IsHydrogen())
////      return(DeleteHydrogen(atom));

    BeginModify();
    //don't need to do anything with coordinates b/c
    //BeginModify() blows away coordinates

    //find bonds to delete
    OBAtom *nbr;
    vector<OBBond*> vdb;
    vector<OBBond*>::iterator j;
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      vdb.push_back(*j);

    for (j = vdb.begin();j != vdb.end();++j)
      RemoveBond((OBBond *)*j); //delete bonds

    _vatom.erase(_vatom.begin()+(atom->GetIdx()-1));
    _natoms--;

    //reset all the indices to the atoms
    int idx;
    vector<OBAtom*>::iterator i;
    OBAtom *atomi;
    for (idx=1,atomi = BeginAtom(i);atomi;atomi = NextAtom(i),++idx)
      atomi->SetIdx(idx);

    EndModify();

 //   DestroyAtom(atom);

    return(true);
}
/*
bool Molecule::RemoveHydrogen (Atom *atom) {
  //deletes the hydrogen atom passed to the function
    if (!atom->IsHydrogen())
      return false;

    //find bonds to delete
    OBAtom *nbr;
    vector<OBBond*> vdb;
    vector<OBBond*>::iterator j;
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      vdb.push_back(*j);

    IncrementMod();
    for (j = vdb.begin();j != vdb.end();++j)
      DeleteBond((OBBond*)*j); //delete bonds
    DecrementMod();

    int idx;
    if (atom->GetIdx() != NumAtoms())
      {
        idx = atom->GetCIdx();
        int size = NumAtoms()-atom->GetIdx();
        vector<double*>::iterator k;
        for (k = _vconf.begin();k != _vconf.end();++k)
          memmove((char*)&(*k)[idx],(char*)&(*k)[idx+3],sizeof(double)*3*size);

      }

    _vatom.erase(_vatom.begin()+(atom->GetIdx()-1));
    _natoms--;

    //reset all the indices to the atoms
    vector<OBAtom*>::iterator i;
    OBAtom *atomi;
    for (idx=1,atomi = BeginAtom(i);atomi;atomi = NextAtom(i),++idx)
      atomi->SetIdx(idx);
    
    UnsetHydrogensAdded();

//    DestroyAtom(atom);
    return(true);

}
*/
bool Molecule::RemoveBond (Bond *bond) {
    BeginModify();

    (bond->GetBeginAtom())->DeleteBond(bond);
    (bond->GetEndAtom())->DeleteBond(bond);
    _vbond.erase(_vbond.begin() + bond->GetIdx()); // bond index starts at 0!!!
    _nbonds--;

    vector<OBBond*>::iterator i;
    int j;
    OBBond *bondi;
    for (bondi = BeginBond(i),j=0;bondi;bondi = NextBond(i),++j)
      bondi->SetIdx(j);

    EndModify();

  //  DestroyBond(bond);

    return(true);
}



color Molecule::get_color_mw (Atom *at) {
//etable should have color data as well...
    vector <double> colors;
    colors = etab.GetRGB (at -> GetAtomicNum ());
  /*  if (at -> GetAtomicNum () == 6) {
        colors [0] = 0.;
        colors [1] = 0.;
        colors [2] = 0.;
    }*/
    return color ((float) colors[0], (float)colors[1],(float) colors[2]);
/*
    float r, g, b;
	switch (at -> GetAtomicNumber) {
		case 1: r=1.0f; g=1.0f; b=1.0f;
			break;
		case 6: r=0.f; g=0.f;  b=0.f;

			break;	
		case 7: r=0.0f; g=0.0f;  b=1.0f;
			break;
		case 8: r=1.0f; g=0.0f;  b=0.0f;
			break;
        case 15: r=0.f; g=1.0f; b=0.5f;
            break;
		case 16: r=1.0f; g=1.0f;  b=0.0f;
			break;	
			
		default: r=0.7f; g= 0.7f;  b=0.7f;
	}
    return color (r, g, b);
*/
}

int Molecule::get_ds_from_neighbour (Atom *at) 
{
	int out = 1;
	Atom *n;
	OBBondIterator i;
	n =  at -> BeginNbrAtom	(i);
	if (n) return get_ds (n);
	return out;	 	
} 


int Molecule::get_ds_from_neighbour (Bond *bo) 
{
	int out = 1;
	Atom *at = bo -> GetBeginAtom ();
	FOR_BONDS_OF_ATOM (b, at) {
		if (&*b != &*bo) return get_ds (&*b);
	}
	at = bo -> GetEndAtom ();
	FOR_BONDS_OF_ATOM (b, at) {
		if (&*b != &*bo) return get_ds (&*b);
	}
	return out;
	return out;
} 
	

void Molecule::ZNAddHydrogens () {
	AddHydrogens ();
	ZNinit ();
	FOR_ATOMS_OF_MOL (at, this) {
		if (at -> IsHydrogen ()) set_ds (&*at, get_ds_from_neighbour (&*at));
	}
	FOR_BONDS_OF_MOL (bo, this) {
		if (bo -> GetBeginAtom () -> IsHydrogen () || bo -> GetEndAtom () -> IsHydrogen ()) set_ds (&*bo, get_ds_from_neighbour (&*bo));
	}
}


void Molecule::set_color_mw (Atom *at){
    set_color (at, get_color_mw (at));
}


void Molecule::add_atom_bonded_to (vect coords, int atomnum, Atom *at) {
    Atom *at2 = new Atom;
    ZNinit_atom (at);
    at2 -> SetVector (coords);
    at2 -> SetAtomicNum (atomnum);
    Bond *b = new Bond;
    ZNinit_bond (b);
    b -> SetBondOrder (1);
    add_atom_bonded_to (at2, b, at);
}



void Molecule::add_atom_bonded_to (Atom *to_add, Bond *bond, Atom *partner) {
 //   to_add -> getVdw ();
    if (!bond -> GetBondOrder ()) bond -> SetBondOrder (1);
    assert (partner -> GetParent () == this);


    ZNAddAtom (to_add);
    set_ds (to_add, get_ds (partner));
    set_color_mw (to_add);

    bond -> SetBegin (to_add);
    bond -> SetEnd (partner);
    ZNAddBond (bond);

   // AddBond (to_add -> GetIdx (), partner -> GetIdx (), 1);
//    bond = GetBond (to_add -> GetIdx (), partner -> GetIdx ());
    ZNinit_bond (bond);
    if (NumBonds ()) {
        //change to get ds from first bond.
        Bond *first_b = GetBond (0); //not safe
        set_ds (bond, get_ds (first_b));
    }
    else set_ds (bond, 1);
}















/*
using namespace std;

class Mol2Par;


*/
Selection::Selection () : Molecule () {
    selection = true;
    SetTitle ("Selection");
}



void Selection::deselect () {
    FOR_ATOMS_OF_MOL (a, this) {
        set_selected (&*a, false);
    }
}


void Selection::add_mol (Molecule *mol) {
    bool add = true;
    for (unsigned int i=0; i<molecules.size (); i++) {
        if (molecules[i] == mol) {
            add = false;
            break;
        }
    }
    if (add) molecules.push_back (mol);
}




Bond * Selection::ZNAddBondToSelection (Bond *bond) {
    int first  = bond -> GetBeginAtomIdx();
    int second = bond -> GetEndAtomIdx();
    int order  = bond -> GetBO();
    int flags  = bond -> GetFlags();
    int insertpos = -1;




    if (first == second)
      return(NULL);

    //    BeginModify();

    if ((unsigned)first <= NumAtoms() && (unsigned)second <= NumAtoms()
        && !GetBond(first, second))
      //atoms exist and bond doesn't
      {
        if (!bond)
          {
            EndModify();
            return(NULL);
          }

        OBAtom *bgn,*end;
        bgn = GetAtom(first);
        end = GetAtom(second);
        if (!bgn || !end)
          {
            obErrorLog.ThrowError(__FUNCTION__, "Unable to add bond - invalid atom index", obDebug);
            return(NULL);
          }
   //     bond->Set(_nbonds,bgn,end,order,flags);
    //    bond->SetParent(this);

        //set aromatic flags if it has the appropriate order
        if (order == 5)
          {
            bond->SetAromatic();
            bgn->SetAromatic();
            end->SetAromatic();
          }

#define OBBondIncrement 100
        if (_vbond.empty() || _nbonds+1 >= _vbond.size())
          {
            _vbond.resize(_nbonds+OBBondIncrement);
            vector<OBBond*>::iterator i;
            for (i = _vbond.begin(),i+=(_nbonds+1);i != _vbond.end();++i)
              *i = (OBBond*)NULL;
          }
#undef  OBBondIncrement

        _vbond[_nbonds] = (OBBond*)bond;
        _nbonds++;

        if (insertpos == -1)
          {
      //      bgn->AddBond(bond);
     //       end->AddBond(bond);
          }
        else
          {
            if (insertpos >= static_cast<int>(bgn->GetValence()))
              bgn->AddBond(bond);
            else //need to insert the bond for the connectivity order to be preserved
              {    //otherwise stereochemistry gets screwed up
                vector<OBBond*>::iterator bi;
                bgn->BeginNbrAtom(bi);
                bi += insertpos;
       //         bgn->InsertBond(bi,bond);
              }
      //      end->AddBond(bond);
          }
      }
    else //at least one atom doesn't exist yet - add to bond_q
      SetData(new OBVirtualBond(first,second,order,flags));

    //    EndModify();
    return(bond);

}





Atom * Selection::ZNAddAtomToSelection (Atom *atom) {

    //    BeginModify();
//    ZNinit_atom (atom);
    OBAtom *obatom = atom;
  //  cerr << "color" << obatom -> HasData ("color") << endl;
  //  obatom->SetIdx(_natoms+1);
//    obatom->SetParent(this);


#define OBAtomIncrement 100

    if (_vatom.empty() || _natoms+1 >= _vatom.size())
      {
        _vatom.resize(_natoms+OBAtomIncrement);
        vector<OBAtom*>::iterator j;
        for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();++j)
          *j = (OBAtom*)NULL;
      }
#undef OBAtomIncrement

    _vatom[_natoms] = (OBAtom*)obatom;
    _natoms++;

//    if (HasData(OBGenericDataType::VirtualBondData))
//      {
//        /*add bonds that have been queued*/
//        OBVirtualBond *vb;
//        vector<OBGenericData*> verase;
//        vector<OBGenericData*>::iterator i;
//        for (i = BeginData();i != EndData();++i)
//          if ((*i)->GetDataType() == OBGenericDataType::VirtualBondData)
//            {
//              vb = (OBVirtualBond*)*i;
//              if (static_cast<unsigned int>(vb->GetBgn()) > _natoms ||
//                  static_cast<unsigned int>(vb->GetEnd()) > _natoms)
//                continue;
//              if (obatom->GetIdx() == static_cast<unsigned int>(vb->GetBgn())
//                  || obatom->GetIdx() == static_cast<unsigned int>(vb->GetEnd()))
//                {
//                  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
//                  verase.push_back(*i);
//                }
//            }
//
//        if (!verase.empty())
//          DeleteData(verase);
//      }

    //    EndModify();

    return(obatom);

}















/*
Molecule::Molecule(Molecule const& mol) {

    cerr << "copy";
    valid = mol.valid;
    selection = mol.selection;
    multi = mol.multi;
    valid = mol.valid;
   
    center = mol.center;
    min_corner = mol.min_corner;
    max_corner = mol.max_corner;
    name = mol.name;
    description = mol.description; 
    
    line_list = mol.line_list;
    stick_list = mol.stick_list;
    surface_list = mol.surface_list;
    backbone_list1 = mol.backbone_list1;
    backbone_list2 = mol.backbone_list2;



    for (unsigned int i=0; i<mol.atoms.size (); i++) atoms.push_back (mol.atoms[i]);

    for (unsigned int i=0; i<mol.bonds.size (); i++) bonds.push_back (mol.bonds[i]);
    cerr <<"bonds "<<bonds.size ()<<endl;
    for (unsigned int i=0; i<mol.residues.size (); i++) residues.push_back (mol.residues[i]);
    for (unsigned int i=0; i<mol.rings.size (); i++) rings.push_back (mol.rings[i]);
    for (unsigned int i=0; i<mol.fragments.size (); i++) fragments.push_back (mol.fragments[i]);
    for (unsigned int i=0; i<mol.res_strings.size (); i++) res_strings.push_back (mol.res_strings[i]);


}
*/

/*

Molecule& Molecule::operator=(const Molecule& mol)
{

  //  cerr << "copy=";
    valid = mol.valid;
    selection = mol.selection;
    multi = mol.multi;
    valid = mol.valid;
   
    center = mol.center;
    min_corner = mol.min_corner;
    max_corner = mol.max_corner;
    name = mol.name;
    description = mol.description; 
    
    line_list = mol.line_list;
    stick_list = mol.stick_list;
    surface_list = mol.surface_list;
    backbone_list1 = mol.backbone_list1;
    backbone_list2 = mol.backbone_list2;

 //   cerr <<"bonds "<<bonds.size ()<<endl;


    for (unsigned int i=0; i<mol.atoms.size (); i++) atoms.push_back (mol.atoms[i]);

    for (unsigned int i=0; i<mol.bonds.size (); i++) bonds.push_back (mol.bonds[i]);

    for (unsigned int i=0; i<mol.residues.size (); i++) residues.push_back (mol.residues[i]);
    for (unsigned int i=0; i<mol.rings.size (); i++) rings.push_back (mol.rings[i]);
    for (unsigned int i=0; i<mol.fragments.size (); i++) fragments.push_back (mol.fragments[i]);
    for (unsigned int i=0; i<mol.res_strings.size (); i++) res_strings.push_back (mol.res_strings[i]);

  //  cerr <<"bond "<<mol.bonds[0]->number <<endl;
    return *this;
}





Molecule::Molecule () 
{
    needs_redraw = false; needs_recolor = false;
    valid=false; selection=false; multi=false;
  //  atoms.clear (); bonds.clear (); rings.clear (); residues.clear ();
}


Molecule::~Molecule() {
    if (!selection) {
        for (unsigned int i=0; i<atoms.size (); i++) {
            delete atoms[i];
        }
        for (unsigned int i=0; i<bonds.size (); i++) {
            delete bonds[i];
        }
        for (unsigned int i=0; i<residues.size (); i++) {
            delete residues[i];
        }
        for (unsigned int i=0; i<rings.size (); i++) {
            delete rings[i];
        }
        for (unsigned int i=0; i<fragments.size (); i++) {
            delete fragments[i];

        }
    }
}




void Molecule::copy_from (Molecule *mol)
{

    cerr << "copy=";
    valid = mol -> valid;
    selection = mol -> selection;
    multi = mol -> multi;

   
    center = mol -> center;
    min_corner = mol -> min_corner;
    max_corner = mol -> max_corner;
    name = mol -> name;
    description = mol -> description; 
    
    line_list = mol -> line_list;
    stick_list = mol -> stick_list;
    surface_list = mol -> surface_list;
    backbone_list1 = mol -> backbone_list1;
    backbone_list2 = mol -> backbone_list2;

    cerr <<"bonds "<<bonds.size ()<<endl;


    for (unsigned int i=0; i<mol -> atoms.size (); i++) atoms.push_back (mol -> atoms[i]);

    for (unsigned int i=0; i<mol -> bonds.size (); i++) bonds.push_back (mol -> bonds[i]);

    for (unsigned int i=0; i<mol -> residues.size (); i++) residues.push_back (mol -> residues[i]);
    for (unsigned int i=0; i<mol -> rings.size (); i++) rings.push_back (mol -> rings[i]);
    for (unsigned int i=0; i<mol -> fragments.size (); i++) fragments.push_back (mol -> fragments[i]);
    for (unsigned int i=0; i<mol -> res_strings.size (); i++) res_strings.push_back (mol -> res_strings[i]);

    cerr <<"bond "<<mol -> bonds[0]->number <<endl;

}












void Atom::find_mol2_type () {
    switch (atomicNumber) {
        case 1:
            atomType = "H";            
        break;
        case 6:
            if (is_aromatic ()) atomType == "C.ar";
            else if (bound.size ()== 2) atomType == "C.1";
            else if (bound.size ()== 3) atomType == "C.2";
            else if (bound.size ()== 4) atomType == "C.3";
        break;
        case 7:
            if (is_aromatic ()) atomType == "N.ar";
            else if (bound.size ()== 2) atomType == "N.1";
            else if (bound.size ()== 3) atomType == "N.2";
            else if (bound.size ()== 4) atomType == "N.3";
        case 8:
            if (is_aromatic ()) atomType == "O.ar";
            else if (bound.size ()== 1) atomType == "O.2";
            else if (bound.size ()== 2) atomType == "O.3";
        break;

    }
}

unsigned int Atom::bonded_to (int bondtype, int atnumb) {
    unsigned int out = 0;
    if (bondtype!=5) {
        for (unsigned int n=0; n<bonds.size(); n++) {
            if (bonds[n]->kekule == bondtype || bondtype == -2) {
                if (bound[n]->atomicNumber == atnumb || atnumb == -2) {
                    out += 1;
                }
            }
        }
    }
    else {
        for (unsigned int n=0; n<bonds.size(); n++) {
            if (bonds[n]->mol2Type == bondtype || bondtype == -2) {
                if (bound[n]->atomicNumber == atnumb || atnumb == -2) {
                    out += 1;
                }
            }
        }
    }
    return out;
}

unsigned int Atom::bonded_to (string MMFFstring) {
    unsigned int out = 0;
    for (unsigned int n=0; n<bound.size(); n++) {
        if (bound[n]->MMFFstring == MMFFstring) {
            out += 1;
            
        }
    }
    return out;
}

bool Atom::bonded_to (Atom *at) {
    for (unsigned int n=0; n<bound.size(); n++) {
        if (bound[n]== at) {
            return true;
            
        }
    }
    return false;
}



bool Atom::is_in_ring (Ring *rin) {
    if (!in_ring.size ()) return false;
    else {
        for (unsigned int r=0; r<in_ring.size (); r++) {
            if (in_ring[r]==rin) return true;
        }
    }
    return false;
}

bool Atom::is_aromatic () {
    bool ar = false;
    for (unsigned r=0; r<in_ring.size (); r++) {
        if (in_ring[r]->aromatic) ar =true;
    }
    return ar;
}


bool Atom::isSp2 (){
    if (mol2Type==2 || mol2Type==3 || mol2Type==6 || mol2Type==9 || mol2Type==11 || mol2Type==18)
            return true; 
    else return false;
}


color Atom::get_color_mw () {
    float r, g, b;
	switch (atomicNumber) {
		case 1: r=1.0f; g=1.0f; b=1.0f;
			break;
		case 6: r=0.f; g=0.f;  b=0.f;

			break;	
		case 7: r=0.0f; g=0.0f;  b=1.0f;
			break;
		case 8: r=1.0f; g=0.0f;  b=0.0f;
			break;
        case 15: r=0.f; g=1.0f; b=0.5f;
            break;
		case 16: r=1.0f; g=1.0f;  b=0.0f;
			break;	
			
		default: r=0.7f; g= 0.7f;  b=0.7f;
	}
    return color (r, g, b);
}

void Atom::set_color_mw (){
    col = get_color_mw ();
}






Bond::Bond () {
    fr_visited = false;
    fr_parent = NULL;
}

bool Bond::is_aromatic () {
    bool ar = false;
    for (unsigned r=0; r<in_ring.size (); r++) {
        if (in_ring[r]->aromatic) ar =true;
    }
    return ar;
}


bool Bond::is_rotatable () {
    if (mol2Type != 1) return false;       //multiple bond
    if (atomPTR[0]->bound.size () == 1 || atomPTR[1]->bound.size () == 1) return false;     //to terminal atom
    if (in_ring.size ()) return false;    //ring bond
    if (atomPTR[0]->atomicNumber == 6 && atomPTR[0]->bonded_to (1, 1) == 3) return false; //to terminal CH3
    if (atomPTR[1]->atomicNumber == 6 && atomPTR[1]->bonded_to (1, 1) == 3) return false; //to terminal CH3
    return true;    
}


bool Bond::conjugated_to (unsigned int order, Ring *ring) {
    for (unsigned int b=0; b<atomPTR[0]->bonds.size (); b++) {
        if (atomPTR[0]->bonds[b]->kekule == order && atomPTR[0]->bonds[b]!=this && atomPTR[0]->bonds[b]->is_in_ring (ring)) return true;
    }
    for (unsigned int b=0; b<atomPTR[1]->bonds.size (); b++) {
        if (atomPTR[1]->bonds[b]->kekule == order && atomPTR[1]->bonds[b]!=this && atomPTR[1]->bonds[b]->is_in_ring (ring)) return true;
    }
    return false;
}

bool Bond::is_in_ring (Ring *rin) {
    if (!in_ring.size ()) return false;
    else {
        for (unsigned int r=0; r<in_ring.size (); r++) {
            if (in_ring[r]==rin) return true;
        }
    }
    return false;
}




void Bond::set_order (int ord) {
    mol2Type = ord;
    kekule = ord;
}






Atom::Atom () {
    col = color (0, 0, 0, 255);
    visible = true;
    selected = false;
    score = 0.f;
    no2bond = false;
	inBackBone = false;
    sphere_already_drawn = false;
    displayStyle = 0;
    fragment = NULL;
    residue = NULL;
    ID = 0;
    atomicNumber = 6;
    substructureNumber = 0;
    substructureName = "";

    vdw = 1.0;
				

}

void Atom::getVdw () {
    vdw = mol2par->getVdw(atomicNumber);
}

/*

int Molecule::readPDB (string filename, ifstream&) {
    string PDBAtom = "HETATM";
    string PDBConnect = "CONECT";
    string PDBEnd = "END";

	string buffer ="";
	unsigned int lineNumber = 0;
	name = filename; 
	bool go_on = TRUE;
    go_on = getline(file, buffer); 
	while (go_on) {
		lineNumber++;// cout <<lineNumber;
		istringstream issGlobal(buffer);
		string token;
		issGlobal >> token;		

		if (token == PDBAtom) {}
        else if (token == PDBConnect) {}
        else if (token == PDBEnd) {
            go_on = false;
        }
    }


}

*/
/*


int Molecule::readMultiMOL2(string filename, ifstream& file, bool& continueRead)
{
	string triposMolecule = "@<TRIPOS>MOLECULE";
	string triposAtom = "@<TRIPOS>ATOM";
	string triposSubStructure = "@<TRIPOS>SUBSTRUCTURE";
	string triposBACKBONE = "BACKBONE";
	string triposDICT = "DICT";
	string triposESSENTIAL = "ESSENTIAL";
	string triposDIRECT = "DIRECT";
	string triposINTERRES = "INTERRES";
	string triposBond = "@<TRIPOS>BOND";

		
	
	string buffer ="";
	unsigned int lineNumber = 0;
	name = filename; 
	
	
	bool newMoleculeFound = false;
	bool nextMoleculeFound = false;
	
	bool res = true;
	if (!continueRead) {
		res = getline(file, buffer);  


	} 
	while (res && !nextMoleculeFound) {

		lineNumber++;// cout <<lineNumber;
		istringstream issGlobal(buffer);
		string token;
		issGlobal >> token;		
        int atomCount, bondCount, substructCount;
   //     cout << token <<" "<<lineNumber<<endl;
		if (token == triposMolecule || continueRead) {
			if (!newMoleculeFound) {
	            newMoleculeFound = true;
				continueRead = false;
                getline(file, buffer);
                lineNumber++;
				description = buffer;
		//		cout << "Description: " << description << endl;
				getline(file, buffer);
				lineNumber++;
				istringstream iss(buffer);
				iss >> atomCount;
				if (!checkInputStream(filename, "Atom count expected.", lineNumber, iss)) {
					return false;
				}
				iss >> bondCount;
				if (!checkInputStream(filename, "Bond count expected.", lineNumber, iss)) {
					return false;
				}
				iss >> substructCount;
	///			/*if (!checkInputStream(filename, "Substructure count expected.", lineNumber, iss)) {
	///				return false;
	///			} 
			}
			else {
				nextMoleculeFound = true;
	//			continueRead = true;
			}
		}
		else if (token == triposAtom) {
			for (unsigned int i = 0; i < atomCount; i++) {
				getline(file, buffer);
				lineNumber++;
				istringstream iss(buffer);
				
				unsigned int atomNumber = 0;
				string atomName = "";
				string atomType = "";
				unsigned int atomSubstructureNumber = 1;
				string atomSubstructureName = "";
				float atomCharge = 0.0;
                float x, y, z;
				iss >> atomNumber;
				if (!checkInputStream(filename, "Atom number expected.", lineNumber, iss)) {
					return false;
				}
				iss >> atomName;
				if (!checkInputStream(filename, "Atom name expected.", lineNumber, iss)) {
					return false;
				}
				iss >> x;
				if (!checkInputStream(filename, "Atom coordinate X expected.", lineNumber, iss)) {
					return false;
				}
				iss >> y; 
				if (!checkInputStream(filename, "Atom coordinate Y expected.", lineNumber, iss)) {
					return false;
				}
	        	iss >> z;
				if (!checkInputStream(filename, "Atom coordinate Z expected.", lineNumber, iss)) {
					return false;
				}
				iss >> atomType;
				if (!checkInputStream(filename, "Atom type expected.", lineNumber, iss)) {
					return false;
				}
				iss >> atomSubstructureNumber; 
				if (!checkInputStream(filename, "Atom substructure ID expected.", lineNumber, iss)) {
					return false;
				}
				iss >> atomSubstructureName;
				if (!checkInputStream(filename, "Atom substructure name expected.", lineNumber, iss)) {
					return false;
				}
				iss >> atomCharge;
				if (!checkInputStream(filename, "Atom charge expected.", lineNumber, iss)) {
					return false;
				}
				
				unsigned int mol2Type = mol2par->getmol2AtomTypeID(atomType);
		///		if (mol2Type == 20) {
		////			cerr << "REMOVED LONE PAIR: " << atomNumber << endl;
		/////			cout << "removed lone pair " << setw(5) << setfill (' ') << atomName <<  "(" << setw(4) << setfill (' ') << atomNumber << ") " << endl;
		////			continue;
		///		} 
				Atom* atom = new Atom();
                atom->atomType = atomType;
				atom->mol2Type = mol2Type; 
				atom->atomicNumber = mol2par->getOrdinalNumber(atomType);

                atom->getVdw ();

				atom->number = atomNumber;
				atom->name = atomName;
				
				atom->substructureNumber = atomSubstructureNumber;
				atom->substructureName = atomSubstructureName;
				atom->charge = atomCharge;
			    atom->inBackBone = false;

                atom->sphere_already_drawn = false;
                atom->displayStyle = 0;
				
				atom-> GetVector () = vect (x, y, z);
				

				atom->ID = atoms.size();
				
		///////		atom->isProtein = isProtein;
	/////			atom->isWater = isWater;
//////				atom->isLigand = !isProtein && !isWater;
				
                atom->set_color_mw ();
   				atoms.push_back(atom);				
				string atomProperty;
				while (getline(iss, atomProperty, '|')) {
					// skip white spaces
					istringstream iss2(atomProperty);
					iss2 >> atomProperty;
					//cout << bondProperty << "|";
					// check for backbone
					if (atomProperty == triposBACKBONE) {
						atom->inBackBone = true;
						atom->mol2BACKBONE = true;
					}
					if (atomProperty == triposDICT) {
						atom->mol2DICT = true;
					}
					if (atomProperty == triposESSENTIAL) {
						atom->mol2ESSENTIAL = true;
					}
					if (atomProperty == triposDIRECT) {
						atom->mol2DIRECT = true;
					}
        



				}
			}	
		}
		else if (token == triposBond) {
			for (unsigned int i = 0; i < bondCount; i++) {
				getline(file, buffer);
				lineNumber++;
				istringstream iss(buffer);
				
				unsigned int bondNumber;
				unsigned int fromAtom;
				unsigned int toAtom;
				string bondType;

				iss >> bondNumber;
				if (!checkInputStream(filename, "Bond number expected.", lineNumber, iss)) {
					return 0;
				}
				iss >> fromAtom;
				if (!checkInputStream(filename, "Bond atom 1 expected.", lineNumber, iss)) {
					return 0;
				}
				iss >> toAtom;
				if (!checkInputStream(filename, "Bond atom 2 expected.", lineNumber, iss)) {
					return 0;
				}
				iss >> bondType;
				if (!checkInputStream(filename, "Bond type expected.", lineNumber, iss)) {
					return 0;
				}
	
				
				Bond* bond = new Bond();
				bond->mol2Type = mol2par->getmol2BondTypeID(bondType);
                bond->kekule = bond->mol2Type;
				bond->number = bondNumber;
       //         bond->hasCoplanar0 = false;
        //        bond->hasCoplanar1 = false;
                bond->displayStyle = 1;
				
				bool atomFound1 = false;
				bool atomFound2 = false;
				for (unsigned int j = 0; j < atoms.size() && !atomFound1; j++) {
					if (atoms[j]->number == fromAtom) { 
						bond->atomID[0] = j;
						bond->atomPTR[0] = atoms[j];
						atomFound1 = true;
					}
				}
				for (unsigned int j = 0; j < atoms.size() && !atomFound2; j++) {
					if (atoms[j]->number == toAtom) { 
						bond->atomID[1] = j;
						bond->atomPTR[1] = atoms[j];
						atomFound2 = true;
					}
				}
               

				if (!atomFound1 || !atomFound2) {
					//cerr << "parse error in file " << filename << " (line " << lineNumber << "): Could not find atom " << fromAtom << ". Exit." << endl;
					//exit(1);
					//cerr << "REMOVED BOND: " << bondNumber << endl;
					cout << "removed bond " << setw(5) << setfill (' ') << bondNumber << endl;
					delete bond;
                }
				bonds.push_back(bond);
            }
        }
       if (!nextMoleculeFound) res = getline(file, buffer);
    }


    find_bound ();
    find_residues ();
    find_fragments ();
    find_rings ();
    find_kekule ();
    find_limits ();
    find_center ();

    cout<<fragments.size()<<" fragments"<<endl;
    if (nextMoleculeFound) {
        return 2;
    }
    return 1;
}

Bond * Molecule::find_bond (Atom *at1, Atom *at2) {
    for (unsigned int i=0; i<at1->bonds.size (); i++) {
        if (at1->bonds[i]->atomPTR[0]==at1 && at1->bonds[i]->atomPTR[1]==at2) return at1->bonds[i];
        if (at1->bonds[i]->atomPTR[1]==at1 && at1->bonds[i]->atomPTR[0]==at2) return at1->bonds[i];
    }
    return NULL;
} 




void Molecule::remove (Bond* bo) {
    for (unsigned int i=0; i<bonds.size (); i++) {
        if (bonds[i] == bo) {
            bonds.erase (bonds.begin ()+i);
            return;
        }
    }
}

void Molecule::remove (Atom* at) {
    for (unsigned int i=0; i<atoms.size (); i++) {
        if (atoms[i] == at) {
            atoms.erase (atoms.begin ()+i);
            return;
        }
    }
}




void Molecule::del (Bond* bo) {
    for (unsigned int i=0; i<bonds.size (); i++) {
        if (bonds[i] == bo) {
            bonds.erase (bonds.begin ()+i);
            delete bo;
            return;
        }
    }
}

void Molecule::del (Atom* at) {
    for (unsigned int i=0; i<atoms.size (); i++) {
        if (atoms[i] == at) {
            atoms.erase (atoms.begin ()+i);
           delete at;
            return;
        }
    }
}

void Molecule::add (Atom *at) {
    atoms.push_back (at);
}

void Molecule::add (Bond *bo) {
    bonds.push_back (bo);
}





void Molecule::add_atom_bonded_to (Atom *to_add, Atom *partner) {
    Bond *bond = new Bond;
    add_atom_bonded_to (to_add, bond, partner);
}




Atom *Molecule::add_atom_bonded_to (vect coords, int atomnum, Atom *at) {
    Molecule *mol = at->residue->molecule;
    Atom *new_at = new Atom;
    new_at->atomicNumber = atomnum;
    new_at-> GetVector () = coords;  
    new_at->substructureNumber = at->substructureNumber;
    new_at->substructureName = at->substructureName;

    new_at->getVdw ();
    new_at->set_color_mw ();
    new_at->displayStyle = at->displayStyle;
    mol->atoms.push_back (new_at);
    Bond *bond = new Bond;
    bond->atomPTR[0] = at;
    bond->atomPTR[1] = new_at;
    bond->mol2Type = 1;
    bond->kekule = 1;
    if (mol->bonds.size ()) bond->displayStyle = mol->bonds[0]->displayStyle;
    else bond->displayStyle = 1;
    mol->bonds.push_back (bond);
    return new_at;
}

void Molecule::number_atoms () {
    for (unsigned int i=0; i<atoms.size (); i++) atoms[i]->ID = i;
}

void Molecule::number_bonds () {
    for (unsigned int j=0; j<bonds.size (); j++)         {
        bonds[j]->number = j;
        bonds[j]->atomID[0] = bonds[j]->atomPTR[0]->ID;
        bonds[j]->atomID[1] = bonds[j]->atomPTR[1]->ID;
    }
}

void Molecule::find_bound () {
    for (unsigned int i=0; i<atoms.size (); i++) {

        atoms[i]->no2bond = false;
        atoms[i]->bonds.clear ();
        atoms[i]->bound.clear ();
    }
    for (unsigned int j=0; j<bonds.size (); j++) {
        Bond *bond = bonds[j];
        bond->kekule = bond->mol2Type;
        bond->atomPTR[0]->bound.push_back (bond->atomPTR[1]);
        bond->atomPTR[1]->bound.push_back (bond->atomPTR[0]);
        bond->atomPTR[0]->bonds.push_back (bond); 
        bond->atomPTR[1]->bonds.push_back (bond);                   
                 
    }
}



void Molecule::find_center () {
    center = find_mass_center (atoms);
}


void Molecule::find_kekule () {
    bool go_on = true;
    int nn=0;
    while (go_on) {
        go_on = false;
        for (unsigned int r=0; r<rings.size (); r++) { 
            Ring *rin = rings[r];
            if (rin->aromatic) {
                if (rin->kekule_type>nn) go_on=true;
                if (rin->kekule_type==nn) rin->find_kekule ();
            }
        }
        for (unsigned int r=0; r<rings.size (); r++) { 
            if (!rings[r]->aromatic) {
                for (unsigned int b=0; b<rings[r]->bonds.size (); b++) {             
                    if (rings[r]->bonds[b]->mol2Type==5) {
                        rings[r]->find_kekule_pseudoaromatic ();
                        break;
                    }
                }
            }
        }
        nn++;       
    }
    for (unsigned int r=0; r<rings.size (); r++) {
        Ring *rin = rings[r];
        if (rin->aromatic) rin->test_kekule_aromaticity ();
    }
}

void Molecule::find_limits () {
    if (atoms.size ()) {
        float xm, ym, zm, xM, yM, zM;
        xm = xM = atoms[0]-> GetVector ().x();
        ym = yM = atoms[0]-> GetVector ().y();
        zm = zM = atoms[0]-> GetVector ().z();
        for (unsigned a=0; a<atoms.size(); a++){
            if (atoms[a]-> GetVector ().x()<xm) xm= atoms[a]-> GetVector ().x();
            if (atoms[a]-> GetVector ().x()>xM) xM= atoms[a]-> GetVector ().x();
            if (atoms[a]-> GetVector ().y()<ym) ym= atoms[a]-> GetVector ().y();
            if (atoms[a]-> GetVector ().y()>yM) yM= atoms[a]-> GetVector ().y();
            if (atoms[a]-> GetVector ().z()<zm) zm= atoms[a]-> GetVector ().z();
            if (atoms[a]-> GetVector ().z()>zM) zM= atoms[a]-> GetVector ().z();
        }
        min_corner = vect (xm, ym, zm);
        max_corner = vect (xM, yM, zM);
    }
}

void Molecule::find_fragments () {
    fragments.clear ();
    for (unsigned int a=0; a<atoms.size(); a++) {
        atoms[a]->visited = false;
    }
    for (unsigned int a=0; a<atoms.size(); a++) {
        if (atoms[a]->visited) continue;
        Fragment *fragment = new Fragment ();
        queue <Atom*> queue;
        queue.push (atoms[a]);
        atoms[a]->visited = true;
        while (queue.size () ) {
            Atom *at = queue.front();
            queue.pop ();
            fragment->atoms.push_back (at);
            for (unsigned int i=0; i<at->bonds.size (); i++) {
                if (!at->bonds[i]->is_rotatable () && !at->bound[i]->visited) {
                    queue.push (at->bound[i]);       
                    at->bound[i]->visited = true;
                }
            }
        }
        fragments.push_back (fragment);

    }

}

void Molecule::find_rings () {
    for (unsigned int i=0; i<rings.size (); i++) {
        delete rings[i];
    }
    rings.clear ();
    for (unsigned int a=0; a<atoms.size (); a++) {
        atoms[a]->in_ring.clear ();
    }
    for (unsigned int b=0; b<bonds.size (); b++) {
        bonds[b]->in_ring.clear ();
    }


    for (unsigned int b=0; b<bonds.size(); b++){
    //    cout << "b ";
        if ((bonds[b]->atomPTR[0]->bonds.size ()==1) || (bonds[b]->atomPTR[1]->bonds.size ()==1) || bonds[b]->in_ring.size ()) continue;
        queue <Bond*> queue;
   //     cout <<"ok"<<endl;
        for (unsigned int b2=0; b2<bonds.size(); b2++) {
            bonds[b2]->fr_visited = false;
            bonds[b2]->fr_parent = NULL;
        
        }

//    int nr = 0;
  //      cout <<"direction "<<b0dir<<endl;        

        bonds[b]->fr_dir = 0;
        bonds[b]->fr_visited = true;
        bonds[b]->fr_parent = NULL;
        queue.push (bonds[b]);
        bool ring_found = false;
     //   cout <<"starting from bond "<<bonds[b]->number<<endl;
        while (queue.size () && !ring_found) {
            Bond *bo = queue.front();
            queue.pop ();
       //     cout << "popping "<<bo->number<<endl;
            for (unsigned int bb = 0; bb<bo->atomPTR[bo->fr_dir]->bonds.size (); bb++) {
                Bond *nextbo = bo->atomPTR[bo->fr_dir]->bonds[bb];
                if (nextbo!=bo && nextbo->atomPTR[0]->bound.size()!=1 && nextbo->atomPTR[1]->bound.size()!=1 ) {
    
                    if (nextbo->fr_visited) {
         //               cout <<nextbo->number<<" already visited"<<endl;
                        if (nextbo==bonds[b]) 
                                ring_found=close_ring (bo);
                    }
                    else{
                        int dir = 0;
                        for (unsigned int d=0; d<nextbo->atomPTR[0]->bonds.size(); d++){
                            if (nextbo->atomPTR[0]->bonds[d]==bo) {
                            dir = 1;
                            break;
                            }
                        }
                        nextbo->fr_dir = dir;
                        nextbo->fr_parent = bo;

                        nextbo->fr_visited = true;
           //             cout << "pushing "<<nextbo->number<<" "<<nextbo->atomPTR[0]->ID<<" "<< nextbo->atomPTR[1]->ID<<endl;
                        queue.push (nextbo);
                    }
                }
            }  

        }

    }
 //   cout<<"rings "<<nr<<endl;
}


bool Molecule::close_ring (Bond *bo) {
    vector <Bond*> out_bonds;
    Bond *curr_bo;
 //   cout <<"close rings"<<endl;
   curr_bo = bo;
   out_bonds.push_back (bo); 
    while (curr_bo->fr_parent) {
        out_bonds.push_back (curr_bo->fr_parent);
        curr_bo = curr_bo->fr_parent;
    }

 //   for (unsigned int q=0; q<out_bonds.size(); q++) cout<<out_bonds[q]<<" "; cout <<endl;


    if (true) {


        bool already_found = false;
        if (!already_found) {
            Ring *ring = new Ring;
  //          bool aromatic = true;
            for  (unsigned int i=0; i<out_bonds.size(); i++) {
  //              if (out_bonds[i]->mol2Type!=5) aromatic = false;
                out_bonds[i]->in_ring.push_back (ring);
                bool not_put0 = true;
                bool not_put1 = true;
                for(unsigned int at=0; at<ring->atoms.size();at++){
                    if (ring->atoms[at]==out_bonds[i]->atomPTR[0]) {
                        not_put0 = false;
                    }
                    if (ring->atoms[at]==out_bonds[i]->atomPTR[1]) {
                        not_put1 = false;
                    }
    
                }
                if (not_put0) {
                    ring->atoms.push_back (out_bonds[i]->atomPTR[0]);
                    out_bonds[i]->atomPTR[0]->in_ring.push_back (ring);
                }
                if (not_put1) {
                    ring->atoms.push_back (out_bonds[i]->atomPTR[1]);
                    out_bonds[i]->atomPTR[1]->in_ring.push_back (ring);
                }
                ring->bonds.push_back (out_bonds[i]);       
            }
            ring->center = find_mass_center (ring->atoms);

            ring->test_aromaticity ();
            rings.push_back (ring);
            return 1;
        }
    }
    else return 0;
}

void Molecule::findBonds (Atom *atom, vector <Bond*>& out){
    unsigned int ID=atom->ID;
    for (unsigned int i=0; i<bonds.size ();i++) {
        for (unsigned int j=0; j<2;j++) {
            if (bonds[i]->atomID[j]==ID) {
   //             cout <<"bond found"<<endl;
                out.push_back (bonds[i]);
            }
        }
    }
  //  cout <<"atom "<<atom->atomType<<" bound to "<<out.size()<<" other atoms"<<endl;


}

Atom* Molecule::findOtherAtom (Bond *bond, Atom *atom){
    if (bond->atomID[0]!=atom->ID) return bond->atomPTR[0];
    else if   (bond->atomID[1]!=atom->ID) return bond->atomPTR[1];
    else {
        cerr << "Atom "<< atom->ID << "not in bond "<< bond->number<<endl;
        return bond->atomPTR[0];
    }
}


void Molecule::find_residues (){
    vector<string> vec;
    vector<int> vecn;
    int lastn = -1;
    bool resf = false;
    for (unsigned int i=0; i<residues.size (); i++) {
        delete residues[i];
    }
    residues.clear ();

    for (unsigned int i=0; i<atoms.size(); i++) {
        resf = false;
        int ssnum =atoms[i]->substructureNumber;

        if (ssnum != lastn) {
            for (unsigned int j=0; j<vecn.size(); j++) {
                if (ssnum == vecn[j]) {
                    resf = true;
                    break;
                }
            }
            if (!resf){
                lastn = ssnum;
                vecn.push_back (ssnum);
                assert (i>=0 && i<atoms.size ());
                string ssname =atoms[i]->substructureName;
                stringstream  ssn;
                ssn << ssnum ;
                vec.push_back (ssn.str ()+"  "+ssname+ssn.str());
                Residue *res = new Residue ();
                res->number = ssnum;
                res->name = ssname;

                res->molecule = this;
                residues.push_back (res);
            }
        }
    }
    for (unsigned int i = 0; i<atoms.size (); i++) {
        for (unsigned int j = 0; j<residues.size (); j++) {
            if (atoms[i]->substructureNumber==residues[j]->number) {
                residues[j]->atoms.push_back (atoms[i]);
                atoms[i]->residue=residues[j];
                break;
            }
        }
    }

    for (unsigned int j = 0; j<residues.size (); j++) {
        for (unsigned int i = 0; i<residues[j]->atoms.size (); i++) {
            if (residues[j]->atoms[i]->inBackBone) {
                if (residues[j]->atoms[i]->atomType == "N.am" ) {
                    residues[j]->NB = residues[j]->atoms[i];
                    residues[j]->hasNB = true;
                }
                else if (residues[j]->atoms[i]->atomType == "C.2" ) {residues[j]->CB = residues[j]->atoms[i];
                    residues[j]->hasCB = true;
                }
                else if (residues[j]->atoms[i]->atomType == "C.3" ) {residues[j]->CA = residues[j]->atoms[i];
                    residues[j]->hasCA = true;
                }
            }
        }

        if (residues[j]->hasNB && residues[j]->hasCB && residues[j]->hasCA){
            vect point;

            for (unsigned int n=0; n<residues[j]->NB->bound.size ();n++) {
                if (residues[j]->NB->bound[n]->atomType == "C.2") {
                    vect pre;
                    pre = residues[j]->NB->bound[n]-> GetVector ();
                    point = mean (pre, residues[j]->NB-> GetVector ());
                    residues[j]->backbone_list.push_back (point);


                }
            }

            point = mean (residues[j]->NB-> GetVector (), residues[j]->CA-> GetVector ());
            residues[j]->backbone_list.push_back (point);


            point = mean (residues[j]->CB-> GetVector (), residues[j]->CA-> GetVector ());
            residues[j]->backbone_list.push_back (point);

            for (unsigned int n=0; n<residues[j]->CB->bound.size ();n++) {
                if (residues[j]->CB->bound[n]->atomType == "N.am") {
                    vect fol;
                    fol = residues[j]->CB->bound[n]-> GetVector ();
                    point = mean (fol, residues[j]->CB-> GetVector ());
                    residues[j]->backbone_list.push_back (point);
            
                }
            }

        }
    if (residues[j]->backbone_list.size()) {
        for (unsigned int n =0; n<8; n++) {
            residues[j]->refine_backbone ();
        }

    }

    }
    res_strings = vec;

}


Residue::Residue (){
    backbone_style = 0;
    backbone_visible = true;
    hasCA = false;
    hasNB = false;
    hasCB = false;
    backbone_color = color (0.8, 0.8, 0.f, 1.f);

    backbone_list.clear ();
};

void Residue::refine_backbone (){
    vector <vect> nb;
    nb.push_back (backbone_list[0]);
        for (unsigned int i=0; i<backbone_list.size ()-1; i++){
        vect point;
        point = mean (backbone_list[i], backbone_list [i+1]);
        nb.push_back (point);
        } 
    nb.push_back (backbone_list[backbone_list.size ()-1]);
    backbone_list = nb;
}


Ring::Ring () {
    alpha_atom = NULL;
    IM_CAT = false;
    N5ANION = false;
}

void Ring::complete_kekule () {
    bool go_on = true;
    while (go_on) {
        go_on = false;
        for (unsigned int i=0; i<bonds.size (); i++) {
            if (bonds[i]->kekule==5) {
                if (bonds[i]->conjugated_to (2, this) && bonds[i]->conjugated_to (1, this)) {
                    bonds[i]->kekule =1;
                    cout<<"bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<" set to 1"<<endl;
                }   
                else if (bonds[i]->conjugated_to (2, this)) {
                    bonds[i]->kekule =1;
                    cout<<"bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<" set to 1"<<endl;
                }   
                else if (bonds[i]->conjugated_to (1, this) && !bonds[i]->atomPTR[0]->no2bond && !bonds[i]->atomPTR[1]->no2bond) {
                    bonds[i]->kekule =2;
                    cout<<"bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<" set to 2"<<endl;
                    bonds[i]->atomPTR[0]->no2bond = true;
                    bonds[i]->atomPTR[1]->no2bond = true;
                }
             
                else if (!bonds[i]->conjugated_to (2, this) && !bonds[i]->conjugated_to (1, this)) {
                    go_on = true;
                }
            }
        }
    }
}

bool Ring::all_aromatic () {
    bool all_arom = true;
    for (unsigned int i=0; i<bonds.size (); i++) {
        if (bonds[i]->kekule!=5) all_arom=false;
    }
    return all_arom;
}

void Ring::find_kekule_pseudoaromatic () {
    cout <<"pseudo"<<endl;
     for (unsigned int i=0; i<bonds.size (); i++) {   
        if ((bonds[i]->atomPTR[0]->no2bond || bonds[i]->atomPTR[1]->no2bond) && bonds[i]->kekule==5) {
            cout<<"bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<"set as 1"<<endl;
            bonds[i]->kekule=1;
        }    
    }
    bool go_on = true;
    while (go_on) {
        go_on = false;    
        for (unsigned int i=0; i<bonds.size (); i++) {
            if (((bonds[i]->atomPTR[0]->atomicNumber==7 && bonds[i]->atomPTR[0]->bound.size ()==2) || (bonds[i]->atomPTR[0]->atomicNumber==6 && bonds[i]->atomPTR[0]->bound.size ()==3))  &&
((bonds[i]->atomPTR[1]->atomicNumber==7 && bonds[i]->atomPTR[1]->bound.size ()==2) || (bonds[i]->atomPTR[1]->atomicNumber==6 && bonds[i]->atomPTR[1]->bound.size ()==3))  && 
bonds[i]->conjugated_to (1, this) &&  
!bonds[i]->conjugated_to (2, this) && 
bonds[i]->kekule==5 &&
!bonds[i]->atomPTR[0]->no2bond && 
!bonds[i]->atomPTR[1]->no2bond) {
                
                cout<<"bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<" set as 2"<<endl;
                bonds[i]->atomPTR[0]->no2bond = true;
                bonds[i]->atomPTR[1]->no2bond = true;


                bonds[i]->kekule=2;
                go_on = true;
                
            }    
            if (bonds[i]->conjugated_to (2, this) && bonds[i]->kekule==5) {
                cout<<"bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<" set as 1"<<endl;
                bonds[i]->kekule=1;
                go_on = true;
            }
        }
    }

    go_on = true;
    while (go_on) {
        go_on = false;    
        for (unsigned int i=0; i<bonds.size (); i++) {
            if (((bonds[i]->atomPTR[0]->atomicNumber==7) || (bonds[i]->atomPTR[0]->atomicNumber==6 && bonds[i]->atomPTR[0]->bound.size ()==3))  &&
((bonds[i]->atomPTR[1]->atomicNumber==7) || (bonds[i]->atomPTR[1]->atomicNumber==6 && bonds[i]->atomPTR[1]->bound.size ()==3))  && 
bonds[i]->conjugated_to (1, this) &&  
!bonds[i]->conjugated_to (2, this) && 
bonds[i]->kekule==5 &&
!bonds[i]->atomPTR[0]->no2bond && 
!bonds[i]->atomPTR[1]->no2bond) {
                
                cout<<"bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<" set as 2"<<endl;
                bonds[i]->atomPTR[0]->no2bond = true;
                bonds[i]->atomPTR[1]->no2bond = true;


                bonds[i]->kekule=2;
                go_on = true;
                
            }    
            if (bonds[i]->conjugated_to (2, this) && bonds[i]->kekule==5) {
                cout<<"bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<" set as 1"<<endl;
                bonds[i]->kekule=1;
                go_on = true;
            }
        }



        for (unsigned int i=0; i<bonds.size (); i++) {
            if (bonds[i]->kekule == 5) {
                bonds[i]->kekule=1;
                cout<<"bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<" set as 1"<<endl;
            }
        }
    }
}

void Ring::find_kekule () {

     for (unsigned int i=0; i<bonds.size (); i++) {   
        if ((bonds[i]->atomPTR[0]->no2bond || bonds[i]->atomPTR[1]->no2bond) && bonds[i]->kekule==5) {
            cout<<"find_kekule: bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<"set as 1"<<endl;
            bonds[i]->kekule=1;
        }    
    }

    if (atoms.size ()==5) {
        cout<<"find kekule 5 membered ring"<<endl;
        if (alpha_atom) {
            for (unsigned int i=0; i<bonds.size (); i++) {
                if ((bonds[i]->atomPTR[0]==alpha_atom || bonds[i]->atomPTR[1]==alpha_atom) && bonds[i]->kekule==5) {
                    bonds[i]->kekule=1;

                } 
            }
            complete_kekule ();
        }
        else {
            if (all_aromatic ()) bonds[0]->kekule =1;
            complete_kekule ();
        }

    }
    else if (atoms.size ()==6) {
        cout<<"find kekule 6 membered ring"<<endl;


        if (all_aromatic ()) {
            cout <<"all aromatic"<<endl;
            for (unsigned int i=0; i<bonds.size (); i++) {
                if (bonds[i]->in_ring.size()>1 && (!bonds[i]->atomPTR[0]->no2bond && !bonds[i]->atomPTR[1]->no2bond)) {
                    bonds[i]->kekule = 2;
                    cout<<"ring fusion bond "<<bonds[i]->atomPTR[0]->ID<<" "<<bonds[i]->atomPTR[1]->ID<<"set as 2"<<endl;
                    bonds[i]->atomPTR[0]->no2bond = true;
                    bonds[i]->atomPTR[1]->no2bond = true;

                }
            }
        }

        if (all_aromatic ()) {
            bonds[0]->kekule =1;
            cout<<"random bond "<<bonds[0]->atomPTR[0]->ID<<" "<<bonds[0]->atomPTR[1]->ID<<"set as 1"<<endl;
        }
        complete_kekule ();
    
    }

}
unsigned int Ring::count (int an) {
    unsigned int out = 0;
    for (unsigned int i=0; i<atoms.size (); i++) {
        if (atoms[i]->atomicNumber == an) out+=1;
    }
    return out;
}

void Ring::test_kekule_aromaticity () {
    for (unsigned int i=0; i<atoms.size (); i++) {
        bool not_ar = true;
        for (unsigned int j=0; j<atoms[i]->bonds.size (); j++) {
            if (atoms[i]->bonds[j]->kekule==2 && atoms[i]->bonds[j]->mol2Type==5) not_ar=false;
        }
        if (atoms[i]->atomicNumber==8 || atoms[i]->atomicNumber==16) not_ar=false;
        if (atoms[i]->atomicNumber==7 && atoms[i]->bound.size ()==3 && atoms.size ()!=6) not_ar=false;
        if (not_ar) {
            aromatic = false;
            break;
        } 
    }
}

void Ring::test_aromaticity () {

    bool ar = false;
    bool quit = false;
    unsigned int size = atoms.size ();
    for (unsigned int i=0; i<atoms.size (); i++) {
        if (atoms[i]->bound.size ()>3) {
            aromatic =false;
            quit = true;
            atoms[i]->no2bond = true;
            cout <<atoms[i]->ID<<" set as no2bond"<<endl;
  //          return;
        }
        for (unsigned int j=0; j<atoms[i]->bound.size (); j++) {
            if (!atoms[i]->bound[j]->is_in_ring (this) && atoms[i]->bonds[j]->mol2Type==2) {
                aromatic=false;
                quit = true;
                atoms[i]->no2bond = true;
            cout <<atoms[i]->ID<<" set as no2bond"<<endl;
 //               return;
            }
        }
    }
    if (quit) return;
    int pyrrole_n = 0;
    int pyridine_n = 0;
    int os = 0;
    Atom *lasto = NULL;
    Atom *lastn = NULL;
    for (unsigned int i=0; i<atoms.size (); i++) {
        if (atoms[i]->atomicNumber==8 || atoms[i]->atomicNumber==16) {
            os ++;
            lasto = atoms[i];
        }
        if (atoms[i]->atomicNumber==7 && atoms[i]->bound.size ()==3) {
            pyrrole_n++;
            lastn = atoms[i];
        }
        if (atoms[i]->atomicNumber==7 && atoms[i]->bound.size ()==2) {
            pyridine_n++;
        }
    }


    if (size<5) ar=false;

    else if (size ==5) {

        if (os) { //furane thiophene
            if (os>1) ar=false;
            else {
                ar=true;
                kekule_type = 0;
                alpha_atom = lasto;
                cout <<lasto->ID<<" set as no2bond"<<endl;
                lasto->no2bond = true;
                if (pyrrole_n) IM_CAT = true;
            }
        }
        else if (pyrrole_n==1) { //pyrrole , pyrazole, imidazole, triazole && tetrazole
            ar=true;
            alpha_atom = lastn;
            lastn->no2bond = true;  
            kekule_type = 0;        

        }
        else if (pyrrole_n==2 && pyridine_n) {
            cout <<"false"<<endl;
            ar = false;
        }
        else if (pyrrole_n==2) {//imidaziolium cation 

            int nCC =0, NH=0;
            Bond *lastCC;
            for (unsigned int i=0; i<bonds.size (); i++) {
                if (bonds[i]->atomPTR[0]->atomicNumber==6 && bonds[i]->atomPTR[1]->atomicNumber==6) {
                    nCC++; 
                    lastCC = bonds[i];
                }
            }
            for (unsigned int i=0; i<atoms.size (); i++) {
                if (atoms[i]->atomicNumber==7 && atoms[i]->bonded_to (1, 1)) {
                    NH++;
                }
            }
            if (nCC==1 && NH==2) {//imidazolium
                cout <<"imidazolium"<<endl;
                lastCC->kekule=2;
                lastCC->atomPTR[0]->no2bond = true;
                lastCC->atomPTR[0]->no2bond = true;
                kekule_type =4;
                ar = true;
                alpha_atom = NULL;
                IM_CAT = true;
            }
            else {
                ar = false;
            }
        
        }
        else if (pyridine_n ==4 && !pyrrole_n) { //tetrazole anion
            cout<<"tetrazole anion"<<endl;
            ar = true;
            kekule_type = 3;
            N5ANION = true;
            alpha_atom = NULL;
        }

        else ar=false;
    }
    else if (size ==6) {
        kekule_type = 3;
        if (os) ar=false;
        else ar=true;
    }
    else ar=false;

    aromatic = ar;

}


Fragment::Fragment () {
}

//Fragment::~Fragment () {}




vect find_mass_center (vector<Atom*>& invec){

    unsigned int i;
    vect midCoo;

	
	unsigned int numMid = 0;
	for (i = 0; i < invec.size(); i++) {
	    vect a = sum (midCoo, invec[i]-> GetVector ());
        midCoo = a;
		numMid++;
		
	}
	
    midCoo.multiply (1.0f/(float) numMid);
	
/*
	float radius = 0.0;
	for (i = 0; i < invec.size(); i++) {

			float r2 = (invec[i]-> GetVector ()[0] - midCoo[0]) *(invec[i]-> GetVector ()[0] - midCoo[0]) +
						(invec[i]-> GetVector ()[1] - midCoo[1]) *(invec[i]-> GetVector ()[1] - midCoo[1]) +
						(invec[i]-> GetVector ()[2] - midCoo[2]) *(invec[i]-> GetVector ()[2] - midCoo[2]);
			if (r2 > radius) {
				radius = r2;
			}
	}
	


    return midCoo;

}

*/

