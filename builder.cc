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





#include "builder.h"


Builder::Builder (DDWin *ddw) {
    ddwin = ddw;
    start_magic_pencil_atom = NULL;
    last_magic_pencil_atom = NULL;
    mcounter = 1;

}


void Builder::save_Hs (Atom *at, vector <Atom *> &prev_ats, vector <ZNBond *> &prev_bonds) {
    prev_ats.clear ();
    prev_bonds.clear ();
    FOR_NBORS_OF_ATOM (n, at) {
        if (n -> IsHydrogen ()) {
            ZNBond *bo = at -> GetBond (&*n);
            if (bo) {
                prev_ats.push_back (&*n);
                prev_bonds.push_back (bo);
            }
        }
    }
}


void Builder::add_Hs (Atom *at, vector <Atom *> &atoms, vector <ZNBond *> &bonds) {
    ZNMolecule *mol = (ZNMolecule *) at -> GetParent ();

    if (atoms.size () && bonds.size ()) {
        assert (atoms.size() == bonds.size ());
        for (unsigned int i = 0; i< atoms.size (); i++) {
            mol -> add_atom_bonded_to (atoms[i], bonds[i], at);
        }
    }
}


void Builder::delete_Hs (ZNMolecule *mol) {
//    cerr << "removing_Hs" << endl;
    vector <Atom *>to_del;
    FOR_ATOMS_OF_MOL (n, mol) {
        if (n) {
            if (n -> IsHydrogen ()) {
                to_del.push_back (&*n);
            }
        }
    }  
    for (unsigned int i = 0; i < to_del.size (); i++) { 
        mol -> RemoveAtom (to_del[i]);
    }
 //   cerr << "removed_Hs" << endl;
}





void Builder::delete_Hs (Atom *at) {
  //  cerr << "removing_Hs" << endl;
    ZNMolecule *mol = (ZNMolecule *) at -> GetParent ();
    vector <Atom *>to_del;
    FOR_NBORS_OF_ATOM (n, at) {
        if (n) {
            if (n -> IsHydrogen ()) {
                to_del.push_back (&*n);
            }
        }
    }  
    for (unsigned int i = 0; i < to_del.size (); i++) { 
        mol -> RemoveAtom (to_del[i]);
    }
 //   cerr << "removed_Hs" << endl;
}



void Builder::add_atom (int atomnum) {
	vector <Atom *> to_select;
    bool continue_mol = find_target ();
    if (!continue_mol) {
        vect coord;
        Atom * new_at = new_atom (coord, atomnum);
		FOR_NBORS_OF_ATOM (n, new_at) {
			if (n->IsHydrogen ()) {
				to_select.push_back (&*n);
				break;
			}
		}
    }
    else {
		FOR_ATOMS_OF_MOL (at, ddwin ->target_molecule) {
			
			if (!at->IsHydrogen ()) {
				mutate_atom_to (&*at, atomnum);
				FOR_NBORS_OF_ATOM (n, &*at) {
					if (n->IsHydrogen ()) {
						to_select.push_back (&*n);
						break;
					}
				}
			}
			
			else { //adding to an H
				OBBondIterator i;
				Atom *parent = at -> BeginNbrAtom (i);
				if (parent) {
				//				vect v = parent -> GetNewBondVector ()
				Atom *new_atom =add_atom_bonded_to (get_coordinates (&*at), atomnum, parent);
				FOR_NBORS_OF_ATOM (n, new_atom) {
					if (n->IsHydrogen ()) {
						to_select.push_back (&*n);
						break;
					}
				}
				}
			}
		}
	}
	ddwin ->deselect ();
	ddwin ->select (to_select);
	ddwin ->set_current_target (-1);
	ddwin->gl->draw_molecule (ddwin->target_molecule);

}

////////////////////////////////////////////////////////////////////////////////
void Builder::add_fragment (string str) {
	OBConversion conv;
	ZNMolecule *mol = new ZNMolecule ();
	conv.SetInFormat ("SMI");
	conv.ReadString (mol, str);
    bool continue_mol = find_target ();
    if (!continue_mol) {
		add_mol (str);
    }
    else {
		if (ddwin ->target_molecule ->NumAtoms () ==1) {
			Atom *to_add = ddwin ->target_molecule ->GetAtom (1);
			ddwin ->deselect ();
			ZNMolecule *old_mol = (ZNMolecule *) to_add ->GetParent ();
			old_mol ->DeleteHydrogens (to_add);
			old_mol ->ZNinit ();
			int oldat_idx = to_add ->GetIdx ();

			OBConversion conv;
			conv.SetInFormat ("SMI");
			conv.ReadString (mol, str);

	//		mol -> AddHydrogens ();
			mend_coordinates (mol);

			int newat_idx = old_mol ->NumAtoms () +1;
			(*old_mol) += (*mol);
			OBBuilder obbuild;
			obbuild.Connect (*old_mol, oldat_idx, newat_idx);
			mend_coordinates (old_mol);
			
			old_mol ->AddHydrogens ();
			old_mol ->ZNinit ();


		}
	/*
		FOR_ATOMS_OF_MOL (at, ddwin ->target_molecule) {
			if (!at->IsHydrogen ()) {
				mutate_atom_to (&*at, atomnum);
				FOR_NBORS_OF_ATOM (n, &*at) {
					if (n->IsHydrogen ()) {
						to_select.push_back (&*n);
						break;
					}
				}
			}
			
			else { //adding to an H
				OBBondIterator i;
				Atom *parent = at -> BeginNbrAtom (i);
				if (parent) {
				//				vect v = parent -> GetNewBondVector ()
				Atom *new_atom =add_atom_bonded_to (get_coordinates (&*at), atomnum, parent);
				FOR_NBORS_OF_ATOM (n, new_atom) {
					if (n->IsHydrogen ()) {
						to_select.push_back (&*n);
						break;
					}
				}
				}
			}
		}
		*/
	}
//	ddwin ->deselect ();
//	ddwin ->select (to_select);
//	ddwin ->set_current_target (-1);
//	ddwin->gl->draw_molecule (ddwin->target_molecule);
}


void Builder::add_mol (string str) {
	vector <Atom *> to_select;
	bool continue_mol = find_target ();
	if (!continue_mol) {
		OBConversion conv;
		ZNMolecule *mol = new ZNMolecule ();
		conv.SetInFormat ("SMI");
		conv.ReadString (mol, str);
		mol -> AddHydrogens ();
		mend_coordinates (mol);
		CreateZNMoleculeCommand *command = new CreateZNMoleculeCommand (mol, ddwin);
		ddwin -> execute (command);
	}
	else {
//cerr<< "add_mol: 3"<<endl;
/*
		OBConversion conv;
		ZNMolecule *mol = new ZNMolecule ();
		conv.SetInFormat ("SMI");
		conv.ReadString (mol, str);



		mutate_atom_to (&*at, atomnum);
		FOR_NBORS_OF_ATOM (n, &*at) {
			if (n->IsHydrogen ()) {
				to_select.push_back (&*n);
				break;
			}
		}
*/

/*
		FOR_ATOMS_OF_MOL (at, ddwin ->target_molecule) {
			
			if (!at->IsHydrogen ()) {
				mutate_atom_to (&*at, atomnum);
				FOR_NBORS_OF_ATOM (n, &*at) {
					if (n->IsHydrogen ()) {
						to_select.push_back (&*n);
						break;
					}
				}
			}
			
			else { //adding to an H
				OBBondIterator i;
				Atom *parent = at -> BeginNbrAtom (i);
				if (parent) {
				//				vect v = parent -> GetNewBondVector ()
				Atom *new_atom =add_atom_bonded_to (get_coordinates (&*at), atomnum, parent);
				FOR_NBORS_OF_ATOM (n, new_atom) {
					if (n->IsHydrogen ()) {
						to_select.push_back (&*n);
						break;
					}
				}
				}
			}
		}
*/

	}
	ddwin ->deselect ();
/*
	ddwin ->select (to_select);
	ddwin ->set_current_target (-1);
	ddwin->gl->draw_molecule (ddwin->target_molecule);
*/
}
///////////////////////////////////////////////////////////////////////////////

void Builder::mutate_atom_to (Atom *at, unsigned int atomnum) {
    if (at -> GetAtomicNum () != atomnum) {
        MutateAtomCommand *command = new MutateAtomCommand (at, atomnum, ddwin);
        ddwin -> execute (command);
    }

}



Atom *Builder::add_atom_bonded_to (vect coord, int atmnum, Atom *at) {
 //   cerr << "add" << endl;
    Atom *new_at = new Atom ();
    new_at -> SetAtomicNum (atmnum);
	ZNMolecule *mol = (ZNMolecule *) at -> GetParent ();
	mol->ZNinit_atom (new_at);
    set_coordinates (new_at, coord);  
    ZNBond *bo = new ZNBond ();
    AddAtomCommand *command = new AddAtomCommand (new_at, bo, at, ddwin);
    ddwin -> execute (command);
    return command -> atom;
}








Atom *Builder::new_atom (vect coord, int atomnum) {
    Atom *at = new Atom ();
    at->SetAtomicNum (atomnum);


    ZNMolecule *mol = new ZNMolecule ();
   // mol -> ZNinit_atom (at);


///    mol -> set_color_mw (at);



    stringstream ssname;
    ssname << "New Molecule " << mcounter;
    string s = ssname.str ();
    mol -> SetTitle (s);
    mcounter++;
    mol -> ZNAddAtom (at);
	set_coordinates (at, coord);
    mol -> ZNSetConformers (); 
	finalise_molecule(mol);
    mol -> ZNAddHydrogens (at);

    CreateZNMoleculeCommand *command = new CreateZNMoleculeCommand (mol, ddwin);
    ddwin -> execute (command);
    return at;

}

bool Builder::find_target () {
    if (ddwin->target_molecule->selection) return true;
    else {
        return false;
    }
}





ZNBond *Builder::new_bond (Atom *at1, Atom *at2, int order) {
    if (at1 -> GetParent () == at2 -> GetParent ()) {
        ZNMolecule *mol = (ZNMolecule *) at1 -> GetParent ();
        ZNBond *bo = mol -> GetBond (at1, at2);
        if (!bo) {
            bo = new ZNBond;
            mol -> ZNinit_bond (bo);
            bo -> SetBegin (at1);
            bo -> SetEnd (at2);
            AddBondCommand *command = new AddBondCommand (bo, ddwin);
            ddwin -> execute (command);
        } 
        bo -> SetBondOrder (order);
        ddwin->gl->draw_molecule (mol);

        return bo;
    }
	else {
			ZNMolecule *mol1 = (ZNMolecule *) at1 -> GetParent ();
			ZNMolecule *mol2 = (ZNMolecule *) at2 -> GetParent ();
			ZNMolecule *mol3 = sum (mol1, mol2);
			Atom *new_at1 = mol3 ->GetAtom (at1 ->GetIdx ());
		//	cerr << new_at1 -> GetVector () << at1 ->GetVector ()<< endl;
			Atom *new_at2 = mol3 ->GetAtom (at2 ->GetIdx () + mol1 -> NumAtoms ());
		//	cerr << new_at2 -> GetVector () << at2 ->GetVector ()<< endl;			
			ZNBond *bo = new ZNBond;
			bo -> SetBegin (new_at1);
            bo -> SetEnd (new_at2);
            mol3 -> ZNinit_bond (bo);
			bo -> SetBondOrder (1);
			
			AddBondCommand *bond_command = new AddBondCommand (bo, ddwin);
			bond_command -> redo ();
			delete bond_command;
			
			CreateZNMoleculeCommand *command = new CreateZNMoleculeCommand (mol3, ddwin);
			DeleteZNMoleculeCommand *command1 = new DeleteZNMoleculeCommand (mol1, ddwin);
			DeleteZNMoleculeCommand *command2 = new DeleteZNMoleculeCommand (mol2, ddwin);
			ddwin ->data -> undo_stack -> beginMacro ("Merge ZNMolecules");
			ddwin ->execute (command2);
			ddwin ->execute (command1);
			ddwin ->execute (command);

			ddwin ->data -> undo_stack -> endMacro ();			
			
	}
    return NULL;
}

void Builder::set_bond (ZNBond *bo, int order) {
    ModifyBondCommand *command = new ModifyBondCommand (bo, order, ddwin);    
    ddwin -> execute (command);
}

void Builder::delete_bond (ZNBond *bo) {
    DeleteBondCommand *command = new DeleteBondCommand (bo, ddwin);
    ddwin -> execute (command);
}


void Builder::set_aromatic (Ring *ring) {
  /*  ddwin -> data -> undo_stack -> beginMacro ("Set Aromatic");
    RedefineZNMoleculeCommand *redefine1 = new RedefineZNMoleculeCommand (ring->atoms[0]->residue->molecule, this, 1);
    ddwin -> execute (redefine1);
    for (unsigned int i=0; i<ring->bonds.size (); i++) {
        ModifyBondCommand *command = new ModifyBondCommand (ring -> bonds[i], 5, ddwin, false);    
        ddwin -> execute (command);
    }
    RedefineZNMoleculeCommand *redefine = new RedefineZNMoleculeCommand (ring->atoms[0]->residue->molecule, this, 0);
    ddwin -> execute (redefine);
    ddwin -> data ->undo_stack -> endMacro ();
*/
}




void Builder::set_non_aromatic (Ring *ring) {
  /*  ddwin -> data -> undo_stack -> beginMacro ("Set Non Aromatic");
    RedefineZNMoleculeCommand *redefine1 = new RedefineZNMoleculeCommand (ring->atoms[0]->residue->molecule, this, 1);
    ddwin -> execute (redefine1);
    for (unsigned int i=0; i<ring->bonds.size (); i++) {
        ModifyBondCommand *command = new ModifyBondCommand (ring -> bonds[i], 1, ddwin, false);    
        ddwin -> execute (command);
    }
    RedefineZNMoleculeCommand *redefine = new RedefineZNMoleculeCommand (ring->atoms[0]->residue->molecule, this, 0);
    ddwin -> execute (redefine);
    ddwin -> data ->undo_stack -> endMacro ();
*/
}






///////////////////////////////////////////////////////////////








void Builder::add_H (Atom *at, vector <Atom *> &atoms, vector <ZNBond *> &bonds) {

    ZNMolecule *mol = (ZNMolecule *) at -> GetParent ();
    bool planar = at -> HasDoubleBond ();
    bool linear = at -> HasBondOfOrder (3);


    if (atoms.size () && bonds.size ()) {
        assert (atoms.size() == bonds.size ());
        for (unsigned int i = 0; i< atoms.size (); i++) {
            mol -> add_atom_bonded_to (atoms[i], bonds[i], at);
        }
    }
    else {

        int ntoadd = at -> GetImplicitValence() - at -> GetValence();
        int nbrs = at -> GetValence ();
        if (ntoadd == 4) {
            Atom *dummy = new Atom;
            vect c = get_coordinates (at);
            c.x () -= 1.5;
            set_coordinates (dummy, c);

            vect coord2, coord3, coord4;
            four_coordinates (dummy, at, coord2, coord3, coord4);
            mol->add_atom_bonded_to (c, 1, at);
            mol->add_atom_bonded_to (coord2, 1, at);
            mol->add_atom_bonded_to (coord3, 1, at);
            mol->add_atom_bonded_to (coord4, 1, at);
            delete dummy;
        }


        else if (ntoadd == 3) {
            if (!nbrs) {
                Atom *dummy = new Atom;
                vect c = get_coordinates (dummy);
                c.x() -= 1.5;
                set_coordinates (dummy, c);
                vect coord2, coord3, coord4;
                four_coordinates (dummy, at, coord2, coord3, coord4);

                mol->add_atom_bonded_to (coord2, 1, at);
                mol->add_atom_bonded_to (coord3, 1, at);
                mol->add_atom_bonded_to (coord4, 1, at);

                delete dummy;

            }
            else {
                vect coord2, coord3, coord4;
                OBBondIterator i;
                Atom *neigh1 = at -> BeginNbrAtom (i);


                four_coordinates (neigh1, at, coord2, coord3, coord4);
                mol->add_atom_bonded_to (coord2, 1, at);
                mol->add_atom_bonded_to (coord3, 1, at);
                mol->add_atom_bonded_to (coord4, 1, at);

            }
        }
        else if (ntoadd == 2) {
            if (nbrs == 2) {
                vect coord1, coord2;

                OBBondIterator i;
                Atom *neigh1 = at -> BeginNbrAtom (i);
                Atom *neigh2 = at ->  NextNbrAtom (i);
                assert (neigh2 != at);
                assert (neigh1 != at);
                assert (neigh2 != neigh1);

                two_tetrahedral_coordinates (neigh1, neigh2, at, coord1, coord2);


                mol->add_atom_bonded_to (coord1, 1, at);
                mol->add_atom_bonded_to (coord2, 1, at);

            }                
            else if (planar) {
                vect coord1, coord2, coord3;
                OBBondIterator i;
                Atom *neigh1 = at -> BeginNbrAtom (i);
                three_coordinates (neigh1, at, coord1, coord2);

                mol->add_atom_bonded_to (coord1, 1, at);
                mol->add_atom_bonded_to (coord2, 1, at);

            }
            else if (!CountBonds (at)) {
                Atom *dummy = new Atom;
                vect v (-1.5, 0, 0);
                set_coordinates (dummy, v);
                vect coord1, coord2, coord3;
                four_coordinates (dummy, at, coord1, coord2, coord3);

                mol->add_atom_bonded_to (coord1, 1, at);
                mol->add_atom_bonded_to (coord2, 1, at);
                delete dummy;
            }
            else {
                vect coord1, coord2, coord3;
                OBBondIterator i;
                Atom *neigh1 = at -> BeginNbrAtom (i);
                four_coordinates (neigh1, at, coord1, coord2, coord3);

                mol->add_atom_bonded_to (coord1, 1, at);
                mol->add_atom_bonded_to (coord2, 1, at);
            }
        }
        else if (ntoadd == 1) {
            if (planar || linear || nbrs == 3) {
                vect coord1;
                one_H (at, coord1);

               mol->add_atom_bonded_to (coord1, 1, at);



            }
                else if (!nbrs) {
                vect coord;
                coord = get_coordinates (at);
                coord.x() -= 1.5;

               mol->add_atom_bonded_to (coord, 1, at);

            }
            else {
                vect coord1, coord2, coord3, coord4;
                OBBondIterator i;
                Atom *neigh1 = at -> BeginNbrAtom (i);
                four_coordinates (neigh1, at, coord2, coord3, coord4);

               mol->add_atom_bonded_to (coord2, 1, at);

    
            }
        }

//    mol->find_bound ();

    }
}






void Builder::delete_atom (Atom *at) {
    ZNMolecule *mol = (ZNMolecule *) at->GetParent ();
    int heavy_atoms = 0;
	FOR_ATOMS_OF_MOL (a, mol){
		if (!a -> IsHydrogen ()) {
			heavy_atoms++;
			if (heavy_atoms > 1) break;
		}
	}

    if (heavy_atoms <= 1) {
        DeleteZNMoleculeCommand *command = new DeleteZNMoleculeCommand (mol, ddwin);
        ddwin -> execute (command);
    }
    else {
        DeleteAtomCommand *command = new DeleteAtomCommand (at, ddwin);
        ddwin -> execute (command);
    }

}




void Builder::redefine_mol (ZNMolecule *mol) {
/* //   cout <<"find_bound"<<endl;
    mol->find_bound ();
    mol->number_atoms ();
    mol->number_bonds ();
//    cout <<"find_res"<<endl;
    mol->find_residues ();
//    cout <<"find_rings"<<endl;
    mol->find_rings ();
    mol->find_kekule ();
    mol->find_limits ();
    mol->find_center ();
    for (unsigned int i=0; i<mol->atoms.size (); i++) {
        mol->atoms[i]->find_mol2_type ();
    }
//    cout <<"find_strings"<<endl;
//    cout <<"find_init"<<endl;
    ddwin->data->mmff->initialize_mol (mol);
//    cout <<"done"<<endl;
*/
}

void Builder::add_bond (int order) {
	ZNMolecule *mol = ddwin ->target_molecule;
    if (mol ->selection && mol ->NumAtoms () ==2) {
		Atom *at1 = mol ->GetAtom (1); //atoms numeration migh change in next releases of openbabel to 0 and 1
		Atom *at2 = mol ->GetAtom (2);
		new_bond (at1, at2, order);
    }

}





void Builder::delete_Hs (Atom *at, vector <Atom *> &del_atoms, vector <ZNBond *> &del_bonds) {
    del_atoms.clear ();
    del_bonds.clear ();
    ZNMolecule* mol = (ZNMolecule *) at -> GetParent (); 

    FOR_NBORS_OF_ATOM (n, at) {
        if (n -> IsHydrogen ()) {
            del_atoms.push_back (&*n);
            ZNBond * b;
            b = mol -> GetBond (&*n, at); 
            del_bonds.push_back (b);

      //      mol -> RemoveBond (b);

        }
    }
    for (int i = (del_atoms.size () -1); i >= 0; i--) {

         mol -> RemoveAtom (del_atoms[i]);  

    }
 //   mol->find_bound ();
}

void Builder::one_H (Atom *center, vect &coord) {
    float dis = 1.5f;
    if (CountBonds (center)) {
        coord.null ();
        
        FOR_NBORS_OF_ATOM (n, center) {
            coord = sum (coord, get_coordinates (&*n));   
        }
        coord.multiply (1.f/(float) CountBonds (center));
        vect cv = get_coordinates(center);
        vect sumv = cv;
        sumv.multiply (2);
        coord = subtract (sumv, coord);
        vect coord_cent = subtract (coord, cv);
        float coor_mod = coord_cent.module ();
        coord = subtract (coord, cv);
        coord.multiply (dis / coor_mod);
        coord = sum (cv, coord);
        
    }
    else {
        coord = get_coordinates(center);
        coord.x() += 1.5f;
    }
}



void Builder::two_tetrahedral_coordinates (Atom *root1, Atom *root2, Atom *center, vect &coord2, vect &coord3) {
    float dis = 1.5f;
    vect r1c, r2c, cc;
    r1c = get_coordinates(root1);
    r2c = get_coordinates(root2);
    cc  = get_coordinates(center);
    vect middle_point = mean (r1c, r2c);
    coord2 = middle_point;
    coord3 = middle_point;
    vect root_axis1 = subtract (r1c, cc);
    vect root_axis2 = subtract (r2c, cc);
    vect rotation_axis = cross_product (root_axis1, root_axis2);
    rotate_around_vector (middle_point, rotation_axis, cc, PI/2);
    vect rot_axis2;
    rot_axis2 = subtract (middle_point, cc);

    double mid_mod = rot_axis2.module ();
    if (mid_mod < 0.00002) mid_mod = 0.00002;
    vect v = subtract (coord3, cc);
    v.multiply (dis / mid_mod);
    coord2 = sum (cc, v); 
    coord3 = coord2;
    rotate_around_vector (coord2, rot_axis2, cc, PI*(125.25)/180);
    rotate_around_vector (coord3, rot_axis2, cc, -PI*(125.25)/180); 

}

void Builder::three_coordinates (Atom* root, Atom *center, vect &coord2, vect &coord3) {
    float d = 1.5f;
    vect root_bond, rotation_axis, plane_bond, reference, rv, cv;
    rv = get_coordinates(root);
    cv = get_coordinates(center);
    if (CountBonds (root)) {
        OBBondIterator i;
        Atom *neigh1 = root -> BeginNbrAtom (i);
        reference = get_coordinates(neigh1);
    }        
    else {
        reference = subtract (rv, vect (1., 1., 1.));
    }
    root_bond = subtract (rv, cv);
    plane_bond = subtract (reference, rv);
    root_bond.scale_to (d);
    coord2 = sum (cv, root_bond);
    coord3 = coord2;
    rotation_axis = cross_product (root_bond, plane_bond);
    
    rotate_around_vector (coord2, rotation_axis, cv, 2*PI/3);
    rotate_around_vector (coord3, rotation_axis, cv, -2*PI/3);

}

void Builder::set_magic_pencil_atomic_number (int atomn) {
    magic_pencil_atomic_number = atomn;
    switch (atomn){
        case -2:
			ddwin ->gl ->setCursor(QCursor (QPixmap (":icons/pointer_rubber.png"), 2, 27));
            break;
        case 6:
			ddwin ->gl ->setCursor(QCursor (QPixmap (":icons/pencil_C.png"), 0, 29));
            break;
        case 7:
			ddwin ->gl ->setCursor(QCursor (QPixmap (":icons/pencil_N.png"), 0, 29));
            break;
        case 8:
			ddwin ->gl ->setCursor(QCursor (QPixmap (":icons/pencil_O.png"), 0, 29));
            break;
		case 16:
			ddwin ->gl ->setCursor(QCursor (QPixmap (":icons/pencil_S.png"), 0, 29));
            break;
    }

}

void Builder::four_coordinates (Atom* root, Atom *center, vect &coord2, vect &coord3, vect &coord4) {
    vect root_bond, rotation_axis, plane_bond, reference, rv, cv;
    rv = get_coordinates(root);
    cv = get_coordinates(center);
    float dist = 1.5f;
//    if (root -> NumAtoms ()) {
 //       reference = root -> BeginAtom () -> GetVector ();
 //   }        
    if (false) {}
    else {
        reference = subtract (rv, vect (1., 1., 1.));
    }
        root_bond = subtract (cv, rv);
        plane_bond = subtract (rv, reference);

    rotation_axis = cross_product (root_bond, plane_bond);
    root_bond.scale_to (dist);
    coord2 = subtract (cv, root_bond);
    coord3 = subtract (cv, root_bond);
    coord4 = subtract (cv, root_bond);

    if (!rotation_axis.square_module ()) rotation_axis.z () += 1;
    rotate_around_vector (coord2, rotation_axis, cv, PI*109.5/180);
    rotate_around_vector (coord3, rotation_axis, cv,  PI*109.5/180);
    rotate_around_vector (coord4, rotation_axis, cv,  PI*109.5/180);

    rotate_around_vector (coord3, root_bond, cv,  2*PI/3);
    rotate_around_vector (coord4, root_bond, cv,  -2*PI/3);

}
