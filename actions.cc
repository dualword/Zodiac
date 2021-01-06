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

#include "actions.h"
class Data;

Actions::Actions (Data *dat) {
	data = dat;
}


double Actions::compute_total_energy (ZNMolecule *mol) {
Chemscore *chem =  new Chemscore;
	chem -> initialize_interaction (mol, data ->ddwin ->molecules);
	chem -> load_nonbonded_interactions ();
	float e = chem ->compute_total_energy ();
	delete chem;
	return e;
//	data->mmff->initialize (mol, data -> ddwin ->molecules);
//	return data -> mmff->compute_total_energy ();
}

vector <double> Actions::compute_total_energy (Database *dat) {
	vector <double> v;
	//	cerr << "entries "<< dat->count_entries ()<< endl;
	for (unsigned int i = 0; i < dat ->count_entries (); i ++) {
		v.push_back (compute_total_energy (dat -> get_molecule (i)));
	}
	return v;
}


double Actions::compute_logP (ZNMolecule *mol) {
	
	OBDescriptor* pDesc = OBDescriptor::FindType("logP");
	if(pDesc)
        return pDesc->Predict(mol);
	else return 0. ;
	
}

vector <double> Actions::compute_logP (Database *dat) {
	vector <double> v;
	for (unsigned int i = 0; i < dat -> count_entries (); i ++) {
		v.push_back (compute_logP (dat -> get_molecule (i)));
	}
	return v;
}



void Actions::minimise (ZNMolecule *mol) {
	
    data -> undo_stack -> beginMacro ("Energy Minimisation");
    MoveAtomsCommand *command = new MoveAtomsCommand (data -> ddwin -> gl, 1);
    FOR_ATOMS_OF_MOL(a, mol) {
        command -> add (&*a, get_coordinates(&*a));
    }
    data -> ddwin -> execute (command);
	
	//	cerr << "create" << endl;
    MinimiseThread *thread = new MinimiseThread (0, data -> ddwin);
	//	cerr << "set" << endl;
    thread -> set_molecule (mol);
	//	cerr << "go" << endl;
	data ->ddwin ->run_thread (thread);
	//	    thread -> start ();
    data -> ddwin -> connect (thread, SIGNAL (finished ()), data -> ddwin, SLOT (end_minimisation ()));
}




void Actions::minimise (Database *dat) {
	DatabaseMinimiseThread *thread = new DatabaseMinimiseThread (0, dat, data -> ddwin);
	data ->ddwin ->run_thread (thread);
	//thread -> start ();	
}





void Actions::change_display_style (ZNMolecule *mol, int ats, int bds, int backds) {
	ChangeDisplayStyleCommand *command = new ChangeDisplayStyleCommand (data -> ddwin -> gl);
	FOR_ATOMS_OF_MOL(a, mol) {
		//		int new_style = ats;
		command -> add (&*a, ats);
	}
	FOR_BONDS_OF_MOL (b, mol) {
		//	int new_style = bds;
		command -> add (&*b, bds);
	}
	if (mol -> selection) { //selected bonds are not in selection... we must find them checking the atoms. 
		FOR_ATOMS_OF_MOL(a, mol) {
			//		int new_style = ats;
			ZNMolecule * parent_mol = (ZNMolecule *) a -> GetParent ();
			FOR_NBORS_OF_ATOM (n, &*a) {
				
				if (get_selected (&*n)) {
					if (a ->GetIdx () < n -> GetIdx ()) {
						ZNBond *bo = parent_mol ->GetBond (&*n, &*a);
						command -> add (bo, bds);
					}
				}
			}
		}	
	}
	command -> add (mol, ats, bds, backds); 
	command -> set_name ();
	data ->ddwin ->execute (command);
	data ->ddwin ->gl ->draw_molecule (mol);
}




void Actions::change_display_style (Database *dat, int ats, int bds, int backds) {
	data -> undo_stack -> beginMacro ("Change Database Display Style");
	for (unsigned int i = 0; i < dat -> count_entries (); i ++) {
		change_display_style (dat ->get_molecule (i), ats, bds, backds);
	}
	data ->undo_stack ->endMacro ();
	
}



void Actions::hide_hydrogens (ZNMolecule *mol) {
	data ->ddwin ->gl ->hide_hydrogens (mol); 
}

void Actions::hide_hydrogens (Database *dat) {
	for (unsigned int i = 0; i < dat -> count_entries (); i ++) {
		hide_hydrogens (dat ->get_molecule (i));
	}
	
}

void Actions::hide_nonpolar_hydrogens (ZNMolecule *mol) {
	data ->ddwin ->gl ->hide_nonpolar_hydrogens (mol); 
}

void Actions::hide_nonpolar_hydrogens (Database *dat) {
	for (unsigned int i = 0; i < dat -> count_entries (); i ++) {
		hide_nonpolar_hydrogens (dat ->get_molecule (i));
	}
	
}

void Actions::hide_all_atoms (ZNMolecule *mol) {
	data ->ddwin ->gl ->hide_all_atoms (mol); 
}

void Actions::hide_all_atoms (Database *dat) {
	for (unsigned int i = 0; i < dat -> count_entries (); i ++) {
		hide_all_atoms (dat ->get_molecule (i));
	}
	
}



void Actions::show_all_atoms (ZNMolecule *mol) {
	data ->ddwin ->gl ->show_all_atoms (mol); 
}

void Actions::show_all_atoms (Database *dat) {
	for (unsigned int i = 0; i < dat -> count_entries (); i ++) {
		show_all_atoms (dat ->get_molecule (i));
	}
	
}

void Actions::apply_color_masks (vector <color_mask> masks, ZNMolecule *mol) {
	data ->ddwin ->gl ->apply_color_masks (masks, mol); 
}

void Actions::apply_color_masks (vector <color_mask> masks, Database *dat) {
	data ->undo_stack -> beginMacro ("Color Database");
	for (unsigned int i = 0; i < dat -> count_entries (); i ++) {
		apply_color_masks (masks, dat ->get_molecule (i));
	}
	data ->undo_stack -> endMacro ();
	
}

void Actions::save_as (ZNMolecule *mol, string filename) {
	
	ofstream ofs(filename.c_str ());
	OBConversion conv;
	conv.SetOutStream(&ofs);
	OBFormat* outFormat = conv.FormatFromExt(filename.c_str ());
	if (outFormat) {
		conv.SetOutFormat (outFormat);
		conv.Write (mol);
	}
	ofs.close ();
	
}

void Actions::save_session_as (string filename) {
	ofstream ofs(filename.c_str ());
	ofs << "#ROTATION_MATRIX "<<data ->ddwin ->gl ->Transform.M [0]<<" "<<data ->ddwin ->gl ->Transform.M [1]<<" "<<data ->ddwin ->gl ->Transform.M [2]<<" "<<
	data ->ddwin ->gl ->Transform.M [4]<<" "<<data ->ddwin ->gl ->Transform.M [5]<<" "<<data ->ddwin ->gl ->Transform.M [6]<<" "<<
	data ->ddwin ->gl ->Transform.M [8]<<" "<<data ->ddwin ->gl ->Transform.M [9]<<" "<<data ->ddwin ->gl ->Transform.M [10]<<endl;
	ofs << "#TRANSLATION "<<data ->ddwin ->gl ->view_translations.x()<<" "<<data ->ddwin ->gl ->view_translations.y()<<" "<<data ->ddwin ->gl ->view_translations.z()<<endl;
	ofs << "#CENTER_OF_ROTATION "<<data ->ddwin ->gl ->center_of_rotation.x()<<" "<<data ->ddwin ->gl ->center_of_rotation.y()<<" "<<data ->ddwin ->gl ->center_of_rotation.z()<<endl;
	for (unsigned int i = 1; i < data ->ddwin ->molecules.size (); i++) {
		ZNMolecule *mol = data ->ddwin ->molecules[i];
		if (!mol ->selection && !mol ->multi) {
			ofs << "#BEGIN_MOLECULE"<<endl;
			ofs << "#BEGIN_MOLECULAR_DATA"<<endl;
			OBConversion conv;
			conv.SetOutStream(&ofs);
			conv.SetOutFormat ("MOL2");
			conv.Write (mol);
			ofs << "#END_MOLECULAR_DATA"<<endl;
			int ads =  get_atoms_display_style (mol);
			int bds =  get_bonds_display_style (mol);			
			ofs << "#ATOMS_DISPLAY_STYLE "<< ads<< endl;
			ofs << "#BONDS_DISPLAY_STYLE "<< bds<< endl;
			FOR_ATOMS_OF_MOL (a, mol) {
				color c = get_color (&*a);
				int ds = get_ds (&*a);
				ofs << "#ATOMIC_COLOR "<<a->GetIdx ()<<" "<<c.redF ()<<" "<<c.greenF ()<<" "<<c.blueF ()<<" "<<c.alphaF ()<<endl;
				if (ds != ads) ofs << "#ATOMIC_DISPLAY_STYLE "<< a ->GetIdx ()<<" "<<ds<<endl;
			}
			FOR_BONDS_OF_MOL (b, mol) {
				int ds = get_ds (&*b);
				if (ds != bds) ofs << "#BOND_DISPLAY_STYLE "<< b ->GetIdx ()<<" "<<ds<<endl;
			}
			ofs << "#END_MOLECULE"<<endl;
		}		
	}
	for (unsigned int j = 0; j < data ->ddwin ->graphical_objects.size (); j++) {
		if (data ->ddwin ->graphical_objects[j] -> is_surface ()) {
			ofs <<"#BEGIN_SURFACE"<<endl;
			Surface *surf = (Surface *) data ->ddwin ->graphical_objects[j];
			for (unsigned int f = 0; f < surf ->faces.size (); f++) {
			 ofs <<"#FACE "<<surf ->faces[f] ->v1 -> coordinates.x ()<<" "<<surf ->faces[f] ->v1 -> coordinates.y ()<<" "<<surf ->faces[f] ->v1 -> coordinates.z ()<<" ";
			 ofs <<surf ->faces[f] ->v1 -> normal.x ()<<" "<<surf ->faces[f] ->v1 -> normal.y ()<<" "<<surf ->faces[f] ->v1 -> normal.z ()<<" ";
			 ofs <<surf ->faces[f] ->v1 -> col.redF ()<<" "<<surf ->faces[f] ->v1 -> col.greenF ()<<" "<<surf ->faces[f] ->v1 -> col.blueF ()<<" "<<surf ->faces[f] ->v1 -> col.alphaF ()<<" ";

			 ofs <<surf ->faces[f] ->v2 -> coordinates.x ()<<" "<<surf ->faces[f] ->v2 -> coordinates.y ()<<" "<<surf ->faces[f] ->v2 -> coordinates.z ()<<" ";
			 ofs <<surf ->faces[f] ->v2 -> normal.x ()<<" "<<surf ->faces[f] ->v2 -> normal.y ()<<" "<<surf ->faces[f] ->v2 -> normal.z ()<<" ";
			 ofs <<surf ->faces[f] ->v2 -> col.redF ()<<" "<<surf ->faces[f] ->v2 -> col.greenF ()<<" "<<surf ->faces[f] ->v2 -> col.blueF ()<<" "<<surf ->faces[f] ->v2 -> col.alphaF ()<<" ";

			 ofs <<surf ->faces[f] ->v3 -> coordinates.x ()<<" "<<surf ->faces[f] ->v3 -> coordinates.y ()<<" "<<surf ->faces[f] ->v3 -> coordinates.z ()<<" ";
			 ofs <<surf ->faces[f] ->v3 -> normal.x ()<<" "<<surf ->faces[f] ->v3 -> normal.y ()<<" "<<surf ->faces[f] ->v3 -> normal.z ()<<" ";
			 ofs <<surf ->faces[f] ->v3 -> col.redF ()<<" "<<surf ->faces[f] ->v3 -> col.greenF ()<<" "<<surf ->faces[f] ->v3 -> col.blueF ()<<" "<<surf ->faces[f] ->v3 -> col.alphaF ()<<endl;
			}
			ofs <<"#END_SURFACE"<<endl;
		} 
	}
	
	
	ofs.close ();
}

void Actions::load_session (string filename) {
	bool reading_molecule = false, reading_surf = false;
	//	bool matrix = false;
    data -> undo_stack -> beginMacro ("Load session ");
	ifstream ifs (filename.c_str ());
	string buffer;
	stringstream ss;
	ZNMolecule *last_mol = NULL;
	Surface *last_surf = NULL;
	while (getline(ifs, buffer)) {
		istringstream line(buffer);
		string token;
		line >> token;
		if (token == "#END_MOLECULAR_DATA") {
			ZNMolecule *mol = new ZNMolecule ();
			last_mol = mol;
			reading_molecule = false;
			OBConversion conv;
			stringstream ss1 (ss.str ());
			//	cerr << ss.str ()<<endl;
			conv.SetInStream(&ss1);
			conv.SetInFormat ("MOL2");
			conv.Read (mol);
			finalise_molecule (mol);
			CreateZNMoleculeCommand *command = new CreateZNMoleculeCommand (mol,  data ->ddwin);
			data -> ddwin -> execute (command);
			ss.str(" ");
			//		cerr << ss.str ()<<endl;
		}	
		else if (token == "#END_SURFACE" && reading_surf) {
			Surface *surf = last_surf;
			surf -> lst = data ->ddwin -> gl -> new_list ();
			surf ->render ();
			CreateGraphicalObjectCommand *command = new CreateGraphicalObjectCommand (surf,  data ->ddwin);
			data -> ddwin -> execute (command);		
			last_surf = NULL;
			reading_surf = false;
		}
		else if (reading_molecule) {
			ss << buffer<<endl;
		}
		else if (token == "#ATOMS_DISPLAY_STYLE" && last_mol) {
			int ds;
			line >> ds;
			set_atoms_display_style (last_mol, ds);
		}
		else if (token == "#BONDS_DISPLAY_STYLE" && last_mol) {
			int ds;
			line >> ds;
			set_bonds_display_style (last_mol, ds);
		}
		else if (token == "#ATOMIC_COLOR" && last_mol) {
			int n;
			float r, g, b, a;
			line >> n;
			line >> r;
			line >> g;
			line >> b;
			line >> a;			
			set_color (last_mol ->GetAtom (n), color (r, g, b, a));
		}
		else if (token == "#ATOMIC_DISPLAY_STYLE" && last_mol) {
			int n, ds;
			line >> n;
			line >> ds;
			set_ds (last_mol ->GetAtom (n), ds);
		}
		else if (token == "#BOND_DISPLAY_STYLE" && last_mol) {
			int n, ds;
			line >> n;
			line >> ds;	
			set_ds (last_mol ->GetBond (n), ds);
		}
		else if (token == "#FACE" and reading_surf) {
			SurfFace *face = new SurfFace;
			float x, y, z;
			float r, g, b, a;
			SurfVertex *v1 = new SurfVertex ();
			line >> x;
			line >> y;
			line >> z;
			v1 ->coordinates = vect (x, y, z);
			line >> x;
			line >> y;
			line >> z;
			v1 ->normal = vect (x, y, z);
			line >> r;
			line >> g;
			line >> b;
			line >> a;	
			v1 ->col = color (r, g, b, a);
			
			SurfVertex *v2 = new SurfVertex;
			line >> x;
			line >> y;
			line >> z;
			v2 ->coordinates = vect (x, y, z);
			line >> x;
			line >> y;
			line >> z;
			v2 ->normal = vect (x, y, z);
			line >> r;
			line >> g;
			line >> b;
			line >> a;	
			v2 ->col = color (r, g, b, a);
			
			SurfVertex *v3 = new SurfVertex;
			line >> x;
			line >> y;
			line >> z;
			v3 ->coordinates = vect (x, y, z);
			line >> x;
			line >> y;
			line >> z;
			v3 ->normal = vect (x, y, z);
			line >> r;
			line >> g;
			line >> b;
			line >> a;	
			v3 ->col = color (r, g, b, a);
			face ->v1 = v1;
			face ->v2 = v2;
			face ->v3 = v3;
			last_surf ->faces.push_back (face);
			last_surf ->vertices.push_back (v1);
			last_surf ->vertices.push_back (v2);
			last_surf ->vertices.push_back (v3);			
			
		}		
		else if (token == "#ROTATION_MATRIX") {
			float m0, m1, m2, m3, m4, m5, m6, m7, m8;
			line >> m0;
			line >> m1;
			line >> m2;
			line >> m3;
			line >> m4;
			line >> m5;
			line >> m6;
			line >> m7;
			line >> m8;	
			ChangeFloatCommand *command0 = new ChangeFloatCommand (data ->ddwin ->gl ->Transform.M [0], m0, data ->ddwin ->gl);
			ChangeFloatCommand *command1 = new ChangeFloatCommand (data ->ddwin ->gl ->Transform.M [1], m1, data ->ddwin ->gl);	
			ChangeFloatCommand *command2 = new ChangeFloatCommand (data ->ddwin ->gl ->Transform.M [2], m2, data ->ddwin ->gl);	
			ChangeFloatCommand *command3 = new ChangeFloatCommand (data ->ddwin ->gl ->Transform.M [4], m3, data ->ddwin ->gl);	
			ChangeFloatCommand *command4 = new ChangeFloatCommand (data ->ddwin ->gl ->Transform.M [5], m4, data ->ddwin ->gl);	
			ChangeFloatCommand *command5 = new ChangeFloatCommand (data ->ddwin ->gl ->Transform.M [6], m5, data ->ddwin ->gl);	
			ChangeFloatCommand *command6 = new ChangeFloatCommand (data ->ddwin ->gl ->Transform.M [8], m6, data ->ddwin ->gl);	
			ChangeFloatCommand *command7 = new ChangeFloatCommand (data ->ddwin ->gl ->Transform.M [9], m7, data ->ddwin ->gl);	
			ChangeFloatCommand *command8 = new ChangeFloatCommand (data ->ddwin ->gl ->Transform.M [10], m8, data ->ddwin ->gl);	
			ChangeFloatCommand *command10 = new ChangeFloatCommand (data ->ddwin ->gl ->LastRot.M [0], m0, data ->ddwin ->gl);
			ChangeFloatCommand *command11 = new ChangeFloatCommand (data ->ddwin ->gl ->LastRot.M [1], m1, data ->ddwin ->gl);	
			ChangeFloatCommand *command12 = new ChangeFloatCommand (data ->ddwin ->gl ->LastRot.M [2], m2, data ->ddwin ->gl);	
			ChangeFloatCommand *command13 = new ChangeFloatCommand (data ->ddwin ->gl ->LastRot.M [3], m3, data ->ddwin ->gl);	
			ChangeFloatCommand *command14 = new ChangeFloatCommand (data ->ddwin ->gl ->LastRot.M [4], m4, data ->ddwin ->gl);	
			ChangeFloatCommand *command15 = new ChangeFloatCommand (data ->ddwin ->gl ->LastRot.M [5], m5, data ->ddwin ->gl);	
			ChangeFloatCommand *command16 = new ChangeFloatCommand (data ->ddwin ->gl ->LastRot.M [6], m6, data ->ddwin ->gl);	
			ChangeFloatCommand *command17 = new ChangeFloatCommand (data ->ddwin ->gl ->LastRot.M [7], m7, data ->ddwin ->gl);	
			ChangeFloatCommand *command18 = new ChangeFloatCommand (data ->ddwin ->gl ->LastRot.M [8], m8, data ->ddwin ->gl);
			data -> ddwin -> execute (command0);
			data -> ddwin -> execute (command1);
			data -> ddwin -> execute (command2);
			data -> ddwin -> execute (command3);
			data -> ddwin -> execute (command4);
			data -> ddwin -> execute (command5);
			data -> ddwin -> execute (command6);
			data -> ddwin -> execute (command7);
			data -> ddwin -> execute (command8);
			data -> ddwin -> execute (command10);
			data -> ddwin -> execute (command11);
			data -> ddwin -> execute (command12);
			data -> ddwin -> execute (command13);
			data -> ddwin -> execute (command14);
			data -> ddwin -> execute (command15);
			data -> ddwin -> execute (command16);
			data -> ddwin -> execute (command17);
			data -> ddwin -> execute (command18);																																															
		}
		else if (token == "#BEGIN_MOLECULAR_DATA") {
			reading_molecule = true;
		}
		else if (token == "#BEGIN_SURFACE") {
			last_surf = new Surface ();
			reading_surf = true;
		}
		else if (token == "#TRANSLATION") {
			float trvx, trvy, trvz;
			line >> trvx;
			line >> trvy;
			line >> trvz;
			ChangeVectorCommand *command = new ChangeVectorCommand (data ->ddwin ->gl ->view_translations, vect (trvx, trvy, trvz), data ->ddwin ->gl);
			data -> ddwin -> execute (command);
		}
		else if (token == "#CENTER_OF_ROTATION") {
			float crx, cry, crz;
			line >> crx;
			line >> cry;
			line >> crz;
			ChangeVectorCommand *command = new ChangeVectorCommand (data ->ddwin ->gl ->center_of_rotation, vect (crx, cry, crz), data ->ddwin ->gl);
			data -> ddwin -> execute (command);
		}
		
	}
	data ->undo_stack ->endMacro ();
}

void Actions::save_as (Database *dat, string filename) {
	ofstream ofs(filename.c_str ());
	OBConversion conv;
	conv.SetOutStream(&ofs);
	OBFormat* outFormat = conv.FormatFromExt(filename.c_str ());
	conv.SetOutFormat (outFormat);
	for (unsigned int i = 0; i < dat -> count_entries (); i ++) {
		conv.Write (dat ->get_molecule (i));
	}
	string csv_name = filename+".csv"; 
	save_csv (dat, csv_name);
}

void Actions::save_csv (Database *dat, string filename) {
	ofstream f (filename.c_str ());
	for (unsigned int i = 0; i < dat -> field_names.size (); i ++) {
		f << dat ->field_names[i];
		if (i != dat -> count_fields ()-1) f<<",";
	}
	f<<endl;
	for (unsigned int j = 0; j < dat -> count_entries (); j ++) {
		for (unsigned int i = 0; i < dat -> count_fields (); i ++) {
			f << dat ->entries[j]->cells[i]->get_string ();
			if (i != dat -> count_fields ()-1) f<<",";
		}
		f<<endl;
	}		
}

void Actions::reprotonate (ZNMolecule *mol, double ph) {
	ReprotonateCommand *command = new ReprotonateCommand (mol, data -> ddwin, ph);
	data -> ddwin -> execute (command);
}



void Actions::reprotonate (Database *dat, double ph) {
	data ->undo_stack -> beginMacro ("Reprotonate Database");
	for (unsigned int i = 0; i < dat -> count_entries (); i ++) {
		ReprotonateCommand *command = new ReprotonateCommand (dat ->get_molecule (i), data -> ddwin, ph);
		data -> ddwin -> execute (command);
	}
	data ->undo_stack -> endMacro ();
	
}



void Actions::set_scores_from_charges (ZNMolecule *mol) {
	
}

void Actions::set_scores_from_charges (Database *dat) {
	
}
