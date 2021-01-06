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


double Actions::compute_total_energy (Molecule *mol) {
	data->mmff->initialize (mol, data -> ddwin ->molecules);
	return data -> mmff->compute_total_energy ();
}

vector <double> Actions::compute_total_energy (Database *dat) {
	vector <double> v;
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
		v.push_back (compute_total_energy (dat -> molecules[i]));
	}
	return v;
}


double Actions::compute_logP (Molecule *mol) {
        OBLogP logp;
        return logp.Predict (*mol);
}

vector <double> Actions::compute_logP (Database *dat) {
	vector <double> v;
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
		v.push_back (compute_logP (dat -> molecules[i]));
	}
	return v;
}



void Actions::minimise (Molecule *mol) {

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
	    thread -> start ();
    data -> ddwin -> connect (thread, SIGNAL (finished ()), data -> ddwin, SLOT (end_minimisation ()));
}




void Actions::minimise (Database *dat) {
	DatabaseMinimiseThread *thread = new DatabaseMinimiseThread (0, dat, data -> ddwin);
	thread -> start ();	
}





void Actions::change_display_style (Molecule *mol, int ats, int bds) {
	ChangeDisplayStyleCommand *command = new ChangeDisplayStyleCommand (data -> ddwin -> gl);
	FOR_ATOMS_OF_MOL(a, mol) {
//		int new_style = ats;
		command -> add (&*a, ats);
	}
	FOR_BONDS_OF_MOL (b, mol) {
	//	int new_style = bds;
		command -> add (&*b, bds);
	}
	command -> add (mol, ats, bds); 
	command -> set_name ();
	data ->ddwin ->execute (command);
	data ->ddwin ->gl ->draw_molecule (mol);
}




void Actions::change_display_style (Database *dat, int ats, int bds) {
	data -> undo_stack -> beginMacro ("Change Database Display Style");
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
		change_display_style (dat -> molecules[i], ats, bds);
	}
	data ->undo_stack ->endMacro ();
	
}



void Actions::hide_hydrogens (Molecule *mol) {
	data ->ddwin ->gl ->hide_hydrogens (mol); 
}

void Actions::hide_hydrogens (Database *dat) {
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
		hide_hydrogens (dat -> molecules[i]);
	}
	
}

void Actions::hide_nonpolar_hydrogens (Molecule *mol) {
	data ->ddwin ->gl ->hide_nonpolar_hydrogens (mol); 
}

void Actions::hide_nonpolar_hydrogens (Database *dat) {
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
		hide_nonpolar_hydrogens (dat -> molecules[i]);
	}
	
}

void Actions::hide_all_atoms (Molecule *mol) {
	data ->ddwin ->gl ->hide_all_atoms (mol); 
}

void Actions::hide_all_atoms (Database *dat) {
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
		hide_all_atoms (dat -> molecules[i]);
	}
	
}



void Actions::show_all_atoms (Molecule *mol) {
	data ->ddwin ->gl ->show_all_atoms (mol); 
}

void Actions::show_all_atoms (Database *dat) {
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
		show_all_atoms (dat -> molecules[i]);
	}
	
}

void Actions::apply_color_masks (vector <color_mask> masks, Molecule *mol) {
	data ->ddwin ->gl ->apply_color_masks (masks, mol); 
}

void Actions::apply_color_masks (vector <color_mask> masks, Database *dat) {
	data ->undo_stack -> beginMacro ("Color Database");
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
		apply_color_masks (masks, dat -> molecules[i]);
	}
	data ->undo_stack -> endMacro ();
	
}

void Actions::save_as (Molecule *mol, string filename) {
	
	ofstream ofs(filename.c_str ());
	OBConversion conv;
	conv.SetOutStream(&ofs);
	OBFormat* outFormat = conv.FormatFromExt(filename.c_str ());
	if (outFormat) {
		conv.SetOutFormat (outFormat);
		conv.Write (mol);
	}
	
}

void Actions::save_as (Database *dat, string filename) {
	ofstream ofs(filename.c_str ());
	OBConversion conv;
	conv.SetOutStream(&ofs);
	OBFormat* outFormat = conv.FormatFromExt(filename.c_str ());
	conv.SetOutFormat (outFormat);
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
		conv.Write (dat ->molecules[dat -> grid -> real_index_of_line (i)]);
	}
	string csv_name = filename+".csv"; 
	save_csv (dat, csv_name);
}

void Actions::save_csv (Database *dat, string filename) {
	ofstream f (filename.c_str ());
	for (unsigned int i = 0; i < dat -> fields.size (); i ++) {
		f << dat ->fields[i].name;
		if (i != dat -> fields.size ()-1) f<<",";
	}
	f<<endl;
	for (unsigned int j = 0; j < dat -> molecules.size (); j ++) {
		for (unsigned int i = 0; i < dat -> fields.size (); i ++) {
			f << dat ->fields[i].data[dat -> grid -> real_index_of_line (j)];
			if (i != dat -> fields.size ()-1) f<<",";
		}
		f<<endl;
	}		
}

void Actions::reprotonate (Molecule *mol) {
	ReprotonateCommand *command = new ReprotonateCommand (mol, data -> ddwin);
	data -> ddwin -> execute (command);
}



void Actions::reprotonate (Database *dat) {
	data ->undo_stack -> beginMacro ("Reprotonate Database");
	for (unsigned int i = 0; i < dat -> molecules.size (); i ++) {
	ReprotonateCommand *command = new ReprotonateCommand (dat ->molecules[i], data -> ddwin);
	data -> ddwin -> execute (command);
	}
	data ->undo_stack -> endMacro ();

}



void Actions::set_scores_from_charges (Molecule *mol) {

}

void Actions::set_scores_from_charges (Database *dat) {

}
