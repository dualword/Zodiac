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

#include "database.h"


Database_molecule::Database_molecule () : ZNMolecule (), database (NULL) {
	multi = true;

}




Database_molecule::Database_molecule (const ZNMolecule &mol)
//: ZNMolecule( *(ZNMolecule::ZNMolecule *)&mol ) 
: ZNMolecule( mol ) 
{

	// OBMol::OBMol (*this);
	// (*(OBMol::OBMol *)this)(mol);
    multi = true;
}

Database_molecule &Database_molecule::operator =(const Database_molecule &mol) 
{
//	*(ZNMolecule::ZNMolecule *)this = *(ZNMolecule::ZNMolecule *)&mol;
    ZNMolecule::operator =(mol);
    multi = true;
	return *this;
}




Database::Database () : _needs_redraw (false), _hidden (true), _extend_enabled (false), dummy_mol (new Database_molecule ()){
	dummy_mol = new Database_molecule ();
	dummy_mol -> database = this;
	dummy_mol -> number = -1;
	finalise_molecule (dummy_mol);
	mutex = new QReadWriteLock ();
	field_names.push_back ("Molecule");
	sorting_column = new int (0);
}

void Database::safe_add_field (string name, vector <double> data) {
	mutex ->lockForWrite ();
	add_field (name, data);
	mutex ->unlock ();
}
void Database::add_field (string name, vector <double> data) {
	if (data.size () == count_entries ()) {
		field_names.push_back (name);
		for (unsigned int i = 0; i < data.size (); i++) {
		DoubleDatabaseCell *cell = new DoubleDatabaseCell (data[i]);
		entries[i] ->cells.push_back (cell);

		}

	}
	set_needs_redraw (true);
}

void Database::safe_add_field (string name, vector <string> data, char type) {
	mutex ->lockForWrite ();
	add_field (name, data, type);
	mutex ->unlock ();
}

void Database::add_field (string name, vector <string> data, char type) {

	if (data.size () == count_entries ()) {
		field_names.push_back (name);
		for (unsigned int i = 0; i < data.size (); i++) {
			if (type == 'd') {
				double doub = string_to_double (data[i]);
				DoubleDatabaseCell *cell = new DoubleDatabaseCell (doub);
				entries[i]->cells.push_back (cell);
			}
		}

	}
	set_needs_redraw (true);

}

ZNMolecule *Database::get_molecule (int n) {
	if (count_entries () > n) {
		return ((ZNMoleculeDatabaseCell *) entries[n]->cells[0])->get_molecule ();
	}
	else return NULL;
}

void Database::safe_add_mol (Database_molecule *mol) {
	mutex ->lockForWrite ();
	add_mol (mol);
	mutex ->unlock ();
}

void Database::add_mol (Database_molecule *mol) {

	DatabaseEntry *entry = new DatabaseEntry (mol);
	ZNMoleculeDatabaseCell *cell = new ZNMoleculeDatabaseCell (mol);
	entry ->cells.clear ();

    entry ->cells.push_back (cell);
	for (unsigned int i = 1; i < field_names.size (); i++) {
		DoubleDatabaseCell *cell2 = new DoubleDatabaseCell (0);
		entry ->cells.push_back (cell2);
	}
	entry ->sorting_column = sorting_column;
	//cerr << "assign "<<entry ->sorting_column <<" "<< &sorting_column;
	entries.push_back (entry);
    mol -> database = this;
    mol -> number = count_entries ()-1;
	set_needs_redraw (true);

}

void Database::delete_entry (int i) {
	if (0<=i<count_entries()) {
		delete entries [i];
		entries.erase (entries.begin ()+i);
	}
}

void Database::safe_sort_up (int column) {
		mutex ->lockForWrite ();
		sort_up (column);
		mutex ->unlock ();
}

void Database::safe_sort_down (int column) {
		mutex ->lockForWrite ();
		sort_down (column);
		mutex ->unlock ();
}

void Database::sort_up (int column) {
//	cerr <<"sorting up" << column <<endl;
	if (0 <= column < count_fields ()) {
		*sorting_column = column;
		sort (entries.begin (), entries.end (), compare_cell ());
		set_needs_redraw (true);	
	}
	//	cerr <<"sorting up end" << endl;

}


void Database::sort_down (int column) {
	//	cerr <<"sorting down" << column <<endl;
	if (0 <= column < count_fields ()) {	
		sort (entries.begin (), entries.end (), compare_cell ());
		reverse (entries.begin (), entries.end ());		
		set_needs_redraw (true);	
	}
	//	cerr <<"sorting down end" << endl;
}

void Database::safe_merge_with (Database *dat) {
		mutex ->lockForWrite ();
		merge_with (dat);
		mutex ->unlock ();
}


void Database::merge_with (Database *dat) {
	for (unsigned int i = 0; i < dat ->count_entries (); i++) {
		Database_molecule *new_mol = new Database_molecule (*dat ->get_molecule (i));
		add_mol (new_mol);
	}
}


bool Database::has_extend_enabled () {
	return _extend_enabled;
}

/*
void Database::renumber () {
	for (unsigned int i = 0; i < ref_molecules.size (); i++) {
		ref_molecules[i] -> number = i;
	}
}
*/

void Database::import_csv (string filename) {
	ifstream ifs (filename.c_str ());
	string line, token;
	vector < vector < string > > dats;
	vector <string> names;
	vector <string> value_line;
	getline (ifs, line);
	istringstream iss (line);
	while (getline (iss, token, ',')) {
		names.push_back (token);
	}
	for (unsigned int i = 0; i < names.size (); i++) {
		vector <string> dat;
		dats.push_back (dat);
	}
	while (getline(ifs, line)) {
		value_line.clear ();
		istringstream iss(line);
		string token;
		while (getline(iss, token, ',')){
			value_line.push_back (token);
		}
		if (value_line.size () == names.size ()) {
			for (unsigned int i = 0; i < value_line.size (); i ++) {
				dats[i].push_back (value_line[i]);
			}
		}
	}
	ifs.close();
	if (dats.size ()) {
		cerr <<dats[0].size ()<<" "<<count_entries ()<<endl;
		if (dats[0].size () == count_entries ()) {
			for (unsigned int i = 1; i < dats.size (); i++) {
				add_field (names[i], dats[i]);
			}
		}		
	}
}
