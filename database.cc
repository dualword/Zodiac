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


DatabaseField::DatabaseField (string nam, vector <string> dat) {
	name = nam;
	data = dat;
}


Database_molecule::Database_molecule () : Molecule () {
    multi = true;
    database = NULL;

}

Database::Database (Data *dat) {
	data = dat;
	grid = NULL;
}


void Database::add_field (string name, vector <double> data) {
	if (data.size () == molecules.size ()) {
		vector <string> data_s;
		for (unsigned int i = 0; i < data.size (); i++) {
			stringstream ss;
			ss << data [i];
			data_s.push_back (ss.str ());
		}
		DatabaseField df (name, data_s);
		fields.push_back (df);
		grid -> add_field (&df);
		
	}
}


void Database::add_mol (Database_molecule *mol) {
//	ref_molecules.push_back (mol);
    molecules.push_back (mol);
    mol -> database = this;
    mol -> number = molecules.size ()-1;
}


void Database::set_graphics () {
	if (!grid) {
		grid = new DatabaseGrid (0, this, data);
	}
	grid -> show ();
	grid -> raise ();
}

bool Database::has_extend_enabled () {
	return grid -> has_extend_enabled ();
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
		if (dats[0].size () == molecules.size ()) {
			for (unsigned int i = 0; i < dats.size (); i++) {
				DatabaseField df (names[i], dats[i]);
				fields.push_back (df);
				grid -> add_field (&df);
			}
		}		
	}
}
