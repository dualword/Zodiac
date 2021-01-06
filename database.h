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

#ifndef DATABASE_H
#define DATABASE_H

#include "menu.h"

class DatabaseGrid;

class Database_molecule : public Molecule {
    public:
    Database_molecule ();
    Database *database;
    int number;

};

class DatabaseField {
	public:
	DatabaseField (string, vector <string>);
	vector <string> data;
	string name;
};

class Database {
public:
	Data *data;
    Database (Data *);
	DatabaseGrid *grid;
	string name;
    void add_mol (Database_molecule *mol);
    vector <Database_molecule *> molecules;
//	vector <Database_molecule *> ref_molecules;    
	void set_graphics ();
	bool has_extend_enabled ();
	vector <DatabaseField> fields;
	void add_field (string, vector <double>);
	void import_csv (string filename);
	//void renumber ();
	
};


#endif