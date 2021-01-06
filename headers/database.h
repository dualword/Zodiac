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


class Database_molecule : public ZNMolecule {
    public:
	Database_molecule (const ZNMolecule &mol);
	Database_molecule &operator =(const Database_molecule &mol);
    Database_molecule ();
    Database *database;
    int number;

};


class DatabaseCell {
public: 
	DatabaseCell () {};
	virtual void set_value (string s) {};
	virtual string get_string () {return "no data";};
	virtual double get_double () {return 0.;};
	virtual void set_double (double d) {};
};




class DatabaseEntry {
	public:
	DatabaseEntry (ZNMolecule *mol) {_mol = mol;};
	vector <DatabaseCell *> cells;
	ZNMolecule *_mol;
	ZNMolecule *get_molecule () {return _mol;};
	int *sorting_column;

};

class Database {
public:

    Database ();
	Database_molecule *dummy_mol;
	string name;
    void add_mol (Database_molecule *mol);
	void delete_entry (int i);
	void safe_add_mol (Database_molecule *mol);
//    vector <Database_molecule *> molecules; 
	bool has_extend_enabled ();
	vector <DatabaseEntry *> entries;
	vector <string> field_names;
	void add_field (string, vector <double>);
	void add_field (string, vector <string>, char type = 'd');
	void safe_add_field (string, vector <double>);
	void safe_add_field (string, vector <string>, char type = 'd');
	void safe_merge_with (Database *dat);
	void merge_with (Database *dat);
	ZNMolecule *get_molecule (int);
	void import_csv (string filename);
	bool &get_extend_enabled () {return _extend_enabled;}
	bool &get_is_hidden () {return _hidden;}
	bool get_needs_redraw () {return _needs_redraw;}
	void set_needs_redraw (bool b) {_needs_redraw = b;}
	int count_fields () {return field_names.size ();}
	int count_entries () {return entries.size ();};
	int *sorting_column;
	//void renumber ();
	QReadWriteLock *mutex;
	void safe_sort_up (int column);
	void safe_sort_down (int column);
	void sort_up (int column);
	void sort_down (int column);
	void set_double (int field, int entry, double value) {if (entry < count_entries () && (field < count_fields())) entries [entry] ->cells [field] ->set_double (value);};


private:
	bool _extend_enabled, _hidden;
	bool _needs_redraw;

	
};


struct compare_cell
{
    bool operator()( const DatabaseEntry* a, const DatabaseEntry* b)
    {
		cerr << *a ->sorting_column << " "<<*b ->sorting_column<<" " <<a ->sorting_column << " "<<b ->sorting_column<<endl;
           return (a->cells[*a ->sorting_column]) ->get_double () < (b->cells[*b ->sorting_column]) ->get_double ();
    }
};

class ZNMoleculeDatabaseCell: public DatabaseCell {
	public:
	ZNMoleculeDatabaseCell (ZNMolecule *mol) {_molecule = mol;};
	string get_string () {return _molecule ->GetTitle ();};
	ZNMolecule *get_molecule () {return _molecule;};
	double get_double () {if (_molecule ->multi) return ((Database_molecule *) _molecule) -> number; else return _molecule->NumAtoms ();}
	void set_value (string s) {_molecule ->SetTitle (s);};
private:
	ZNMolecule *_molecule;
};

class DoubleDatabaseCell: public DatabaseCell {
public:
	DoubleDatabaseCell (double d) {_value = d;};
	string get_string () {stringstream ss; ss << _value; return ss.str ();}
	double get_double () {return _value;}
	void set_value (string s) {_value = string_to_double (s);};
	void set_double (double d) {_value = d;};
private:
	double _value;
	
} ;

#endif