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

#ifndef ACTIONS_H
#define ACTIONS_H


#include "ZNdata.h"
#include "obabel_includes.h"

class Actions {
	public:
	Actions (Data *dat);
	
	
	double compute_total_energy (ZNMolecule *mol);
	vector <double> compute_total_energy (Database *dat);
	
	double compute_logP (ZNMolecule *mol);
	vector <double> compute_logP (Database *dat);	
	
	
	void reprotonate (ZNMolecule *mol, double ph = 7.4);
	void reprotonate (Database *dat, double ph = 7.4);
	
	void minimise (ZNMolecule *mol);
	void minimise (Database *dat);	
	
	void change_display_style (Database *dat, int, int, int);	
	void change_display_style (ZNMolecule *mol, int, int, int);

	void hide_hydrogens (ZNMolecule *mol);
	void hide_hydrogens (Database *dat);
	
	void hide_nonpolar_hydrogens (ZNMolecule *mol);
	void hide_nonpolar_hydrogens (Database *dat);
	
	void hide_all_atoms (ZNMolecule *mol);
	void hide_all_atoms (Database *dat);
	
	void show_all_atoms (ZNMolecule *mol);
	void show_all_atoms (Database *dat);
			
			
	void apply_color_masks (vector <color_mask> masks, ZNMolecule *mol);		
	void apply_color_masks (vector <color_mask> masks, Database *dat);
	
	void save_session_as (string filename);
	void load_session (string filename);
	void save_as (ZNMolecule *mol, string filename);
	void save_as (Database *dat, string filename);
	void save_csv (Database *dat, string filename);
	
	void set_scores_from_charges (ZNMolecule *mol);
	void set_scores_from_charges (Database *dat);
	
	
	private:	
	Data *data;
};









#endif

