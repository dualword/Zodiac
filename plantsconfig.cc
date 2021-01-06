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

#include <iostream>
#include <sstream>
#include "plantsconfig.h"

#include <string>




Confline::Confline (char ctype,  string cfilestr, string cdefval, MyLine* ctarget)
{
	defval = cdefval;filestr = cfilestr; targetl = ctarget; type = ctype;

	Default ();
}

Confline::Confline (char ctype,  string cfilestr,string cdefval, QPushButton* ctarget)
{
	defval = cdefval;filestr = cfilestr; targetb = ctarget;type = ctype;
	Default ();
}

Confline::Confline (char ctype,  string cfilestr,string cdefval, QComboBox* ctarget)
{
	defval = cdefval;filestr = cfilestr; targetcb = ctarget;type = ctype;
	Default ();
}

Confline::Confline (char ctype,  string cfilestr, string cdefval,QLineEdit* ctarget)
{
	defval = cdefval;filestr = cfilestr; targetle = ctarget;type = ctype;
	Default ();
}

Confline::Confline (char ctype,  string cfilestr, string cdefval,MyLineF* ctarget)
{
	defval = cdefval;filestr = cfilestr; targetlf = ctarget; type = ctype;
	Default ();
}

Confline::Confline (char ctype,  string cfilestr, string cdefval,Q3ListBox* ctarget)
{
	defval = cdefval;filestr = cfilestr; targetlb = ctarget; type = ctype;
	Default ();
}

Confline::Confline (char ctype,  string cfilestr, string cdefval,PlantsMainWin* ctarget)
{
	defval = cdefval;filestr = cfilestr; targetmw = ctarget; type = ctype;
	Default ();
}

int Confline::getInt (string mystr){
	istringstream iss (mystr);
	int i;
	iss>>i;
	return i;
}


void Confline::Default () {
	current = defval;
	if (type=='m' || type=='M' ) {vect.clear();};
}



void Confline::Apply () {
	if (type=='l') {
		targetl->ins ( QString(current.c_str()) );

	}

	if (type=='f') {
		targetlf->ins ( QString(current.c_str()) );

	}
	else if (type=='b') {

		targetb->setChecked (getInt (current));

	}
	else if (type=='c') {
        int i = targetcb -> findText ( QString(current.c_str()) );
        if (i > -1) targetcb ->setCurrentIndex (i);


	}

	else if (type=='m') {
		targetlb->clear ();
		for (unsigned int j=0; j<vect.size(); j++)
		{
			targetlb->insertItem( QString(vect[j].c_str()) );

		}
	}

	else if (type=='M') {
		targetlb->clear ();
		for (unsigned int j=0; j<vect.size(); j++){
			targetlb->insertItem(new MyListBoxItem (vect[j], targetlb));

		}
	}



	else if (type=='s') {
		istringstream line(current);
		string nbsx, nbsy, nbsz;
		line >> nbsx >>nbsy >> nbsz;
		targetmw->bsx->ins ("");
		targetmw->bsy->ins ("");
		targetmw->bsz->ins ("");
		targetmw->bsx->ins ( QString(nbsx.c_str()) );
		targetmw->bsy->ins ( QString(nbsy.c_str()) );
		targetmw->bsz->ins ( QString(nbsz.c_str()) );

	}


}


void Confline::Get () {
	if (type=='l') {
		current = targetl->val ();

	}

	if (type=='f') {
		current = targetlf->val ();

	}
	else if (type=='b') {
		current = boolstr (targetb->isChecked ());      


	}
	else if (type=='c') {

		QString t = targetcb->currentText ();    current = t.toStdString().c_str ();


	}

	else if (type=='m') {

		if (targetlb->count ()) {
			vect.clear ();
			for (unsigned int j=0; j<targetlb->count() ; j++){
				vect.push_back (targetlb->text(j).toStdString());
			}
		}

	}

	else if (type=='M') {

		if (targetlb->count ()) {
			vect.clear ();
			for (unsigned int j=0; j<targetlb->count() ; j++){
				vect.push_back (targetlb->item (j)->text ().toStdString());
				cout << "check this" << targetlb->item (j)->text ().toStdString() << endl;
			}
		}


	}


	else if (type=='s') {
		istringstream line(current);
		string nbsx, nbsy, nbsz;
		line >> nbsx >>nbsy >> nbsz;
		current = targetmw->bsx->val ();
		current.append (" ");
		current.append (targetmw->bsy->val ());
		current.append (" ");
		current.append (targetmw->bsz->val ());
	}


}



string Confline ::boolstr (bool boo){
	if (boo) return "1";
	else return "0";
}


////////////////////////////////////CONFIGURATION//////////////////////////////////////////////////////////////////

Configuration::Configuration (Plants *plant) {
	plants = plant;
	//  setdata (dat);


	
	
	conf.push_back (new Confline ('l',"aco_ants","20", plants->mainwin->aants));
    conf.push_back (new Confline ('l',"aco_evap", "0.15", plants->mainwin->aevap));
    conf.push_back (new Confline ('l',"aco_sigma", "1.0", plants->mainwin->asigma)); 
    conf.push_back (new Confline ('b',"flip_amide_bonds", "0", plants->mainwin->flipab));  
    conf.push_back (new Confline ('b',"flip_planar_n", "1", plants->mainwin->flipn));
    conf.push_back (new Confline ('b',"force_flipped_bonds_planarity", "0", plants->mainwin->ffbp));
    conf.push_back (new Confline ('c',"scoring_function", "chemplp", plants->mainwin->sfchoice));
    conf.push_back (new Confline ('b',"chemplp_weak_cho", "1", plants->mainwin->cho));
    conf.push_back (new Confline ('c',"ligand_intra_score", "clash2", plants->mainwin->lintra));
    conf.push_back (new Confline ('l',"chemplp_charged_hb_weight", "1.75", plants->mainwin->chbw));
    conf.push_back (new Confline ('l',"chemplp_charged_metal_weight", "2.0", plants->mainwin->cmw));
    conf.push_back (new Confline ('l',"chemplp_hbond_weight", "-4.0", plants->mainwin->hbw));
    conf.push_back (new Confline ('l',"chemplp_hbond_cho_weight", "-2.0", plants->mainwin->hbchow));
    conf.push_back (new Confline ('l',"chemplp_metal_weight", "-6.0", plants->mainwin->mw));
    conf.push_back (new Confline ('l',"chemplp_plp_weight", "1.0", plants->mainwin->plpw));
    conf.push_back (new Confline ('l',"chemplp_intercept_weight", "-20.0", plants->mainwin->iw));
    conf.push_back (new Confline ('l',"cluster_rmsd", "2.0", plants->mainwin->crmsd));
    conf.push_back (new Confline ('l',"cluster_structures", "5", plants->mainwin->cs));
	
	
	

	conf.push_back (new Confline ('M',"ifile", "", plants->mainwin->iflb));
	conf.push_back (new Confline ('f',"protein_file", "", plants->mainwin->pfile));
	conf.push_back (new Confline ('l',"output_dir", "", plants->mainwin->odir));
	conf.push_back (new Confline ('s',"bindingsite_center", "", plants->mainwin));
	conf.push_back (new Confline ('l',"bindingsite_radius", "", plants->mainwin->bsr));
	conf.push_back (new Confline ('b',"write_protein_conformations", "1", plants->mainwin->wpc));
	conf.push_back (new Confline ('b',"write_rescored_structures", "0", plants->mainwin->wrs));
	conf.push_back (new Confline ('b',"write_multi_mol2", "1", plants->mainwin->wmm2));
	conf.push_back (new Confline ('b',"write_ranking_links", "0", plants->mainwin->wrl));
	conf.push_back (new Confline ('b',"write_ranking_multi_mol2", "0", plants->mainwin->wrmm2));
	conf.push_back (new Confline ('b',"write_protein_bindingsite", "1", plants->mainwin->wpb));
	conf.push_back (new Confline ('b',"write_protein_splitted", "1", plants->mainwin->wps));
	conf.push_back (new Confline ('b',"write_per_atom_scores", "0", plants->mainwin->wpas));

	conf.push_back (new Confline ('m',"constr", "", plants->mainwin->constrlb));

	conf.push_back (new Confline ('m',"flex", "", plants->mainwin->flexlb));
	conf.push_back (new Confline ('b',"rigid_ligand", "0", plants->mainwin->rigidl));
	conf.push_back (new Confline ('b',"rigid_all", "0", plants->mainwin->rigida));
	conf.push_back (new Confline ('l',"intra_protein_score_weight", "0.4", plants->mainwin->ipsw));

	conf.push_back (new Confline ('f',"water_molecule_definition", "", plants->mainwin->wref));
	conf.push_back (new Confline ('m',"water", "", plants->mainwin->waterlb));
	conf.push_back (new Confline ('l',"water_protein_hb_weight", "1.0", plants->mainwin->wphb));
	conf.push_back (new Confline ('l',"water_ligand_hb_weight", "1.0", plants->mainwin->wlhb));
	conf.push_back (new Confline ('l',"water_water_hb_weight", "1.0", plants->mainwin->wwhb));
	conf.push_back (new Confline ('l',"no_water_ligand_hb_penalty", "0.0", plants->mainwin->nwhbp));









	Default ();
}

//void Configuration::setdata (Data *datas){
//    data = datas;
//    data->config = this;
//   win->conf = this;

//}


void Configuration::Default () {

	for (unsigned i=0; i<conf.size (); i++){
		conf[i]->Default ();

	}


}


void Configuration::Apply (){

	for (unsigned i=0; i<conf.size (); i++){
		conf[i]->Apply ();
		plants->mainwin->check_ligands ();
		//       data->ddwin->load_waters ();
	}


}


void Configuration::Get (){

	for (unsigned i=0; i<conf.size (); i++){
		conf[i]->Get ();
	}
}



void Configuration::Read (const char * filename){
	if (filename) {
		ifstream ifs (filename);
		string buffer;
		vector<string> ifiles, vflex, vconstr, vwater;
		ifiles.clear(); vflex.clear (); vconstr.clear (); vwater.clear ();
		while (getline(ifs, buffer)) {
			istringstream line2(buffer);
			string token;
			line2 >> token;


			if (token == "ligand_list" || token == "ligand_file") {
				string s = token;
				while (!line2.eof ()) {
					string q;
					line2 >> q;
					s.append (" ");s.append (q);
				}
				ifiles.push_back (s);
			} 




			else if (token == "chemplp_protein_hb_constraint" || token == "shape_constraint" || token == "surface_distance_constraint" || token == "ligand_intra_distance_constraint" || token == "protein_ligand_distance_constraint") {
				string s = token;
				while (!line2.eof ()) {
					string q;
					line2 >> q;
					s.append (" ");s.append (q);
				}
				vconstr.push_back (s);
			} 

			else if (token == "flexible_protein_side_chain_number" || token == "fixed_protein_bond" || token == "flexible_protein_side_chain_string") {
				string s = token;
				while (!line2.eof ()) {
					string q;
					line2 >> q;
					s.append (" ");s.append (q);
				}
				vflex.push_back (s);
			} 

			else if (token == "water_molecule") {
				string s = token;
				vector <float> wat;
				wat.clear ();
				while (!line2.eof ()) {

					float f;
					string q;

					line2 >> q;
					istringstream lin (q);
					lin >> f;   
					wat.push_back (f);
					s.append (" ");s.append (q);
				}
				vwater.push_back (s);

				//           data->waters.push_back(wat);
			} 



			for (unsigned i=0; i<conf.size (); i++){
				if (token == conf[i]->filestr && conf[i]->type!='s') {
					line2>>conf[i]->current;
					break;  
				}
				if  (token == conf[i]->filestr && conf[i]->type=='s') { 
					string x, y, z;  
					line2>>x>>y>>z;
					conf[i]->current=x;conf[i]->current.append (" ");
					conf[i]->current.append(y);conf[i]->current.append (" "); conf[i]->current.append(z);
				}


			}





			for (unsigned i=0; i<conf.size (); i++){
				if (conf[i]->filestr == "ifile") {
					conf[i]->vect = ifiles;

				}         
				if (conf[i]->filestr == "constr") {
					conf[i]->vect = vconstr;

				}    
				if (conf[i]->filestr == "flex") {
					conf[i]->vect = vflex;

				}   
				if (conf[i]->filestr == "water") {
					conf[i]->vect = vwater;

				}    

			}
		}

	}
}






void Configuration::Write (const char * filename){

	ofstream *file = new ofstream(filename);
	*file << "#Generated by PLANTS FE version " <<VERSION<< endl<<endl;

	for (unsigned int i=0;i<conf.size(); i++){
		if (i == 0){*file << endl << endl << "#ALGORITHM" << endl;};
		if (i == 17){*file << endl << endl << "#INPUT-OUTPUT" << endl;};
		if (i == 29){*file << endl << endl << "#CONSTRAINTS && FLEXIBILITY" << endl;};
		if (i == 33){*file << endl << endl << "#WATER" << endl;};
		if (conf[i]->type=='m' || conf[i]->type=='M' ){
			for (unsigned int j=0;j<conf[i]->vect.size(); j++){ 
				*file << conf[i]->vect[j] << endl;  
			}
		}
		else{
			*file << conf[i]->filestr << " " << conf[i]->current << endl;  
		}
	}







	file->close(); 

}


