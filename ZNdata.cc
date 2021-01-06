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

#include "ZNdata.h"




////////////////////////////////////DATA//////////////////////////////////////////////////////////////////

//class MyListBoxItem;

Data::Data (QApplication *master) : QObject (){
	haptic_position_lock = new QReadWriteLock ();
    qapp = master;
	ncounter =0.f;

    twodwin = new MainWindow ();
	current_force_x = current_force_y = current_force_z = 0.;	
	current_position_x = current_position_y = current_position_z = 0.;	
	current_yaw = current_pitch = current_roll = 0.;
	last_yaw = last_pitch = last_roll = 0.;
	
	total_energy_haptic = 0.;
    mmff = new MMFF ();
	actions = new Actions (this); 


    minimize = new Minimize (this);


    charge_begin = -80.f;
    charge_end   = +80.f;

    score_begin = -80.f;
    score_mid   = +0.f ;
    score_end   = +80.f;


    constant_color = color (0, 0, 0);
    charge_begin_color = color (0, 0, 255);
    charge_end_color = color (255, 0, 0);
    score_begin_color = color (0, 255, 0);
    score_mid_color = color (200, 200, 200);
    score_end_color = color (255, 0, 0);


	undo_stack = new QUndoStack ();
	load_defaults ();
	load_preferences ();


}

void Data::load_haptic_device () {
#ifdef HAPI

   // Create a new haptics device, using any device connected.

	haptic_device = new (AnyHapticsDevice);

   // initialize the device
	if (haptic_device ->initDevice () != HAPIHapticsDevice::SUCCESS) {


       // initilization failed, print error message and quit
		cerr << haptic_device ->getLastErrorMsg () << endl;
   	}
	haptic_device ->enableDevice ();
//	H3DUtil::SimpleThread a_simple_thread(haptic_callback, (void *) this);

	HapticForceThread *hapt_thread = new HapticForceThread (this, ddwin);

	hapt_thread ->start ();

#endif //HAPI
}

void Data::set_ddwin (DDWin *ddw){
    ddwin =ddw;
    ddwin->data = this;
   // minimize -> haptic_thread = new HapticThread (0, ddwin);
	
	  connect(twodwin->to3DAct,SIGNAL(triggered()), this, SLOT(from_2D_to_3D()));

}

void Data::from_2D_to_3D () {
	vector <OpenBabel::OBMol *>mols = twodwin ->get_OB_mols ();
	for (unsigned int i = 0; i < mols.size (); i++) {
		ZNMolecule *add_mol = (ZNMolecule *) (mols[i]);
		add_mol ->multi = false;
		add_mol ->selection = false;
		finalise_molecule (add_mol);
		actions -> reprotonate (add_mol);
				finalise_molecule (add_mol);
		ddwin ->add_molecule (add_mol);
		ddwin ->set_current_target (add_mol);
	}
}

void Data::load_defaults () {
	for (unsigned int i = 0; i < vars.size (); i++) {
		delete vars[i];
	}
	vars.clear ();
	vars.push_back (new StringDataVar ("plants_exe", ""));
	ColorDataVar *bgvar = new ColorDataVar ("background_color", color (255, 255, 255, 255));
	background_color = bgvar ->value_ptr ();
	vars.push_back (bgvar);
	DoubleDataVar *qsvar = new DoubleDataVar ("quality_scale", 5.0);
	quality_scale = qsvar ->value_ptr ();
	vars.push_back (qsvar);

	DoubleDataVar *vdwsvar = new DoubleDataVar ("vdw_scale", 0.2);
	vdw_scale = vdwsvar ->value_ptr ();
	vars.push_back (vdwsvar);
	DoubleDataVar *srvar = new DoubleDataVar ("stick_radius", 0.12);
	stick_radius = srvar ->value_ptr ();
	vars.push_back (srvar);		
	DoubleDataVar *sprvar = new DoubleDataVar ("sphere_radius", 0.12);
	sphere_radius = sprvar ->value_ptr ();
	vars.push_back (sprvar);	
	
	DoubleDataVar *dbssvar = new DoubleDataVar ("double_bond_stick_scale", 0.7);
	double_bond_stick_scale = dbssvar ->value_ptr ();
	vars.push_back (dbssvar);
	
	DoubleDataVar *dbsvar = new DoubleDataVar ("double_bond_separation", 0.15);
	double_bond_separation = dbsvar ->value_ptr ();
	vars.push_back (dbsvar);
	
	DoubleDataVar *lwvar = new DoubleDataVar ("line_width", 1.2);
	line_width = lwvar ->value_ptr ();
	vars.push_back (lwvar);		
	
	DoubleDataVar *backbone_tube_helix_a_var = new DoubleDataVar ("backbone_tube_helix_a", 1.0);
	backbone_tube_helix_a = backbone_tube_helix_a_var ->value_ptr ();
	vars.push_back (backbone_tube_helix_a_var);	
	
	DoubleDataVar *backbone_tube_helix_b_var = new DoubleDataVar ("backbone_tube_helix_b", .3);
	backbone_tube_helix_b = backbone_tube_helix_b_var ->value_ptr ();
	vars.push_back (backbone_tube_helix_b_var);	
	
	DoubleDataVar *backbone_tube_helix_c_var = new DoubleDataVar ("backbone_tube_helix_c", 1.0);
	backbone_tube_helix_c = backbone_tube_helix_c_var ->value_ptr ();
	vars.push_back (backbone_tube_helix_c_var);	
	
	DoubleDataVar *backbone_tube_sheet_a_var = new DoubleDataVar ("backbone_tube_sheet_a", 1.0);
	backbone_tube_sheet_a = backbone_tube_sheet_a_var ->value_ptr ();
	vars.push_back (backbone_tube_sheet_a_var);	
	
	DoubleDataVar *backbone_tube_sheet_b_var = new DoubleDataVar ("backbone_tube_sheet_b", .3);
	backbone_tube_sheet_b = backbone_tube_sheet_b_var ->value_ptr ();
	vars.push_back (backbone_tube_sheet_b_var);	
	
	DoubleDataVar *backbone_tube_sheet_c_var = new DoubleDataVar ("backbone_tube_sheet_c", 1.0);
	backbone_tube_sheet_c = backbone_tube_sheet_c_var ->value_ptr ();
	vars.push_back (backbone_tube_sheet_c_var);	
	
	DoubleDataVar *backbone_tube_random_a_var = new DoubleDataVar ("backbone_tube_random_a", .2);
	backbone_tube_random_a = backbone_tube_random_a_var ->value_ptr ();
	vars.push_back (backbone_tube_random_a_var);	
	
	DoubleDataVar *backbone_tube_random_b_var = new DoubleDataVar ("backbone_tube_random_b", .2);
	backbone_tube_random_b = backbone_tube_random_b_var ->value_ptr ();
	vars.push_back (backbone_tube_random_b_var);	
	
	DoubleDataVar *backbone_tube_random_c_var = new DoubleDataVar ("backbone_tube_random_c", 0.0);
	backbone_tube_random_c = backbone_tube_random_c_var ->value_ptr ();
	vars.push_back (backbone_tube_random_c_var);	
	
}

void Data::write_preferences () {
	string pref_name = "Zodiac.ini";

	ofstream ofs(pref_name.c_str ());
	for (unsigned int i=0; i<vars.size (); i++) {
		ofs << vars[i] ->write_string ()<<endl;
	}

	ofs.close ();

}

void Data::load_preferences () {
	string pref_name = "Zodiac.ini";

	ifstream ifs (pref_name.c_str ());
	if (ifs.good ()) {
		string buffer;
		while (getline(ifs, buffer)) {
			for (unsigned int i=0; i<vars.size (); i++) {
					vars[i] ->load (buffer);
			}
		}	
	}

}

void *haptic_callback (void * data) {

#ifdef HAPI
     Data *dat = (Data *) data; 



Vec3 position = dat ->haptic_device->getRawPosition();
Vec3 angles (0.f, 0.f, 0.f);
//Rotation haptic_device ->getOrientation();

	dat ->haptic_position_lock ->lockForWrite ();
	dat ->current_position_x = position.x*2000;
	dat ->current_position_y = position.y*2000;
	dat ->current_position_z = position.z*2000;
	dat ->haptic_position_lock ->unlock ();


//	dat ->current_position_y = 0.;
//	dat ->current_position_z = 0.;

//	dat ->current_position_x = 0;
//	dat ->current_position_y = 0;
//	dat ->current_position_z = 0;

	//cerr << position << endl;
	//cerr << position.x << " position"<<endl;

	dat ->current_pitch = angles[2];
	dat ->current_roll = angles[0];
	dat ->current_yaw = angles[1];

	dat ->current_pitch = 0.;
	dat ->current_roll = 0.;
	dat ->current_yaw = 0.;

	//Vec3 force_to_render (dat ->current_force_x, dat ->current_force_y, dat ->current_force_z);
	//dat ->haptic_device ->sendForce (force_to_render);

#endif //HAPI
	return 0;

}


void ColorDataVar::load (string buffer) {
	istringstream line(buffer);
	string token;
	line >> token;
	
	if (token == name) {
		int count = 0;
		vector <int> vals;
		string q;
		while (!line.eof ()) {
			count++;
			line >> q;
			vals.push_back (string_to_int (q));
		}
		if (count < 3) cerr <<name<<" not enough values";
		else if (count > 4 ) cerr <<name<< "too many values";
		else if (count == 4) {_value = color (vals [0], vals[1], vals [2], vals[3]); cerr <<vals [0]<<" "<<vals[1]<<" "<<vals[2]<<" "<<vals[3]<<endl; }
		else {_value = color (vals [1], vals[2], vals [3]);}
	}

}
	

void StringDataVar::load (string buffer) {
	istringstream line(buffer);
	string token;
	line >> token;

	if (token == name) {
		int count = 0;
		string q;
		while (!line.eof ()) {
			count++;

			line >> q;
		}
		if (!count) cerr <<name<<" no value";
		else if (count > 1 ) cerr <<name<< "multiple values";
		else _value = q;
	}

}


	
void DoubleDataVar::load (string buffer) {
	istringstream line(buffer);
	string token;
	line >> token;

	if (token == name) {
		int count = 0;
		string q;
		while (!line.eof ()) {
			count++;

			line >> q;
		}
		if (!count) cerr <<name<<" no value";
		else if (count > 1 ) cerr <<name<< "multiple values";
		else _value = string_to_double (q);
	}

}


/*
void Data::set_FF (MMFF *mmf, TriposFF *triposf, PLP *pl){
    mmff =mmf;
    triposff =triposf;
    plp = pl;


}
*/

ZNMolecule* Data::check_mol2 (string filename){
    return 0;

 /*   cout << "reading file "<<filename<<endl; 
    ZNMolecule *mol = new ZNMolecule ();
    ifstream ifs (filename.c_str() );
    bool b = false;
    bool multi = (mol->readMultiMOL2 (filename, ifs, b)==2);
    unsigned int asize =mol->atoms.size ();
    if (asize>0) {
        mol->valid=true;
    }
    if (multi) {
        cout <<" multimol2"<<endl;
        delete mol;
        bool go_on = true;
        Database* database = new Database ();
        ifstream mfs (filename.c_str() );
        b = false;
        while (go_on) {
            Database_molecule *new_mol = new Database_molecule ();
            cout <<"b="<<b<<endl;
            new_mol->readMultiMOL2 (filename, mfs, b); 
            b = true;
            if (new_mol->atoms.size () < 1) {
                go_on = false;
                delete new_mol;
            }
            else {
                new_mol->database = database;
                new_mol->number = database->molecules.size ()+1;
                new_mol->valid = true;
                database->molecules.push_back (new_mol);
           //     cerr << database->molecules.size ()<<endl;
            }
        }
        mol = (ZNMolecule *) database->molecules[0];
    }    
    return mol;
*/
}



/*

////////////////////////////////////CONFLINE//////////////////////////////////////////////////////////////////



#include <iostream>
#include <sstream>


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

Confline::Confline (char ctype,  string cfilestr, string cdefval,MainWin* ctarget)
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
        if (type=='m' || type=='M' ) {vector.clear();};
        }



void Confline::Apply () {
    if (type=='l') {
        targetl->ins (current);
 
    }

    if (type=='f') {
        targetlf->ins (current);
 
    }
    else if (type=='b') {
        
        targetb->setOn (getInt (current));
 
    }
    else if (type=='c') {
        
        targetcb->setCurrentText (current);
 
    }

    else if (type=='m') {
        targetlb->clear ();
    for (unsigned int j=0; j<vector.size(); j++){
        targetlb->insertItem(vector[j]);

        }
       }

    else if (type=='M') {
        targetlb->clear ();
    for (unsigned int j=0; j<vector.size(); j++){
        targetlb->insertItem(new MyListBoxItem (vector[j], targetlb));

        }
       }


 
    else if (type=='s') {
    istringstream line(current);
    string nbsx, nbsy, nbsz;
    line >> nbsx >>nbsy >> nbsz;
    targetmw->bsx->ins ("");
    targetmw->bsy->ins ("");
    targetmw->bsz->ins ("");
    targetmw->bsx->ins (nbsx);
    targetmw->bsy->ins (nbsy);
    targetmw->bsz->ins (nbsz);
 
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
      current = boolstr (targetb->isOn ());      

 
    }
    else if (type=='c') {
        
        QString t = targetcb->currentText ();    current = t.latin1();

 
    }

    else if (type=='m') {

        if (targetlb->count ()) {
            vector.clear ();
            for (unsigned int j=0; j<targetlb->count() ; j++){
                vector.push_back (targetlb->text (j));
            }
        }

    }

    else if (type=='M') {

        if (targetlb->count ()) {
            vector.clear ();
            for (unsigned int j=0; j<targetlb->count() ; j++){
                vector.push_back (targetlb->item (j)->text ());
                cout << "check this" << targetlb->item (j)->text () << endl;
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

Configuration::Configuration (Data *dat) {
    setdata (dat);



    conf.push_back (new Confline ('l',"aco_ants","20", data->window->mainwin->aants));
    conf.push_back (new Confline ('l',"aco_evap", "0.15", data->window->mainwin->aevap));
    conf.push_back (new Confline ('l',"aco_sigma", "1.0", data->window->mainwin->asigma)); 
    conf.push_back (new Confline ('b',"flip_amide_bonds", "0", data->window->mainwin->flipab));  
    conf.push_back (new Confline ('b',"flip_planar_n", "1", data->window->mainwin->flipn));
    conf.push_back (new Confline ('b',"force_flipped_bonds_planarity", "0", data->window->mainwin->ffbp));
    conf.push_back (new Confline ('c',"scoring_function", "chemplp", data->window->mainwin->sfchoice));
    conf.push_back (new Confline ('b',"chemplp_weak_cho", "0", data->window->mainwin->cho));
    conf.push_back (new Confline ('c',"ligand_intra_score", "clash", data->window->mainwin->lintra));
    conf.push_back (new Confline ('l',"chemplp_charged_hb_weight", "1.5", data->window->mainwin->chbw));
    conf.push_back (new Confline ('l',"chemplp_charged_metal_weight", "1.5", data->window->mainwin->cmw));
    conf.push_back (new Confline ('l',"chemplp_hbond_weight", "-4.0", data->window->mainwin->hbw));
    conf.push_back (new Confline ('l',"chemplp_hbond_cho_weight", "-3.0", data->window->mainwin->hbchow));
    conf.push_back (new Confline ('l',"chemplp_metal_weight", "-9.0", data->window->mainwin->mw));
    conf.push_back (new Confline ('l',"chemplp_plp_weight", "0.7", data->window->mainwin->plpw));
    conf.push_back (new Confline ('l',"chemplp_intercept_weight", "-20.0", data->window->mainwin->iw));
    conf.push_back (new Confline ('l',"cluster_rmsd", "2.0", data->window->mainwin->crmsd));
    conf.push_back (new Confline ('l',"cluster_structures", "5", data->window->mainwin->cs));


    conf.push_back (new Confline ('M',"ifile", "", data->window->mainwin->iflb));
    conf.push_back (new Confline ('f',"protein_file", "", data->window->mainwin->pfile));
    conf.push_back (new Confline ('l',"output_dir", "", data->window->mainwin->odir));
    conf.push_back (new Confline ('s',"bindingsite_center", "", data->window->mainwin));
    conf.push_back (new Confline ('l',"bindingsite_radius", "", data->window->mainwin->bsr));
    conf.push_back (new Confline ('b',"write_protein_conformations", "1", data->window->mainwin->wpc));
    conf.push_back (new Confline ('b',"write_rescored_structures", "0", data->window->mainwin->wrs));
    conf.push_back (new Confline ('b',"write_multi_mol2", "1", data->window->mainwin->wmm2));
    conf.push_back (new Confline ('b',"write_ranking_links", "0", data->window->mainwin->wrl));
    conf.push_back (new Confline ('b',"write_ranking_multi_mol2", "0", data->window->mainwin->wrmm2));
    conf.push_back (new Confline ('b',"write_protein_bindingsite", "1", data->window->mainwin->wpb));
    conf.push_back (new Confline ('b',"write_protein_splitted", "1", data->window->mainwin->wps));
    conf.push_back (new Confline ('b',"write_per_atom_scores", "0", data->window->mainwin->wpas));

    conf.push_back (new Confline ('m',"constr", "", data->window->mainwin->constrlb));

    conf.push_back (new Confline ('m',"flex", "", data->window->mainwin->flexlb));
    conf.push_back (new Confline ('b',"rigid_ligand", "0", data->window->mainwin->rigidl));
    conf.push_back (new Confline ('b',"rigid_all", "0", data->window->mainwin->rigida));
    conf.push_back (new Confline ('l',"intra_protein_score_weight", "0.4", data->window->mainwin->ipsw));

    conf.push_back (new Confline ('f',"water_molecule_definition", "", data->window->mainwin->wref));
    conf.push_back (new Confline ('m',"water", "", data->window->mainwin->waterlb));
    conf.push_back (new Confline ('l',"water_protein_hb_weight", "1.0", data->window->mainwin->wphb));
    conf.push_back (new Confline ('l',"water_ligand_hb_weight", "1.0", data->window->mainwin->wlhb));
    conf.push_back (new Confline ('l',"water_water_hb_weight", "1.0", data->window->mainwin->wwhb));
    conf.push_back (new Confline ('l',"no_water_ligand_hb_penalty", "0.0", data->window->mainwin->nwhbp));









    Default ();
}

void Configuration::setdata (Data *datas){
    data = datas;
    data->config = this;
 //   win->conf = this;

}


void Configuration::Default () {

    for (unsigned i=0; i<conf.size (); i++){
        conf[i]->Default ();
       
    }


}


void Configuration::Apply (){

    for (unsigned i=0; i<conf.size (); i++){
        conf[i]->Apply ();
        data->window->mainwin->check_ligands ();
        data->ddwin->load_waters ();
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
                while (not line2.eof ()) {
                    string q;
				    line2 >> q;
                    s.append (" ");s.append (q);
                }
                ifiles.push_back (s);
            } 
                


			
			else if (token == "chemplp_protein_hb_constraint" || token == "shape_constraint" || token == "surface_distance_constraint" || token == "ligand_intra_distance_constraint" || token == "protein_ligand_distance_constraint") {
                string s = token;
                while (not line2.eof ()) {
                    string q;
				    line2 >> q;
                    s.append (" ");s.append (q);
                }
                vconstr.push_back (s);
            } 
                
			else if (token == "flexible_protein_side_chain_number" || token == "fixed_protein_bond" || token == "flexible_protein_side_chain_string") {
                string s = token;
                while (not line2.eof ()) {
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
                while (not line2.eof ()) {

                    float f;
                    string q;

				    line2 >> q;
                    istringstream lin (q);
                    lin >> f;   
                    wat.push_back (f);
                     s.append (" ");s.append (q);
                }
                vwater.push_back (s);
  
                data->waters.push_back(wat);
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
                    conf[i]->vector = ifiles;
                       
                }         
                if (conf[i]->filestr == "constr") {
                    conf[i]->vector = vconstr;
                       
                }    
                if (conf[i]->filestr == "flex") {
                    conf[i]->vector = vflex;
                       
                }   
                if (conf[i]->filestr == "water") {
                    conf[i]->vector = vwater;
                       
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
           for (unsigned int j=0;j<conf[i]->vector.size(); j++){ 
                  *file << conf[i]->vector[j] << endl;  
            }
        }
        else{
            *file << conf[i]->filestr << " " << conf[i]->current << endl;  
        }
    }





    


    file->close(); 

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

*/

bool is_db_extended (ZNMolecule *mol) {
			bool ext = false;
			if (mol -> multi) {
				Database_molecule *dm = (Database_molecule *) mol;
				if (dm -> database -> has_extend_enabled ()) ext = true;
			}
			return ext;
}
