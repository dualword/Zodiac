
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

#include "menu.h"
#include "ddwin.h"

#define WIREFRAME 0
#define STICK 1
#define CPK 2
#define BALLANDSTICK 3
#define BALLANDLINE 4

#define NO_ATOMS 0
#define SPHERES 1
#define CPK_SPHERES 2
#define SCALED_CPK_SPHERES 3

#define NO_BONDS 0
#define LINES 1
#define STICKS 2

#define BACKBONE_LINE 0
#define BACKBONE_STICK 1


#define AROMATIC_RINGS 0
#define KEKULE 1
#define AROMATIC_BONDS 2

#define ELEMENT 0
#define CHARGE  1
#define SCORE   2
#define COLOR   3

#define FLOAT  32




ZNWidget::ZNWidget (QWidget *parent) : QWidget (parent) {
	QVBoxLayout *layout = new QVBoxLayout;
	setLayout (layout);
 	layout -> setAlignment (Qt::AlignTop);
}




ZNAdvancedWidget::ZNAdvancedWidget (const char *advance, bool show) : ZNWidget () {

	_private = new ZNWidget ();
	_public = new ZNWidget ();
	_private ->hide();


	layout () -> addWidget (_public);
 	layout () -> setAlignment (Qt::AlignTop);

	basic_layout () ->setContentsMargins (0, 0, 0, 0);
	basic_layout () -> setAlignment (Qt::AlignTop);


	if (show == true) {
		_group_box = new MyGroupBox (layout (), advance, true, false);
		_group_box -> layout () -> addWidget (_private);
		_private ->layout () -> setAlignment (Qt::AlignTop);
		advanced_layout () ->setContentsMargins (0, 0, 0, 0);

		connect (_group_box, SIGNAL (clicked (bool)), this, SLOT (hide_advanced_slot () ) );
	}
}

void ZNAdvancedWidget::hide_advanced_slot () {
	if (_group_box ->isChecked () == true) {
		_private ->show();
	}
	else {
		_private ->hide();
	}

}



ZNMenu::ZNMenu (QWidget *parent, Data *dat, bool tabs, bool advanced) : QMainWindow (parent) {
	if (advanced) _main_widget = new ZNAdvancedWidget ();
	else _main_widget = new ZNWidget ();
	setWindowTitle ("Zodiac");
	_about_action = new QAction (tr ("&About"), this);

	 _about_action ->setShortcut (tr ("Control+A"));
	connect (_about_action, SIGNAL (triggered ()), this, SLOT (_show_about ()));


	_data = dat;

	_tabs = new QTabWidget ();

	if (tabs)	setCentralWidget (_tabs);
	else setCentralWidget (_main_widget);

//	_tabs ->addTab (_main_tab, "Main");
//	ZNAdvancedWidget *_test_tab = new ZNAdvancedWidget ();
//  	_tabs ->addTab (_test_tab, "Test");

}



void ZNMenu::add_help () {
    	QMenu *about = new QMenu(tr("&Help"), this );
    	Q_CHECK_PTR( about );
    	about -> addAction (_about_action);
	menuBar () ->addMenu (about);
}



bool ZNMenu::check_for_a_set_target () {
	return _data ->ddwin ->current_target;
}


void ZNMenu::_load_bot_lf (string& buffer, const char *reference, MyLineFile *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		string s = token;
		while (!line2.eof ()) {
			string q;
			line2 >> q;
			s.append (" ");s.append (q);
			target ->linedit ->setText(q.c_str ());
		}
	}
}


void ZNMenu::_load_bot_fel (string& buffer, const char *reference, MyFloatEditLine *target, const char *sep) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		string s = token;
		while (!line2.eof ()) {
			string q;
			line2 >> q;
//			s.append (sep);s.append (q);
			s.append (" ");s.append (q);
			target ->spinbox ->setValue(QString::fromStdString (q).toDouble());
		}
	}
}


void ZNMenu::_load_bot_iel (string& buffer, const char *reference, MyIntegerEditLine *target, const char *sep) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		string s = token;
		while (!line2.eof ()) {
			string q;
			line2 >> q;
			s.append (" ");s.append (q);
			target ->spinbox ->setValue(QString::fromStdString (q).toInt());
		}
	}
}


void ZNMenu::_load_bot_chb (string& buffer, const char *reference, MyCheckBox *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		string s = token;
		while (!line2.eof ()) {
			string q;
			line2 >> q;
			s.append (" ");s.append (q);
			if (q == "0") target ->setCheckState (Qt::Unchecked);
			if (q == "1") target ->setCheckState (Qt::Checked);
		}
	}
}


void ZNMenu::_load_bot_box (string& buffer, const char *reference, MyGroupBox *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		string s = token;
		while (!line2.eof ()) {
			string q;
			line2 >> q;
			s.append (" ");s.append (q);
			if (q == "0") target ->setChecked (false);
			if (q == "1") target ->setChecked (true);
		}
	}
}


void ZNMenu::_load_bot_cb (string& buffer, const char *reference, MyComboBox *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		string s = token;
		while (!line2.eof ()) {
			string q;
			line2 >> q;
			s.append (" ");s.append (q);
// le combo box non funzionano!?!?!?!?!?!?!?
		}
	}
}

void ZNMenu::_load_bot_le (string& buffer, const char *reference, MyLineEdit *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		string s = token;
		while (!line2.eof ()) {
			string q;
			line2 >> q;
			s.append (" ");s.append (q);
			target ->linedit ->setText(q.c_str ());
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	CONFORMATIONAL SEARCHES
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


 ConformersMenu::ConformersMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, false, false) {
 	setWindowTitle ("Conformational Search");
	add_help ();
	
	_type_hcb = new MyHideComboBox (main_widget () ->layout (), "Search type");
	_stochastic_widget = new ZNWidget ();	
	_systematic_widget = new ZNWidget ();
	main_widget ()->layout () -> addWidget (_stochastic_widget);
	main_widget ()->layout () -> addWidget (_systematic_widget);
	

	_type_hcb ->insertItem (_stochastic_widget, 0, "Stochastic", "Stochastic");
	_type_hcb ->insertItem (_systematic_widget, 1, "Systematic", "Systematic");
	
	//stochastic search
	_results = 20;
	_results_iel = new MyIntegerEditLine (_stochastic_widget ->layout (), "Number of results:", _results, 0);
	
	
	_timlim = 60;
	_timlim_iel = new MyIntegerEditLine (_stochastic_widget ->layout (), "Time limit (seconds):", _timlim, 0);
	
	
	
	
	
	QPushButton *_ok_btn = new QPushButton("Start");
		main_widget ()->layout () -> addWidget (_ok_btn);
	connect (_ok_btn, SIGNAL( clicked() ), this, SLOT( start_thread () ) );
 }

void ConformersMenu::start_thread () {
	if (_type_hcb ->currentData ().toString ().toStdString () == "Stochastic") {
		ZNMolecule *mol = _data -> ddwin ->target_molecule;
		Database *database = new Database;
		_data ->ddwin ->add_database(database);
		StochasticConformationThread *thread = new StochasticConformationThread (0, mol, database, _data ->ddwin) ;
		thread ->set_results_number (_results);
		thread ->set_time_limit (_timlim);
		_data ->ddwin ->run_thread (thread);
//		thread ->start ();
	}
	else if (_type_hcb ->currentData ().toString ().toStdString () == "Systematic") {
		ZNMolecule *mol = _data -> ddwin ->target_molecule;
		Database *database = new Database;
		_data ->ddwin ->add_database(database);
		SystematicConformationThread *thread = new SystematicConformationThread (0, mol, database, _data ->ddwin) ;
		_data ->ddwin ->run_thread (thread);
//		thread ->start ();
	} 
}



DockingMenu::DockingMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, true, true) {
	ZNAdvancedWidget *main_tab = new ZNAdvancedWidget ("", false);
   	_tabs ->addTab (main_tab, "Input");
	add_menu ();
	add_help ();
	_binding_site_box = new MyGroupBox (main_tab ->layout (), "Binding site options");
	_site_x_el = new MyFloatEditLine (_binding_site_box ->layout (), "Binding site center X:", _x, -1000, 1000);
	_site_y_el = new MyFloatEditLine (_binding_site_box ->layout (), "Binding site center Y:", _y, -1000, 1000);
	_site_z_el = new MyFloatEditLine (_binding_site_box ->layout (), "Binding site center Z:", _z, -1000, 1000);
	_site_r_el = new MyFloatEditLine (_binding_site_box ->layout (), "Binding site radius:", _r, -1000, 1000);
	_start_docking_bt = new MyPushButton (main_tab ->layout (), "Start Docking");
	connect (_start_docking_bt ->fbutton, SIGNAL( clicked() ), SLOT( start_docking_slot() ) );
}

void DockingMenu::start_docking_slot () {
		Database *database = new Database;
		ZNMolecule *mol = _data ->ddwin ->target_molecule;
		_data ->ddwin ->add_database(database);
		DockingThread *thread = new DockingThread (0, mol, database, _data ->ddwin) ;
		thread ->set_bindingsite_radius (_r);
		thread ->set_bindingsite_center (vect (_x, _y, _z));
		thread ->set_results_number (10);
		thread ->set_time_limit (60);
		_data ->ddwin ->run_thread (thread);
	//	thread ->start ();
	
}


bool DockingMenu::display_requirements_met () {
	bool b = check_for_a_set_target ();
	if (b) {
		ZNMolecule *mol = _data ->ddwin ->target_molecule;
		vect c = get_center (mol);
		_site_x_el ->set_value (c.x ());
		_site_y_el ->set_value (c.y ());
		_site_z_el ->set_value (c.z ());
	}
	return b;
};

ThreadWidget::ThreadWidget (QWidget *parent) : ZNWidget (parent) {
	_name = new QLabel ("");
	layout () -> addWidget (_name);
	_progress_bar = new QProgressBar ();
	layout () ->addWidget (_progress_bar);
	_action = new QLabel ("");
	layout () -> addWidget (_action);

}

ThreadMenu::ThreadMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, false, false) {
int MAX_THREADS = 10;
	for (unsigned int i = 0; i < MAX_THREADS; i++) {
		ThreadWidget *widget = new ThreadWidget ();
		widgets.push_back (widget);
		_main_widget ->layout () -> addWidget (widget);
		widget ->hide ();
	}
}

void ThreadMenu::clear (int n) {
	if (n < widgets.size ()) {
		for (unsigned int i = n; i < widgets.size (); i++) {
			widgets [i] -> hide ();
		}
	}
}

void ThreadMenu::display_thread (int n, Thread *thread) {
	if (n < widgets.size ()) {
		widgets[n] ->show ();
		widgets[n] -> _name ->setText (tr(thread ->_name.c_str ()));
		widgets[n] -> _progress_bar ->setRange (thread ->_min_step, thread ->_max_step);
		widgets[n] -> _progress_bar ->setValue (thread ->_current_step);
		widgets[n] -> _action ->setText (tr(thread ->_doing.c_str ()));	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Preferences menu
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

PrefMenu::PrefMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, true, true) {
	setWindowTitle ("Preferences");


	ZNAdvancedWidget *main_tab = new ZNAdvancedWidget ("", false);
	ZNAdvancedWidget *rendering_tab = new ZNAdvancedWidget ("", false);	
	_tabs ->addTab (main_tab, "External programs");
	_tabs ->addTab (rendering_tab, "rendering");
	

//	_input_box = new MyGroupBox (main_tab ->basic_layout (), "Input options");
//	_input_box -> setMaximumHeight( 250 );

	label = new QLabel( "Set the absolute path for each executable.", this ); 
	main_tab ->basic_layout () -> addWidget (label);

	plants_file = new MyLineFile (main_tab ->basic_layout (), "PLANTS", 4);

//	gamess_file = new MyLineFile (main_tab ->basic_layout (), "GAMESS", 1);


	_save_preferences_b = new MyPushButton (main_tab ->basic_layout (), "Apply and Save Zodiac.ini", 0, 200);
	connect (_save_preferences_b ->fbutton, SIGNAL( clicked() ), SLOT( _save_preferences () ) );

// backgorund color
//	ZNAdvancedWidget *background_tab = new ZNAdvancedWidget ("", false);
//	_tabs ->addTab (background_tab, "Background color");
	MyFloatEditLine *qsel = new MyFloatEditLine (rendering_tab ->basic_layout (),"Quality scale", *_data ->quality_scale, 0, 1000);


	new MyFloatEditLine (rendering_tab ->basic_layout (), "Sphere radius", *_data ->sphere_radius);
	new MyFloatEditLine (rendering_tab ->basic_layout (), "VdW scale", *_data ->vdw_scale);
	new MyFloatEditLine (rendering_tab ->basic_layout (), "Stick radius", *_data ->stick_radius);
	new MyFloatEditLine (rendering_tab ->basic_layout (), "Double Bond inter Distance", *_data ->double_bond_separation);
//	new MyFloatEditLine (rendering_tab ->basic_layout (), "Aromatic Bond inter Distance", _data ->ddwin ->gl->aromatic_bond_inter_distance);
	new MyFloatEditLine (rendering_tab ->basic_layout (), "Double bond scale", *_data ->double_bond_stick_scale);
	new MyFloatEditLine (rendering_tab ->basic_layout (), "Line Width", *_data ->line_width);
	
	new MyFloatEditLine (rendering_tab ->basic_layout (), "Helix tube a", *_data ->backbone_tube_helix_a);
	new MyFloatEditLine (rendering_tab ->basic_layout (), "Helix tube b", *_data ->backbone_tube_helix_b);
	new MyFloatEditLine (rendering_tab ->basic_layout (), "Helix tube c", *_data ->backbone_tube_helix_c);
		new MyFloatEditLine (rendering_tab ->basic_layout (), "Sheet tube a", *_data ->backbone_tube_sheet_a);
		new MyFloatEditLine (rendering_tab ->basic_layout (), "Sheet tube b", *_data ->backbone_tube_sheet_b);
		new MyFloatEditLine (rendering_tab ->basic_layout (), "Sheet tube c", *_data ->backbone_tube_sheet_c);
	new MyFloatEditLine (rendering_tab ->basic_layout (), "Random tube a", *_data ->backbone_tube_random_a);
		new MyFloatEditLine (rendering_tab ->basic_layout (), "Random tube b", *_data ->backbone_tube_random_b);
		new MyFloatEditLine (rendering_tab ->basic_layout (), "Random tube c", *_data ->backbone_tube_random_c);
	
	
	//rendering_tab ->basic_layout ()->addWidget (qsel);

}


void PrefMenu::_save_preferences () {
		cerr<< "Save Zodiac.ini"<<endl;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	PLANTS menu FE
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


PLANTSMenu::PLANTSMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, true, true) {
	setWindowTitle ("PLANTS 1.1 input");

	QString s = "Zodiac.ini";

	if (!s.isNull()) {
		ifstream ifs (s.toStdString ().c_str ());
		string buffer;

		while (getline(ifs, buffer)) {
			
			istringstream line2(buffer);
			string token;
			line2 >> token;

			if (token == "PLANTS") {
				string s = token;
				while (!line2.eof ()) {
					string q;
					line2 >> q;
					s.append (" ");s.append (q);
		
					plants_exe = q.c_str ();
					if (plants_exe != "") runPlants = true;
				}
			}
		}	
	}

// Main tab
	ZNAdvancedWidget *main_tab = new ZNAdvancedWidget ("", false);
   	_tabs ->addTab (main_tab, "Input");

	_input_box = new MyGroupBox (main_tab ->basic_layout (), "Input options");
	_input_box -> setMaximumHeight( 250 );

//	_input_tc = new MyTwoColumn (_input_box ->layout () );


	load_protein_file = new MyLineFile (_input_box ->layout (), "Load protein file", 1);
	connect(load_protein_file->linedit, SIGNAL (textChanged(const QString &)), SLOT (load_protein () ));

	load_ligand_file = new MyLineFile (_input_box ->layout (), "Load ligand file", 1);
	connect(load_ligand_file->linedit, SIGNAL (textChanged(const QString &)), SLOT (load_ligand () ));

	input_file = new MyListView (_input_box ->layout (), 74, 0, false);
input_file ->hide();

	load_ligand_list_file = new MyLineFile (_input_box ->layout (), "Load ligand list file", 3);
	connect(load_ligand_list_file->linedit, SIGNAL (textChanged(const QString &)), SLOT (load_ligand_list () ));

	input_file_list = new MyListView (_input_box ->layout (), 74, 0, false);
input_file_list ->hide();

	_binding_site_box = new MyGroupBox (_input_box ->layout (), "Binding site options");

	_binding_site_from_ligand_b = new MyPushButton (_binding_site_box ->layout (), "Add binding site from ligand");
	connect (_binding_site_from_ligand_b ->fbutton, SIGNAL( clicked() ), SLOT( set_binding_site() ) );

	_binding_site_tc = new MyTwoColumn (_binding_site_box ->layout () );

	_zero_value = 0;
	_min_radius = 10;
	_site_x_el = new MyFloatEditLine (_binding_site_tc ->left_layout (), "Binding site center X:", _zero_value, -1000, 1000);
	_site_y_el = new MyFloatEditLine (_binding_site_tc ->right_layout (), "Binding site center Y:", _zero_value, -1000, 1000);
	_site_z_el = new MyFloatEditLine (_binding_site_tc ->left_layout (), "Binding site center Z:", _zero_value, -1000, 1000);
	_site_radius_el = new MyFloatEditLine (_binding_site_tc ->right_layout (), "Binding site radius:", _min_radius, 0, 100);
	
// output tab
	ZNAdvancedWidget *output_tab = new ZNAdvancedWidget ();
    	_tabs ->addTab (output_tab, "Output");

	_output_box = new MyGroupBox (output_tab ->basic_layout (), "Output options");

	_output_dir = new MyLineEdit (_output_box ->layout (), "Output dir");
//_output_dir ->linedit ->setPalette(QColor(255,0,0));

	_output_tc = new My3Column (output_tab ->advanced_layout () );

	_write_protein_conformations_value = false;
	_write_protein_conformations_chb = new MyCheckBox (_output_tc ->left_layout (), _write_protein_conformations_value, "Write protein conformations");

	_write_protein_bindingsite_value = false;
	_write_protein_bindingsite_chb = new MyCheckBox (_output_tc ->center_layout (), _write_protein_bindingsite_value, "Write protein binding site");

	_write_protein_splitted_value = false;
	_write_protein_splitted_chb = new MyCheckBox (_output_tc ->right_layout (), _write_protein_splitted_value, "Write protein splitted");

	_write_rescored_structures_value = false;
	_write_rescored_structures_chb = new MyCheckBox (_output_tc ->left_layout (), _write_rescored_structures_value, "Write rescored structures");

	_write_multi_mol2_value = true;
	_write_multi_mol2_chb = new MyCheckBox (_output_tc ->center_layout (), _write_multi_mol2_value, "Write multi mol2");

	_write_ranking_links_value = false;
	_write_ranking_links_chb = new MyCheckBox (_output_tc ->left_layout (), _write_ranking_links_value, "Write ranking links");

	_write_ranking_multi_mol2_value = false;
	_write_ranking_multi_mol2_chb = new MyCheckBox (_output_tc ->right_layout (), _write_ranking_multi_mol2_value, "Write ranking multi mol2");

	_write_per_atom_scores_value = false;
	_write_per_atom_scores_chb = new MyCheckBox (_output_tc ->center_layout (), _write_per_atom_scores_value, "Write score per atoms");

	_write_merged_ligand_value = false;
	_write_merged_ligand_chb = new MyCheckBox (_output_tc ->left_layout (), _write_merged_ligand_value, "Write merged ligand");

	_write_merged_protein_value = false;
	_write_merged_protein_chb = new MyCheckBox (_output_tc ->right_layout (), _write_merged_protein_value, "Write merged protein");

	_write_merged_water_value = false;
	_write_merged_water_chb = new MyCheckBox (_output_tc ->center_layout (), _write_merged_water_value, "Write merged water");

	_keep_original_mol2_description_value = true;
	_keep_original_mol2_description_chb = new MyCheckBox (_output_tc ->right_layout (), _keep_original_mol2_description_value, "Keep original mol2 description");

	_merge_multi_conf_output_box = new MyGroupBox (output_tab ->advanced_layout (), "Merge multiconformer output", true, false);

	_merge_multi_conf_output_tc = new MyTwoColumn (_merge_multi_conf_output_box ->layout ()); 

	_merge_multi_conf_character_value = "_";
	_merge_multi_conf_character = new MyLineEdit (_merge_multi_conf_output_tc ->left_layout (), "Character used for merging:");
//	_merge_multi_conf_character ->linedit -> setMaximumWidth( 50 );
	_merge_multi_conf_character ->linedit -> setText ("_");


	_merge_multi_conf_after_characters_value = 1;
	_merge_multi_conf_after_characters_el = new MyIntegerEditLine (_merge_multi_conf_output_tc ->right_layout (), "Xxx:", _merge_multi_conf_after_characters_value, 1, 1000);


// Algorithm tab
	ZNAdvancedWidget *_algorithm_tab = new ZNAdvancedWidget ("Show scoring function options");
    	_tabs ->addTab (_algorithm_tab, "Algorithm");

	_search_box = new MyGroupBox (_algorithm_tab ->basic_layout (), "Search algorithm");

	_search_tc = new My3Column (_search_box ->layout () );

	_search_speed_cb = new MyComboBox (_search_tc ->left_layout (), "Search speed setting:");
	_search_speed_cb ->combo_box () ->insertItem (0, "normal (x1)", "speed1");
	_search_speed_cb ->combo_box () ->insertItem (1, "fast (x2)", "speed2");
	_search_speed_cb ->combo_box () ->insertItem (2, "very fast (x4)", "speed4");

	_rescore_mode_cb = new MyComboBox (_search_tc ->right_layout (), "Rescore mode:");
	_rescore_mode_cb ->combo_box () ->insertItem (0, "simplex", "simplex");
	_rescore_mode_cb ->combo_box () ->insertItem (1, "none", "no_simplex");

	_aco_ants_value = 20;
	_aco_ants_el = new MyIntegerEditLine (_search_tc ->center_layout (), "Number of ants:", _aco_ants_value, 1);

	_search_adv_box = new MyGroupBox (_search_box ->layout (), "Advanced options", true, false);
	connect (_search_adv_box, SIGNAL (clicked (bool)), this, SLOT (hide_group_box () ) );

	_search_adv_tc = new MyTwoColumn (_search_adv_box ->layout () );
	_search_adv_tc ->hide();

	_aco_evap_value = 1;
	_aco_evap_el = new MyFloatEditLine (_search_adv_tc ->left_layout (), "Evaporation factor:", _aco_evap_value, 0, 1);

	_aco_sigma_value = 1;
	_aco_sigma_el = new MyFloatEditLine (_search_adv_tc ->left_layout (), "Sigma scaling factor:", _aco_sigma_value, 0);

	_flip_amide_bonds_value = false;
	_flip_amide_bonds_chb = new MyCheckBox (_search_adv_tc ->left_layout (), _flip_amide_bonds_value, "Flip amide bonds");

	_flip_planar_n_value = false;
	_flip_planar_n_chb = new MyCheckBox (_search_adv_tc ->right_layout (), _flip_planar_n_value, "Flip planar N");

	_force_flipped_bonds_planarity_value = false;
	_force_flipped_bonds_planarity_chb = new MyCheckBox (_search_adv_tc ->right_layout (), _force_flipped_bonds_planarity_value, "Force flipped bond planary");

	_force_planar_bond_rotation_value = true;
	_force_planar_bond_rotation_chb = new MyCheckBox (_search_adv_tc ->right_layout (), _force_planar_bond_rotation_value, "Force planar bond rotation");


	_cluster_docking_tc = new MyTwoColumn (_algorithm_tab ->basic_layout () );

	_cluster_box = new MyGroupBox (_cluster_docking_tc ->left_layout (), "Cluster algorithm");

	_cluster_rmsd_value = 2.0;
	_cluster_rmsd_el = new MyFloatEditLine (_cluster_box ->layout (), "Cluster RMSD:", _cluster_rmsd_value, 0);

	_cluster_structures_value = 10;
	_cluster_structures_el = new MyIntegerEditLine (_cluster_box ->layout (), "Cluster structures:", _cluster_structures_value, 1);


	_docking_box = new MyGroupBox (_cluster_docking_tc ->right_layout (), "Docking typology");

	_rigid_ligand_value = false;
	_rigid_ligand_chb = new MyCheckBox (_docking_box ->layout (), _rigid_ligand_value, "Rigid ligand docking");

	_rigid_all_value = false;
	_rigid_all_chb = new MyCheckBox (_docking_box ->layout (), _rigid_all_value, "Rigid protein and ligand docking");


	_algorithm_type_cb = new MyComboBox (_algorithm_tab ->advanced_layout (), "Scoring function type:");
	_algorithm_type_cb ->combo_box () ->insertItem (0, "chemplp", "chemplp");
	_algorithm_type_cb ->combo_box () ->insertItem (1, "plp", "plp");
	_algorithm_type_cb ->combo_box () ->insertItem (2, "plp95", "plp95");
	connect (_algorithm_type_cb ->combo_box (), SIGNAL (currentIndexChanged ( int )), this, SLOT (plp_weight ( int )) );

	_scoring_inter_box = new MyGroupBox (_algorithm_tab ->advanced_layout (), "Intermolecular scoring");
	_score_inter_tc = new MyTwoColumn (_scoring_inter_box ->layout () );

	_outside_binding_site_penalty_value = 50.0;
	_outside_binding_site_penalty_el = new MyFloatEditLine (_score_inter_tc ->left_layout (), "Outside binding site penalty:", _outside_binding_site_penalty_value, 0, 100);

	_enable_sulphur_acceptors_value = false;
	_enable_sulphur_acceptors_chb = new MyCheckBox (_score_inter_tc ->right_layout (), _enable_sulphur_acceptors_value, "Enable sulphur acceptors");


	_scoring_intra_box = new MyGroupBox (_algorithm_tab ->advanced_layout (), "Intramolecular scoring");
	_score_intra_tc = new My3Column (_scoring_intra_box ->layout () );

	_ligand_intra_score_cb = new MyComboBox (_score_intra_tc ->left_layout (), "Ligand scoring:");
	_ligand_intra_score_cb ->combo_box () ->insertItem (0, "clash2", "clash2");
	_ligand_intra_score_cb ->combo_box () ->insertItem (1, "clash", "clash");
	_ligand_intra_score_cb ->combo_box () ->insertItem (2, "all atom Lennard-Jones", "lj");

	_chemplp_clash_include_14_value = 0.25;
	_chemplp_clash_include_14_el = new MyFloatEditLine (_score_intra_tc ->center_layout (), "Scoring 1-4 interaction:", _chemplp_clash_include_14_value, 0, 1);

	_chemplp_clash_include_HH_value = false;
	_chemplp_clash_include_HH_chb = new MyCheckBox (_score_intra_tc ->right_layout (), _chemplp_clash_include_HH_value, "Scoring H-H interaction");


	_chemplp_weights_box = new MyGroupBox (_scoring_intra_box ->layout (), "Chemplp interaction weights");


	_chemplp_weights_tc = new My3Column (_chemplp_weights_box ->layout () );

	_chemplp_charged_hb_weight_value = 2.0;
	_chemplp_charged_hb_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->left_layout (), "Charged HB:", _chemplp_charged_hb_weight_value, -100);

	_chemplp_hbond_cho_weight_value = -3.0;
	_chemplp_hbond_cho_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->center_layout (), "HB CHO:", _chemplp_hbond_cho_weight_value, -100);

	_chemplp_hbond_weight_value = -3.0;
	_chemplp_hbond_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->right_layout (), "Neutral HB:", _chemplp_hbond_weight_value, -100);

	_chemplp_metal_weight_value = -3.0;
	_chemplp_metal_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->left_layout (), "Metal:", _chemplp_metal_weight_value, -100);

	_chemplp_charged_metal_weight_value = 2.0;
	_chemplp_charged_metal_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->center_layout (), "Charged metal:", _chemplp_charged_metal_weight_value, -100);

	_chemplp_plp_weight_value = 1.0;
	_chemplp_plp_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->left_layout (), "PLP:", _chemplp_plp_weight_value, -100);

	_chemplp_plp_steric_e_value = -0.4;
	_chemplp_plp_steric_e_el = new MyFloatEditLine (_chemplp_weights_tc ->center_layout (), "PLP steric:", _chemplp_plp_steric_e_value, -100);

	_chemplp_plp_burpolar_e_value = -0.1;
	_chemplp_plp_burpolar_e_el = new MyFloatEditLine (_chemplp_weights_tc ->right_layout (), "Polar PLP:", _chemplp_plp_burpolar_e_value, -100);

	_chemplp_plp_hbond_e_value = -1.0;
	_chemplp_plp_hbond_e_el = new MyFloatEditLine (_chemplp_weights_tc ->left_layout (), "HB PLP:", _chemplp_plp_hbond_e_value, -100);

	_chemplp_plp_metal_e_value = -1.0;
	_chemplp_plp_metal_e_el = new MyFloatEditLine (_chemplp_weights_tc ->center_layout (), "Metal PLP:", _chemplp_plp_metal_e_value, -100);

	_chemplp_plp_repulsive_weight_value = 1.0;
	_chemplp_plp_repulsive_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->right_layout (), "Repulsive PLP:", _chemplp_plp_repulsive_weight_value, -100);

	_chemplp_lipo_weight_value = 0.0;
	_chemplp_lipo_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->left_layout (), "Lipophilic:", _chemplp_lipo_weight_value, -100);

	_chemplp_tors_weight_value = 2.0;
	_chemplp_tors_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->center_layout (), "Ligand torsional:", _chemplp_tors_weight_value, -100);

	_chemplp_intercept_weight_value = -20.0;
	_chemplp_intercept_weight_el = new MyFloatEditLine (_chemplp_weights_tc ->right_layout (), "Intercept:", _chemplp_intercept_weight_value, -100);

	_chemplp_weak_cho_value = true;
	_chemplp_weak_cho_chb = new MyCheckBox (_chemplp_weights_tc ->right_layout (), _chemplp_weak_cho_value, "Weak CHO");


	_plp_weights_box = new MyGroupBox (_scoring_intra_box ->layout (), "Plp/plp95 interaction weights");
	_plp_weights_box ->hide ();

	_plp_weights_tc = new My3Column (_plp_weights_box ->layout () );

	_plp_steric_e_value = -0.4;
	_plp_steric_e_el = new MyFloatEditLine (_plp_weights_tc ->left_layout (), "Steric:", _plp_steric_e_value, -100);

	_plp_burpolar_e_value = -0.05;
	_plp_burpolar_e_el = new MyFloatEditLine (_plp_weights_tc ->left_layout (), "Polar:", _plp_burpolar_e_value, -100);

	_plp_hbond_e_value = -2.0;
	_plp_hbond_e_el = new MyFloatEditLine (_plp_weights_tc ->center_layout (), "HB:", _plp_hbond_e_value, -100);

	_plp_metal_e_value = -4.0;
	_plp_metal_e_el = new MyFloatEditLine (_plp_weights_tc ->center_layout (), "Metal:", _plp_metal_e_value, -100);

	_plp_repulsive_weight_value = 0.5;
	_plp_repulsive_weight_el = new MyFloatEditLine (_plp_weights_tc ->right_layout (), "Repulsive:", _plp_repulsive_weight_value, -100);

	_plp_tors_weight_value = 1.0;
	_plp_tors_weight_el = new MyFloatEditLine (_plp_weights_tc ->right_layout (), "Ligand torsional:", _plp_tors_weight_value, -100);


// Costraints tab
	ZNAdvancedWidget *_costraints_tab = new ZNAdvancedWidget ("", false);
    	_tabs ->addTab (_costraints_tab, "Costraints");

	_costraints_box = new MyGroupBox (_costraints_tab ->basic_layout (), "Costraints");
//	_costraints_box -> setMaximumHeight( 361 );

	_costraints_tc = new MyTwoColumn (_costraints_box ->layout () );


	_protein_hb_costraint_box = new MyGroupBox (_costraints_tc ->right_layout (), "Protein HB costraints");

	_protein_hb_costraint_tc = new MyTwoColumn (_protein_hb_costraint_box ->layout ());

	_protein_hb_costraint_cb = new MyComboBox (_protein_hb_costraint_tc ->left_layout (),  "Atom number");
	_protein_hb_costraint_fel = new MyFloatEditLine (_protein_hb_costraint_tc ->right_layout (), "Weight", _zero_value, 0);
	_protein_hb_costraint_b = new MyPushButton (_protein_hb_costraint_box ->layout (), "Set constraint");
	connect (_protein_hb_costraint_b ->fbutton, SIGNAL( clicked() ), SLOT( set_protein_hb_constraint_slot() ) );


	_shape_costraint_box = new MyGroupBox (_costraints_tc ->right_layout (), "Shape constraints");

	load_shape_file = new MyLineFile (_shape_costraint_box ->layout (), "Load shape file", 1);
	connect(load_shape_file->linedit, SIGNAL (textChanged(const QString &)), SLOT (load_shape () ));

	_shape_costraint_tc = new MyTwoColumn (_shape_costraint_box ->layout ());

	_shape_weight = new MyFloatEditLine (_shape_costraint_tc ->left_layout (), "Weight: ", _zero_value, -10, 0);

	_shape_weight_b = new MyPushButton (_shape_costraint_tc ->right_layout (), "Set constraint");
	connect (_shape_weight_b ->fbutton, SIGNAL( clicked() ), SLOT( set_shape_constraint_slot() ) );


	_surface_distance_costraint_box = new MyGroupBox (_costraints_tc ->left_layout (), "Surface distance constraints");

	_surface_distance_costraint_tc = new MyTwoColumn (_surface_distance_costraint_box ->layout ());

	_from_surface_distance = new MyFloatEditLine (_surface_distance_costraint_tc ->left_layout (), "From: ", _zero_value, -10, 10);
	_to_surface_distance = new MyFloatEditLine (_surface_distance_costraint_tc ->right_layout (), "To: ", _zero_value, -10, 10);

	_ligand_surface_distance_cb = new MyComboBox (_surface_distance_costraint_tc ->left_layout (),  "Ligand atom number");

	_weight_surface_distance_fel = new MyFloatEditLine (_surface_distance_costraint_tc ->right_layout (), "Weight: ", _zero_value, -10, 0);

	_surface_distance_b = new MyPushButton (_surface_distance_costraint_box ->layout (), "Set constraint");
	connect (_surface_distance_b ->fbutton, SIGNAL( clicked() ), SLOT( set_surface_distance_slot() ) );


	_ligand_intra_distance_costraint_box = new MyGroupBox (_costraints_tc ->left_layout (), "Ligand intra distance constraints");

	_ligand_intra_distance_costraint_tc = new MyTwoColumn (_ligand_intra_distance_costraint_box ->layout ());

	_from_ligand_intra_distance = new MyFloatEditLine (_ligand_intra_distance_costraint_tc ->left_layout (), "From: ", _zero_value, -10, 10);
	_to_ligand_intra_distance = new MyFloatEditLine (_ligand_intra_distance_costraint_tc ->right_layout (), "To: ", _zero_value, -10, 10);

	_number_a_ligand_intra_distance_cb = new MyComboBox (_ligand_intra_distance_costraint_tc ->left_layout (),  "Ligand atom number");

	_number_b_ligand_intra_distance_cb = new MyComboBox (_ligand_intra_distance_costraint_tc ->right_layout (),  "Ligand atom number");

	_weight_ligand_intra_distance = new MyFloatEditLine (_ligand_intra_distance_costraint_tc ->left_layout (), "Weight: ", _zero_value, -10, 0);

	_ligand_intra_distance_b = new MyPushButton (_ligand_intra_distance_costraint_tc ->right_layout (), "Set constraint");
	connect (_ligand_intra_distance_b ->fbutton, SIGNAL( clicked() ), SLOT( set_ligand_intra_distance_slot () ) );


	_protein_ligand_distance_costraint_box = new MyGroupBox (_costraints_tc ->left_layout (), "Protein ligand distance constraints");

	_protein_ligand_distance_costraint_tc = new MyTwoColumn (_protein_ligand_distance_costraint_box ->layout ());

	_from_protein_ligand_distance = new MyFloatEditLine (_protein_ligand_distance_costraint_tc ->left_layout (), "From: ", _zero_value, -10, 10);
	_to_protein_ligand_distance = new MyFloatEditLine (_protein_ligand_distance_costraint_tc ->right_layout (), "To: ", _zero_value, -10, 10);

	_number_a_protein_ligand_distance_cb = new MyComboBox (_protein_ligand_distance_costraint_tc ->left_layout (),  "Protein atom number");

	_number_b_protein_ligand_distance_cb = new MyComboBox (_protein_ligand_distance_costraint_tc ->right_layout (),  "Ligand atom number");

	_weight_protein_ligand_distance = new MyFloatEditLine (_protein_ligand_distance_costraint_tc ->left_layout (), "Weight: ", _zero_value, -10, 10);

	_protein_ligand_distance_b = new MyPushButton (_protein_ligand_distance_costraint_tc ->right_layout (), "Set constraint");
	connect (_protein_ligand_distance_b ->fbutton, SIGNAL( clicked() ), SLOT( set_protein_ligand_distance_slot () ) );

	costraints_file_list = new MyListView (_costraints_tc ->right_layout (), 138 );


// Flexibility tab
	ZNAdvancedWidget *_flexibility_tab = new ZNAdvancedWidget ("", false);
    	_tabs ->addTab (_flexibility_tab, "Flexibility");

	_flexibility_box = new MyGroupBox (_flexibility_tab ->basic_layout (), "Flexibility");
//	_flexibility_box -> setMaximumHeight( 145 );

	_flexibility_tc = new MyTwoColumn (_flexibility_box ->layout () );
	_flexibility_tc ->_left -> setMinimumWidth ( 350 );

	_flexible_protein_side_chain_string = new MyListButton (_flexibility_tc ->left_layout (), "Add flexible side chain");
	connect (_flexible_protein_side_chain_string->fbutton, SIGNAL( clicked() ), SLOT( add_fscs_slot() ) );

	_fix_protein_bond = new MyListButton (_flexibility_tc ->left_layout (), "Add fixed protein bond");
	connect (_fix_protein_bond->fbutton, SIGNAL( clicked() ), SLOT( add_fpb_slot() ) );

	_flexible_protein_side_chain_number = new MyListButton (_flexibility_tc ->left_layout (), "Add flexible side chain");
	connect (_flexible_protein_side_chain_number->fbutton, SIGNAL( clicked() ), SLOT( add_fscn_slot() ) );

	_intra_protein_score_weight_value = 0.3;
	_intra_protein_score_weight_el = new MyFloatEditLine (_flexibility_tc ->left_layout (), "Intra protein score weight:", _intra_protein_score_weight_value, 0, 1);


	flexibility_file_list = new MyListView (_flexibility_tc ->right_layout (), 85 );


// Water tab
	ZNAdvancedWidget *_water_tab = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (_water_tab, "Water");

	_water_box = new MyGroupBox (_water_tab ->basic_layout (), "Water");

	_water_site_b = new MyPushButton (_water_box ->layout (), "Add water molecule");
	connect (_water_site_b ->fbutton, SIGNAL( clicked() ), SLOT( set_water_site() ) );

	_water_site_tc = new MyTwoColumn (_water_box ->layout () );

	_water_min_radius = 2;
	_water_site_x_el = new MyFloatEditLine (_water_site_tc ->left_layout (), "Water X:", _zero_value, -1000, 1000);
	_water_site_y_el = new MyFloatEditLine (_water_site_tc ->right_layout (), "Water Y:", _zero_value, -1000, 1000);
	_water_site_z_el = new MyFloatEditLine (_water_site_tc ->left_layout (), "Water Z:", _zero_value, -1000, 1000);
	_water_site_radius_el = new MyFloatEditLine (_water_site_tc ->right_layout (), "Radius:", _water_min_radius, 0, 10);


	load_water_file = new MyLineFile (_water_box ->layout (), "Water reference", 1);
//	connect(load_water_file->linedit, SIGNAL (textChanged(const QString &)), SLOT (load_water () ));


	_water_weights_box = new MyGroupBox (_water_box ->layout (), "Water molecule weights", true, false);
	connect (_water_weights_box, SIGNAL (clicked (bool)), this, SLOT (hide_group_box2 () ) );

	_water_weights_tc = new MyTwoColumn (_water_weights_box ->layout () );
	_water_weights_tc -> hide();

	_no_water_ligand_hb_penalty_value = 0.0;
	_no_water_ligand_hb_penalty_el = new MyFloatEditLine (_water_weights_tc ->left_layout (), "No HB penalty:", _no_water_ligand_hb_penalty_value, -100);

	_water_enable_penalty_value = 8.0;
	_water_enable_penalty_el = new MyFloatEditLine (_water_weights_tc ->right_layout (), "No water penalty:", _water_enable_penalty_value, -100);

	_water_protein_hb_weight_value = 1.0;
	_water_protein_hb_weight_el = new MyFloatEditLine (_water_weights_tc ->left_layout (), "Water-protein HB:", _water_protein_hb_weight_value, -100);

	_water_ligand_hb_weight_value = 1.0;
	_water_ligand_hb_weight_el = new MyFloatEditLine (_water_weights_tc ->right_layout (), "Water-ligand HB:", _water_ligand_hb_weight_value, -100);

	_water_water_hb_weight_value = 1.0;
	_water_water_hb_weight_el = new MyFloatEditLine (_water_weights_tc ->left_layout (), "Water-water HB:", _water_water_hb_weight_value, -100);


//
	add_menu ();
	add_help ();
}


void PLANTSMenu::update_boxes (MyListView *parent, const char *mol_name, int type){
//	parent ->_lw ->insertItem(type, tr( mol_name ));
	parent ->_lw ->addItem(tr( mol_name ));

}


////////////////////////////////////////////////////////////////////////////////
// not very good
// By Tosh
void PLANTSMenu::hide_group_box () {
	if (_search_adv_box ->isChecked () == true) {
		_search_adv_tc ->show();
	}
	else {
		_search_adv_tc ->hide();
	}
}


void PLANTSMenu::hide_group_box2 () {
	if (_water_weights_box ->isChecked () == true) {
		_water_weights_tc ->show();
	}
	else {
		_water_weights_tc ->hide();
	}
}
//
//////////////////////////////////////////////////////////////////////////////////////


void PLANTSMenu::load_protein () {
	string p_name = load_protein_file->val ();

	load_protein_file ->control_no -> hide();
	load_protein_file ->control_yes -> hide();

	if (p_name.find (".mol2")!=string::npos) {

		ifstream ifs(p_name.c_str ());
		OBConversion conv(&ifs);
		OBFormat* inFormat = conv.FormatFromExt(p_name.c_str ());

		OBMol *mol = new OBMol ();

		if(conv.SetInFormat(inFormat) && conv.Read(mol)) {
			OBMol *mol2 = new OBMol ();

			if (conv.Read(mol2)) {
				load_protein_file ->linedit ->clear ();
				load_protein_file ->control_no ->show();
			}
			else {
				load_protein_file ->control_yes ->show();
				char chain;
				int res_number = 0;
				int atom = 0;
				int bond = 0;
				string res = "";
				string residue_plus = "";
				FOR_BONDS_OF_MOL(b, mol)
				{
					bond += b;
					stringstream ss;
					ss << bond;
					string bond_number = ss.str();
					_fix_protein_bond ->_combo_box ->addItem (bond_number.c_str(), bond_number.c_str());
				}
				FOR_ATOMS_OF_MOL(a, mol)
				{
					atom += a;
					stringstream ss;
					ss << atom;
					string atom_number = ss.str();
					_protein_hb_costraint_cb ->_combo_box ->addItem (atom_number.c_str(), atom_number.c_str());
					_number_a_protein_ligand_distance_cb ->_combo_box ->addItem (atom_number.c_str(), atom_number.c_str());
				}
				FOR_RESIDUES_OF_MOL(r, mol)
				{
					res_number =  r->GetNum();
					res =  r->GetName();
					chain =  r->GetChain();
					QString residue = res.c_str();
					stringstream ss;
					ss << res_number;
					string residue_number = ss.str();
					_flexible_protein_side_chain_string ->_combo_box ->addItem (residue, residue);
					_flexible_protein_side_chain_number ->_combo_box ->addItem (residue_number.c_str(), residue_number.c_str());
				}

			}
		}
        	else {
			load_protein_file ->control_no -> show();
			delete mol;
		}
	}
	else {
//		_err_mol ();
	}
}


void PLANTSMenu::load_ligand () {
	string l_name = load_ligand_file->val ();

	load_ligand_file ->control_no -> hide();
	load_ligand_file ->control_yes -> hide();

	if (l_name.find (".mol2")!=string::npos) {

		ifstream ifs(l_name.c_str ());
		OBConversion conv(&ifs);
		OBFormat* inFormat = conv.FormatFromExt(l_name.c_str ());

		OBMol *mol = new OBMol ();

		if(conv.SetInFormat(inFormat) && conv.Read(mol)) {
			OBMol *mol2 = new OBMol ();

			const char *title = "";

			if (conv.Read(mol2)) {
				load_ligand_file ->control_yes -> show();

				title =  mol->GetTitle();
				update_boxes (input_file, title, 2);
			}
			else {
				load_ligand_file ->control_yes -> show();
				int atom = 0;
				FOR_ATOMS_OF_MOL(a, mol)
				{
					atom += a;
					stringstream ss;
					ss << atom;
					string atom_number = ss.str();
					_ligand_surface_distance_cb ->_combo_box ->addItem (atom_number.c_str(), atom_number.c_str());
					_number_a_ligand_intra_distance_cb ->_combo_box ->addItem (atom_number.c_str(), atom_number.c_str());
					_number_b_ligand_intra_distance_cb ->_combo_box ->addItem (atom_number.c_str(), atom_number.c_str());
					_number_b_protein_ligand_distance_cb ->_combo_box ->addItem (atom_number.c_str(), atom_number.c_str());
				}
				title =  mol->GetTitle();
				update_boxes (input_file, title, 2);
			}

		}
        	else {
			load_ligand_file ->control_no -> hide();
			delete mol;
		}
	}
	else {
//		_err_mol ();
//		load_protein_file ->linedit ->setBackgroundRole (QPalette::Highlight);
	}
}


void PLANTSMenu::load_shape () {
	string s_name = load_shape_file->val ();

	if (s_name.find (".mol2")!=string::npos) {

		ifstream ifs(s_name.c_str ());
		OBConversion conv(&ifs);
		OBFormat* inFormat = conv.FormatFromExt(s_name.c_str ());

		OBMol *mol = new OBMol ();

		if(conv.SetInFormat(inFormat) && conv.Read(mol)) {
			OBMol *mol2 = new OBMol ();

			if (conv.Read(mol2)) {
				_err_multi ();
			}
			else {
//				const char *title_shape = "";
				title_shape =  mol->GetTitle();
cerr << "i have read this file "<<title_shape<<endl;
			}
		}
        	else {
			_err_mol ();
			delete mol;
		}
	}
	else {
		_err_mol ();
//		load_shape_file ->linedit ->setBackgroundRole (QPalette::Highlight);
	}
}


void PLANTSMenu::set_shape_constraint_slot() {
	string mol = load_shape_file->val ();
	stringstream ss3;
	ss3 << _shape_weight ->spinbox -> value ();
	string weight = ss3.str();
	if (mol != ""){
		string st = "shape_costraint "+mol+" "+weight;
		update_boxes (costraints_file_list, st.c_str(), 2);
	}
}


void PLANTSMenu::load_ligand_list (){
	string ss = load_ligand_list_file->val ();
//	QString s = QFileDialog::getOpenFileName(this, 
//		tr ("Open file"), "",tr("PLANTS config File (*.pcfg);;All files (*)"));

	QString s(ss.c_str ());

	if (!s.isNull()) {
		ifstream ifs (s.toStdString ().c_str ());
		string buffer;

		input_file_list ->_lw ->clear();

		while (getline(ifs, buffer)) {
			istringstream line2(buffer);
			string token;
			line2 >> token;

			update_boxes (input_file_list, token.c_str(), 3);
		}
	}
}


void PLANTSMenu::set_binding_site () {
	QString s = QFileDialog::getOpenFileName(this, tr ("Open file"), "",tr("Tripos Mol2 File (*.mol2)"));
	string bs_name = s.toStdString ();

	if (bs_name.find (".mol2")!=string::npos) {

		ifstream ifs(bs_name.c_str ());
		OBConversion conv(&ifs);
		OBFormat* inFormat = conv.FormatFromExt(bs_name.c_str ());

		OBMol *mol = new OBMol ();

		if(conv.SetInFormat(inFormat) && conv.Read(mol)) {
			OBMol *mol2 = new OBMol ();

			if (conv.Read(mol2)) {
				_err_multi ();
			}
			else {
				coordinates =  mol->Center(1);
				site_x = coordinates.x ();
				site_y = coordinates.y ();
				site_z = coordinates.z ();

				_site_x_el ->spinbox ->setValue (site_x);
				_site_y_el ->spinbox ->setValue (site_y);
				_site_z_el ->spinbox ->setValue (site_z);
			}
		}
        	else {
			_err_mol ();
			delete mol;
		}
	}
	else {
//		_err_mol ();
//		load_protein_file ->linedit ->setBackgroundRole (QPalette::Highlight);
	}
}


void PLANTSMenu::set_water_site () {
	QString s = QFileDialog::getOpenFileName(this, tr ("Open file"), "",tr("Tripos Mol2 File (*.mol2)"));
	string w_name = s.toStdString ();

	if (w_name.find (".mol2")!=string::npos) {

		ifstream ifs(w_name.c_str ());
		OBConversion conv(&ifs);
		OBFormat* inFormat = conv.FormatFromExt(w_name.c_str ());

		OBMol *mol = new OBMol ();

		if(conv.SetInFormat(inFormat) && conv.Read(mol)) {
			OBMol *mol2 = new OBMol ();

			if (conv.Read(mol2)) {
				_err_multi ();
			}
			else {
				water_coordinates =  mol->Center(1);
				water_site_x = water_coordinates.x ();
				water_site_y = water_coordinates.y ();
				water_site_z = water_coordinates.z ();

				_water_site_x_el ->spinbox ->setValue (water_site_x);
				_water_site_y_el ->spinbox ->setValue (water_site_y);
				_water_site_z_el ->spinbox ->setValue (water_site_z);
			}
		}
        	else {
			_err_mol ();
			delete mol;
		}
	}
	else {
//		_err_mol ();
//		load_protein_file ->linedit ->setBackgroundRole (QPalette::Highlight);
	}
}


void PLANTSMenu::add_fscs_slot() {
	string res = _flexible_protein_side_chain_string->_combo_box->currentText ().toStdString ();
	if (res != ""){
		string st = "flexible_protein_side_chain_string "+res;
		update_boxes (flexibility_file_list, st.c_str(), 1);
	}
}


void PLANTSMenu::add_fscn_slot() {
	string res_number = _flexible_protein_side_chain_number->_combo_box->currentText ().toStdString ();
	if (res_number != ""){
		string st = "flexible_protein_side_chain_number "+res_number;
		update_boxes (flexibility_file_list, st.c_str(), 2);
	}
}

void PLANTSMenu::add_fpb_slot() {
	string bond = _fix_protein_bond->_combo_box->currentText ().toStdString ();
	if (bond != ""){
		string st = "fix_protein_bond "+bond;
		update_boxes (flexibility_file_list, st.c_str(), 3);
	}
}


void PLANTSMenu::set_protein_hb_constraint_slot() {
	string atom = _protein_hb_costraint_cb ->_combo_box ->currentText ().toStdString ();
	stringstream ss;
	ss << _protein_hb_costraint_fel ->spinbox -> value ();
	string weight = ss.str();
	if (atom != ""){
		string st = "chemplp_protein_hb_costraint "+atom+" "+weight;
		update_boxes (costraints_file_list, st.c_str(), 1);
	}
}

void PLANTSMenu::set_surface_distance_slot () {
	string atom = _ligand_surface_distance_cb ->_combo_box ->currentText ().toStdString ();
	stringstream ss;
	ss << _from_surface_distance ->spinbox -> value ();
	string from = ss.str();
	stringstream ss2;
	ss2 << _to_surface_distance ->spinbox -> value ();
	string to = ss2.str();
	stringstream ss3;
	ss3 << _weight_surface_distance_fel ->spinbox -> value ();
	string weight = ss3.str();
	if (atom != ""){
		string st = "surface_distance_costraint "+from+to+weight+"("+atom+")";
		update_boxes (costraints_file_list, st.c_str(), 3);
	}
}

void PLANTSMenu::set_ligand_intra_distance_slot () {
	string atom_a = _number_a_ligand_intra_distance_cb ->_combo_box ->currentText ().toStdString ();
	string atom_b = _number_b_ligand_intra_distance_cb ->_combo_box ->currentText ().toStdString ();
	stringstream ss;
	ss << _from_ligand_intra_distance ->spinbox -> value ();
	string from = ss.str();
	stringstream ss2;
	ss2 << _to_ligand_intra_distance ->spinbox -> value ();
	string to = ss2.str();
	stringstream ss3;
	ss3 << _weight_ligand_intra_distance ->spinbox -> value ();
	string weight = ss3.str();
	if (atom_a != ""){
		string st = "ligand_intra_distance_costraint "+from+" "+to+" "+weight+" "+atom_a+" "+atom_b;
		update_boxes (costraints_file_list, st.c_str(), 4);
	}
}


void PLANTSMenu::set_protein_ligand_distance_slot () {
	string atom_a = _number_a_protein_ligand_distance_cb ->_combo_box ->currentText ().toStdString ();
	string atom_b = _number_b_protein_ligand_distance_cb ->_combo_box ->currentText ().toStdString ();
	stringstream ss;
	ss << _from_protein_ligand_distance ->spinbox -> value ();
	string from = ss.str();
	stringstream ss2;
	ss2 << _to_protein_ligand_distance ->spinbox -> value ();
	string to = ss2.str();
	stringstream ss3;
	ss3 << _weight_protein_ligand_distance ->spinbox -> value ();
	string weight = ss3.str();
	if (atom_a != "" & atom_b != ""){
		string st = "protein_ligand_distance_costraint "+from+" "+to+" "+weight+" "+atom_a+" "+atom_b;
		update_boxes (costraints_file_list, st.c_str(), 5);
	}

}

void PLANTSMenu::plp_weight (int i) {
	if (i == 0) {
		_plp_weights_box ->hide ();
		_chemplp_weights_box ->show ();
	} 
	else {
		_chemplp_weights_box ->hide ();
		_plp_weights_box ->show ();
	} 
}

void PLANTSMenu::add_menu () {
    	QMenu *file = new QMenu(tr("&File"), this );
	_load_action = new QAction (tr ("&Open"), this);
	connect (_load_action, SIGNAL (triggered ()), this, SLOT (_load_slot ()));
    	file -> addAction (_load_action);

	_save_action = new QAction (tr ("&Save"), this);
	connect (_save_action, SIGNAL (triggered ()), this, SLOT (_save_slot ()));
    	file -> addAction (_save_action);

	if (runPlants == true) {
		_run_action = new QAction (tr ("Save and Run"), this);
		connect (_run_action, SIGNAL (triggered ()), this, SLOT (_run_slot ()));
		file -> addAction (_run_action);
	}

	addMenu (file);

    	QMenu *settings = new QMenu(tr("&Settings"), this );
	_restore_action = new QAction (tr ("&Restore to default"), this);
	connect (_restore_action, SIGNAL (triggered ()), this, SLOT (_restore_slot ()));
    	settings -> addAction (_restore_action);
	addMenu (settings);

}


void PLANTSMenu::_err_mol () {
	QMessageBox::about( _data ->ddwin, "PLANTS 1.1 config file front end" ,
			QString(("ZODIAC. Version "+VERSION+ "\n\n"
			"Could not read this file\n").c_str()) );
}

void PLANTSMenu::_err_multi () {
	QMessageBox::about( _data ->ddwin, "PLANTS 1.1 config file front end" ,
			QString(("ZODIAC. Version "+VERSION+ "\n\n"
			"Could not use a multy mol2 file\n").c_str()) );
}

void PLANTSMenu::_show_about () {
	QMessageBox::about( _data ->ddwin, "About PLANTS 1.1 config file front end" ,
			QString(("ZODIAC. Version "+VERSION+ "\n"
			"PLANTS 1.1 config file front end. Code by Alberto Massarotti. All rights reserved.\n"
			"For bug reports, suggestions or any other feedback please visit our website at www.zeden.org\n\n").c_str()) );
}


void PLANTSMenu::_load_slot () {
	QString s = QFileDialog::getOpenFileName(this, tr ("Open file"), "",tr("PLANTS config File (*.pcfg);;All files (*)"));

	if (!s.isNull()) {
		ifstream ifs (s.toStdString ().c_str ());
		string buffer;

		_restore_slot ();

		while (getline(ifs, buffer)) {
//			1.1
			_load_bot_cb (buffer, "search_speed", _search_speed_cb);
			_load_bot_iel (buffer, "aco_ants", _aco_ants_el);
			_load_bot_fel (buffer, "aco_evap", _aco_evap_el);
			_load_bot_fel (buffer, "aco_sigma", _aco_sigma_el);
			_load_bot_chb (buffer, "flip_amide_bonds", _flip_amide_bonds_chb);
			_load_bot_chb (buffer, "flip_planar_n", _flip_planar_n_chb);
			_load_bot_chb (buffer, "force_flipped_bonds_planarity", _force_flipped_bonds_planarity_chb);
			_load_bot_chb (buffer, "force_planar_bond_rotation", _force_planar_bond_rotation_chb);
			_load_bot_cb (buffer, "rescore_mode", _rescore_mode_cb);
//			1.2
			_load_bot_bindingsite (buffer, "bindingsite_center");
			_load_bot_fel (buffer, "bindingsite_radius", _site_radius_el);
//			1.3
			_load_bot_fel (buffer, "cluster_rmsd", _cluster_rmsd_el);
			_load_bot_iel (buffer, "cluster_structures", _cluster_structures_el);
//			1.4
			_load_bot_cb (buffer, "scoring_function", _algorithm_type_cb);
			_load_bot_fel (buffer, "outside_binding_site_penalty", _outside_binding_site_penalty_el);
			_load_bot_chb (buffer, "enable_sulphur_acceptors", _enable_sulphur_acceptors_chb);
			_load_bot_cb (buffer, "ligand_intra_score", _ligand_intra_score_cb);
			_load_bot_fel (buffer, "chemplp_clash_include_14", _chemplp_clash_include_14_el);
			_load_bot_chb (buffer, "chemplp_clash_include_HH", _chemplp_clash_include_HH_chb);
			_load_bot_fel (buffer, "plp_steric_e", _plp_steric_e_el);
			_load_bot_fel (buffer, "plp_burpolar_e", _plp_burpolar_e_el);
			_load_bot_fel (buffer, "plp_hbond_e", _plp_hbond_e_el);
			_load_bot_fel (buffer, "plp_metal_e", _plp_metal_e_el);
			_load_bot_fel (buffer, "plp_repulsive_weight", _plp_repulsive_weight_el);
			_load_bot_fel (buffer, "plp_tors_weight", _plp_tors_weight_el);
			_load_bot_chb (buffer, "chemplp_weak_cho", _chemplp_weak_cho_chb);
			_load_bot_fel (buffer, "chemplp_charged_hb_weight", _chemplp_charged_hb_weight_el);
			_load_bot_fel (buffer, "chemplp_charged_metal_weight", _chemplp_charged_metal_weight_el);
			_load_bot_fel (buffer, "chemplp_hbond_weight", _chemplp_hbond_weight_el);
			_load_bot_fel (buffer, "chemplp_hbond_cho_weight", _chemplp_hbond_cho_weight_el);
			_load_bot_fel (buffer, "chemplp_metal_weight", _chemplp_metal_weight_el);
			_load_bot_fel (buffer, "chemplp_plp_weight", _chemplp_plp_weight_el);
			_load_bot_fel (buffer, "chemplp_plp_steric_e", _chemplp_plp_steric_e_el);
			_load_bot_fel (buffer, "chemplp_plp_burpolar_e", _chemplp_plp_burpolar_e_el);
			_load_bot_fel (buffer, "chemplp_plp_hbond_e", _chemplp_plp_hbond_e_el);
			_load_bot_fel (buffer, "chemplp_plp_metal_e", _chemplp_plp_metal_e_el);
			_load_bot_fel (buffer, "chemplp_plp_repulsive_weight", _chemplp_plp_repulsive_weight_el);
			_load_bot_fel (buffer, "chemplp_tors_weight", _chemplp_tors_weight_el);
			_load_bot_fel (buffer, "chemplp_lipo_weight", _chemplp_lipo_weight_el);
			_load_bot_fel (buffer, "chemplp_intercept_weight", _chemplp_intercept_weight_el);
//			1.5
			_load_bot_lf (buffer, "protein_file", load_protein_file);
			_load_bot_lf (buffer, "ligand_file", load_ligand_file);
			_load_bot_lf (buffer, "ligand_list", load_ligand_list_file);
//			1.6
			_load_bot_le (buffer, "output_dir", _output_dir);
			_load_bot_chb (buffer, "write_protein_conformations", _write_protein_conformations_chb);
			_load_bot_chb (buffer, "write_protein_bindingsite", _write_protein_bindingsite_chb);
			_load_bot_chb (buffer, "write_protein_splitted", _write_protein_splitted_chb);
			_load_bot_chb (buffer, "write_rescored_structures", _write_rescored_structures_chb);
			_load_bot_chb (buffer, "write_multi_mol2", _write_multi_mol2_chb);
			_load_bot_chb (buffer, "write_ranking_links", _write_ranking_links_chb);
			_load_bot_chb (buffer, "write_ranking_multi_mol2", _write_ranking_multi_mol2_chb);
			_load_bot_chb (buffer, "write_per_atom_scores", _write_per_atom_scores_chb);
			_load_bot_chb (buffer, "write_merged_ligand", _write_merged_ligand_chb);
			_load_bot_chb (buffer, "write_merged_protein", _write_merged_protein_chb);
			_load_bot_chb (buffer, "keep_original_mol2_description", _keep_original_mol2_description_chb);
			_load_bot_box (buffer, "merge_multi_conf_output", _merge_multi_conf_output_box);
			_load_bot_le (buffer, "merge_multi_conf_character", _merge_multi_conf_character);
			_load_bot_iel (buffer, "merge_multi_conf_after_characters", _merge_multi_conf_after_characters_el);
//			1.7
			_load_bot_list_view (buffer, "chemplp_protein_hb_costraint", costraints_file_list);
			_load_bot_list_view (buffer, "shape_costraint", costraints_file_list);
			_load_bot_list_view (buffer, "surface_distance_costraint", costraints_file_list);
			_load_bot_list_view (buffer, "ligand_intra_distance_costraint", costraints_file_list);
			_load_bot_list_view (buffer, "protein_ligand_distance_costraint", costraints_file_list);
//			1.8
			_load_bot_list_view (buffer, "flexible_protein_side_chain_string", flexibility_file_list);
			_load_bot_list_view (buffer, "flexible_protein_side_chain_number", flexibility_file_list);
			_load_bot_fel (buffer, "intra_protein_score_weight", _intra_protein_score_weight_el);
			_load_bot_list_view (buffer, "fix_protein_bond", flexibility_file_list);
//			1.9
			_load_bot_chb (buffer, "rigid_ligand", _rigid_ligand_chb);
			_load_bot_chb (buffer, "rigid_all", _rigid_all_chb);
//			1.10
			_load_bot_water (buffer, "water_molecule");
			_load_bot_lf (buffer, "water_molecule_definition", load_water_file);
			_load_bot_fel (buffer, "water_protein_hb_weight", _water_protein_hb_weight_el);
			_load_bot_fel (buffer, "water_ligand_hb_weight", _water_ligand_hb_weight_el);
			_load_bot_fel (buffer, "water_water_hb_weight", _water_water_hb_weight_el);
			_load_bot_fel (buffer, "no_water_ligand_hb_penalty", _no_water_ligand_hb_penalty_el);
			_load_bot_fel (buffer, "water_enable_penalty", _water_enable_penalty_el);
		}
	}
}


void PLANTSMenu::_load_bot_bindingsite (string& buffer, const char *reference) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		int n = 1;
		string x = "";
		string y = "";
		string z = ""; 
		string r = "";
		while (!line2.eof ()) {
			string q;
			line2 >> q;
			if (n == 1) x = q;
			else if (n == 2) y = q;
			else if (n == 3) z = q;
			n ++;
		}
		_site_x_el ->spinbox ->setValue(QString::fromStdString (x).toDouble());
		_site_y_el ->spinbox ->setValue(QString::fromStdString (y).toDouble());
		_site_z_el ->spinbox ->setValue(QString::fromStdString (z).toDouble());
	}
}


void PLANTSMenu::_load_bot_water (string& buffer, const char *reference) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		int n = 1;
		string x = "";
		string y = "";
		string z = ""; 
		string r = "";
		while (!line2.eof ()) {
			string q;
			line2 >> q;
			if (n == 1) x = q;
			if (n == 2) y = q;
			if (n == 3) z = q;
			if (n == 4) r = q;
			n = n + 1;
		}
		_water_site_x_el ->spinbox ->setValue(QString::fromStdString (x).toDouble());
		_water_site_y_el ->spinbox ->setValue(QString::fromStdString (y).toDouble());
		_water_site_z_el ->spinbox ->setValue(QString::fromStdString (z).toDouble());
		_water_site_radius_el ->spinbox ->setValue(QString::fromStdString (r).toDouble());
	}
}


void PLANTSMenu::_load_bot_list_view (string& buffer, const char *reference, MyListView *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	if (token == reference) {
		string s = token;
		int n = 1;
		string a = "";
		string b = "";
		string c = ""; 
		string d = "";
		string e = "";
		while (!line2.eof ()) {
			string q;
			line2 >> q;
			s.append (" ");s.append (q);
//			target ->linedit ->setText(q.c_str ());
			if (n == 1) a = q;
			if (n == 2) b = q;
			if (n == 3) c = q;
			if (n == 4) d = q;
			if (n == 5) e = q;
			n = n + 1;
		}
		string st = token+" "+a+" "+b+" "+c+" "+d+" "+e;
		update_boxes (target, st.c_str(), 1);
	}
}


void PLANTSMenu::_run_slot () {
	_save_slot ();

	cerr << plants_exe << " --mode screen "<< saved_file <<endl;

}


void PLANTSMenu::_save_slot () {
	QFileDialog save_file(this);

	QStringList filters;
	filters << "PLANTS config file (*.pcfg)";
	
	save_file.setNameFilters (filters);
	save_file.setAcceptMode (QFileDialog::AcceptSave);
	save_file.exec();

	QString ext = save_file.selectedNameFilter();

	QString filter = save_file.selectedNameFilter();
	filter.truncate (filter.lastIndexOf (')'));
	filter.remove (0, filter.indexOf ('*')+1);
	QString s = save_file.selectedFile();
	if (!s.contains ('.') ) s+= filter;


	saved_file = s.toStdString ().c_str ();

	if (!s.isNull()) {
		ofstream ofs(s.toStdString ().c_str ());
		ofs << "#################################################################################\n"
			"#\n"
			"#        PLANTS 1.1 config file written by Zodiac "<<VERSION<<" (www.zeden.org)\n"
			"#\n"
			"##################################################################################\n"
			"\n"
			"# 1.1 search algorithm\n"
			"search_speed ";
		ofs << _search_speed_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << "aco_ants ";
		ofs << _aco_ants_el ->spinbox ->value()<<endl;
		ofs << "aco_evap ";
		ofs << _aco_evap_el ->spinbox ->value()<<endl;
		ofs << "aco_sigma ";
		ofs << _aco_sigma_el ->spinbox ->value()<<endl;
		ofs << "flip_amide_bonds ";
		ofs << _flip_amide_bonds_chb ->isChecked()<<endl;
		ofs << "flip_planar_n ";
		ofs << _flip_planar_n_chb ->isChecked()<<endl;
		ofs << "force_flipped_bonds_planarity ";
		ofs << _force_flipped_bonds_planarity_chb ->isChecked()<<endl;
		ofs << "force_planar_bond_rotation ";
		ofs << _force_planar_bond_rotation_chb ->isChecked()<<endl;
		ofs << "rescore_mode ";
		ofs << _rescore_mode_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << "\n"
			"# 1.2 binding site\n"
			"bindingsite_center ";
		ofs << _site_x_el ->spinbox ->value();
		ofs << 	" ";
		ofs << _site_y_el ->spinbox ->value();
		ofs << 	" ";
		ofs << _site_z_el ->spinbox ->value()<<endl;
		ofs << 	"bindingsite_radius ";
		ofs << _site_radius_el ->spinbox ->value()<<endl;
		ofs << "\n"
			"# 1.3 cluster algorithm\n"
			"cluster_rmsd ";
		ofs << _cluster_rmsd_el ->spinbox ->value()<<endl;
		ofs << "cluster_structures ";
		ofs << _cluster_structures_el ->spinbox ->value()<<endl;
		ofs << "\n"
			"# 1.4 scoring functions\n"
			"scoring_function ";
		ofs << _algorithm_type_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << "outside_binding_site_penalty ";
		ofs << _outside_binding_site_penalty_el ->spinbox ->value()<<endl;
		ofs << "enable_sulphur_acceptors ";
		ofs << _enable_sulphur_acceptors_chb ->isChecked()<<endl;
		ofs << "ligand_intra_score ";
		ofs << _ligand_intra_score_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << "chemplp_clash_include_14 ";
		ofs << _chemplp_clash_include_14_el ->spinbox ->value()<<endl;
		ofs << "chemplp_clash_include_HH ";
		ofs << _chemplp_clash_include_HH_chb ->isChecked()<<endl;

		if (_algorithm_type_cb ->currentData ().toString ().toStdString () == "chemplp") {
			ofs << "chemplp_weak_cho ";
			ofs << _chemplp_weak_cho_chb ->isChecked()<<endl;
			ofs << "chemplp_charged_hb_weight ";
			ofs << _chemplp_charged_hb_weight_el ->spinbox ->value()<<endl;
			ofs << "chemplp_charged_metal_weight ";
			ofs << _chemplp_charged_metal_weight_el ->spinbox ->value()<<endl;
			ofs << "chemplp_hbond_weight ";
			ofs << _chemplp_hbond_weight_el ->spinbox ->value()<<endl;
			ofs << "chemplp_hbond_cho_weight ";
			ofs << _chemplp_hbond_cho_weight_el ->spinbox ->value()<<endl;
			ofs << "chemplp_metal_weight ";
			ofs << _chemplp_metal_weight_el ->spinbox ->value()<<endl;
			ofs << "chemplp_plp_weight ";
			ofs << _chemplp_plp_weight_el ->spinbox ->value()<<endl;
			ofs << "chemplp_plp_steric_e ";
			ofs << _chemplp_plp_steric_e_el ->spinbox ->value()<<endl;
			ofs << "chemplp_plp_burpolar_e ";
			ofs << _chemplp_plp_burpolar_e_el ->spinbox ->value()<<endl;
			ofs << "chemplp_plp_hbond_e ";
			ofs << _chemplp_plp_hbond_e_el ->spinbox ->value()<<endl;
			ofs << "chemplp_plp_metal_e ";
			ofs << _chemplp_plp_metal_e_el ->spinbox ->value()<<endl;
			ofs << "chemplp_plp_repulsive_weight ";
			ofs << _chemplp_plp_repulsive_weight_el ->spinbox ->value()<<endl;
			ofs << "chemplp_tors_weight ";
			ofs << _chemplp_tors_weight_el ->spinbox ->value()<<endl;
			ofs << "chemplp_lipo_weight ";
			ofs << _chemplp_lipo_weight_el ->spinbox ->value()<<endl;
			ofs << "chemplp_intercept_weight ";
			ofs << _chemplp_intercept_weight_el ->spinbox ->value()<<endl;
		}
		else {
			ofs << "plp_steric_e ";
			ofs << _plp_steric_e_el ->spinbox ->value()<<endl;
			ofs << "plp_burpolar_e ";
			ofs << _plp_burpolar_e_el ->spinbox ->value()<<endl;
			ofs << "plp_hbond_e ";
			ofs << _plp_hbond_e_el ->spinbox ->value()<<endl;
			ofs << "plp_metal_e ";
			ofs << _plp_metal_e_el ->spinbox ->value()<<endl;
			ofs << "plp_repulsive_weight ";
			ofs << _plp_repulsive_weight_el ->spinbox ->value()<<endl;
			ofs << "plp_tors_weight ";
			ofs << _plp_tors_weight_el ->spinbox ->value()<<endl;
		}

		ofs << "\n"
			"# 1.5 input option\n"
			"protein_file ";
		ofs << load_protein_file ->val ()<<endl;

		if (load_ligand_file ->val () != "") {
			ofs << "ligand_file ";
			ofs << load_ligand_file ->val ()<<endl;
		}

		if (load_ligand_list_file ->val () != "") {
			ofs << "ligand_list ";
			ofs << load_ligand_list_file ->val ()<<endl;
		}

		ofs << "\n"
			"# 1.6 output option\n"
			"output_dir ";
		ofs << _output_dir ->linedit ->text ().toStdString ()<<endl;
		ofs << "write_protein_conformations ";
		ofs << _write_protein_conformations_chb ->isChecked()<<endl;
		ofs << "write_protein_bindingsite ";
		ofs << _write_protein_bindingsite_chb ->isChecked()<<endl;
		ofs << "write_protein_splitted ";
		ofs << _write_protein_splitted_chb ->isChecked()<<endl;
		ofs << "write_rescored_structures ";
		ofs << _write_rescored_structures_chb ->isChecked()<<endl;
		ofs << "write_multi_mol2 ";
		ofs << _write_multi_mol2_chb ->isChecked()<<endl;
		ofs << "write_ranking_links ";
		ofs << _write_ranking_links_chb ->isChecked()<<endl;
		ofs << "write_ranking_multi_mol2 ";
		ofs << _write_ranking_multi_mol2_chb ->isChecked()<<endl;
		ofs << "write_per_atom_scores ";
		ofs << _write_per_atom_scores_chb ->isChecked()<<endl;
		ofs << "write_merged_ligand ";
		ofs << _write_merged_ligand_chb ->isChecked()<<endl;
		ofs << "write_merged_protein ";
		ofs << _write_merged_protein_chb ->isChecked()<<endl;
		ofs << "write_merged_water ";
		ofs << _write_merged_water_chb ->isChecked()<<endl;
		ofs << "keep_original_mol2_description ";
		ofs << _keep_original_mol2_description_chb ->isChecked()<<endl;
		ofs << "merge_multi_conf_output ";
		ofs << _merge_multi_conf_output_box ->isChecked()<<endl;
		ofs << "merge_multi_conf_character ";
		ofs << _merge_multi_conf_character ->linedit ->text ().toStdString ()<<endl;
		ofs << "merge_multi_conf_after_characters ";
		ofs << _merge_multi_conf_after_characters_el ->spinbox ->value()<<endl;
		ofs << "\n"
			"# 1.7 costraints\n";

		for (unsigned int i = 0; i <costraints_file_list ->_lw ->count (); i++) {
			QString q = costraints_file_list ->_lw ->item(i) ->text ();
			ofs << q.toStdString ()<<endl;
		}

		ofs << "\n"
			"# 1.8 flexible side-chains\n";
		for (unsigned int i = 0; i <flexibility_file_list ->_lw ->count (); i++) {
			QString q = flexibility_file_list ->_lw ->item(i) ->text ();
			ofs << q.toStdString ()<<endl;
		}
		ofs << "\n"
			"# 1.9 multiconformer docking\n"
			"rigid_ligand ";
		ofs << _rigid_ligand_chb ->isChecked()<<endl;
		ofs << "rigid_all ";
		ofs << _rigid_all_chb ->isChecked()<<endl;

		if (_water_site_x_el ->spinbox ->value() ) {
			ofs << "\n"
				"# 1.10 water\n"
				"water_molecule ";
			ofs << _water_site_x_el ->spinbox ->value();
			ofs << 	" ";
			ofs << _water_site_y_el ->spinbox ->value();
			ofs << 	" ";
			ofs << _water_site_z_el ->spinbox ->value();
			ofs << 	" ";
			ofs << _water_site_radius_el ->spinbox ->value()<<endl;
			ofs << "water_molecule_definition ";
			ofs << load_water_file ->val ()<<endl;
			ofs << "water_protein_hb_weight ";
			ofs << _water_protein_hb_weight_el ->spinbox ->value()<<endl;
			ofs << "water_ligand_hb_weight ";
			ofs << _water_ligand_hb_weight_el ->spinbox ->value()<<endl;
			ofs << "water_water_hb_weight ";
			ofs << _water_water_hb_weight_el ->spinbox ->value()<<endl;
			ofs << "no_water_ligand_hb_penalty ";
			ofs << _no_water_ligand_hb_penalty_el ->spinbox ->value()<<endl;
			ofs << "water_enable_penalty ";
			ofs << _water_enable_penalty_el ->spinbox ->value()<<endl;
		}

		ofs << "\n"
			"##################################################################################\n" ;
	}
}


void PLANTSMenu::_restore_slot () {
//			1.1
			_search_speed_cb->_combo_box ->setCurrentIndex(0);
			_aco_ants_el ->spinbox ->setValue(20);
			_aco_evap_el ->spinbox ->setValue(1);
			_aco_sigma_el ->spinbox ->setValue(1);
			_flip_amide_bonds_chb ->setCheckState (Qt::Unchecked);
			_flip_planar_n_chb ->setCheckState (Qt::Unchecked);
			_force_flipped_bonds_planarity_chb ->setCheckState (Qt::Unchecked);
			_force_planar_bond_rotation_chb ->setCheckState (Qt::Checked);
			_rescore_mode_cb->_combo_box ->setCurrentIndex(0);
//			1.2
			_site_x_el ->spinbox ->setValue(0);
			_site_y_el ->spinbox ->setValue(0);
			_site_z_el ->spinbox ->setValue(0);
			_site_radius_el ->spinbox ->setValue(10);
//			1.3
			_cluster_rmsd_el ->spinbox ->setValue(2);
			_cluster_structures_el ->spinbox ->setValue(10);
//			1.4
			_algorithm_type_cb->_combo_box ->setCurrentIndex(0);
			_outside_binding_site_penalty_el ->spinbox ->setValue(50);
			_enable_sulphur_acceptors_chb ->setCheckState (Qt::Unchecked);
			_ligand_intra_score_cb->_combo_box ->setCurrentIndex(0);
			_chemplp_clash_include_14_el ->spinbox ->setValue(0.25);
			_chemplp_clash_include_HH_chb ->setCheckState (Qt::Unchecked);
			_plp_steric_e_el ->spinbox ->setValue(-0.4);
			_plp_burpolar_e_el ->spinbox ->setValue(-0.05);
			_plp_hbond_e_el ->spinbox ->setValue(-2);
			_plp_metal_e_el ->spinbox ->setValue(-4);
			_plp_repulsive_weight_el ->spinbox ->setValue(0.5);
			_plp_tors_weight_el ->spinbox ->setValue(1);
			_chemplp_weak_cho_chb ->setCheckState (Qt::Checked);
			_chemplp_charged_hb_weight_el ->spinbox ->setValue(2);
			_chemplp_charged_metal_weight_el ->spinbox ->setValue(2);
			_chemplp_hbond_weight_el ->spinbox ->setValue(-3);
			_chemplp_hbond_cho_weight_el ->spinbox ->setValue(-3);
			_chemplp_metal_weight_el ->spinbox ->setValue(-6);
			_chemplp_plp_weight_el ->spinbox ->setValue(1);
			_chemplp_plp_steric_e_el ->spinbox ->setValue(-0.4);
			_chemplp_plp_burpolar_e_el ->spinbox ->setValue(-0.1);
			_chemplp_plp_hbond_e_el ->spinbox ->setValue(-1);
			_chemplp_plp_metal_e_el ->spinbox ->setValue(-1);
			_chemplp_plp_repulsive_weight_el ->spinbox ->setValue(1);
			_chemplp_tors_weight_el ->spinbox ->setValue(2);
			_chemplp_lipo_weight_el ->spinbox ->setValue(0);
			_chemplp_intercept_weight_el ->spinbox ->setValue(-20);
//			1.5
			load_protein_file ->linedit ->clear();
			load_ligand_file ->linedit ->clear();
			input_file ->_lw ->clear();
			load_ligand_list_file ->linedit ->clear();
			input_file_list ->_lw ->clear();
//			1.6
			_output_dir ->linedit ->clear();
			_write_protein_conformations_chb ->setCheckState (Qt::Unchecked);
			_write_protein_bindingsite_chb ->setCheckState (Qt::Unchecked);
			_write_protein_splitted_chb ->setCheckState (Qt::Unchecked);
			_write_rescored_structures_chb ->setCheckState (Qt::Unchecked);
			_write_multi_mol2_chb ->setCheckState (Qt::Checked);
			_write_ranking_links_chb ->setCheckState (Qt::Unchecked);
			_write_ranking_multi_mol2_chb ->setCheckState (Qt::Unchecked);
			_write_per_atom_scores_chb ->setCheckState (Qt::Unchecked);
			_write_merged_ligand_chb ->setCheckState (Qt::Unchecked);
			_write_merged_protein_chb ->setCheckState (Qt::Unchecked);
			_keep_original_mol2_description_chb ->setCheckState (Qt::Checked);
			_merge_multi_conf_output_box ->setChecked (Qt::Unchecked);
			_merge_multi_conf_character ->linedit ->setText("_");
			_merge_multi_conf_after_characters_el ->spinbox ->setValue(1);
//			1.7
			costraints_file_list ->_lw ->clear();
			_protein_hb_costraint_cb ->_combo_box ->clear();
			_ligand_surface_distance_cb ->_combo_box ->clear();
			_number_a_ligand_intra_distance_cb ->_combo_box ->clear();
			_number_b_ligand_intra_distance_cb ->_combo_box ->clear();
			_number_a_protein_ligand_distance_cb ->_combo_box ->clear();
			_number_b_protein_ligand_distance_cb ->_combo_box ->clear();
//			1.8
			flexibility_file_list ->_lw ->clear();
			_flexible_protein_side_chain_string ->_combo_box ->clear();
			_fix_protein_bond ->_combo_box ->clear();
			_flexible_protein_side_chain_number ->_combo_box ->clear();
//			1.9
			_rigid_ligand_chb ->setCheckState (Qt::Unchecked);
			_rigid_all_chb ->setCheckState (Qt::Unchecked);
//			1.10
			_water_site_x_el ->spinbox ->setValue(0);
			_water_site_y_el ->spinbox ->setValue(0);
			_water_site_z_el ->spinbox ->setValue(0);
			_water_site_radius_el ->spinbox ->setValue(2);
			load_water_file ->linedit ->clear();
			_water_protein_hb_weight_el ->spinbox ->setValue(1);
			_water_ligand_hb_weight_el ->spinbox ->setValue(1);
			_water_water_hb_weight_el ->spinbox ->setValue(1);
			_no_water_ligand_hb_penalty_el ->spinbox ->setValue(0);
			_water_enable_penalty_el ->spinbox ->setValue(8);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	GAMESS menu FE
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


GamessMenu::GamessMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, true, true) {
	setWindowTitle ("GAMESS input");

// Main tab
	ZNAdvancedWidget *main_tab = new ZNAdvancedWidget ();
    	_tabs ->addTab (main_tab, "Main");

	_base_box = new MyGroupBox (main_tab ->basic_layout (), "Base settings");
	
	_run_type_cb = new MyComboBox (_base_box ->layout (), "Calculation type:");
	_run_type_cb ->combo_box () ->insertItem (0, "Equilibrium geometry", "optimize");
	_run_type_cb ->combo_box () ->insertItem (1, "Single point energy", "energy");
	_run_type_cb ->combo_box () ->insertItem (2, "Gradient", "gradient");
	_run_type_cb ->combo_box () ->insertItem (3, "Transition state", "sadpoint");
	_run_type_cb ->combo_box () ->insertItem (4, "Vibrational frequencies", "hessian");
//	_run_type_cb ->combo_box () ->insertItem (5, "Surface", "surf");
//	_run_type_cb ->combo_box () ->insertItem (6, "Properties", "prop");
//	_run_type_cb ->combo_box () ->insertItem (7, "IRC calculation", "irc");

	_semi_emp_cb = new QComboBox ();
	_semi_emp_cb ->insertItem (0, "AM1", "AM1");
	_semi_emp_cb ->insertItem (1, "PM3", "PM3");
	_semi_emp_cb ->insertItem (2, "MNDO", "MNDO");
	_semi_emp_cb ->insertItem (3, "MIDI", "MIDI");

	_H_F_cb = new QComboBox ();
	_H_F_cb ->insertItem (0, "STO-3G", "130000");
	_H_F_cb ->insertItem (1, "STO-4G", "140000");
	_H_F_cb ->insertItem (2, "3-21G", "230000");
	_H_F_cb ->insertItem (3, "3-21G*", "231000");
	_H_F_cb ->insertItem (4, "3-21G**", "231100");
	_H_F_cb ->insertItem (5, "3-21+G*", "231010");
	_H_F_cb ->insertItem (6, "3-21+G**", "231110");
	_H_F_cb ->insertItem (7, "4-31G", "340000");
	_H_F_cb ->insertItem (8, "4-31G*", "341000");
	_H_F_cb ->insertItem (9, "4-31G**", "341100");
	_H_F_cb ->insertItem (10, "4-31+G*", "341010");
	_H_F_cb ->insertItem (11, "4-31+G**", "341110");
	_H_F_cb ->insertItem (12, "6-31G", "360000");
	_H_F_cb ->insertItem (13, "6-31G*", "361000");
	_H_F_cb ->insertItem (14, "6-31G**", "361100");
	_H_F_cb ->insertItem (15, "6-31+G*", "361010");
	_H_F_cb ->insertItem (16, "6-31+G**", "361110");
	_H_F_cb ->insertItem (17, "6-311G", "460000");
	_H_F_cb ->insertItem (18, "6-311G*", "461000");
	_H_F_cb ->insertItem (19, "6-311G**", "461100");
	_H_F_cb ->insertItem (20, "6-311+G*", "461010");
	_H_F_cb ->insertItem (21, "6-311+G**", "461110");
//	_H_F_cb ->insertItem (22, "6-311++G(2d,p)", "462111");


	_gbasis_hcb = new MyHideComboBox (_base_box ->layout (), "Basis set:");

	_gbasis_hcb -> layout () -> addWidget (_semi_emp_cb);
	_gbasis_hcb -> layout () -> addWidget (_H_F_cb);

	_gbasis_hcb ->insertItem (_semi_emp_cb, 0, "Semi-empirical", "S_e");
	_gbasis_hcb ->insertItem (_H_F_cb, 1, "Harteee-Fock", "H_F");

	connect (_gbasis_hcb ->combo_box (), SIGNAL (currentIndexChanged ( int )), this, SLOT (hide_solv ( int )) );


// Solvation group
	_solv_box = new MyGroupBox (main_tab ->basic_layout (), "Solvation", true, false);
	_solv_box ->hide ();

	_solvent_cb = new MyComboBox (_solv_box ->layout (), "Solvent:");
	_solvent_cb ->combo_box () ->insertItem (0, "No solvent", "gas");
	_solvent_cb ->combo_box () ->insertItem (1, "Water", "water");
	_solvent_cb ->combo_box () ->insertItem (2, "1,2-dichloroethane", "12dclet");
	_solvent_cb ->combo_box () ->insertItem (3, "Acetone", "acetone");
	_solvent_cb ->combo_box () ->insertItem (4, "Aniline", "aniline");
	_solvent_cb ->combo_box () ->insertItem (5, "Benzene", "benzene");
	_solvent_cb ->combo_box () ->insertItem (6, "Chlorobenzene", "clbenz");
	_solvent_cb ->combo_box () ->insertItem (7, "Chloroform", "clform");
	_solvent_cb ->combo_box () ->insertItem (8, "CCl4", "ctcl");
	_solvent_cb ->combo_box () ->insertItem (9, "Cyclohexane", "cychex");
	_solvent_cb ->combo_box () ->insertItem (10, "DMSO", "dmso");
	_solvent_cb ->combo_box () ->insertItem (11, "Ethanol", "c2h5oh");
	_solvent_cb ->combo_box () ->insertItem (12, "Methanol", "ch3oh");
	_solvent_cb ->combo_box () ->insertItem (13, "Methyl chloride", "methycl");
	_solvent_cb ->combo_box () ->insertItem (14, "Neptane", "neptane");
	_solvent_cb ->combo_box () ->insertItem (15, "Nitromethane", "nitmet");
	_solvent_cb ->combo_box () ->insertItem (16, "THF", "thf");
	_solvent_cb ->combo_box () ->insertItem (17, "Toluene", "toluene");

	connect (_solvent_cb ->combo_box (), SIGNAL (currentIndexChanged ( int )), this, SLOT (hide_temp ( int )) );

	_temp = 298.0;
	_tabs_el = new MyFloatEditLine (_solv_box ->layout (), "Absolute temperature (K):", _temp, 0);
	_tabs_el ->hide ();


// Main tab - advanced
	_scftyp_cb = new MyComboBox (main_tab  ->advanced_layout (), "Wavefunction:");
	_scftyp_cb ->combo_box () ->insertItem (0, "RHF", "RHF");
	_scftyp_cb ->combo_box () ->insertItem (1, "ROHF", "ROHF");
	_scftyp_cb ->combo_box () ->insertItem (2, "UHF", "UHF");

	_mult_cb = new MyComboBox (main_tab  ->advanced_layout (), "Moltiplicity:");
	_mult_cb ->combo_box () ->insertItem (0, "Singlet/Closed", "1");
	_mult_cb ->combo_box () ->insertItem (1, "Doublet/Radical", "2");
	_mult_cb ->combo_box () ->insertItem (2, "Triplet/Biradical", "3");
	_mult_cb ->combo_box () ->insertItem (3, "Quartet", "4");
	_mult_cb ->combo_box () ->insertItem (4, "Quintet", "5");
	_mult_cb ->combo_box () ->insertItem (5, "Sextet", "6");
	_mult_cb ->combo_box () ->insertItem (6, "Septet", "7");

	_guess_cb = new MyComboBox (main_tab  ->advanced_layout (), "Type of initial orbtal guess:");
	_guess_cb ->combo_box () ->insertItem (0, "HUCKEL", "huckel");
	_guess_cb ->combo_box () ->insertItem (1, "HCORE", "hcore");
	_guess_cb ->combo_box () ->insertItem (2, "MOREAD", "moread");
	_guess_cb ->combo_box () ->insertItem (3, "RDMINI", "rdmini");
	_guess_cb ->combo_box () ->insertItem (4, "MOSAVED", "mosaved");
	_guess_cb ->combo_box () ->insertItem (5, "SKIP", "skip");
	_guess_cb ->combo_box () ->insertItem (6, "FMO", "fmo");


// ADD ICHARG, but can'I obtain this data from zodiac????
// ADD UNITSICHARG, but can'I obtain this data from zodiac????

// elpot tab
	ZNAdvancedWidget *elpot_tab = new ZNAdvancedWidget ();
    	_tabs ->addTab (elpot_tab, "ELPOT");

	_elpot_box = new MyGroupBox (elpot_tab ->basic_layout (), "Settings", true, false);

	_where_cb = new MyComboBox (_elpot_box ->layout (), "Where:");
	_where_cb ->combo_box () ->insertItem (0, "At each nucleus", "nuclei");
	_where_cb ->combo_box () ->insertItem (1, "Center of mass", "comass");
	_where_cb ->combo_box () ->insertItem (2, "At points givent in $points", "points");
	_where_cb ->combo_box () ->insertItem (3, "At grid givent in $grid", "grid");
	_where_cb ->combo_box () ->insertItem (4, "At points controlled by $pdc", "pdc");

	connect (_where_cb ->combo_box (), SIGNAL (currentIndexChanged ( int )), this, SLOT (hide_elpot ( int )) );

	_output_cb = new MyComboBox (_elpot_box ->layout (), "Output:");
	_output_cb ->combo_box () ->insertItem (0, "punch and paper", "both");
	_output_cb ->combo_box () ->insertItem (1, "punch", "punch");
	_output_cb ->combo_box () ->insertItem (2, "paper", "paper");

	_points_box = new MyGroupBox (elpot_tab ->basic_layout (), "$POINTS");
	_points_box ->hide ();

	_grid_box = new MyGroupBox (elpot_tab ->basic_layout (), "$GRID");
	_grid_box ->hide ();

	_pdc_box = new MyGroupBox (elpot_tab ->basic_layout (), "$PDC");
	_pdc_box ->hide ();

	_ptsel_cb = new MyComboBox (_pdc_box ->layout (), "PTSEL:");
	_ptsel_cb ->combo_box () ->insertItem (0, "Geodesic", "geodesic");
	_ptsel_cb ->combo_box () ->insertItem (1, "Connolly", "connolly");
	_ptsel_cb ->combo_box () ->insertItem (2, "Chelpg", "chelpg");

	_constr_cb = new MyComboBox (_pdc_box ->layout (), "CONSTR:");
	_constr_cb ->combo_box () ->insertItem (0, "None", "none");
	_constr_cb ->combo_box () ->insertItem (1, "Charge", "charge");
	_constr_cb ->combo_box () ->insertItem (2, "Dipole", "dipole");
	_constr_cb ->combo_box () ->insertItem (3, "Qupole", "qupole");

	_vdwscl = 1.4;
	_vdwscl_fel = new MyFloatEditLine (_pdc_box ->layout (), "First shell of VDW spheres:", _vdwscl, 0);

	_vdwinc = 0.2;
	_vdwinc_fel = new MyFloatEditLine (_pdc_box ->layout (), "Increment for successive shells:", _vdwinc, 0);

	_layer = 4;
	_layer_iel = new MyIntegerEditLine (_pdc_box ->layout (), "Number of layers:", _layer, 0);

	_maxpdc = 4;
	_maxpdc_iel = new MyIntegerEditLine (_pdc_box ->layout (), "Total number of points:", _maxpdc, 0);

// System tab
	ZNWidget *_system_tab = new ZNWidget ();
    	_tabs ->addTab (_system_tab, "System");

	_run_box = new MyGroupBox (_system_tab ->layout (), "Run settings");

	_exetyp_cb = new MyComboBox (_run_box ->layout (), "Type:");
	_exetyp_cb ->combo_box () ->insertItem (0, "Run calculation", "run");
	_exetyp_cb ->combo_box () ->insertItem (1, "Check errors", "check");

	_timlim = 600;
	_timlim_iel = new MyIntegerEditLine (_run_box ->layout (), "Time limit (minutes):", _timlim, 0);

	_memory = 2000;
	_memory_iel = new MyIntegerEditLine (_run_box ->layout (), "Memory limit (MB):", _memory, 0);


// test tab
/*	ZNAdvancedWidget *_test_tab = new ZNAdvancedWidget ();
    	_tabs ->addTab (_test_tab, "Test");

	MyComboBox *_sopra_cb;
	_sopra_cb = new MyComboBox (_test_tab ->basic_layout (), "Run type:");

	MyComboBox *_sotto_cb;
	_sotto_cb = new MyComboBox (_test_tab ->advanced_layout (), "Run type2:");
*/

// menu
	add_menu ();
	add_help ();
}


void GamessMenu::_show_about ()
{
	QMessageBox::about( _data ->ddwin, "About GAMESS input front end" ,
			QString(("ZODIAC. Version "+VERSION+ "\n"
			"GAMESS input front end. Code by Alberto Massarotti. All rights reserved.\n"
			"For bug reports, suggestions or any other feedback please visit our website at www.zeden.org\n\n").c_str()) );
}


void GamessMenu::hide_temp (int i) {
	if (i == 0) {
		_tabs_el ->hide ();
		_solv = 0;
	}
	else {
		_tabs_el ->show ();
		_solv = 1;
	} 
}


void GamessMenu::hide_solv (int i) {
	if (i == 0) {
		_solv_box ->hide ();
	} 
	else {
		_solv_box ->show ();
	} 
}

void GamessMenu::hide_elpot (int i) {
	_points_box ->hide ();
	_grid_box ->hide ();
	_pdc_box ->hide ();
	if (i == 2) {
		_points_box ->show ();
	} 
	else if (i == 3) {
		_grid_box ->show ();
	}
	else if (i == 4) {
		_pdc_box ->show ();
	}
}


void GamessMenu::add_menu () {
    	QMenu *file = new QMenu(tr("&File"), this );
	_load_action = new QAction (tr ("&Load"), this);
	connect (_load_action, SIGNAL (triggered ()), this, SLOT (_load_slot ()));
	file -> addAction (_load_action);

	_save_action = new QAction (tr ("&Save"), this);
	connect (_save_action, SIGNAL (triggered ()), this, SLOT (_save_slot ()));
    	file -> addAction (_save_action);
	addMenu (file);

    	QMenu *settings = new QMenu(tr("&Settings"), this );
	_restore_action = new QAction (tr ("&Restore to default"), this);
	connect (_restore_action, SIGNAL (triggered ()), this, SLOT (_restore_slot ()));
    	settings -> addAction (_restore_action);
	addMenu (settings);

}


void GamessMenu::_restore_slot () {
_run_type_cb ->_combo_box ->setCurrentIndex(0);
_gbasis_hcb ->_combo_box ->setCurrentIndex(0);

_semi_emp_cb ->setCurrentIndex(0);

_H_F_cb ->setCurrentIndex(0);
_mult_cb ->_combo_box ->setCurrentIndex(0);

_scftyp_cb ->_combo_box ->setCurrentIndex(0);
_solvent_cb ->_combo_box ->setCurrentIndex(0);

_tabs_el ->spinbox ->setValue(298);
_exetyp_cb ->_combo_box ->setCurrentIndex(0);
_timlim_iel ->spinbox ->setValue(600);
_memory_iel ->spinbox ->setValue(2000);

_guess_cb ->_combo_box ->setCurrentIndex(0);


_where_cb ->_combo_box ->setCurrentIndex(0);
_output_cb ->_combo_box ->setCurrentIndex(0);
_ptsel_cb ->_combo_box ->setCurrentIndex(0);
_constr_cb ->_combo_box ->setCurrentIndex(0);

_elpot_box -> setChecked (false);
_vdwscl_fel ->spinbox ->setValue(1.4);
_vdwinc_fel ->spinbox ->setValue(0.2);
_layer_iel ->spinbox ->setValue(4);
_maxpdc_iel ->spinbox ->setValue(4);

}


void GamessMenu::_load_bot_iel (string& buffer, const char *reference, MyIntegerEditLine *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == reference) {
		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);

		if (end.toStdString () == "MEMORY") {
			target ->spinbox ->setValue(filter.toInt()/1000);
		}
		else {
			target ->spinbox ->setValue(filter.toInt());
		}
	}
}

void GamessMenu::_load_bot_fel (string& buffer, const char *reference, MyFloatEditLine *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == reference) {
		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);
		target ->spinbox ->setValue(filter.toDouble());
	}
}

void GamessMenu::_load_bot_box (string& buffer, const char *reference, MyGroupBox *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == reference) {
		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);
		if (filter == "0") target ->setChecked (false);
		if (filter == "1") target ->setChecked (true);
	}
}

void GamessMenu::_load_bot_cb (string& buffer, const char *reference, MyComboBox *target) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == reference) {

		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);
		target ->_combo_box ->setCurrentIndex (target ->_combo_box ->findData (filter));
	}
}

void GamessMenu::_load_bot_gbasis (string& buffer) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == "GBASIS") {

		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);

		if (filter.toStdString() == "AM1")  {
			_gbasis_hcb ->_combo_box ->setCurrentIndex (0);
			_semi_emp_cb ->setCurrentIndex (_semi_emp_cb ->findData (filter));
		}
		else if (filter.toStdString() == "PM3") {
			_gbasis_hcb ->_combo_box ->setCurrentIndex (0);
			_semi_emp_cb ->setCurrentIndex (_semi_emp_cb ->findData (filter));
		}
		else if (filter.toStdString() == "MNDO") {
			_gbasis_hcb ->_combo_box ->setCurrentIndex (0);
			_semi_emp_cb ->setCurrentIndex (_semi_emp_cb ->findData (filter));
		}
		else if (filter.toStdString() == "MIDI") {
			_gbasis_hcb ->_combo_box ->setCurrentIndex (0);
			_semi_emp_cb ->setCurrentIndex (_semi_emp_cb ->findData (filter));
		}
		else {
			_gbasis_hcb ->_combo_box ->setCurrentIndex (1);

			if (filter.toStdString() == "STO") {
				_gbasis_load = 100000;
			}
			else if (filter.toStdString() == "N21") {
				_gbasis_load = 200000;
			}
			else if (filter.toStdString() == "N31") {
				_gbasis_load = 300000;
			}
			else if (filter.toStdString() == "N311") {
				_gbasis_load = 400000;
			}
		}
	}
}

void GamessMenu::_load_bot_ngauss (string& buffer) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == "NGAUSS") {

		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);

		if (filter.toStdString() == "3")  {
			_gbasis_load = _gbasis_load + 30000;
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
		else if (filter.toStdString() == "4") {
			_gbasis_load = _gbasis_load + 40000;
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
		else if (filter.toStdString() == "6") {
			_gbasis_load = _gbasis_load + 60000;
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
	}
}

void GamessMenu::_load_bot_ndfunc (string& buffer) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == "NDFUNC") {

		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);

		if (filter.toStdString() == "0")  {
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
		else if (filter.toStdString() == "1") {
			_gbasis_load = _gbasis_load + 1000;
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
		else if (filter.toStdString() == "2") {
			_gbasis_load = _gbasis_load + 2000;
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
	}
}

void GamessMenu::_load_bot_npfunc (string& buffer) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == "NPFUNC") {

		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);

		if (filter.toStdString() == "0")  {
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
		else if (filter.toStdString() == "1")  {
			_gbasis_load = _gbasis_load + 100;
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
	}
}

void GamessMenu::_load_bot_diffsp (string& buffer) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == "DIFFSP") {

		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);

		if (filter.toStdString() == ".true.")  {
			_gbasis_load = _gbasis_load + 10;
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
	}
}

void GamessMenu::_load_bot_diffs (string& buffer) {
	istringstream line2(buffer);
	string token;
	line2 >> token;

	QString end = QString::fromStdString (token);
	end.truncate (end.lastIndexOf ('='));

	if (end.toStdString () == "DIFFS") {

		QString filter = QString::fromStdString (token);
		filter.remove (0, filter.indexOf ('=')+1);

		if (filter.toStdString() == ".true.")  {
			_gbasis_load = _gbasis_load + 1;
			_H_F_cb ->setCurrentIndex (_H_F_cb ->findData (_gbasis_load));
		}
	}
}

void GamessMenu::_load_slot () {
	QString s = QFileDialog::getOpenFileName(this, tr ("Open file"), "",tr("GAMESS input file (*.inp)"));

	_gbasis_load = 0;

	if (!s.isNull()) {
		ifstream ifs (s.toStdString ().c_str ());
		string buffer;

		_restore_slot ();

		while (getline(ifs, buffer)) {

			_load_bot_iel (buffer, "TIMLIM", _timlim_iel);
			_load_bot_iel (buffer, "MEMORY", _memory_iel);
			_load_bot_iel (buffer, "LAYER", _layer_iel);
			_load_bot_fel (buffer, "TABS", _tabs_el);
			_load_bot_iel (buffer, "MAXPDC", _maxpdc_iel);
			_load_bot_fel (buffer, "VDWSCL", _vdwscl_fel);
			_load_bot_fel (buffer, "VDWINC", _vdwinc_fel);

			_load_bot_cb (buffer, "RUNTYP", _run_type_cb);
			_load_bot_cb (buffer, "MULT", _mult_cb);
			_load_bot_cb (buffer, "SCFTYP", _scftyp_cb);
			_load_bot_cb (buffer, "SOLVENT", _solvent_cb);
			_load_bot_cb (buffer, "EXETYP", _exetyp_cb);
			_load_bot_cb (buffer, "GUESS", _guess_cb);
			_load_bot_cb (buffer, "WHERE", _where_cb);
			_load_bot_cb (buffer, "OUTPUT", _output_cb);
			_load_bot_cb (buffer, "PTSEL", _ptsel_cb);
			_load_bot_cb (buffer, "CONSTR", _constr_cb);

			_load_bot_box (buffer, "ELPOT", _elpot_box);

			_load_bot_gbasis (buffer);
			_load_bot_ngauss (buffer);
			_load_bot_ndfunc (buffer);
			_load_bot_npfunc (buffer);
			_load_bot_diffsp (buffer);
			_load_bot_diffs (buffer);

		}
	}
}


void GamessMenu::_save_slot () {
	QFileDialog save_file(this);

	QStringList filters;
	filters << "GAMESS input file (*.inp)";
	
	save_file.setNameFilters (filters);
	save_file.setAcceptMode (QFileDialog::AcceptSave);
	save_file.exec();

	QString ext = save_file.selectedNameFilter();

	QString filter = save_file.selectedNameFilter();
	filter.truncate (filter.lastIndexOf (')'));
	filter.remove (0, filter.indexOf ('*')+1);
	QString s = save_file.selectedFile();
	if (!s.contains ('.') ) s+= filter;

	
    if (!s.isNull()) {
        if (_data -> ddwin ->current_target) {
		ofstream ofs(s.toStdString ().c_str ());
		ofs << "!##################################################################################\n"
			"!#\n"
			"!#        GAMESS input file written by Zodiac "<<VERSION<<" (www.zeden.org)\n"
			"!#\n"
			"!##################################################################################\n"
			"\n"
			" $BASIS\n"
			"GBASIS=" ;
		if (_gbasis_hcb ->currentData ().toString ().toStdString () == "S_e") {
			ofs << _semi_emp_cb ->currentText ().toStdString () <<endl;
		}
		else {
			_gbasis = _H_F_cb ->itemData (_H_F_cb ->currentIndex ()).toInt ();
			if (_gbasis < 200000) {
				_gbasis = _gbasis - 100000;
				ofs << "STO" <<endl;
			}
			else if (_gbasis < 300000) {
					_gbasis = _gbasis - 200000;
					ofs << "N21" <<endl;
				}
				else if (_gbasis < 400000) {
						_gbasis = _gbasis - 300000;
						ofs << "N31" <<endl;
					}
					else {
						_gbasis = _gbasis - 400000;
						ofs << "N311" <<endl;
			}

			ofs << "NGAUSS=" ;
			if (_gbasis < 40000) {
				_gbasis = _gbasis - 30000;
				ofs << "3" <<endl;
			}
			else if (_gbasis < 60000) {
				_gbasis = _gbasis - 40000;
				ofs << "4" << endl;
				}
				else {
				_gbasis = _gbasis - 60000;
				ofs << "6" << endl;
			}

			ofs << "NDFUNC=" ;
			if (_gbasis < 1000) {
				ofs << "0" << endl;
			}
			else if (_gbasis < 2000) {
				_gbasis = _gbasis - 1000;
				ofs << "1" << endl;
			}
			else {
				_gbasis = _gbasis - 2000;
				ofs << "2" << endl;
			}

			ofs << "NPFUNC=" ;
			if (_gbasis < 100) {
				ofs << "0" << endl;
			}
			else {
				_gbasis = _gbasis - 100;
				ofs << "1" << endl;
			}

			if (_gbasis < 10) {
				
			}
			else {
				ofs << "DIFFSP=.true." << endl;
				_gbasis = _gbasis - 10;
			}

			if (_gbasis < 1) {
				
			}
			else {
				ofs << "DIFFS=.true." << endl;
			}
		}
		ofs << " $END\n\n!##################################################################################\n\n"
			" $CONTRL\n"
			"SCFTYP=" ;
			ofs << _scftyp_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << "RUNTYP=" ;
			ofs << _run_type_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << "MULT=" ;
			ofs << _mult_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << "EXETYP=" ;
			ofs << _exetyp_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << " $END\n\n!##################################################################################\n\n";
		if (_solv  == 1) {
			ofs << " $PCM\n"
				"SOLVNT=" ;
				ofs<< _solvent_cb ->currentData ().toString ().toStdString ()<<endl;
			ofs << "TABS=" ;
				ofs<< _tabs_el ->currentData ()<<endl;
//			ofs << "RSOLV=" ;
//				ofs<< _rsolv_el ->currentData ().toString ().toStdString ()<<endl;
			ofs << " $END\n\n!##################################################################################\n\n";
		}
		ofs << " $SYSTEM\n"
			"TIMLIM=" ;
			ofs << _timlim_iel ->spinbox -> value ()<<endl;
//		ofs << "MWORDS=" ;
//			ofs << _mwords_cb ->currentData ().toString ().toStdString ()<<endl;
//		ofs << "MEMDD=" ;
//			ofs << _medd_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << "MEMORY=" ;
			_memory2 = _memory_iel ->spinbox -> value ();
			_memory2 = 1000 * _memory2;
			ofs << _memory2<<endl;
		ofs << " $END\n\n!##################################################################################\n\n" ;
		if (_elpot_box -> isChecked () == true) {
			ofs << " $ELPOT\n"
				"IEPOT=1\n"
				"WHERE=";
				ofs<< _where_cb ->currentData ().toString ().toStdString ()<<endl;
			ofs <<	"OUTPUT=";
				ofs<< _output_cb ->currentData ().toString ().toStdString ()<<endl;
			ofs << " $END\n\n!##################################################################################\n\n";
			if (_where_cb ->currentData () == "points") {
				ofs << " $POINTS\n";

				ofs << " $END\n\n!##################################################################################\n\n";
			}
			else if (_where_cb ->currentData () == "grid") {
				ofs << " $GRID\n";

				ofs << " $END\n\n!##################################################################################\n\n";
			}
			else if (_where_cb ->currentData () == "pdc") {
				ofs << " $PDC\n"
					"PTSEL=";
					ofs<< _ptsel_cb ->currentData ().toString ().toStdString ()<<endl;
				ofs <<	"CONSTR=";
					ofs<< _constr_cb ->currentData ().toString ().toStdString ()<<endl;
				ofs << "VSWSCL=" ;
					ofs<< _vdwscl_fel ->spinbox -> value ()<<endl;
				ofs << "VSWINC=" ;
					ofs<< _vdwinc_fel ->spinbox -> value ()<<endl;
				ofs << "LAYER=" ;
					ofs<< _layer_iel ->spinbox -> value ()<<endl;
				ofs << "MAXPDC=" ;
					ofs<< _maxpdc_iel ->spinbox -> value ()<<endl;
				ofs << " $END\n\n!##################################################################################\n\n";
			}
		}
		ofs << " $GUESS\n"
			"GUESS=" ;
			ofs<< _guess_cb ->currentData ().toString ().toStdString ()<<endl;
		ofs << " $END\n\n!##################################################################################\n\n" ;

		ZNMolecule *mol = _data -> ddwin ->target_molecule;

		OBConversion conv;
		conv.SetOutStream(&ofs);
		OBFormat* outFormat = conv.FormatFromExt(".inp");
		if (outFormat) {
			conv.SetOutFormat (outFormat);
			conv.Write (mol);
		}
	}
		
   }
	
}



///////////////////////////////////////////////////////////////////////////////////////////////////////


IODeviceMenu::IODeviceMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, true) {
//	setWindowTitle ("Input/Output devices");
//	_input_listview  = new QListView ();
//	_output_listview = new QListView ();
//	QWidget *hbox = new QWidget ();
//	hbox ->setLayout (new QHBoxLayout ());
//	main_tab ->layout ()->addWidget (hbox);
//	hbox ->layout () ->addWidget (_input_listview);
//	hbox ->layout () ->addWidget (_output_listview);
}


SequenceMenu::SequenceMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, false, false) {
	setWindowTitle ("Sequence");

	ZNAdvancedWidget *main_tab = new ZNAdvancedWidget ("", false);;
	setCentralWidget (main_tab);

	tab = new QTableWidget();
//	connect (tab, SIGNAL (cellClicked ( int , int )) ,SLOT (clicked_cell (int, int)));
	connect (tab, SIGNAL (itemSelectionChanged ()) ,SLOT (selected_cells ()));

	main_tab ->basic_layout () ->addWidget (tab);

	tab -> setSortingEnabled (false);
//	set_database(db);



//	display_mode_box = new MyGroupBox (main_tab ->basic_layout(), "Display mode");

//	_dsetbutts_gc = new MyGridColumn (display_mode_box ->layout(), 1, 5);

	add_menu ();
//	add_help ();
	connect (_data -> ddwin, SIGNAL (non_selection_molecules_updated ()), this, SLOT (update_tab ()));

}

void SequenceMenu::update_tab () {
	tab ->clear ();
	int c = 0;
	for (unsigned int m = 1; m < _data ->ddwin ->molecules.size () ; m ++) {
		if (!_data ->ddwin ->molecules[m]-> selection) c++;
	}
	tab ->setRowCount(c);
	int res_number = 0;
	int res_number_total = 0;
	int num = 0;
	string res = "";
	char chain;
	for (unsigned int m = 1; m < _data ->ddwin ->molecules.size () ; m ++) {
		ZNMolecule *mol = _data ->ddwin ->molecules [m];
		if (mol -> selection) continue;
		QTableWidgetItem *mol_name_row = new QTableWidgetItem (mol ->GetTitle ());
		tab ->setVerticalHeaderItem (m-1, mol_name_row);
		FOR_RESIDUES_OF_MOL(r, mol)
		{
			res_number_total++;
		}
		column = tab ->columnCount ();
		if (res_number_total > column)	{
			tab ->setColumnCount (res_number_total);
		}
		
		FOR_RESIDUES_OF_MOL(r, mol)
		{
			res_number =  r->GetNum();
			res =  r->GetName();
			chain =  r->GetChain();
			QString residue = res.c_str();
			stringstream ss;
			ss << res_number;
			string residue_number = ss.str();
			
			QTableWidgetItem *field_name = new QTableWidgetItem (residue);
			tab -> setItem (m -1, num++, field_name);
			
		}
		
	}
	connect (tab, SIGNAL (itemClicked ( QTableWidgetItem *)) ,SLOT (pressed_cell (QTableWidgetItem *)));

	for (unsigned int j = 0; j < tab ->columnCount (); j ++) {
		tab ->resizeColumnToContents (j);
	}
	/*
	char chain;
	int res_number = 0;
	int res_number_total = 0;
	string res = "";
	string residue_plus = "";
	FOR_RESIDUES_OF_MOL(r, mol)
	{
		res_number_total++;
	}

	column = tab ->columnCount ();
	if (res_number_total > column)	{
		tab ->setColumnCount (res_number_total);
	}

	row = tab ->rowCount ();
	row = row + 1;
	tab ->setRowCount (row);

	QTableWidgetItem *mol_name_row = new QTableWidgetItem (mol ->GetTitle ());
	tab ->setVerticalHeaderItem (row - 1, mol_name_row);
//	connect (tab, SIGNAL (itemClicked ( QTableWidgetItem *)) ,SLOT (pressed_cell (QTableWidgetItem *)));


	int num =0;
	FOR_RESIDUES_OF_MOL(r, mol)
	{
		res_number =  r->GetNum();
		res =  r->GetName();
		chain =  r->GetChain();
		QString residue = res.c_str();
		stringstream ss;
		ss << res_number;
		string residue_number = ss.str();

		QTableWidgetItem *field_name = new QTableWidgetItem (residue);
		tab -> setItem (row -1, num++, field_name);

	}

	for (unsigned int j = 0; j < tab ->columnCount (); j ++) {
		tab ->resizeColumnToContents (j);
	}
	 */
}


void SequenceMenu::del_from_tab (ZNMolecule *mol) {

}


void SequenceMenu::selected_cells () {

	QList<QTableWidgetItem *> list;
	list = tab ->selectedItems();
	cerr << "size "<< list.size();
	cerr <<endl;

	if (_data ->ddwin ->current_target) {
		ZNMolecule *target = 0;
			_data ->ddwin ->deselect ();
		Selection *sel = new Selection;
		for (unsigned int i = 0; i < list.size(); i++) {
			int column = tab ->column (list.at(i));
			int row = tab ->row (list.at (i));

			target = _data ->ddwin -> molecules[row+1];
			Resid *res = target ->GetResidue (column);
			FOR_ATOMS_OF_RESIDUE (at, res) {
				sel ->select_atom (&*at);
			}
		}
		find_center (sel);
		find_limits (sel);
		_data ->ddwin ->add_molecule (sel);
		_data ->ddwin ->set_current_target (-1);

	}
}


void SequenceMenu::clicked_cell (int x, int y) {

//cerr << "Clicked. x: "<< x << "; y: " << y << endl;
//cerr << mol_name_row -> text().toStdString() << endl;

//cerr << tab ->verticalHeaderItem (x) ->text().toStdString() << endl;

	if (_data ->ddwin ->current_target) {
		ZNMolecule *target = 0;
		target = _data ->ddwin -> molecules[x+1];
		_data ->ddwin ->deselect ();
		Selection *sel = new Selection;
		Resid *res = target ->GetResidue (y);
		FOR_ATOMS_OF_RESIDUE (at, res) {
			sel ->select_atom (&*at);
//			to_add.push_back (&*at);
		}
		find_center (sel);
		find_limits (sel);
		_data ->ddwin ->add_molecule (sel);
		_data ->ddwin ->set_current_target (-1);
	}
}

void SequenceMenu::pressed_cell (QTableWidgetItem *item) {

cerr << "Pressed. item: "<< endl;
/*
	if (_data ->ddwin ->current_target) {
		ZNMolecule *target = 0;
		if (!_data ->ddwin ->target_molecule -> selection)  target = _data ->ddwin ->target_molecule;
		else target = ((Selection *) _data ->ddwin ->target_molecule) ->get_molecules ()[0];
		_data ->ddwin ->deselect ();
		Selection *sel = new Selection;
		FOR_ATOMS_OF_MOL (a, target){
			sel ->select_atom (&*a);
		}
		find_center (sel);
		find_limits (sel);
		_data ->ddwin ->add_molecule (sel);
		_data ->ddwin ->set_current_target (-1);
	}
*/
}


void SequenceMenu::update_graphics (int row, int r, QString residue) {


	QTableWidgetItem *field_name = new QTableWidgetItem (residue);
	tab -> setItem (row -1, r, field_name);

/*
	for (unsigned int j = 0; j < database ->field_names.size (); j ++) {
		QTableWidgetItem *field_name = new QTableWidgetItem (QString (database ->field_names[j].c_str ()));
		tab ->setHorizontalHeaderItem (j, field_name);
	}

	for (unsigned int i=0; i< database ->count_entries (); i++) {
		for (unsigned int j = 0; j < database ->count_fields (); j ++) {
			QTableWidgetItem *value_widget = new QTableWidgetItem ();
			value_widget ->setData (Qt::DisplayRole, QString (database -> entries[i] ->cells[j] ->get_string ().c_str ()));
			tab -> setItem (i, j, value_widget);
//			cerr << "insert " << i << " "<<j<<" "<< database -> entries[i] ->cells[j] ->get_string () << endl;
		}
	}
	database ->set_needs_redraw (false);
	database ->mutex ->unlock ();
*/
}

void SequenceMenu::add_menu () {
    	QMenu *file = new QMenu(tr("&File"), this );
	_save_sequence = new QAction (tr ("&Save sequence"), this);
	connect (_save_sequence, SIGNAL (triggered ()), this, SLOT (_save_sequence_slot ()));
    	file -> addAction (_save_sequence);
	addMenu (file);


    	QMenu *display = new QMenu(tr("&Display residue"), this );
 
	_residue_indice = new QAction (tr ("&Three letters"), this);
	connect (_residue_indice, SIGNAL (triggered ()), this, SLOT (_residue_indice_slot ()));
	_residue_indice ->setCheckable(true);

	_residue_uid = new QAction (tr ("&One letter"), this);
	connect (_residue_uid, SIGNAL (triggered ()), this, SLOT (_residue_uid_slot ()));
	_residue_uid ->setCheckable(true);


	aminoacidGroup = new QActionGroup(this);
	aminoacidGroup ->setExclusive(true);
	aminoacidGroup ->addAction (_residue_indice);
	aminoacidGroup ->addAction (_residue_uid);
	_residue_indice ->setChecked(true);

    	display ->addAction (_residue_indice);
    	display ->addAction (_residue_uid);

	addMenu (display);

}


void SequenceMenu::_save_sequence_slot () {
//    QString s = QFileDialog::getSaveFileName(this, tr ("Save As"), "",tr("Fasta (*.fasta)"));

	QFileDialog save_file(this);

	save_file.setNameFilter (tr("Fasta (*.fasta)"));
	save_file.setAcceptMode (QFileDialog::AcceptSave);
	save_file.exec();

	QString ext = save_file.selectedNameFilter();

	QString filter = save_file.selectedNameFilter();
	filter.truncate (filter.lastIndexOf (')'));
	filter.remove (0, filter.indexOf ('*')+1);
	QString s = save_file.selectedFile();
	if (!s.contains ('.') ) s+= filter;



	
    if (!s.isNull()) {
        if (_data ->ddwin ->current_target) {
		_data -> actions -> save_as (_data ->ddwin ->target_molecule, s.toStdString ());
		}
	}
}


void SequenceMenu::_residue_indice_slot () {

	for (unsigned int i = 0; i < tab ->rowCount (); i++) {

		for (unsigned int j = 0; j < tab ->columnCount (); j++) {

			QTableWidgetItem *item = tab ->item(i, j);
			if (item) {		

			QString label = item ->text();

			if (label == "G")	tab ->item(i, j) ->setText("GLY");
			else if (label == "P")	tab ->item(i, j) ->setText("PRO");
			else if (label == "A")	tab ->item(i, j) ->setText("ALA");
			else if (label == "V")	tab ->item(i, j) ->setText("VAL");
			else if (label == "L")	tab ->item(i, j) ->setText("LEU");
			else if (label == "I")	tab ->item(i, j) ->setText("ILE");
			else if (label == "M")	tab ->item(i, j) ->setText("MET");
			else if (label == "C")	tab ->item(i, j) ->setText("CYS");
			else if (label == "F")	tab ->item(i, j) ->setText("PHE");
			else if (label == "Y")	tab ->item(i, j) ->setText("TYR");
			else if (label == "W")	tab ->item(i, j) ->setText("TRP");
			else if (label == "H")	tab ->item(i, j) ->setText("HIS");
			else if (label == "K")	tab ->item(i, j) ->setText("LYS");
			else if (label == "R")	tab ->item(i, j) ->setText("ARG");
			else if (label == "Q")	tab ->item(i, j) ->setText("GLN");
			else if (label == "N")	tab ->item(i, j) ->setText("ASN");
			else if (label == "E")	tab ->item(i, j) ->setText("GLU");
			else if (label == "D")	tab ->item(i, j) ->setText("ASP");
			else if (label == "S")	tab ->item(i, j) ->setText("SER");
			else if (label == "T")	tab ->item(i, j) ->setText("THR");

			}
		}
	}
	for (unsigned int j = 0; j < tab ->columnCount (); j ++) {
		tab ->resizeColumnToContents (j);
	}
}

void SequenceMenu::_residue_uid_slot () {

	for (unsigned int i = 0; i < tab ->rowCount (); i++) {

		for (unsigned int j = 0; j < tab ->columnCount (); j++) {		

			QTableWidgetItem *item = tab ->item(i, j);
			if (item) {

			QString label = item ->text();

			if (label == "GLY")	tab ->item(i, j) ->setText("G");
			else if (label == "PRO")	tab ->item(i, j) ->setText("P");
			else if (label == "ALA")	tab ->item(i, j) ->setText("A");
			else if (label == "VAL")	tab ->item(i, j) ->setText("V");
			else if (label == "LEU")	tab ->item(i, j) ->setText("L");
			else if (label == "ILE")	tab ->item(i, j) ->setText("I");
			else if (label == "MET")	tab ->item(i, j) ->setText("M");
			else if (label == "CYS")	tab ->item(i, j) ->setText("C");
			else if (label == "PHE")	tab ->item(i, j) ->setText("F");
			else if (label == "TYR")	tab ->item(i, j) ->setText("Y");
			else if (label == "TRP")	tab ->item(i, j) ->setText("W");
			else if (label == "HIS")	tab ->item(i, j) ->setText("H");
			else if (label == "LYS")	tab ->item(i, j) ->setText("K");
			else if (label == "ARG")	tab ->item(i, j) ->setText("R");
			else if (label == "GLN")	tab ->item(i, j) ->setText("Q");
			else if (label == "ASN")	tab ->item(i, j) ->setText("N");
			else if (label == "GLU")	tab ->item(i, j) ->setText("E");
			else if (label == "ASP")	tab ->item(i, j) ->setText("D");
			else if (label == "SER")	tab ->item(i, j) ->setText("S");
			else if (label == "THR")	tab ->item(i, j) ->setText("T");
}
		}
	}
	for (unsigned int j = 0; j < tab ->columnCount (); j ++) {
		tab ->resizeColumnToContents (j);
	}
}

void SequenceMenu::add_help () {
    	QMenu *about = new QMenu(tr("&Help"), this );
    	Q_CHECK_PTR( about );
    	about -> addAction (_about_action);
	menuBar () ->addMenu (about);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////

DisplayMenu::DisplayMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, false, true) {
	setWindowTitle ("Display settings");

	ZNAdvancedWidget *main_tab = (ZNAdvancedWidget *) _main_widget;

	display_mode_box = new MyGroupBox (main_tab ->basic_layout(), "Display mode");

	_dsetbutts_gc = new MyGridColumn (display_mode_box ->layout(), 1, 5);

	dset_wireframe = new QPushButton("Wireframe");
	dset_stick = new QPushButton("Sticks");
	dset_cpk = new QPushButton("CPK");
	dset_ball_and_stick = new QPushButton("Ball and Stick");
	dset_ball_and_line = new QPushButton("Ball and Line");

	_dsetbutts_gc ->gridlayout ->addWidget (dset_wireframe, 1, 1);
	_dsetbutts_gc ->gridlayout ->addWidget (dset_stick, 1, 2);
	_dsetbutts_gc ->gridlayout ->addWidget (dset_cpk, 1, 3);
	_dsetbutts_gc ->gridlayout ->addWidget (dset_ball_and_stick, 1, 4);
	_dsetbutts_gc ->gridlayout ->addWidget (dset_ball_and_line, 1, 5);


	connect (dset_wireframe, SIGNAL( clicked() ), SLOT( set_wireframe() ) );
	connect (dset_stick, SIGNAL( clicked() ), SLOT( set_stick() ) );
	connect (dset_cpk, SIGNAL( clicked() ), SLOT( set_cpk() ) );
	connect (dset_ball_and_stick, SIGNAL( clicked() ), SLOT( set_ballandstick() ) );
	connect (dset_ball_and_line, SIGNAL( clicked() ), SLOT( set_ballandline() ) );


	disp_gc = new MyGridColumn (display_mode_box ->layout(), 1, 2);

	at_disp = new QComboBox ();
	disp_gc ->gridlayout ->addWidget (at_disp, 1, 1);
	at_disp->insertItem(0, "No Atoms" );
	at_disp->insertItem(1, "Spheres" );
	at_disp->insertItem(2, "CPK Spheres" );
	at_disp->insertItem(3, "Scaled CPK Spheres" );



	bo_disp = new QComboBox();
	disp_gc ->gridlayout ->addWidget (bo_disp, 1, 2);
	bo_disp->insertItem(0, "No Bonds" );
	bo_disp->insertItem(1, "Lines" );
	bo_disp->insertItem(2, "Sticks" );
	bo_disp -> setCurrentIndex (1);

	backbone_disp = new MyComboBox (display_mode_box ->layout(), "Backbone mode");
	backbone_disp ->_combo_box ->insertItem(0, "None" );
	backbone_disp ->_combo_box ->insertItem(1, "Line" );
	backbone_disp ->_combo_box ->insertItem(2, "Tube" );
	backbone_disp ->_combo_box ->setCurrentIndex (0);


	QWidget *hwid = new QWidget;
	hwid -> setLayout (new QHBoxLayout ());
	hwid ->layout () ->setAlignment (Qt::AlignHCenter);
	main_tab ->basic_layout() ->addWidget (hwid);
	ok_p = new MyPushButton (hwid ->layout (), "Ok", 0, 100);
	connect (ok_p ->fbutton, SIGNAL( clicked() ), SLOT( disp_ok() ) );







/*
	ar_disp = new QComboBox();
	disp_gc ->gridlayout ->addWidget (ar_disp, 2, 1);
	ar_disp->insertItem(0, "Aromatic rings" );
	ar_disp->insertItem(1, "Kekule resonance structures" );
//	ar_disp->insertItem(2, "Aromatic bonds" );
*/




}

void DisplayMenu::set_wireframe (){
    set_conf (WIREFRAME);
}

void DisplayMenu::set_stick (){
    set_conf (STICK);
}

void DisplayMenu::set_cpk (){
    set_conf (CPK);
}

void DisplayMenu::set_ballandline (){
    set_conf (BALLANDLINE);
}

void DisplayMenu::set_ballandstick (){
    set_conf (BALLANDSTICK);
}

void DisplayMenu::set_conf (int conf){

    int atom_mode, bond_mode;


	switch (conf) {
		case WIREFRAME:
            atom_mode = NO_ATOMS;
            bond_mode = LINES;   
            
			break;
		case STICK:
            atom_mode = NO_ATOMS;
            bond_mode = STICKS;  
           
            break;
		case CPK: 
            atom_mode = CPK_SPHERES;
            bond_mode = NO_BONDS; 
            break;
		case BALLANDSTICK: 
            atom_mode = SCALED_CPK_SPHERES;
            bond_mode = STICKS;            
            break;

		case BALLANDLINE: 
            atom_mode = SPHERES;
            bond_mode = LINES;            
            break;
			
		default:
            atom_mode = NO_ATOMS;
            bond_mode = LINES; 

	}

    at_disp->setCurrentIndex (atom_mode);
    bo_disp->setCurrentIndex (bond_mode);

}

void DisplayMenu::disp_ok () {
    for (unsigned int i=0;i<at_opts.size();i++) at_opts[i]->set ();
	
//	_data -> ddwin ->gl ->aromatic_display_style = style_str_to_i (ar_disp->currentText().toStdString());
	int at_st = style_str_to_i (at_disp->currentText().toStdString());
	int bo_st = style_str_to_i (bo_disp->currentText().toStdString());
	int backbone_st = backbone_disp ->_combo_box ->currentIndex();
	
    if (_data -> ddwin ->current_target) {
			bool ext = false;
			if (_data -> ddwin ->target_molecule -> multi) {
				Database_molecule *dm = (Database_molecule *) _data -> ddwin ->target_molecule;
				if (dm -> database -> has_extend_enabled ()) ext = true;
			}
			if (!ext) {
				_data -> actions -> change_display_style (_data -> ddwin ->target_molecule, at_st, bo_st, backbone_st);
			}
			else {
				Database_molecule *dm = (Database_molecule *) _data -> ddwin ->target_molecule;
				Database *dat = dm -> database;
				_data -> actions -> change_display_style (dat, at_st, bo_st, backbone_st);
			}

		
	}
}


int DisplayMenu::style_str_to_i (string style){


    if (style=="No Atoms") return NO_ATOMS;
    else if (style=="Spheres") return SPHERES;
    else if (style=="CPK Spheres") return CPK_SPHERES;
    else if (style=="Scaled CPK Spheres") return SCALED_CPK_SPHERES;

    else if (style=="No Bonds") return NO_BONDS;
    else if (style=="Lines") return LINES;
    else if (style=="Sticks") return STICKS;

    else if (style=="Aromatic rings") return AROMATIC_RINGS;
    else if (style=="Kekule resonance structures") return KEKULE;
//    else if (style=="Aromatic bonds") return AROMATIC_BONDS;
    return 0;


}


void DDWin::set_popups (){

/*
    dsetpopup = new QWidget ();
    dsetpopup->setWindowTitle ("Display Settings");


//    dsetpopup->setFrameStyle( Q3Frame::WinPanel|Q3Frame::Raised );
    QVBoxLayout *dsetpopupv = new QVBoxLayout (dsetpopup);
     dsetpopup->setMinimumSize (500, 250);   
  //


    QHBoxLayout *dsetbutts = new QHBoxLayout ();
    dsetpopupv -> addLayout (dsetbutts);
    QPushButton *dset_wireframe = new QPushButton("Wireframe");
    QPushButton *dset_stick = new QPushButton("Sticks");
    QPushButton *dset_cpk = new QPushButton("CPK");
    QPushButton *dset_ball_and_stick = new QPushButton("Ball and Stick");
    QPushButton *dset_ball_and_line = new QPushButton("Ball and Line");

    dsetbutts -> addWidget (dset_wireframe);
    dsetbutts -> addWidget (dset_stick);
    dsetbutts -> addWidget (dset_cpk);
    dsetbutts -> addWidget (dset_ball_and_stick);
    dsetbutts -> addWidget (dset_ball_and_line);




    connect (dset_wireframe, SIGNAL( clicked() ), gl, SLOT( set_wireframe() ) );
    connect (dset_stick, SIGNAL( clicked() ), gl, SLOT( set_stick() ) );
    connect (dset_cpk, SIGNAL( clicked() ), gl, SLOT( set_cpk() ) );
    connect (dset_ball_and_stick, SIGNAL( clicked() ), gl, SLOT( set_ballandstick() ) );
    connect (dset_ball_and_line, SIGNAL( clicked() ), gl, SLOT( set_ballandline() ) );

    QHBoxLayout *mainhbox = new QHBoxLayout ();
    dsetpopupv -> addLayout (mainhbox);


    ar_disp = new QComboBox();
    ar_disp->insertItem(0, "Aromatic rings" );
    ar_disp->insertItem(1, "Kekule resonance structures" );

    ar_disp -> setCurrentIndex (gl -> aromatic_display_style);
//    ar_disp->insertItem(2, "Aromatic bonds" );
    
    dsetpopupv -> addWidget (ar_disp);

	
	backbone_disp = new QComboBox();
    backbone_disp->insertItem(0, "None" );
    backbone_disp->insertItem(1, "Line" );
    backbone_disp->insertItem(2, "Tube" );
	
    backbone_disp -> setCurrentIndex (0);
	//    ar_disp->insertItem(2, "Aromatic bonds" );
    
    dsetpopupv -> addWidget (backbone_disp);
	

    QHBoxLayout *butt2   = new QHBoxLayout ();
    dsetpopupv -> addLayout (butt2);
    QVBoxLayout *atomsvb = new QVBoxLayout ();
    mainhbox -> addLayout (atomsvb);
    QVBoxLayout *bondsvb = new QVBoxLayout ();
    mainhbox -> addLayout (bondsvb);

    at_disp = new QComboBox ();
    atomsvb -> addWidget (at_disp);
    at_disp->insertItem(0, "No Atoms" );
    at_disp->insertItem(1, "Spheres" );
    at_disp->insertItem(2, "CPK Spheres" );
    at_disp->insertItem(3, "Scaled CPK Spheres" );



    bo_disp = new QComboBox();
    bondsvb -> addWidget (bo_disp);
    bo_disp->insertItem(0, "No Bonds" );
    bo_disp->insertItem(1, "Lines" );
    bo_disp->insertItem(2, "Sticks" );
    bo_disp -> setCurrentIndex (1);


    at_opts.push_back (new MyFloatEditLine (atomsvb, "Sphere radius", gl->sphere_radius));
    at_opts.push_back (new MyFloatEditLine (bondsvb, "Stick radius", gl->stick_rad));
    at_opts.push_back (new MyFloatEditLine (bondsvb, "Double ZNBond inter Distance", gl->double_bond_inter_distance));
    at_opts.push_back (new MyFloatEditLine (bondsvb, "Aromatic ZNBond inter Distance", gl->aromatic_bond_inter_distance));
    at_opts.push_back (new MyFloatEditLine (atomsvb, "VdW scale", gl->vdw_scale));
//    at_opts.push_back (new MyFloatEditLine (atomsvb, "VdW Precision", gl->vdw_precision));
    at_opts.push_back (new MyFloatEditLine (bondsvb, "Double bond scale", gl->double_bond_stick_radius_scale));
    
    QPushButton *ok = new QPushButton("Ok");
    butt2 -> addWidget (ok);
    connect (ok, SIGNAL( clicked() ), SLOT( disp_ok() ) );

*/


	b_color_menu = new BackboneColorMenu (0, data);
	go_color_menu = new GraphicalObjectsColorMenu (0, data);
    color_menu = new ColorMenu (0, this);
	color_settings_menu = new ColorSettingsMenu (0, this);
    builder_menu = new BuilderMenu (0, data, builder);

    pref_menu = new PrefMenu (0, data);

    haptic_menu = new HapticMenu (0, data);
    gamess_menu = new GamessMenu (0, data);
    plants_menu = new PLANTSMenu (0, data);
    display_menu = new DisplayMenu (0, data);
    sequence_menu = new SequenceMenu (0, data);
	docking_menu = new DockingMenu (0, data);
	thread_menu = new ThreadMenu (0, data);
	iodevice_menu = new IODeviceMenu (0, data);
    clicked_atom_menu = new Clicked_atomMenu (0, data);
    surface_menu = new SurfaceMenu (0, data);
    sphere_menu = new SphereMenu (0, data);
	conformers_menu = new ConformersMenu (0, data);
    map_menu = new MapMenu (0, data);
    graphical_objects_menu = new GraphicalObjectsMenu (0, this);
    DDsettings_menu = new DDSettingsMenu (0, this);


}



void SurfaceMenu::draw_surface () {

    if (data -> ddwin -> current_target) {
        ZNMolecule *mol = data -> ddwin -> target_molecule;
        mesh = (gtype->currentIndex ()==1);
        float a = ((float) alpha)/100;
        surface = new Surface ();
		surface -> near_to_dist = near_to_f;
		ZNMolecule *neartm = NULL;
		if (near_to ->currentIndex ()!= 0) neartm = data ->ddwin ->molecules[near_to ->currentIndex ()];
		else neartm = 0;
        surface -> set_molecule (mol, neartm);
        surface -> lst = data -> ddwin -> gl -> new_list ();
        surface -> name = string ("Surface ") + mol -> GetTitle ();
        surface -> mesh = mesh;


        SurfaceThread *thread = new SurfaceThread (0, surface, data -> ddwin);

        thread -> alph = a;
        thread -> res = res;
		data ->ddwin ->run_thread (thread);
        //thread -> start ();


        connect (thread, SIGNAL (finished ()), this, SLOT (add_surface ()));

    }
}



void SurfaceMenu::add_surface () {
    surface -> render ();
    CreateGraphicalObjectCommand *command = new CreateGraphicalObjectCommand (surface, data -> ddwin);
    data -> ddwin -> execute (command);
	data ->ddwin ->go_color_menu ->set_target(surface);
	data ->ddwin ->go_color_menu ->display ();
}



SurfaceMenu::SurfaceMenu (QWidget *parent, Data *dat )
   :    QWidget(parent)
{
	near_to_f = 4.f;
    data = dat;
  //  molecule = mol;
    res = data->ddwin->gl->surface_resolution;


    QVBoxLayout *vbox = new QVBoxLayout (this);
  //  this->setMinimumSize (500, 250);   
 //   this->setMaximumSize (500, 250);
 //   this->setMinimumSize (500, 250);   
 //   this->setMaximumSize (500, 250);
    this->setWindowTitle("Surfaces");

    resolution = new MyFloatEditLine (vbox, "Resolution", res);

    stype = new QComboBox();
    vbox -> addWidget (stype);
    stype->insertItem(0, "Connolly" );
//    stype->insertItem(1, "Gaussian contact" );


    gtype = new QComboBox();
    vbox -> addWidget (gtype);
    gtype->insertItem(0, "Surface" );
    gtype->insertItem(1, "Mesh" );
	
	near_to_dle = new MyFloatEditLine (vbox, "Near to", near_to_f);
	near_to = new QComboBox ();
	vbox -> addWidget (near_to);

//    QSlider *slider = new QSlider( Horizontal, this, "slider" );
 //   alpha_p = new MyFloatEditLine (vbox, "Opacity", alpha);
//    connect (slider, SIGNAL(valueChanged(int)), alpha_p, SLOT(set_value(int)) );
//    connect ();
    alpha_s = new MySlider (vbox, "Opacity",alpha, 0, 100);
    alpha_s->slider->setValue (100);
    QPushButton *ok = new QPushButton ("Ok");
    vbox -> addWidget (ok);
    connect (ok, SIGNAL (clicked ()) ,SLOT (draw_surface ()));
	connect (data -> ddwin, SIGNAL (targets_updated ()), this, SLOT (update_near_to ()));
	//update_near_to ();

}



void SurfaceMenu::update_near_to () {
	near_to ->clear ();
	for (unsigned int i = 0; i < data -> ddwin -> target ->count (); i++) {
		near_to ->insertItem (i, data -> ddwin -> target ->itemText (i));
	}
}

///////////////////////////////////////////////////////////////////////


MapMenu::MapMenu (QWidget *parent, Data *dat) : ZNMenu (parent, dat, true, true) {
	_r = 10.; _resolution = 1.5; _threshold = 0.;
	
	ZNAdvancedWidget *main_tab = new ZNAdvancedWidget ("", false);
   	_tabs ->addTab (main_tab, "Map");
	add_menu ();
	add_help ();
	_site_box = new MyGroupBox (main_tab ->layout (), "Binding site options");
	_site_x_el = new MyFloatEditLine (_site_box ->layout (), "Site center X:", _x, -1000, 1000);
	_site_y_el = new MyFloatEditLine (_site_box ->layout (), "Site center Y:", _y, -1000, 1000);
	_site_z_el = new MyFloatEditLine (_site_box ->layout (), "site center Z:", _z, -1000, 1000);
	_site_r_el = new MyFloatEditLine (_site_box ->layout (), "Site radius:", _r, 0, 100);
	_threshold_el = new MyFloatEditLine (_site_box ->layout (), "threshold:", _threshold, -1000, 1000);
	_resolution_el = new MyFloatEditLine (_site_box ->layout (), "resolution:", _resolution, 0, 100);
	_name_el = new MyLineEdit (_site_box ->layout (), "Name: ");
	_compute_map_bt = new MyPushButton (main_tab ->layout (), "Compute Map");
	_type_cb = new MyComboBox (_site_box ->layout (), "Potential:");
	_type_cb ->combo_box () ->insertItem (-1, "Chemscore Lipophilic potential", 0);
	_type_cb ->combo_box () ->insertItem (1, "Chemscore HB acceptor potential", 1);
	_type_cb ->combo_box () ->insertItem (2, "Chemscore HB donor potential", 2);
//	_type_cb ->combo_box () ->insertItem (3, "Electrostatic Potential", 3);
	
    MyCompleteColorSettings *_solid_color_button = new MyCompleteColorSettings (main_tab ->layout (), _solid_color);
	connect (_compute_map_bt ->fbutton, SIGNAL( clicked() ), SLOT( compute_map_slot() ) );
}

bool MapMenu::display_requirements_met () {
	bool b = check_for_a_set_target ();
	if (b) {
		ZNMolecule *mol = _data ->ddwin ->target_molecule;
		vect c = get_center (mol);
		_site_x_el ->set_value (c.x ());
		_site_y_el ->set_value (c.y ());
		_site_z_el ->set_value (c.z ());
		string q = string ("Map ") + string (_data ->ddwin ->target_molecule ->GetTitle ());
		_name_el   ->linedit ->setText(q.c_str ());
	}
	return b;
};

void MapMenu::compute_map_slot () {
		string name = _name_el ->linedit ->text ().toStdString ();
		bool found = false;
		for (unsigned int i = 0; i < _data ->ddwin ->graphical_objects.size (); i++) {
			if (_data ->ddwin ->graphical_objects[i] -> is_map () && _data ->ddwin ->graphical_objects[i] -> name == name) {
				map = (Map *) _data ->ddwin ->graphical_objects[i];
				found = true;
				break;
			}
		}
		if (!found) {	
			map = new Map ();
			map ->solid_color = _solid_color;
			map -> lst = _data -> ddwin -> gl -> new_list ();
			map ->molecule = _data ->ddwin ->target_molecule;
			map ->site_center = vect (_x, _y, _z);
			map ->site_radius = _r;	
			map ->resolution = _resolution;
			map ->name = name;	
			map ->type = _type_cb -> currentData ().toInt ();
	
		}

 //       map -> name = string ("Map");		

			map ->threshold = _threshold;


        MapThread *thread = new MapThread (0, map, _data -> ddwin);

		_data ->ddwin ->run_thread (thread);

        connect (thread, SIGNAL (finished ()), this, SLOT (add_map ()));
}


void MapMenu::add_map () {
    map -> render ();
	bool found = false;
	for (unsigned int i = 0; i < _data ->ddwin ->graphical_objects.size (); i ++) {
		if (map == (Map *) _data ->ddwin ->graphical_objects[i]) {
			found = true;
			break;
		}
	}
	if (!found) {

		CreateGraphicalObjectCommand *command = new CreateGraphicalObjectCommand (map, _data -> ddwin);
		_data -> ddwin -> execute (command);
	}

}

///////////////////////////////////////////////////////////////////////


SphereMenu::SphereMenu (QWidget *parent, Data *dat )
   :    QWidget(parent)
{
    data = dat;
  //  molecule = mol;

    x = 0; y=0, z=0;
    QVBoxLayout *vbox = new QVBoxLayout (this);

    setWindowTitle("Spheres");

    cent_x = new MyFloatEditLine (vbox, "Center x", x);
    cent_y = new MyFloatEditLine (vbox, "Center y", y);
    cent_z = new MyFloatEditLine (vbox, "Center z", z);
    rad = new MyFloatEditLine (vbox, "Radius", radius);
    alpha_s = new MySlider (vbox, "Opacity",alpha, 0, 100);
    alpha_s->slider->setValue (100);
    MyColorButton *colorbutt = new MyColorButton (vbox, col);


    QPushButton *ok = new QPushButton ("Ok");
    vbox -> addWidget (ok);
    connect (ok, SIGNAL (clicked ()) ,SLOT (draw_sphere ()));
}



void SphereMenu::draw_sphere () {

    float a = (float) alpha/100;
    col.setAlphaF (a);
    Sphere *sphere = new Sphere ();
    sphere -> lst = data -> ddwin -> gl -> new_list ();
    sphere -> name = string ("Sphere");
    sphere -> set_center (vect (x, y, z));
    sphere -> set_radius (radius);
    sphere -> set_color (col);
    sphere -> render ();
    CreateGraphicalObjectCommand *command = new CreateGraphicalObjectCommand (sphere, data -> ddwin);
    data -> ddwin -> execute (command);
}


///////////////////////////////////////////////////////////////////////




/*
void SphereMenu::draw_sphere () {

    float a = alpha/100;
    col.setAlphaF (a);
    Sphere *sphere = new Sphere ();
    sphere -> list = data -> ddwin -> gl -> new_list ();
    sphere -> name = string ("Sphere");
    sphere -> set_center (center);
    sphere -> set_radius (radius);
    sphere -> set_color (col);
    sphere -> render ();
    CreateGraphicalObjectCommand *command = new CreateGraphicalObjectCommand (sphere, data -> ddwin);
    data -> ddwin -> execute (command);
}
*/

///////////////////////////////////////////////////////////////////////

BuilderMenu::BuilderMenu (QWidget *parent, Data *dat, Builder *build) : ZNMenu (parent, dat, true, true) {
	setWindowTitle ("Builder");

	builder = build;

	ZNAdvancedWidget *main_tab = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (main_tab, "Basic");

	_basic_gc = new MyGridColumn (main_tab ->basic_layout (), 3, 3);

	_atoms_box = new MyGroupBox ("Atoms");
	_basic_gc ->gridlayout ->addWidget (_atoms_box, 1, 1);

	_atoms_gc = new MyGridColumn (_atoms_box ->layout (), 1, 10);

	_C = new MyPushButton ("C", 0, 30);
	connect (_C ->fbutton, SIGNAL (clicked ()) ,SLOT (add_C ()));
	_atoms_gc ->gridlayout ->addWidget (_C, 1, 1);

	_N = new MyPushButton ("N", 0, 30);
	connect (_N ->fbutton, SIGNAL (clicked ()) ,SLOT (add_N ()));
	_atoms_gc ->gridlayout ->addWidget (_N, 1, 2);

	_O = new MyPushButton ("O", 0, 30);
	connect (_O ->fbutton, SIGNAL (clicked ()) ,SLOT (add_O ()));
	_atoms_gc ->gridlayout ->addWidget (_O, 1, 3);

	_F = new MyPushButton ("F", 0, 30);
	connect (_F ->fbutton, SIGNAL (clicked ()) ,SLOT (add_F ()));
	_atoms_gc ->gridlayout ->addWidget (_F, 2, 1);

//	_H = new MyPushButton ("H", 0, 30);
//	connect (_H ->fbutton, SIGNAL (clicked ()) ,SLOT (add_H ()));
//	_atoms_gc ->gridlayout ->addWidget (_H, 1, 5);

	_P = new MyPushButton ("P", 0, 30);
	connect (_P ->fbutton, SIGNAL (clicked ()) ,SLOT (add_P ()));
	_atoms_gc ->gridlayout ->addWidget (_P, 2, 2);

	_S = new MyPushButton ("S", 0, 30);
	connect (_S ->fbutton, SIGNAL (clicked ()) ,SLOT (add_S ()));
	_atoms_gc ->gridlayout ->addWidget (_S, 2, 3);

	_Cl = new MyPushButton ("Cl", 0, 30);
	connect (_Cl ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cl ()));
	_atoms_gc ->gridlayout ->addWidget (_Cl, 3, 1);

	_I = new MyPushButton ("I", 0, 30);
	connect (_I ->fbutton, SIGNAL (clicked ()) ,SLOT (add_S ()));
	_atoms_gc ->gridlayout ->addWidget (_I, 3, 2);

	_Br = new MyPushButton ("Br", 0, 30);
	connect (_Br ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Br ()));
	_atoms_gc ->gridlayout ->addWidget (_Br, 3, 3);


	_bonds_box = new MyGroupBox ("Bonds");
	_basic_gc ->gridlayout ->addWidget (_bonds_box, 1, 2);

	_bonds_gc = new MyGridColumn (_bonds_box ->layout (), 4, 1);

	_single_bond_b = new MyPushButton ("Single bond", 0, 120);
	connect (_single_bond_b ->fbutton, SIGNAL (clicked ()) ,SLOT (single_bond ()));
	_bonds_gc ->gridlayout ->addWidget (_single_bond_b, 1, 1);

	_double_bond_b = new MyPushButton ("Double bond", 0, 120);
	connect (_double_bond_b ->fbutton, SIGNAL (clicked ()) ,SLOT (double_bond ()));
	_bonds_gc ->gridlayout ->addWidget (_double_bond_b, 2, 1);

	_triple_bond_b = new MyPushButton ("Triple bond", 0, 120);
	connect (_triple_bond_b ->fbutton, SIGNAL (clicked ()) ,SLOT (triple_bond ()));
	_bonds_gc ->gridlayout ->addWidget (_triple_bond_b, 3, 1);

//	_no_bond_b = new MyPushButton ("Delete bond", 0, 100);
//	connect (_no_bond_b ->fbutton, SIGNAL (clicked ()) ,SLOT (no_bond ()));
//	_bonds_gc ->gridlayout ->addWidget (_no_bond_b, 4, 1);


	_rings_box = new MyGroupBox ("Rings");
	_basic_gc ->gridlayout ->addWidget (_rings_box, 1, 3);

	_rings_gc = new MyGridColumn (_rings_box ->layout (), 3, 3);

	_ring3_b = new MyPushButton ("", 0, 30);
	connect (_ring3_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_ring3 ()));
	_rings_gc ->gridlayout ->addWidget (_ring3_b, 1, 1);
	_ring3_b ->fbutton ->setIcon (QIcon (":icons/ring3.png"));

	_ring4_b = new MyPushButton ("", 0, 30);
	connect (_ring4_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_ring4 ()));
	_rings_gc ->gridlayout ->addWidget (_ring4_b, 1, 2);
	_ring4_b ->fbutton ->setIcon (QIcon (":icons/ring4.png"));

	_ring5_b = new MyPushButton ("", 0, 30);
	connect (_ring5_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_ring5 ()));
	_rings_gc ->gridlayout ->addWidget (_ring5_b, 1, 3);
	_ring5_b ->fbutton ->setIcon (QIcon (":icons/ring5.png"));

	_ring6_b = new MyPushButton ("", 0, 30);
	connect (_ring6_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_ring6 ()));
	_rings_gc ->gridlayout ->addWidget (_ring6_b, 2, 1);
	_ring6_b ->fbutton ->setIcon (QIcon (":icons/ring6.png"));

	_ring7_b = new MyPushButton ("", 0, 30);
	connect (_ring7_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_ring7 ()));
	_rings_gc ->gridlayout ->addWidget (_ring7_b, 2, 2);
	_ring7_b ->fbutton ->setIcon (QIcon (":icons/ring7.png"));

	_ring8_b = new MyPushButton ("", 0, 30);
	connect (_ring8_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_ring8 ()));
	_rings_gc ->gridlayout ->addWidget (_ring8_b, 2, 3);
	_ring8_b ->fbutton ->setIcon (QIcon (":icons/ring8.png"));

	_furanO_b = new MyPushButton ("", 0, 30);
	connect (_furanO_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_furanO ()));
	_rings_gc ->gridlayout ->addWidget (_furanO_b, 3, 1);
	_furanO_b ->fbutton ->setIcon (QIcon (":icons/furanO.png"));

	_furan_b = new MyPushButton ("", 0, 30);
	connect (_furan_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_furan ()));
	_rings_gc ->gridlayout ->addWidget (_furan_b, 3, 2);
	_furan_b ->fbutton ->setIcon (QIcon (":icons/furan.png"));

	_benzene_b = new MyPushButton ("", 0, 30);
	connect (_benzene_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_benzene ()));
	_rings_gc ->gridlayout ->addWidget (_benzene_b, 3, 3);
	_benzene_b ->fbutton ->setIcon (QIcon (":icons/benzene.png"));


	_fragments_box = new MyGroupBox ("Fragments");
	_basic_gc ->gridlayout ->addWidget (_fragments_box, 1, 4);

	_fragments_gc = new MyGridColumn (_fragments_box ->layout (), 3, 3);

	_CO_b = new MyPushButton ("C=O", 0, 40);
	connect (_CO_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_CO ()));
	_fragments_gc ->gridlayout ->addWidget (_CO_b, 1, 1);
//	_CO_b ->fbutton ->setIcon (QIcon (":icons/CO.png"));

	_NCO_b = new MyPushButton ("NC=O", 0, 40);
	connect (_NCO_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_NCO ()));
	_fragments_gc ->gridlayout ->addWidget (_NCO_b, 1, 2);
//	_NCO_b ->fbutton ->setIcon (QIcon (":icons/NCO.png"));

	_COOH_b = new MyPushButton ("COOH", 0, 40);
	connect (_COOH_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_COOH ()));
	_fragments_gc ->gridlayout ->addWidget (_COOH_b, 1, 3);
//	_COOH_b ->fbutton ->setIcon (QIcon (":icons/COOH.png"));

	_CCd_b = new MyPushButton ("C=C", 0, 40);
	connect (_CCd_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_CCd ()));
	_fragments_gc ->gridlayout ->addWidget (_CCd_b, 2, 1);
//	_CCd_b ->fbutton ->setIcon (QIcon (":icons/CCd.png"));

	_CCt_b = new MyPushButton ("C#C", 0, 40);
	connect (_CCt_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_CCt ()));
	_fragments_gc ->gridlayout ->addWidget (_CCt_b, 2, 2);
//	_CCt_b ->fbutton ->setIcon (QIcon (":icons/CCt.png"));

	_CN_b = new MyPushButton ("C#N", 0, 40);
	connect (_CN_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_CN ()));
	_fragments_gc ->gridlayout ->addWidget (_CN_b, 2, 3);
//	_CN_b ->fbutton ->setIcon (QIcon (":icons/CN.png"));

	_NO2_b = new MyPushButton ("NO2", 0, 40);
	connect (_NO2_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_NO2 ()));
	_fragments_gc ->gridlayout ->addWidget (_NO2_b, 3, 1);
//	_NO2_b ->fbutton ->setIcon (QIcon (":icons/NO2.png"));

	_SO2_b = new MyPushButton ("SO2", 0, 40);
	connect (_SO2_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_SO2 ()));
	_fragments_gc ->gridlayout ->addWidget (_SO2_b, 3, 2);
//	_SO2_b ->fbutton ->setIcon (QIcon (":icons/SO2.png"));

	_PO3_b = new MyPushButton ("PO3", 0, 40);
	connect (_PO3_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_PO3 ()));
	_fragments_gc ->gridlayout ->addWidget (_PO3_b, 3, 3);
//	_PO3_b ->fbutton ->setIcon (QIcon (":icons/PO3.png"));


	_smiles_box = new MyGroupBox (main_tab ->basic_layout (), "Smiles");

	_smiles_tc = new MyTwoColumn (_smiles_box ->layout () );

	smiles = new QLineEdit( this);
	_smiles_tc ->left_layout () -> addWidget (smiles);

	_smiles_b = new MyPushButton (_smiles_tc ->right_layout (), "Add smiles");
	connect (_smiles_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_smiles ()));

	_del_box = new MyGroupBox (main_tab ->basic_layout (), "Delete");


// Groups
	ZNAdvancedWidget *groups_tab = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (groups_tab, "Groups");

	_aminoacids_box = new MyGroupBox (groups_tab ->basic_layout (), "Amino acids");

	_aminoacids_gc = new MyGridColumn (_aminoacids_box ->layout (), 2, 10);

    _Ala_b = new MyPushButton ("Ala", 0, 35);
    connect (_Ala_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ala ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Ala_b, 1, 1);

    _Arg_b = new MyPushButton ("Arg", 0, 35);
    connect (_Arg_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Arg ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Arg_b, 1, 2);

    _Asn_b = new MyPushButton ("Asn", 0, 35);
    connect (_Asn_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Asn ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Asn_b, 1, 3);

    _Asp_b = new MyPushButton ("Asp", 0, 35);
    connect (_Asp_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Asp ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Asp_b, 1, 4);

    _Cys_b = new MyPushButton ("Cys", 0, 35);
    connect (_Cys_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cys ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Cys_b, 1, 5);

    _Glu_b = new MyPushButton ("Glu", 0, 35);
    connect (_Glu_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Glu ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Glu_b, 1, 6);

    _Cln_b = new MyPushButton ("Cln", 0, 35);
    connect (_Cln_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cln ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Cln_b, 1, 7);

    _Gly_b = new MyPushButton ("Gly", 0, 35);
    connect (_Gly_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Gly ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Gly_b, 1, 8);

    _His_b = new MyPushButton ("His", 0, 35);
    connect (_His_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_His ()));
    _aminoacids_gc ->gridlayout ->addWidget (_His_b, 1, 9);

    _Ile_b = new MyPushButton ("Ile", 0, 35);
    connect (_Ile_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ile ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Ile_b, 1, 10);

    _Leu_b = new MyPushButton ("Leu", 0, 35);
    connect (_Leu_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Leu ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Leu_b, 2, 1);

    _Lys_b = new MyPushButton ("Lys", 0, 35);
    connect (_Lys_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Lys ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Lys_b, 2, 2);

    _Met_b = new MyPushButton ("Met", 0, 35);
    connect (_Met_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Met ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Met_b, 2, 3);

    _Phe_b = new MyPushButton ("Phe", 0, 35);
    connect (_Phe_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Phe ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Phe_b, 2, 4);

    _Pro_b = new MyPushButton ("Pro", 0, 35);
    connect (_Pro_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Pro ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Pro_b, 2, 5);

    _Ser_b = new MyPushButton ("Ser", 0, 35);
    connect (_Ser_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ser ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Ser_b, 2, 6);

    _Thr_b = new MyPushButton ("Thr", 0, 35);
    connect (_Thr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Thr ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Thr_b, 2, 7);

    _Trp_b = new MyPushButton ("Trp", 0, 35);
    connect (_Trp_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Trp ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Trp_b, 2, 8);

    _Tyr_b = new MyPushButton ("Tyr", 0, 35);
    connect (_Tyr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Tyr ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Tyr_b, 2, 9);

    _Val_b = new MyPushButton ("Val", 0, 35);
    connect (_Val_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Val ()));
    _aminoacids_gc ->gridlayout ->addWidget (_Val_b, 2, 10);


	_heterocycles_box = new MyGroupBox (groups_tab ->basic_layout (), "Heterocycles");

	_nucleotides_box = new MyGroupBox (groups_tab ->basic_layout (), "Nucleotides");


// Periodic table
	ZNAdvancedWidget *periodictable_tab = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (periodictable_tab, "Periodic table");

	_periodictable_gc = new MyGridColumn (periodictable_tab ->basic_layout (), 9, 18);

	QFont font("Helvetica", 8);

	_H_b = new MyPushButton ("H", 18, 24);
	connect (_H_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_H ()));
	_periodictable_gc ->gridlayout ->addWidget (_H_b, 1, 1);

	_He_b = new MyPushButton ("He", 18, 24);
	connect (_He_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_He ()));
	_periodictable_gc ->gridlayout ->addWidget (_He_b, 1, 18);

	_Li_b = new MyPushButton ("Li", 18, 24);
	connect (_Li_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Li ()));
	_periodictable_gc ->gridlayout ->addWidget (_Li_b, 2, 1);

	_Be_b = new MyPushButton ("Be", 18, 24);
	connect (_Be_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Be ()));
	_periodictable_gc ->gridlayout ->addWidget (_Be_b, 2, 2);

	_B_b = new MyPushButton ("B", 18, 24);
	connect (_B_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_B ()));
	_periodictable_gc ->gridlayout ->addWidget (_B_b, 2, 13);

	_C_b = new MyPushButton ("C", 18, 24);
	connect (_C_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_C ()));
	_periodictable_gc ->gridlayout ->addWidget (_C_b, 2, 14);

	_N_b = new MyPushButton ("N", 18, 24);
	connect (_N_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_N ()));
	_periodictable_gc ->gridlayout ->addWidget (_N_b, 2, 15);

	_O_b = new MyPushButton ("O", 18, 24);
	connect (_O_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_O ()));
	_periodictable_gc ->gridlayout ->addWidget (_O_b, 2, 16);

	_F_b = new MyPushButton ("F", 18, 24);
	connect (_F_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_F ()));
	_periodictable_gc ->gridlayout ->addWidget (_F_b, 2, 17);

	_Ne_b = new MyPushButton ("Ne", 18, 24);
	connect (_Ne_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ne ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ne_b, 2, 18);

	_Na_b = new MyPushButton ("Na", 18, 24);
	connect (_Na_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Na ()));
	_periodictable_gc ->gridlayout ->addWidget (_Na_b, 3, 1);

	_Mg_b = new MyPushButton ("Mg", 18, 24);
	connect (_Mg_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Mg ()));
	_periodictable_gc ->gridlayout ->addWidget (_Mg_b, 3, 2);

	_Al_b = new MyPushButton ("Al", 18, 24);
	connect (_Al_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Al ()));
	_periodictable_gc ->gridlayout ->addWidget (_Al_b, 3, 13);

	_Si_b = new MyPushButton ("Si", 18, 24);
	connect (_Si_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Si ()));
	_periodictable_gc ->gridlayout ->addWidget (_Si_b, 3, 14);

	_P_b = new MyPushButton ("P", 18, 24);
	connect (_P_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_P ()));
	_periodictable_gc ->gridlayout ->addWidget (_P_b, 3, 15);

	_S_b = new MyPushButton ("S", 18, 24);
	connect (_S_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_S ()));
	_periodictable_gc ->gridlayout ->addWidget (_S_b, 3, 16);

	_Cl_b = new MyPushButton ("Cl", 18, 24);
	connect (_Cl_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cl ()));
	_periodictable_gc ->gridlayout ->addWidget (_Cl_b, 3, 17);

	_Ar_b = new MyPushButton ("Ar", 18, 24);
	connect (_Ar_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ar ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ar_b, 3, 18);

	_K_b = new MyPushButton ("K", 18, 24);
	connect (_K_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_K ()));
	_periodictable_gc ->gridlayout ->addWidget (_K_b, 4, 1);

	_Ca_b = new MyPushButton ("Ca", 18, 24);
	connect (_Ca_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ca ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ca_b, 4, 2);

	_Sc_b = new MyPushButton ("Sc", 18, 24);
	connect (_Sc_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Sc ()));
	_periodictable_gc ->gridlayout ->addWidget (_Sc_b, 4, 3);

	_Ti_b = new MyPushButton ("Ti", 18, 24);
	connect (_Ti_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ti ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ti_b, 4, 4);

	_V_b = new MyPushButton ("V", 18, 24);
	connect (_V_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_V ()));
	_periodictable_gc ->gridlayout ->addWidget (_V_b, 4, 5);

	_Cr_b = new MyPushButton ("Cr", 18, 24);
	connect (_Cr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cr ()));
	_periodictable_gc ->gridlayout ->addWidget (_Cr_b, 4, 6);

	_Mn_b = new MyPushButton ("Mn", 18, 24);
	connect (_Mn_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Mn ()));
	_periodictable_gc ->gridlayout ->addWidget (_Mn_b, 4, 7);

	_Fe_b = new MyPushButton ("Fe", 18, 24);
	connect (_Fe_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Fe ()));
	_periodictable_gc ->gridlayout ->addWidget (_Fe_b, 4, 8);

	_Co_b = new MyPushButton ("Co", 18, 24);
	connect (_Co_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Co ()));
	_periodictable_gc ->gridlayout ->addWidget (_Co_b, 4, 9);

	_Ni_b = new MyPushButton ("Ni", 18, 24);
	connect (_Ni_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ni ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ni_b, 4, 10);

	_Cu_b = new MyPushButton ("Cu", 18, 24);
	connect (_Cu_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cu ()));
	_periodictable_gc ->gridlayout ->addWidget (_Cu_b, 4, 11);

	_Zn_b = new MyPushButton ("Zn", 18, 24);
	connect (_Zn_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Zn ()));
	_periodictable_gc ->gridlayout ->addWidget (_Zn_b, 4, 12);

	_Ga_b = new MyPushButton ("Ga", 18, 24);
	connect (_Ga_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ga ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ga_b, 4, 13);

	_Ge_b = new MyPushButton ("Ge", 18, 24);
	connect (_Ge_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ge ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ge_b, 4, 14);

	_As_b = new MyPushButton ("As", 18, 24);
	connect (_As_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_As ()));
	_periodictable_gc ->gridlayout ->addWidget (_As_b, 4, 15);

	_Se_b = new MyPushButton ("Se", 18, 24);
	connect (_Se_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Se ()));
	_periodictable_gc ->gridlayout ->addWidget (_Se_b, 4, 16);

	_Br_b = new MyPushButton ("Br", 18, 24);
	connect (_Br_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Br ()));
	_periodictable_gc ->gridlayout ->addWidget (_Br_b, 4, 17);

	_Kr_b = new MyPushButton ("Kr", 18, 24);
	connect (_Kr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Kr ()));
	_periodictable_gc ->gridlayout ->addWidget (_Kr_b, 4, 18);

	_Rb_b = new MyPushButton ("Rb", 18, 24);
	connect (_Rb_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Rb ()));
	_periodictable_gc ->gridlayout ->addWidget (_Rb_b, 5, 1);

	_Sr_b = new MyPushButton ("Sr", 18, 24);
	connect (_Sr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Sr ()));
	_periodictable_gc ->gridlayout ->addWidget (_Sr_b, 5, 2);

	_Y_b = new MyPushButton ("Y", 18, 24);
	connect (_Y_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Y ()));
	_periodictable_gc ->gridlayout ->addWidget (_Y_b, 5, 3);

	_Zr_b = new MyPushButton ("Zr", 18, 24);
	connect (_Zr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Zr ()));
	_periodictable_gc ->gridlayout ->addWidget (_Zr_b, 5, 4);

	_Nb_b = new MyPushButton ("Nb", 18, 24);
	connect (_Nb_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Nb ()));
	_periodictable_gc ->gridlayout ->addWidget (_Nb_b, 5, 5);

	_Mo_b = new MyPushButton ("Mo", 18, 24);
	connect (_Mo_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Mo ()));
	_periodictable_gc ->gridlayout ->addWidget (_Mo_b, 5, 6);

	_Tc_b = new MyPushButton ("Tc", 18, 24);
	connect (_Tc_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Tc ()));
	_periodictable_gc ->gridlayout ->addWidget (_Tc_b, 5, 7);

	_Ru_b = new MyPushButton ("Ru", 18, 24);
	connect (_Ru_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ru ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ru_b, 5, 8);

	_Rh_b = new MyPushButton ("Rh", 18, 24);
	connect (_Rh_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Rh ()));
	_periodictable_gc ->gridlayout ->addWidget (_Rh_b, 5, 9);

	_Pd_b = new MyPushButton ("Pd", 18, 24);
	connect (_Pd_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Pd ()));
	_periodictable_gc ->gridlayout ->addWidget (_Pd_b, 5, 10);

	_Ag_b = new MyPushButton ("Ag", 18, 24);
	connect (_Ag_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ag ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ag_b, 5, 11);

	_Cd_b = new MyPushButton ("Cd", 18, 24);
	connect (_Cd_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cd ()));
	_periodictable_gc ->gridlayout ->addWidget (_Cd_b, 5, 12);

	_In_b = new MyPushButton ("In", 18, 24);
	connect (_In_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_In ()));
	_periodictable_gc ->gridlayout ->addWidget (_In_b, 5, 13);

	_Sn_b = new MyPushButton ("Sn", 18, 24);
	connect (_Sn_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Sn ()));
	_periodictable_gc ->gridlayout ->addWidget (_Sn_b, 5, 14);

	_Sb_b = new MyPushButton ("Sb", 18, 24);
	connect (_Sb_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Sb ()));
	_periodictable_gc ->gridlayout ->addWidget (_Sb_b, 5, 15);

	_Te_b = new MyPushButton ("Te", 18, 24);
	connect (_Te_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Te ()));
	_periodictable_gc ->gridlayout ->addWidget (_Te_b, 5, 16);

	_I_b = new MyPushButton ("I", 18, 24);
	connect (_I_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_I ()));
	_periodictable_gc ->gridlayout ->addWidget (_I_b, 5, 17);

	_Xe_b = new MyPushButton ("Xe", 18, 24);
	connect (_Xe_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Xe ()));
	_periodictable_gc ->gridlayout ->addWidget (_Xe_b, 5, 18);

	_Cs_b = new MyPushButton ("Cs", 18, 24);
	connect (_Cs_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cs ()));
	_periodictable_gc ->gridlayout ->addWidget (_Cs_b, 6, 1);

	_Ba_b = new MyPushButton ("Ba", 18, 24);
	connect (_Ba_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ba ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ba_b, 6, 2);

	_Lu_b = new MyPushButton ("Lu", 18, 24);
	connect (_Lu_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Lu ()));
	_periodictable_gc ->gridlayout ->addWidget (_Lu_b, 8, 17);

	_Hf_b = new MyPushButton ("Hf", 18, 24);
	connect (_Hf_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Hf ()));
	_periodictable_gc ->gridlayout ->addWidget (_Hf_b, 6, 4);

	_Ta_b = new MyPushButton ("Ta", 18, 24);
	connect (_Ta_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ta ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ta_b, 6, 5);

	_W_b = new MyPushButton ("W", 18, 24);
	connect (_W_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_W ()));
	_periodictable_gc ->gridlayout ->addWidget (_W_b, 6, 6);

	_Re_b = new MyPushButton ("Re", 18, 24);
	connect (_Re_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Re ()));
	_periodictable_gc ->gridlayout ->addWidget (_Re_b, 6, 7);

	_Os_b = new MyPushButton ("Os", 18, 24);
	connect (_Os_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Os ()));
	_periodictable_gc ->gridlayout ->addWidget (_Os_b, 6, 8);

	_Ir_b = new MyPushButton ("Ir", 18, 24);
	connect (_Ir_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ir ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ir_b, 6, 9);

	_Pt_b = new MyPushButton ("Pt", 18, 24);
	connect (_Pt_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Pt ()));
	_periodictable_gc ->gridlayout ->addWidget (_Pt_b, 6, 10);

	_Au_b = new MyPushButton ("Au", 18, 24);
	connect (_Au_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Au ()));
	_periodictable_gc ->gridlayout ->addWidget (_Au_b, 6, 11);

	_Hg_b = new MyPushButton ("Hg", 18, 24);
	connect (_Hg_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Hg ()));
	_periodictable_gc ->gridlayout ->addWidget (_Hg_b, 6, 12);

	_Tl_b = new MyPushButton ("Tl", 18, 24);
	connect (_Tl_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Tl ()));
	_periodictable_gc ->gridlayout ->addWidget (_Tl_b, 6, 13);

	_Pb_b = new MyPushButton ("Pb", 18, 24);
	connect (_Pb_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Pb ()));
	_periodictable_gc ->gridlayout ->addWidget (_Pb_b, 6, 14);

	_Bi_b = new MyPushButton ("Bi", 18, 24);
	connect (_Bi_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Bi ()));
	_periodictable_gc ->gridlayout ->addWidget (_Bi_b, 6, 15);

	_Po_b = new MyPushButton ("Po", 18, 24);
	connect (_Po_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Po ()));
	_periodictable_gc ->gridlayout ->addWidget (_Po_b, 6, 16);

	_At_b = new MyPushButton ("At", 18, 24);
	connect (_At_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_At ()));
	_periodictable_gc ->gridlayout ->addWidget (_At_b, 6, 17);

	_Rn_b = new MyPushButton ("Rn", 18, 24);
	connect (_Rn_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Rn ()));
	_periodictable_gc ->gridlayout ->addWidget (_Rn_b, 6, 18);

	_Fr_b = new MyPushButton ("Fr", 18, 24);
	connect (_Fr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Fr ()));
	_periodictable_gc ->gridlayout ->addWidget (_Fr_b, 7, 1);

	_Ra_b = new MyPushButton ("Ra", 18, 24);
	connect (_Ra_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ra ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ra_b, 7, 2);

	_Lr_b = new MyPushButton ("Lr", 18, 24);
	connect (_Lr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Lr ()));
	_periodictable_gc ->gridlayout ->addWidget (_Lr_b, 9, 17);

	_Rf_b = new MyPushButton ("Rf", 18, 24);
	connect (_Rf_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Rf ()));
	_periodictable_gc ->gridlayout ->addWidget (_Rf_b, 7, 4);

	_Db_b = new MyPushButton ("Db", 18, 24);
	connect (_Db_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Db ()));
	_periodictable_gc ->gridlayout ->addWidget (_Db_b, 7, 5);

	_Sg_b = new MyPushButton ("Sg", 18, 24);
	connect (_Sg_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Sg ()));
	_periodictable_gc ->gridlayout ->addWidget (_Sg_b, 7, 6);

	_Bh_b = new MyPushButton ("Bh", 18, 24);
	connect (_Bh_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Bh ()));
	_periodictable_gc ->gridlayout ->addWidget (_Bh_b, 7, 7);

	_Hs_b = new MyPushButton ("Hs", 18, 24);
	connect (_Hs_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Hs ()));
	_periodictable_gc ->gridlayout ->addWidget (_Hs_b, 7, 8);

	_Mt_b = new MyPushButton ("Mt", 18, 24);
	connect (_Mt_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Mt ()));
	_periodictable_gc ->gridlayout ->addWidget (_Mt_b, 7, 9);

	_Ds_b = new MyPushButton ("Ds", 18, 24);
	connect (_Ds_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ds ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ds_b, 7, 10);

	_Rg_b = new MyPushButton ("Rg", 18, 24);
	connect (_Rg_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Rg ()));
	_periodictable_gc ->gridlayout ->addWidget (_Rg_b, 7, 11);

	_Uub_b = new MyPushButton ("Uub", 18, 24);
	connect (_Uub_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Uub ()));
	_periodictable_gc ->gridlayout ->addWidget (_Uub_b, 7, 12);

	_Uut_b = new MyPushButton ("Uut", 18, 24);
	connect (_Uut_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Uut ()));
	_periodictable_gc ->gridlayout ->addWidget (_Uut_b, 7, 13);

	_Uuq_b = new MyPushButton ("Uuq", 18, 24);
	connect (_Uuq_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Uuq ()));
	_periodictable_gc ->gridlayout ->addWidget (_Uuq_b, 7, 14);

	_Uup_b = new MyPushButton ("Uup", 18, 24);
	connect (_Uup_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Uup ()));
	_periodictable_gc ->gridlayout ->addWidget (_Uup_b, 7, 15);

	_Uuh_b = new MyPushButton ("Uuh", 18, 24);
	connect (_Uuh_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Uuh ()));
	_periodictable_gc ->gridlayout ->addWidget (_Uuh_b, 7, 16);

	_Uus_b = new MyPushButton ("Uus", 18, 24);
	connect (_Uus_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Uus ()));
	_periodictable_gc ->gridlayout ->addWidget (_Uus_b, 7, 17);

	_Uuo_b = new MyPushButton ("Uuo", 18, 24);
	connect (_Uuo_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Uuo ()));
	_periodictable_gc ->gridlayout ->addWidget (_Uuo_b, 7, 18);

	_La_b = new MyPushButton ("La", 18, 24);
	connect (_La_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_La ()));
	_periodictable_gc ->gridlayout ->addWidget (_La_b, 8, 3);

	_Ce_b = new MyPushButton ("Ce", 18, 24);
	connect (_Ce_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ce ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ce_b, 8, 4);

	_Pr_b = new MyPushButton ("Pr", 18, 24);
	connect (_Pr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Pr ()));
	_periodictable_gc ->gridlayout ->addWidget (_Pr_b, 8, 5);

	_Nd_b = new MyPushButton ("Nd", 18, 24);
	connect (_Nd_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Nd ()));
	_periodictable_gc ->gridlayout ->addWidget (_Nd_b, 8, 6);

	_Pm_b = new MyPushButton ("Pm", 18, 24);
	connect (_Pm_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Pm ()));
	_periodictable_gc ->gridlayout ->addWidget (_Pm_b, 8, 7);

	_Sm_b = new MyPushButton ("Sm", 18, 24);
	connect (_Sm_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Sm ()));
	_periodictable_gc ->gridlayout ->addWidget (_Sm_b, 8, 8);

	_Eu_b = new MyPushButton ("Eu", 18, 24);
	connect (_Eu_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Eu ()));
	_periodictable_gc ->gridlayout ->addWidget (_Eu_b, 8, 9);

	_Gd_b = new MyPushButton ("Gd", 18, 24);
	connect (_Gd_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Gd ()));
	_periodictable_gc ->gridlayout ->addWidget (_Gd_b, 8, 10);

	_Tb_b = new MyPushButton ("Tb", 18, 24);
	connect (_Tb_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Tb ()));
	_periodictable_gc ->gridlayout ->addWidget (_Tb_b, 8, 11);

	_Dy_b = new MyPushButton ("Dy", 18, 24);
	connect (_Dy_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Dy ()));
	_periodictable_gc ->gridlayout ->addWidget (_Dy_b, 8, 12);

	_Ho_b = new MyPushButton ("Ho", 18, 24);
	connect (_Ho_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ho ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ho_b, 8, 13);

	_Er_b = new MyPushButton ("Er", 18, 24);
	connect (_Er_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Er ()));
	_periodictable_gc ->gridlayout ->addWidget (_Er_b, 8, 14);

	_Tm_b = new MyPushButton ("Tm", 18, 24);
	connect (_Tm_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Tm ()));
	_periodictable_gc ->gridlayout ->addWidget (_Tm_b, 8, 15);

	_Yb_b = new MyPushButton ("Yb", 18, 24);
	connect (_Yb_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Yb ()));
	_periodictable_gc ->gridlayout ->addWidget (_Yb_b, 8, 16);

	_Ac_b = new MyPushButton ("Ac", 18, 24);
	connect (_Ac_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Ac ()));
	_periodictable_gc ->gridlayout ->addWidget (_Ac_b, 9, 3);

	_Th_b = new MyPushButton ("Th", 18, 24);
	connect (_Th_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Th ()));
	_periodictable_gc ->gridlayout ->addWidget (_Th_b, 9, 4);

	_Pa_b = new MyPushButton ("Pa", 18, 24);
	connect (_Pa_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Pa ()));
	_periodictable_gc ->gridlayout ->addWidget (_Pa_b, 9, 5);

	_U_b = new MyPushButton ("U", 18, 24);
	connect (_U_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_U ()));
	_periodictable_gc ->gridlayout ->addWidget (_U_b, 9, 6);

	_Np_b = new MyPushButton ("Np", 18, 24);
	connect (_Np_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Np ()));
	_periodictable_gc ->gridlayout ->addWidget (_Np_b, 9, 7);

	_Pu_b = new MyPushButton ("Pu", 18, 24);
	connect (_Pu_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Pu ()));
	_periodictable_gc ->gridlayout ->addWidget (_Pu_b, 9, 8);

	_Am_b = new MyPushButton ("Am", 18, 24);
	connect (_Am_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Am ()));
	_periodictable_gc ->gridlayout ->addWidget (_Am_b, 9, 9);

	_Cm_b = new MyPushButton ("Cm", 18, 24);
	connect (_Cm_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cm ()));
	_periodictable_gc ->gridlayout ->addWidget (_Cm_b, 9, 10);

	_Bk_b = new MyPushButton ("Bk", 18, 24);
	connect (_Bk_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Bk ()));
	_periodictable_gc ->gridlayout ->addWidget (_Bk_b, 9, 11);

	_Cf_b = new MyPushButton ("Cf", 18, 24);
	connect (_Cf_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Cf ()));
	_periodictable_gc ->gridlayout ->addWidget (_Cf_b, 9, 12);

	_Es_b = new MyPushButton ("Es", 18, 24);
	connect (_Es_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Es ()));
	_periodictable_gc ->gridlayout ->addWidget (_Es_b, 9, 13);

	_Fm_b = new MyPushButton ("Fm", 18, 24);
	connect (_Fm_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Fm ()));
	_periodictable_gc ->gridlayout ->addWidget (_Fm_b, 9, 14);

	_Md_b = new MyPushButton ("Md", 18, 24);
	connect (_Md_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_Md ()));
	_periodictable_gc ->gridlayout ->addWidget (_Md_b, 9, 15);

	_No_b = new MyPushButton ("No", 18, 24);
	connect (_No_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_No ()));
	_periodictable_gc ->gridlayout ->addWidget (_No_b, 9, 16);

	_H_b ->fbutton ->setFont (font);
	_Li_b ->fbutton ->setFont (font);
	_Be_b ->fbutton ->setFont (font);
	_Na_b ->fbutton ->setFont (font);
	_Mg_b ->fbutton ->setFont (font);
	_K_b ->fbutton ->setFont (font);
	_Ca_b ->fbutton ->setFont (font);
	_Rb_b ->fbutton ->setFont (font);
	_Sr_b ->fbutton ->setFont (font);
	_Cs_b ->fbutton ->setFont (font);
	_Ba_b ->fbutton ->setFont (font);
	_Fr_b ->fbutton ->setFont (font);
	_Ra_b ->fbutton ->setFont (font);


	_B_b ->fbutton ->setFont (font);
	_C_b ->fbutton ->setFont (font);
	_N_b ->fbutton ->setFont (font);
	_O_b ->fbutton ->setFont (font);
	_F_b ->fbutton ->setFont (font);
	_Al_b ->fbutton ->setFont (font);
	_Si_b ->fbutton ->setFont (font);
	_P_b ->fbutton ->setFont (font);
	_S_b ->fbutton ->setFont (font);
	_Cl_b ->fbutton ->setFont (font);
	_Ga_b ->fbutton ->setFont (font);
	_Ge_b ->fbutton ->setFont (font);
	_As_b ->fbutton ->setFont (font);
	_Se_b ->fbutton ->setFont (font);
	_Br_b ->fbutton ->setFont (font);
	_In_b ->fbutton ->setFont (font);
	_Sn_b ->fbutton ->setFont (font);
	_Sb_b ->fbutton ->setFont (font);
	_Te_b ->fbutton ->setFont (font);
	_I_b ->fbutton ->setFont (font);
	_Tl_b ->fbutton ->setFont (font);
	_Pb_b ->fbutton ->setFont (font);
	_Bi_b ->fbutton ->setFont (font);
	_Po_b ->fbutton ->setFont (font);
	_At_b ->fbutton ->setFont (font);

	_He_b ->fbutton ->setFont (font);
	_Ne_b ->fbutton ->setFont (font);
	_Ar_b ->fbutton ->setFont (font);
	_Kr_b ->fbutton ->setFont (font);
	_Xe_b ->fbutton ->setFont (font);
	_Rn_b ->fbutton ->setFont (font);

	_Uub_b ->fbutton ->setFont (font);
	_Uut_b ->fbutton ->setFont (font);
	_Uuq_b ->fbutton ->setFont (font);
	_Uup_b ->fbutton ->setFont (font);
	_Uuh_b ->fbutton ->setFont (font);
	_Uus_b ->fbutton ->setFont (font);
	_Uuo_b ->fbutton ->setFont (font);

	_Sc_b ->fbutton ->setFont (font);
	_Ti_b ->fbutton ->setFont (font);
	_V_b ->fbutton ->setFont (font);
	_Cr_b ->fbutton ->setFont (font);
	_Mn_b ->fbutton ->setFont (font);
	_Fe_b ->fbutton ->setFont (font);
	_Co_b ->fbutton ->setFont (font);
	_Ni_b ->fbutton ->setFont (font);
	_Cu_b ->fbutton ->setFont (font);
	_Zn_b ->fbutton ->setFont (font);
	_Y_b ->fbutton ->setFont (font);
	_Zr_b ->fbutton ->setFont (font);
	_Nb_b ->fbutton ->setFont (font);
	_Mo_b ->fbutton ->setFont (font);
	_Tc_b ->fbutton ->setFont (font);
	_Ru_b ->fbutton ->setFont (font);
	_Rh_b ->fbutton ->setFont (font);
	_Pd_b ->fbutton ->setFont (font);
	_Ag_b ->fbutton ->setFont (font);
	_Cd_b ->fbutton ->setFont (font);
	_Lu_b ->fbutton ->setFont (font);
	_Hf_b ->fbutton ->setFont (font);
	_Ta_b ->fbutton ->setFont (font);
	_W_b ->fbutton ->setFont (font);
	_Re_b ->fbutton ->setFont (font);
	_Os_b ->fbutton ->setFont (font);
	_Ir_b ->fbutton ->setFont (font);
	_Pt_b ->fbutton ->setFont (font);
	_Au_b ->fbutton ->setFont (font);
	_Hg_b ->fbutton ->setFont (font);
	_Lr_b ->fbutton ->setFont (font);
	_Rf_b ->fbutton ->setFont (font);
	_Db_b ->fbutton ->setFont (font);
	_Sg_b ->fbutton ->setFont (font);
	_Bh_b ->fbutton ->setFont (font);
	_Hs_b ->fbutton ->setFont (font);
	_Mt_b ->fbutton ->setFont (font);
	_Ds_b ->fbutton ->setFont (font);
	_Rg_b ->fbutton ->setFont (font);

	_Ac_b ->fbutton ->setFont (font);
	_Th_b ->fbutton ->setFont (font);
	_Pa_b ->fbutton ->setFont (font);
	_U_b ->fbutton ->setFont (font);
	_Np_b ->fbutton ->setFont (font);
	_Pu_b ->fbutton ->setFont (font);
	_Am_b ->fbutton ->setFont (font);
	_Cm_b ->fbutton ->setFont (font);
	_Bk_b ->fbutton ->setFont (font);
	_Cf_b ->fbutton ->setFont (font);
	_Es_b ->fbutton ->setFont (font);
	_Fm_b ->fbutton ->setFont (font);
	_Md_b ->fbutton ->setFont (font);
	_No_b ->fbutton ->setFont (font);
	_La_b ->fbutton ->setFont (font);
	_Ce_b ->fbutton ->setFont (font);
	_Pr_b ->fbutton ->setFont (font);
	_Nd_b ->fbutton ->setFont (font);
	_Pm_b ->fbutton ->setFont (font);
	_Sm_b ->fbutton ->setFont (font);
	_Eu_b ->fbutton ->setFont (font);
	_Gd_b ->fbutton ->setFont (font);
	_Tb_b ->fbutton ->setFont (font);
	_Dy_b ->fbutton ->setFont (font);
	_Ho_b ->fbutton ->setFont (font);
	_Er_b ->fbutton ->setFont (font);
	_Tm_b ->fbutton ->setFont (font);
	_Yb_b ->fbutton ->setFont (font);


	_H_b ->fbutton ->setPalette(Qt::yellow);
	_Li_b ->fbutton ->setPalette(Qt::yellow);
	_Be_b ->fbutton ->setPalette(Qt::yellow);
	_Na_b ->fbutton ->setPalette(Qt::yellow);
	_Mg_b ->fbutton ->setPalette(Qt::yellow);
	_K_b ->fbutton ->setPalette(Qt::yellow);
	_Ca_b ->fbutton ->setPalette(Qt::yellow);
	_Rb_b ->fbutton ->setPalette(Qt::yellow);
	_Sr_b ->fbutton ->setPalette(Qt::yellow);
	_Cs_b ->fbutton ->setPalette(Qt::yellow);
	_Ba_b ->fbutton ->setPalette(Qt::yellow);
	_Fr_b ->fbutton ->setPalette(Qt::yellow);
	_Ra_b ->fbutton ->setPalette(Qt::yellow);


	_B_b ->fbutton ->setPalette(Qt::cyan);
	_C_b ->fbutton ->setPalette(Qt::cyan);
	_N_b ->fbutton ->setPalette(Qt::cyan);
	_O_b ->fbutton ->setPalette(Qt::cyan);
	_F_b ->fbutton ->setPalette(Qt::cyan);
	_Al_b ->fbutton ->setPalette(Qt::cyan);
	_Si_b ->fbutton ->setPalette(Qt::cyan);
	_P_b ->fbutton ->setPalette(Qt::cyan);
	_S_b ->fbutton ->setPalette(Qt::cyan);
	_Cl_b ->fbutton ->setPalette(Qt::cyan);
	_Ga_b ->fbutton ->setPalette(Qt::cyan);
	_Ge_b ->fbutton ->setPalette(Qt::cyan);
	_As_b ->fbutton ->setPalette(Qt::cyan);
	_Se_b ->fbutton ->setPalette(Qt::cyan);
	_Br_b ->fbutton ->setPalette(Qt::cyan);
	_In_b ->fbutton ->setPalette(Qt::cyan);
	_Sn_b ->fbutton ->setPalette(Qt::cyan);
	_Sb_b ->fbutton ->setPalette(Qt::cyan);
	_Te_b ->fbutton ->setPalette(Qt::cyan);
	_I_b ->fbutton ->setPalette(Qt::cyan);
	_Tl_b ->fbutton ->setPalette(Qt::cyan);
	_Pb_b ->fbutton ->setPalette(Qt::cyan);
	_Bi_b ->fbutton ->setPalette(Qt::cyan);
	_Po_b ->fbutton ->setPalette(Qt::cyan);
	_At_b ->fbutton ->setPalette(Qt::cyan);

	_He_b ->fbutton ->setPalette(Qt::red);
	_Ne_b ->fbutton ->setPalette(Qt::red);
	_Ar_b ->fbutton ->setPalette(Qt::red);
	_Kr_b ->fbutton ->setPalette(Qt::red);
	_Xe_b ->fbutton ->setPalette(Qt::red);
	_Rn_b ->fbutton ->setPalette(Qt::red);

	_Uub_b ->fbutton ->setPalette(Qt::white);
	_Uut_b ->fbutton ->setPalette(Qt::white);
	_Uuq_b ->fbutton ->setPalette(Qt::white);
	_Uup_b ->fbutton ->setPalette(Qt::white);
	_Uuh_b ->fbutton ->setPalette(Qt::white);
	_Uus_b ->fbutton ->setPalette(Qt::white);
	_Uuo_b ->fbutton ->setPalette(Qt::white);

	_Sc_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Ti_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_V_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Cr_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Mn_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Fe_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Co_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Ni_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Cu_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Zn_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Y_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Zr_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Nb_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Mo_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Tc_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Ru_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Rh_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Pd_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Ag_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Cd_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Hf_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Ta_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_W_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Re_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Os_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Ir_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Pt_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Au_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Hg_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Rf_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Db_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Sg_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Bh_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Hs_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Mt_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Ds_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Rg_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_Ac_b ->fbutton ->setPalette(QColor(250, 153, 0));
	_La_b ->fbutton ->setPalette(QColor(250, 153, 0));

	_Lu_b ->fbutton ->setPalette(Qt::lightGray);
	_Lr_b ->fbutton ->setPalette(Qt::lightGray);
	_Th_b ->fbutton ->setPalette(Qt::lightGray);
	_Pa_b ->fbutton ->setPalette(Qt::lightGray);
	_U_b ->fbutton ->setPalette(Qt::lightGray);
	_Np_b ->fbutton ->setPalette(Qt::lightGray);
	_Pu_b ->fbutton ->setPalette(Qt::lightGray);
	_Am_b ->fbutton ->setPalette(Qt::lightGray);
	_Cm_b ->fbutton ->setPalette(Qt::lightGray);
	_Bk_b ->fbutton ->setPalette(Qt::lightGray);
	_Cf_b ->fbutton ->setPalette(Qt::lightGray);
	_Es_b ->fbutton ->setPalette(Qt::lightGray);
	_Fm_b ->fbutton ->setPalette(Qt::lightGray);
	_Md_b ->fbutton ->setPalette(Qt::lightGray);
	_No_b ->fbutton ->setPalette(Qt::lightGray);
	_Ce_b ->fbutton ->setPalette(Qt::lightGray);
	_Pr_b ->fbutton ->setPalette(Qt::lightGray);
	_Nd_b ->fbutton ->setPalette(Qt::lightGray);
	_Pm_b ->fbutton ->setPalette(Qt::lightGray);
	_Sm_b ->fbutton ->setPalette(Qt::lightGray);
	_Eu_b ->fbutton ->setPalette(Qt::lightGray);
	_Gd_b ->fbutton ->setPalette(Qt::lightGray);
	_Tb_b ->fbutton ->setPalette(Qt::lightGray);
	_Dy_b ->fbutton ->setPalette(Qt::lightGray);
	_Ho_b ->fbutton ->setPalette(Qt::lightGray);
	_Er_b ->fbutton ->setPalette(Qt::lightGray);
	_Tm_b ->fbutton ->setPalette(Qt::lightGray);
	_Yb_b ->fbutton ->setPalette(Qt::lightGray);

}

void BuilderMenu::add_mol (string str) {
	_data ->ddwin ->builder ->add_mol (str);
}

void BuilderMenu::add_fragment (string str) {
	_data ->ddwin ->builder ->add_fragment (str);
	/*
	OBConversion conv;
	ZNMolecule *mol = new ZNMolecule ();
	conv.SetInFormat ("SMI");
	conv.ReadString (mol, str);
	mol -> AddHydrogens ();
	mend_coordinates (mol);
	CreateZNMoleculeCommand *command = new CreateZNMoleculeCommand (mol, builder -> ddwin);
	builder -> ddwin -> execute (command);
	*/
}

void BuilderMenu::add_smiles () {
	string str = smiles ->text ().toStdString ();
	_data ->ddwin ->builder ->add_fragment (str);
}


void BuilderMenu::add_benzene () {
	string str = "c1ccccc1";
	add_fragment (str);
}

void BuilderMenu::add_ring3 () {
	string str = "C1CC1";
	add_fragment (str);
}

void BuilderMenu::add_ring4 () {
	string str = "C1CCC1";
	add_fragment (str);
}

void BuilderMenu::add_ring5 () {
	string str = "C1CCCC1";
	add_fragment (str);
}

void BuilderMenu::add_ring6 () {
	string str = "C1CCCCC1";
	add_fragment (str);
}

void BuilderMenu::add_ring7 () {
	string str = "C1CCCCCC1";
	add_fragment (str);
}

void BuilderMenu::add_ring8 () {
	string str = "C1CCCCCCC1";
	add_fragment (str);
}

void BuilderMenu::add_furan () {
	string str = "c1cccO1";
	add_fragment (str);
}

void BuilderMenu::add_furanO () {
	string str = "c1cccC1";
	add_fragment (str);
}

void BuilderMenu::add_CO () {
	string str = "C(=O)";
	add_fragment (str);
}

void BuilderMenu::add_NCO () {
	string str = "NC(=O)";
	add_fragment (str);
}

void BuilderMenu::add_COOH () {
	string str = "C(=O)O";
	add_fragment (str);
}

void BuilderMenu::add_CCd () {
	string str = "C=C";
	add_fragment (str);
}

void BuilderMenu::add_CCt () {
	string str = "C#C";
	add_fragment (str);
}

void BuilderMenu::add_NO2 () {
	string str = "[H]N(=O)(=O)";
	add_fragment (str);
}

void BuilderMenu::add_PO3 () {
	string str = "[H]P(=O)(=O)O";
	add_fragment (str);
}

void BuilderMenu::add_SO2 () {
	string str = "[H]S(=O)([H])=O";
	add_fragment (str);
}

void BuilderMenu::add_CN () {
	string str = "C#N";
	add_fragment (str);
}

void BuilderMenu::single_bond () {
    builder->add_bond (1);
}

void BuilderMenu::double_bond () {
    builder->add_bond (2);
}

void BuilderMenu::triple_bond () {
    builder->add_bond (3);
}

void BuilderMenu::no_bond () {
  //  builder->remove_bond ();
}


void BuilderMenu::add_H () {
    builder->add_atom (1);
}
void BuilderMenu::add_He () {
    builder->add_atom (2);
}
void BuilderMenu::add_Li () {
    builder->add_atom (3);
}
void BuilderMenu::add_Be () {
    builder->add_atom (4);
}
void BuilderMenu::add_B () {
    builder->add_atom (5);
}
void BuilderMenu::add_C () {
    builder->add_atom (6);
}
void BuilderMenu::add_N () {
    builder->add_atom (7);
}
void BuilderMenu::add_O () {
    builder->add_atom (8);
}
void BuilderMenu::add_F () {
    builder->add_atom (9);
}
void BuilderMenu::add_Ne () {
    builder->add_atom (10);
}
void BuilderMenu::add_Na () {
    builder->add_atom (11);
}
void BuilderMenu::add_Mg () {
    builder->add_atom (12);
}
void BuilderMenu::add_Al () {
    builder->add_atom (13);
}
void BuilderMenu::add_Si () {
    builder->add_atom (14);
}
void BuilderMenu::add_P () {
    builder->add_atom (15);
}
void BuilderMenu::add_S () {
    builder->add_atom (16);
}
void BuilderMenu::add_Cl () {
    builder->add_atom (17);
}
void BuilderMenu::add_Ar () {
    builder->add_atom (18);
}
void BuilderMenu::add_K () {
    builder->add_atom (19);
}
void BuilderMenu::add_Ca () {
    builder->add_atom (20);
}
void BuilderMenu::add_Sc () {
    builder->add_atom (21);
}
void BuilderMenu::add_Ti () {
    builder->add_atom (22);
}
void BuilderMenu::add_V () {
    builder->add_atom (23);
}
void BuilderMenu::add_Cr () {
    builder->add_atom (24);
}
void BuilderMenu::add_Mn () {
    builder->add_atom (25);
}
void BuilderMenu::add_Fe () {
    builder->add_atom (26);
}
void BuilderMenu::add_Co () {
    builder->add_atom (27);
}
void BuilderMenu::add_Ni () {
    builder->add_atom (28);
}
void BuilderMenu::add_Cu () {
    builder->add_atom (29);
}
void BuilderMenu::add_Zn () {
    builder->add_atom (30);
}
void BuilderMenu::add_Ga () {
    builder->add_atom (31);
}
void BuilderMenu::add_Ge () {
    builder->add_atom (32);
}
void BuilderMenu::add_As () {
    builder->add_atom (33);
}
void BuilderMenu::add_Se () {
    builder->add_atom (34);
}
void BuilderMenu::add_Br () {
    builder->add_atom (35);
}
void BuilderMenu::add_Kr () {
    builder->add_atom (36);
}
void BuilderMenu::add_Rb () {
    builder->add_atom (37);
}
void BuilderMenu::add_Sr () {
    builder->add_atom (38);
}
void BuilderMenu::add_Y () {
    builder->add_atom (39);
}
void BuilderMenu::add_Zr () {
    builder->add_atom (40);
}
void BuilderMenu::add_Nb () {
    builder->add_atom (41);
}
void BuilderMenu::add_Mo () {
    builder->add_atom (42);
}
void BuilderMenu::add_Tc () {
    builder->add_atom (43);
}
void BuilderMenu::add_Ru () {
    builder->add_atom (44);
}
void BuilderMenu::add_Rh () {
    builder->add_atom (45);
}
void BuilderMenu::add_Pd () {
    builder->add_atom (46);
}
void BuilderMenu::add_Ag () {
    builder->add_atom (47);
}
void BuilderMenu::add_Cd () {
    builder->add_atom (48);
}
void BuilderMenu::add_In () {
    builder->add_atom (49);
}
void BuilderMenu::add_Sn () {
    builder->add_atom (50);
}
void BuilderMenu::add_Sb () {
    builder->add_atom (51);
}
void BuilderMenu::add_Te () {
    builder->add_atom (52);
}
void BuilderMenu::add_I () {
    builder->add_atom (53);
}
void BuilderMenu::add_Xe () {
    builder->add_atom (54);
}
void BuilderMenu::add_Cs () {
    builder->add_atom (55);
}
void BuilderMenu::add_Ba () {
    builder->add_atom (56);
}
void BuilderMenu::add_La () {
    builder->add_atom (57);
}
void BuilderMenu::add_Ce () {
    builder->add_atom (58);
}
void BuilderMenu::add_Pr () {
    builder->add_atom (59);
}
void BuilderMenu::add_Nd () {
    builder->add_atom (60);
}
void BuilderMenu::add_Pm () {
    builder->add_atom (61);
}
void BuilderMenu::add_Sm () {
    builder->add_atom (62);
}
void BuilderMenu::add_Eu () {
    builder->add_atom (63);
}
void BuilderMenu::add_Gd () {
    builder->add_atom (64);
}
void BuilderMenu::add_Tb () {
    builder->add_atom (65);
}
void BuilderMenu::add_Dy () {
    builder->add_atom (66);
}
void BuilderMenu::add_Ho () {
    builder->add_atom (67);
}
void BuilderMenu::add_Er () {
    builder->add_atom (68);
}
void BuilderMenu::add_Tm () {
    builder->add_atom (69);
}
void BuilderMenu::add_Yb () {
    builder->add_atom (70);
}
void BuilderMenu::add_Lu () {
    builder->add_atom (71);
}
void BuilderMenu::add_Hf () {
    builder->add_atom (72);
}
void BuilderMenu::add_Ta () {
    builder->add_atom (73);
}
void BuilderMenu::add_W () {
    builder->add_atom (74);
}
void BuilderMenu::add_Re () {
    builder->add_atom (75);
}
void BuilderMenu::add_Os () {
    builder->add_atom (76);
}
void BuilderMenu::add_Ir () {
    builder->add_atom (77);
}
void BuilderMenu::add_Pt () {
    builder->add_atom (78);
}
void BuilderMenu::add_Au () {
    builder->add_atom (79);
}
void BuilderMenu::add_Hg () {
    builder->add_atom (80);
}
void BuilderMenu::add_Tl () {
    builder->add_atom (81);
}
void BuilderMenu::add_Pb () {
    builder->add_atom (82);
}
void BuilderMenu::add_Bi () {
    builder->add_atom (83);
}
void BuilderMenu::add_Po () {
    builder->add_atom (84);
}
void BuilderMenu::add_At () {
    builder->add_atom (85);
}
void BuilderMenu::add_Rn () {
    builder->add_atom (86);
}
void BuilderMenu::add_Fr () {
    builder->add_atom (87);
}
void BuilderMenu::add_Ra () {
    builder->add_atom (88);
}
void BuilderMenu::add_Ac () {
    builder->add_atom (89);
}
void BuilderMenu::add_Th () {
    builder->add_atom (90);
}
void BuilderMenu::add_Pa () {
    builder->add_atom (91);
}
void BuilderMenu::add_U () {
    builder->add_atom (92);
}
void BuilderMenu::add_Np () {
    builder->add_atom (93);
}
void BuilderMenu::add_Pu () {
    builder->add_atom (94);
}
void BuilderMenu::add_Am () {
    builder->add_atom (95);
}
void BuilderMenu::add_Cm () {
    builder->add_atom (96);
}
void BuilderMenu::add_Bk () {
    builder->add_atom (97);
}
void BuilderMenu::add_Cf () {
    builder->add_atom (98);
}
void BuilderMenu::add_Es () {
    builder->add_atom (99);
}
void BuilderMenu::add_Fm () {
    builder->add_atom (100);
}
void BuilderMenu::add_Md () {
    builder->add_atom (101);
}
void BuilderMenu::add_No () {
    builder->add_atom (102);
}
void BuilderMenu::add_Lr () {
    builder->add_atom (103);
}
void BuilderMenu::add_Rf () {
    builder->add_atom (104);
}
void BuilderMenu::add_Db () {
    builder->add_atom (105);
}
void BuilderMenu::add_Sg () {
    builder->add_atom (106);
}
void BuilderMenu::add_Bh () {
    builder->add_atom (107);
}
void BuilderMenu::add_Hs () {
    builder->add_atom (108);
}
void BuilderMenu::add_Mt () {
    builder->add_atom (109);
}
void BuilderMenu::add_Ds () {
    builder->add_atom (110);
}
void BuilderMenu::add_Rg () {
    builder->add_atom (111);
}
void BuilderMenu::add_Uub () {
    builder->add_atom (112);
}
void BuilderMenu::add_Uut () {
    builder->add_atom (113);
}
void BuilderMenu::add_Uuq () {
    builder->add_atom (114);
}
void BuilderMenu::add_Uup () {
    builder->add_atom (115);
}
void BuilderMenu::add_Uuh () {
    builder->add_atom (116);
}
void BuilderMenu::add_Uus () {
    builder->add_atom (117);
}
void BuilderMenu::add_Uuo () {
    builder->add_atom (118);
}

void BuilderMenu::add_Ala () {
    string str = "C[C@H](N)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Arg () {
    string str = "N[C@@H](CCCNC(N)=N)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Asn () {
    string str = "N[C@@H](CC(N)=O)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Asp () {
    string str = "N[C@@H](CC(O)=O)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Cys () {
    string str = "C([C@@H](C(=O)O)N)S";
	add_fragment (str);
}
void BuilderMenu::add_Glu () {
    string str = "N[C@@H](CCC(O)=O)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Cln () {
    string str = "N[C@@H](CCC(N)=O)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Gly () {
    string str = "NCC(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_His () {
    string str = "N[C@@H](Cc1[nH]cnc1)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Ile () {
    string str = "CC[C@H](C)[C@H](N)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Leu () {
    string str = "CC(C)C[C@H](N)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Lys () {
    string str = "C(CCN)CC(C(=O)O)N";
	add_fragment (str);
}
void BuilderMenu::add_Met () {
    string str = "CSCC[C@H](N)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Phe () {
    string str = "C1=CC=C(C=C1)CC(C(=O)O)N";
	add_fragment (str);
}
void BuilderMenu::add_Pro () {
    string str = "OC(=O)[C@@H]1CCCN1";
	add_fragment (str);
}
void BuilderMenu::add_Ser () {
    string str = "OCC(N)C(=O)O";
	add_fragment (str);
}
void BuilderMenu::add_Thr () {
    string str = "C[C@@H](O)[C@H](N)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Trp () {
    string str = "N[C@@H](Cc1c2ccccc2n([H])c1)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Tyr () {
    string str = "N[C@@H](Cc1ccc(O)cc1)C(O)=O";
	add_fragment (str);
}
void BuilderMenu::add_Val () {
    string str = "CC(C)C(N)C(=O)O";
	add_fragment (str);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


HapticMenu::HapticMenu (QWidget *parent, Data *dat )
   :    ZNMenu (parent, dat, true, true)
{
	haptic_thread = NULL;
	mult = 1.;
	restrain_lock = new QReadWriteLock;
	last_atom = false;
	last_k = 1.;
	last_dist = 1.5;
	label = new QLabel;
	label ->setText ("");
	
	ZNAdvancedWidget *settings_widget = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (settings_widget, "Settings");
	
	ZNAdvancedWidget *restrains_widget = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (restrains_widget, "Restrains");
	
		ZNAdvancedWidget *energy_widget = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (energy_widget, "Energy");
	
	restrains_widget ->layout() ->addWidget(label);

    minimize = _data ->minimize;
    this->setWindowTitle("Haptic");

   // QWidget *settings_widget = new QWidget (this);
    QLayout *setmainlayout = settings_widget ->layout ();


    dofmode = new QComboBox;
	dofmode -> insertItem(0,  "6" );
    dofmode -> insertItem(1,  "3 N" );
	dofmode -> setCurrentIndex(1);
    setmainlayout -> addWidget (dofmode);

    interff = new QComboBox();
    interff->insertItem(1,  "MMFF" );
//    interff->insertItem(2,  "PLP" );
    interff->insertItem(2,  "Chemscore" );
	interff->insertItem(3,  "MMFF+Chemscore");
    setmainlayout -> addWidget (interff);
	automove_b = true;
	color_by_score_b = false;
	cluster_RMSD = 0.5;
	E_tolerance = 2.f;
	saving = true;
    MyCheckBox *automove = new MyCheckBox (setmainlayout,automove_b, "automove");
    MyCheckBox *color_by_scor = new MyCheckBox (setmainlayout, color_by_score_b, "color by score");
	MyCheckBox *saving_cb = new MyCheckBox (setmainlayout, saving, "auto save results");
	MyFloatEditLine *rmsd_line = new MyFloatEditLine (setmainlayout, "cluster RMSD", cluster_RMSD, 0, 5);
	MyFloatEditLine *Egap_line = new MyFloatEditLine (setmainlayout, "Energy tolerance", E_tolerance, 0, 100);
	MyFloatEditLine *mult_line = new MyFloatEditLine (setmainlayout, "Haptic feedback scaling", mult, 0, 100);
	//MyIntegerEditLine *number_of_threads_edit_line = new MyIntegerEditLine (setmainlayout, "number of threads", minimize ->haptic_number_of_threads);

	MyPushButton *add_restr_b = new MyPushButton (restrains_widget ->layout (), "Add restrain");

	    connect (add_restr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (add_restrain ()));

	
	restrain_list = new MyListView (restrains_widget ->layout (), 0, 0, false);
	
	restrain_dist = new MyFloatEditLine (restrains_widget ->layout(), "rest value",last_dist, 0, 1000);
	restrain_k = new MyFloatEditLine (restrains_widget ->layout(), "Strength",last_k, -1000, 1000); 
	connect (restrain_k ->spinbox, SIGNAL (valueChanged (double )) ,SLOT (update_k (double)));
	connect (restrain_dist ->spinbox, SIGNAL (valueChanged (double)) ,SLOT (update_dist (double)));
	
	MyPushButton *clear_restr_b = new MyPushButton (restrains_widget ->layout (), "Clear restrains");
	connect (restrain_list, SIGNAL (deleting (int)), SLOT (delete_restrain (int)));
	connect (restrain_list ->_lw, SIGNAL (currentRowChanged (int)), SLOT (update_spinboxes (int)));
	connect (clear_restr_b ->fbutton, SIGNAL (clicked ()) ,SLOT (clear_restrains ()));
	
    QHBoxLayout *button_lay = new QHBoxLayout;
	QWidget *butt_widg = new QWidget;
	butt_widg ->setLayout (button_lay);
    setmainlayout -> addWidget (butt_widg);
    QPushButton *Ok_b = new QPushButton ("Start");
    connect (Ok_b, SIGNAL (clicked ()) ,SLOT (Ok ()));
    button_lay -> addWidget (Ok_b);

    QPushButton *end_b = new QPushButton ("Stop");
    connect (end_b, SIGNAL (clicked ()) ,SLOT (end ()));


    button_lay -> addWidget (end_b);

  //  addTab( settings_widget, "Settings" );
    


  //  total_E = new MyLabelf(energy_widget, "Total Energy", &minimize->total_E);
    interaction_E = new MyLabelf(energy_widget ->layout (), "Interaction Energy", _data->total_energy_haptic);
  //  interaction_E = new MyLabelf(energy_widget, "Interaction Energy", &minimize->total_interaction_E);
    QPushButton *savepose_b = new QPushButton ("Save Current Pose");
	energy_widget ->layout() ->addWidget(savepose_b);
	connect (savepose_b, SIGNAL (clicked ()) ,SLOT (user_save_current_pose ()));
}
void HapticMenu::add_restrain_atom (Atom *a) { 
	if (!last_atom) {
		last_atom = a;
		label ->setText ("Choose second Atom");
	}
	else if (last_atom != a) {
		label ->setText ("");
		restrain_lock ->lockForWrite ();
		ElasticRestrain *restrain = new ElasticRestrain;
		restrain ->at2 = a;
		restrain ->at1 = last_atom;
		restrain ->k = last_k;
		restrain ->dist0 = dist (get_coordinates (restrain ->at1),get_coordinates ( restrain->at2));
		restrains.push_back (restrain);
		stringstream ss;
		ss << " "<<restrain ->at1 ->GetType ()<<"\t"<<restrain ->at2 ->GetType ()<<"\t" <<restrain ->dist0<<"\t"<<restrain->k;
		restrain_list ->_lw ->addItem(tr( ss.str ().c_str () ));
		last_atom = 0;
		_data ->ddwin ->adding_restrains = false;
		restrain_lock ->unlock ();

	}
}

void HapticMenu::update_k (double) {
				restrain_lock ->lockForWrite ();

	if (restrain_list ->_lw ->currentRow()>-1) {
		int n = restrain_list ->_lw ->currentRow();
		ElasticRestrain *rest = ((ElasticRestrain *)restrains[n]);
		rest->k = last_k;
		stringstream ss;
		ss << " "<<rest ->at1 ->GetType ()<<"\t"<<rest ->at2 ->GetType ()<<"\t" <<rest ->dist0<<"\t"<<rest->k;
		restrain_list ->_lw ->item (n) -> setText(ss.str().c_str());
		
	}
				restrain_lock ->unlock ();
}

void HapticMenu::update_spinboxes (int i) {
	if (i > -1) {
		ElasticRestrain *rest = ((ElasticRestrain *)restrains[i]);
		restrain_dist ->spinbox ->setValue (rest ->dist0);
		restrain_k ->spinbox ->setValue (rest ->k);
	}
}

void HapticMenu::update_dist (double) {
	restrain_lock ->lockForWrite ();
	
	if (restrain_list ->_lw ->currentRow()>-1) {
		int n = restrain_list ->_lw ->currentRow();
		ElasticRestrain *rest = ((ElasticRestrain *)restrains[n]);
		rest->dist0 = last_dist;
		stringstream ss;
		ss << " "<<rest ->at1 ->GetType ()<<"\t"<<rest ->at2 ->GetType ()<<"\t" <<rest ->dist0<<"\t"<<rest->k;
		restrain_list ->_lw ->item (n) -> setText(ss.str().c_str());
		
	}
	restrain_lock ->unlock ();
}

void HapticMenu::clear_restrains () {
			restrain_lock ->lockForWrite ();
	restrain_list ->_lw ->clear ();
	for (unsigned int i = 0; i < restrains.size (); i ++) {
		delete restrains [i];
	}
	restrains.clear (); 
			restrain_lock ->unlock ();
	
}

void HapticMenu::delete_restrain (int i) {
	restrain_lock ->lockForWrite ();
	delete restrains [i];
	restrains.erase (restrains.begin () + i);
			restrain_lock ->unlock ();
}			 
			 
void HapticMenu::add_restrain () {
	_data ->ddwin ->adding_restrains = true;
			label ->setText ("Choose first Atom");
}

void HapticMenu::user_save_current_pose () {
	haptic_thread ->user_save_result = true;
}

void HapticMenu::Ok () {
	if (_data ->ddwin ->current_target) {
	if (!haptic_thread) {
		haptic_thread = new HapticThread (0, _data ->ddwin);
		haptic_thread ->automove = automove_b;
		haptic_thread ->color_by_score = color_by_score_b;
		for (unsigned int i =0; i< minimize ->interaction_ffs.size (); i++) {
			delete minimize ->interaction_ffs [i];
		}
		 minimize->interaction_ffs.clear ();
		string iff = interff->currentText ().toStdString ();
		if (iff == "MMFF") minimize->interaction_ffs.push_back (new MMFF ());
		else if (iff == "Chemscore") minimize->interaction_ffs.push_back ( new Chemscore ());
		else if (iff == "MMFF+Chemscore") {
			minimize->interaction_ffs.push_back ( new MMFF ());
			minimize->interaction_ffs.push_back ( new Chemscore ());
		}
		else minimize->interaction_ffs.push_back (new MMFF ());

		string dofm =dofmode->currentText ().toStdString ();
		if (dofm == "6") minimize -> haptic_dof_mode = 0;
		else if (dofm == "3 N") minimize -> haptic_dof_mode = 1;
		else minimize -> haptic_dof_mode = 1;
		minimize -> start_haptic_mode ();
	}
	}
}



void HapticMenu::end () {
	if (haptic_thread) {
		ZNMolecule *mol = haptic_thread ->molecule;
		Database *dat = haptic_thread ->results;
		haptic_thread -> stop ();
		haptic_thread -> wait ();
		haptic_thread = 0;
		MoveAtomsCommand *command = new MoveAtomsCommand (_data -> ddwin -> gl, 0);
		FOR_ATOMS_OF_MOL(a, mol) {
			command -> add (&*a, get_coordinates (&*a));
		}
		//  command -> name ("haptic simulation");
		_data -> ddwin -> execute (command);
		_data -> ddwin -> add_database(dat);
		_data -> undo_stack -> endMacro ();
		_data->ddwin->haptic_mode = false;
		//haptic_molecule = NULL;
		_data->ddwin ->data ->current_force_x=0.;
		_data->ddwin ->data ->current_force_y=0.;
		_data->ddwin ->data ->current_force_z=0.;
		_data -> ddwin -> unlock_editing ();
	}
}



void HapticMenu::update_energy () {
  //  total_E->update ();
  //  internal_E->update ();
    interaction_E->update ();
}

void HapticMenu::maybe_save_result () {
	if (haptic_thread) {
		float E = _data ->total_energy_haptic;
		if (E < haptic_thread ->save_E_threshold + E_tolerance) {
			if (E > haptic_thread ->save_E_threshold) 
			{
				haptic_thread ->is_worse = true;
			}
			haptic_thread ->save_E_threshold = E;
			haptic_thread ->save_result = true;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Clicked_atomMenu::Clicked_atomMenu (QWidget *parent, Data* dat) 
   :    ZNMenu (parent, dat, true, true)
    {
	setWindowTitle ("Clicked Atom");
	ZNAdvancedWidget *atom_tab = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (atom_tab, "Atom");
	
		ZNAdvancedWidget *residue_tab = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (residue_tab, "Residue");
	
	
		ZNAdvancedWidget *molecule_tab = new ZNAdvancedWidget ("", false);
	_tabs ->addTab (molecule_tab, "ZNMolecule");
	
	
    MyPushButton *setascenterb = new MyPushButton (atom_tab ->layout (), "Set As Center of view");
    connect (setascenterb ->fbutton, SIGNAL (clicked ()),SLOT (set_clicked_atom_as_center_of_view ()));
    MyPushButton *setasrotcenterb = new MyPushButton (atom_tab ->layout (), "Set As Center of Rotation");
    connect (setasrotcenterb ->fbutton, SIGNAL (clicked ()),SLOT (set_clicked_atom_as_center_of_rotation ()));	
	
	QWidget *coors_widget = new QWidget ();
	coors_widget ->setLayout (new QHBoxLayout ());
	atom_tab ->layout ()->addWidget (coors_widget);
		idl = new QLabel ("");
   //  coors_widget->setMaximumHeight (25);
    aplx = new QLineEdit("");
    aply = new QLineEdit("");
    aplz = new QLineEdit("");
	coors_widget ->layout () -> addWidget (new QLabel("coordinates"));
	coors_widget ->layout ()-> addWidget (aplx);
	coors_widget ->layout ()-> addWidget (aply);
	coors_widget ->layout ()-> addWidget (aplz);


	QWidget *charge_widget = new QWidget ();
	charge_widget ->setLayout (new QHBoxLayout ());
	atom_tab ->layout ()->addWidget (charge_widget);
	atom_tab ->layout () ->addWidget (idl);

	formal_charge_le = new MyIntegerEditLine (charge_widget ->layout (), "Formal Charge", formal_charge, -8, 8);
	
/*   
    Q3Grid *atompropts = new Q3Grid (5,atomselpopupv);
    atompropts->setMaximumSize (350, 60);
    (void) new QLabel("Atom ID", atompropts );
    aplid = new QLabel("", atompropts );
    (void) new QLabel("      ", atompropts );
    (void) new QLabel("Partial charge", atompropts );
    aplq = new QLabel("", atompropts );
    (void) new QLabel("Element", atompropts );
    aplat = new QLabel("", atompropts );
    (void) new QLabel("      ", atompropts );
    (void) new QLabel("Formal charge", atompropts );
    aplfc = new QLineEdit("", atompropts );


    Q3HBox *coors = new Q3HBox (atomselpopupv);
     coors->setMaximumHeight (25);
    (void) new QLabel("coordinates", coors );
    aplx = new QLineEdit("", coors );
    aply = new QLineEdit("", coors );
    aplz = new QLineEdit("", coors );
    aplid->setFrameStyle( Q3Frame::Panel | Q3Frame::Sunken );
    aplat->setFrameStyle( Q3Frame::Panel | Q3Frame::Sunken );
    aplq->setFrameStyle( Q3Frame::Panel | Q3Frame::Sunken );





    Q3VBox *resv = new Q3VBox (this);
    Q3Grid *resgrd = new Q3Grid (5,resv);
    (void) new QLabel("Residue Type", resgrd );
    resna = new QLabel("", resgrd );
    (void) new QLabel("      ", resgrd );
    (void) new QLabel("Residue Number", resgrd );
    resnu = new QLabel("", resgrd );
	(void) new QLabel("Atom role in residue", resgrd );
	aptype = new QLabel("", resgrd );
    resna->setFrameStyle( Q3Frame::Panel | Q3Frame::Sunken );
    resnu->setFrameStyle( Q3Frame::Panel | Q3Frame::Sunken );
    addTab( resv, "Residue" );


    Q3VBox *molv = new Q3VBox (this);
    QPushButton *addh = new QPushButton ("Add hydrogens", molv);
    connect (addh, SIGNAL (clicked ()),SLOT (add_Hs ()));
    addTab( molv, "ZNMolecule" );
	
	*/
}



void Clicked_atomMenu::update (){

	
    Atom *at = clicked_atom;
	Resid *res= at ->GetResidue ();
	
	
	
	formal_charge = at->GetFormalCharge ();
	formal_charge_le ->update ();

//    set_value (aplid, at -> GetIdx ());
//    set_value (aplat, string (etab.GetSymbol (at -> GetAtomicNum ())));
//    set_value (aplq, at -> GetPartialCharge ());
//    set_value (aplfc, at-> GetFormalCharge ());
    set_value (aplx, get_coordinates (at).x());
    set_value (aply, get_coordinates (at).y());
    set_value (aplz, get_coordinates (at).z());
	QString atomID = QString(res->GetAtomID(at).c_str());
	set_value (idl, atomID.trimmed().toStdString());
//	set_value (resnu, at->GetResidue ()->GetNum ());
//	int prop = 0;
//	bool a = at -> GetResidue ()-> GetAtomProperty (at, prop);
//	cerr << a << endl;
//	set_value (aptype, prop);
//	int n = at->GetResidue ()->GetResKey ();
//	set_value (resna, Residue[n]);
//	set_value (resna, at->GetResidue ()->GetName ());
	
}



void Clicked_atomMenu::set (Atom *at){
    clicked_atom = at;
    update ();
}



void Clicked_atomMenu::set_value (QLabel *lab, float val) {
    stringstream ss;
    ss << val;
    lab->setText (QString(ss.str().c_str()));
}



void Clicked_atomMenu::set_value (QLabel *lab, string val) {
    lab->setText (QString(val.c_str()));
}



void Clicked_atomMenu::set_value (QLineEdit *lab, float val) {
    stringstream ss;
    ss << val;
    lab->setText (QString(ss.str().c_str()));
}



void Clicked_atomMenu::set_clicked_atom_as_center_of_view (){

    _data ->ddwin->gl->set_center_of_view (get_coordinates (clicked_atom));
}



void Clicked_atomMenu::set_clicked_atom_as_center_of_rotation () {

    _data ->ddwin->gl->set_center_of_rotation (get_coordinates (clicked_atom));

}



void Clicked_atomMenu::add_Hs () {
    ZNMolecule *mol = (ZNMolecule *) clicked_atom -> GetParent ();
    mol -> DeleteHydrogens ();
    mol -> AddHydrogens ();
    finalise_molecule (mol);
    _data ->ddwin -> gl -> draw_molecule (mol);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

/*
BrowserMenu::BrowserMenu (QWidget *parent, DDWin *ddw)
   :    QWidget (parent)
{
    ddwin = ddw;
	target = NULL;
	QHBoxLayout *hbox = new QHBoxLayout;
	setLayout (hbox);
  //  this->setMinimumSize (500, 250);   
  //  this->setMaximumSize (500, 250);
  //  this->setMinimumSize (500, 250);   
  //  this->setMaximumSize (500, 250);
 //   this->setWindowTitle("Database Browser");

    QPushButton *first = new QPushButton ("<<");
    connect (first, SIGNAL (clicked ()) ,SLOT (first_slot ()));
    QPushButton *prev = new QPushButton ("<");
    connect (prev, SIGNAL (clicked ()) ,SLOT (prev_slot ()));

    QPushButton *next = new QPushButton (">");
    connect (next, SIGNAL (clicked ()) ,SLOT (next_slot ()));

    QPushButton *last = new QPushButton (">>");
    connect (last, SIGNAL (clicked ()) ,SLOT (last_slot ()));

	hbox -> addWidget (first);
	hbox -> addWidget (prev);
	hbox -> addWidget (next);
	hbox -> addWidget (last);
}



void BrowserMenu::first_slot () {
    current_number = 0;
    set_mol ();
}



void BrowserMenu::prev_slot () {
    if (current_number > 0) {
        current_number--;
        set_mol ();
    }
}



void BrowserMenu::set_target (Database *db) {
	target = db;
}



void BrowserMenu::next_slot () {
    if (current_number+1 < target->count_entries ()) {
    current_number++;
    set_mol ();
    }
}



void BrowserMenu::last_slot () {
    current_number = target->count_entries ()-1;
    set_mol ();
}



void BrowserMenu::set_mol () {
    for (unsigned int i=0; i<ddwin->molecules.size (); i++) {
        if (ddwin->molecules[i]->multi) {
            Database_molecule *dm;
            dm = (Database_molecule *) ddwin->molecules[i];
            if (dm->database == target) {
                ddwin->molecules[i] = target->get_molecule (current_number);
                ddwin->set_current_target (i);
                ddwin->gl->draw_molecule (ddwin->target_molecule);
            }
        }
    }
}
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////


DatabaseGrid::DatabaseGrid (QWidget *parent, Database *db, Data *dat) 
:ZNMenu (parent, dat, false, false) {
	hide_cb = new MyCheckBox (_main_widget ->layout (), db -> get_is_hidden (), "Hide");
	connect (hide_cb, SIGNAL (toggled (bool)), SLOT (manage_hide (bool)));
	setup_actions ();
	ext_bool = false;
	tab = new QTableWidget();
    connect (tab, SIGNAL (cellChanged ( int , int )) ,SLOT (update_cell (int, int)));
	_main_widget ->layout () ->addWidget (tab);
	tab -> setSortingEnabled (false);
	set_database(db);
	MyCheckBox *extend = new MyCheckBox (_main_widget ->layout (), database -> get_extend_enabled (), "Extend actions to database");

	_main_widget ->layout ()-> addWidget (tab);
	
	current_number = 0;
	QHBoxLayout *hbox = new QHBoxLayout;
	///addLayout (hbox);
	QPushButton *first = new QPushButton ("<<");
    connect (first, SIGNAL (clicked ()) ,SLOT (first_slot ()));
    QPushButton *prev = new QPushButton ("<");
    connect (prev, SIGNAL (clicked ()) ,SLOT (prev_slot ()));
	
	le = new QLineEdit ();
	connect (le, SIGNAL (textEdited (const QString)), SLOT (manage_number_changed (const QString)));
	
    QPushButton *next = new QPushButton (">");
    connect (next, SIGNAL (clicked ()) ,SLOT (next_slot ()));
	
    QPushButton *last = new QPushButton (">>");
    connect (last, SIGNAL (clicked ()) ,SLOT (last_slot ()));
	
	hbox -> addWidget (first);
	hbox -> addWidget (prev);
	hbox -> addWidget (le);
	hbox -> addWidget (next);
	hbox -> addWidget (last);
	
	QWidget *wid = new QWidget;
	wid ->setLayout (hbox);
	_main_widget ->layout ()-> addWidget (wid);
	
	data = dat;
	connect (tab, SIGNAL (cellDoubleClicked (int, int)), this, SLOT (manage_double_click (int, int)));	
	add_menu ();
	add_help ();
}


void DatabaseGrid::add_menu () {
    	QMenu *add = new QMenu(tr("&Add"), this );

	add -> addAction (newFieldAct);
	add -> addAction (addTargetMolAct);
	add -> addAction (loadCsvAct);
	add -> addAction (mergedatabaseAct);
	
	QMenu *edit = new QMenu (tr("&Edit"), this);
/*
	_save_action = new QAction (tr ("&Save"), this);
	connect (_save_action, SIGNAL (triggered ()), this, SLOT (_save_slot ()));
    	file -> addAction (_save_action);

*/

	edit -> addAction (univocalnamesAct);
	addMenu (add);
	addMenu (edit);

/*    	QMenu *settings = new QMenu(tr("&Settings"), this );
	_restore_action = new QAction (tr ("&Restore to default"), this);
	connect (_restore_action, SIGNAL (triggered ()), this, SLOT (_restore_slot ()));
    	settings -> addAction (_restore_action);
	addMenu (settings);
	*/
}

void DatabaseGrid::setup_actions () {
	newFieldAct = new QAction (tr ("&New Field"), this);
	connect (newFieldAct, SIGNAL (triggered ()), this, SLOT (new_field_slot ()));

	loadCsvAct = new QAction (tr ("Import CSV file"), this);
	connect (loadCsvAct, SIGNAL (triggered ()), this, SLOT (load_csv_slot ()));


    toolbar = new QToolBar ("Main");
	addToolBar (toolbar);
	sortupAct = new QAction (tr("Sort up"), this);
    connect (sortupAct, SIGNAL (triggered ()), this, SLOT (sort_up_slot ()));
	
	sortdownAct = new QAction (tr("Sort down"), this);
    connect (sortdownAct, SIGNAL (triggered ()), this, SLOT (sort_down_slot ()));
	
	
	addTargetMolAct = new QAction (tr("Add target Molecule"), this);
	connect (addTargetMolAct, SIGNAL (triggered ()), this,  SLOT (add_target_molecule_slot ()));
	
	mergedatabaseAct = new QAction (tr("Merge Database"), this);
	connect (mergedatabaseAct, SIGNAL (triggered ()), this,  SLOT (merge_database_slot ()));
	
	FiTAct = new QAction (tr("FiT Consensus"), this);
	connect (FiTAct, SIGNAL (triggered ()), this,  SLOT (FiT_slot ()));

	calcAct = new QAction (tr("Calculate descriptors"), this);
	connect (calcAct, SIGNAL (triggered ()), this,  SLOT (calc_slot ()));
	
	univocalnamesAct = new QAction (tr("Univocal Names"), this);
	connect (univocalnamesAct, SIGNAL (triggered ()), this,  SLOT (univocal_names_slot ()));


	toolbar ->addAction (sortupAct);
	toolbar ->addAction (sortdownAct);	
	toolbar ->addAction (FiTAct);
	toolbar ->addAction (calcAct);
}

void DatabaseGrid::manage_number_changed (const QString str) {
	istringstream ss (str.toStdString ());
	int i;
	ss >> i;
	set_current_molecule (i-1);
}


void DatabaseGrid::update_cell (int r, int c) {
//cerr <<" update "<< r << " " << c << endl;
	database ->mutex ->lockForRead ();
//need to make this undoable
if (0<r<database->count_fields () && 0<c<database->count_entries ()) {
	QString st= tab ->item (r, c) ->text ();
	database ->entries [r] -> cells [c] ->set_value (st.toStdString ());
}

	database ->mutex ->unlock ();
	database ->set_needs_redraw (true);
}


void DatabaseGrid::manage_double_click (int r) {
	set_current_molecule (r);
}

void DatabaseGrid::manage_hide (bool b) {
	database ->get_is_hidden () = b;
	set_mol ();
}


void DatabaseGrid::set_current_molecule (int i) {
	if (i>-1 && i< database -> count_entries ()) {
		deselect_row (current_number);
		stringstream ss;
		ss << i+1;
		le -> setText (QString (ss.str ().c_str ()));
		current_number = i;
		select_row (current_number);
		set_mol ();
	}
}



void DatabaseGrid::set_database (Database *db) {
	database = db;
	update_graphics ();
}



void DatabaseGrid::update_graphics () {
//needs deleting all previous cells
	database ->mutex ->lockForRead ();
	if (!database ->count_entries ()) {
		hide_cb ->setEnabled (false);
		hide_cb ->setChecked (true);
	}
	else hide_cb ->setEnabled (true);
//	cerr << database ->count_fields () << " fields      "<<database ->count_entries ()<<"entries"<<endl;
	tab -> setColumnCount (database -> count_fields ());
	tab ->setRowCount (database -> count_entries ());
//	QTableWidgetItem *name = new QTableWidgetItem (QString ("Name"));
//	tab ->setHorizontalHeaderItem (0, name);

	for (unsigned int j = 0; j < database ->field_names.size (); j ++) {
		QTableWidgetItem *field_name = new QTableWidgetItem (QString (database ->field_names[j].c_str ()));
		tab ->setHorizontalHeaderItem (j, field_name);
	}

	for (unsigned int i=0; i< database ->count_entries (); i++) {
//		Database_molecule *mol = database -> molecules [i];

//		QTableWidgetItem *num_widget = new QTableWidgetItem ();
//		num_widget ->setData (Qt::DisplayRole, mol -> number);
//		tab -> setItem (i, 0, num_widget);
		for (unsigned int j = 0; j < database ->count_fields (); j ++) {
			QTableWidgetItem *value_widget = new QTableWidgetItem ();
			value_widget ->setData (Qt::DisplayRole, QString (database -> entries[i] ->cells[j] ->get_string ().c_str ()));
			tab -> setItem (i, j, value_widget);
//			cerr << "insert " << i << " "<<j<<" "<< database -> entries[i] ->cells[j] ->get_string () << endl;
		}
	}
	database ->set_needs_redraw (false);
	database ->mutex ->unlock ();
}

void DatabaseGrid::sort_up_slot () {
	int column = tab ->currentColumn ();
	database -> safe_sort_up (column);
}

void DatabaseGrid::sort_down_slot () {
	int column = tab ->currentColumn ();
	database -> safe_sort_down (column);
}

void DatabaseGrid::first_slot () {
    set_current_molecule (0);
}

void DatabaseGrid::univocal_names_slot () {
	database ->mutex ->lockForWrite ();
	int nf = (int) log10 (database ->count_entries ())+1;
	for (unsigned int i = 0; i < database ->count_entries (); i++) {
		stringstream ss;
		ss << "Molecule_"<<setw(nf) << setfill('0') <<i+1;
		database ->entries[i] ->cells[0] ->set_value (ss.str ());
	}
		database ->set_needs_redraw (true);
		database ->mutex ->unlock ();
	
}

void DatabaseGrid::prev_slot () {
	set_current_molecule (current_number - 1);
}



void DatabaseGrid::next_slot () {
    set_current_molecule (current_number + 1);
}



void DatabaseGrid::last_slot () {
    set_current_molecule ( database->count_entries ()-1);
}

void DatabaseGrid::new_field_slot () {
	vector  <string> dat (database ->count_entries (), "0");
	database ->safe_add_field ("New Field", dat);
}

void DatabaseGrid::add_target_molecule_slot () {
	if (data ->ddwin ->current_target) {
		if (!data ->ddwin ->target_molecule ->multi && !data ->ddwin ->target_molecule ->selection) {
			Database_molecule *new_mol = new Database_molecule (*data ->ddwin ->target_molecule);
			data ->ddwin ->set_lists (new_mol);
			database ->safe_add_mol (new_mol);
		}
	}
}

void DatabaseGrid::merge_database_slot () {
    QString dat_name = QFileDialog::getOpenFileName(this, 
                    tr ("Database"), _data ->ddwin ->last_visited_dir, tr("All Files (*)"));
	if (!dat_name.isEmpty ()) {
		Database *new_dat =_data ->ddwin ->load_multi_file (dat_name.toStdString ());
		database ->safe_merge_with (new_dat);
	}
}


void DatabaseGrid::load_csv_slot () {
    QString file_name = QFileDialog::getOpenFileName(this, 
                    tr ("CSV files"), _data ->ddwin ->last_visited_dir, tr("CSV (*)"));
	if (!file_name.isEmpty ()) {
		database ->import_csv (file_name.toStdString ());
	}
}

void DatabaseGrid::calc_slot () {
	cerr << "calc param" << endl;

	database ->mutex ->lockForWrite ();
	float perc = 0.1;
	int mols = database ->count_entries ();
	int top = (int) (mols *perc);
	vector  <string> dat (mols, "0");
	database ->add_field ("Param", dat);
/*	
	for (unsigned int i = 0; mols; i++) {
		OBMol *mol = new OBMol ();
				int atom = 0;
				FOR_ATOMS_OF_MOL(a, mol)
				{
					atom += a;
					stringstream ss;
					ss << atom;
					string atom_number = ss.str();
	cerr << "atom number:" << ss.toStdString () << endl;
				}
	}



	vector <int> sel_columns = selected_columns ();
	if (sel_columns.size ()) {
		for (unsigned int i = 0; i < sel_columns.size (); i++) {
			database -> sort_up (sel_columns[i]);
			for (unsigned int j = 0; j < top; j++) {
				double d = database ->entries[j] ->cells[database ->count_fields () -1] ->get_double ();
				database ->entries[j] ->cells[database ->count_fields () -1] ->set_double (d + 1);
			}
		}
	}
*/
	database ->set_needs_redraw (true);
	database ->mutex ->unlock ();

}

void DatabaseGrid::FiT_slot () {
	database ->mutex ->lockForWrite ();
	float perc = 0.1;
	int mols = database ->count_entries ();
	int top = (int) (mols *perc);
	vector  <string> dat (mols, "0");
	database ->add_field ("Fit", dat);
	
	vector <int> sel_columns = selected_columns ();
	if (sel_columns.size ()) {
		for (unsigned int i = 0; i < sel_columns.size (); i++) {
			database -> sort_up (sel_columns[i]);
			for (unsigned int j = 0; j < top; j++) {
				double d = database ->entries[j] ->cells[database ->count_fields () -1] ->get_double ();
				database ->entries[j] ->cells[database ->count_fields () -1] ->set_double (d + 1);
			}
		}
	}
	database ->set_needs_redraw (true);
	database ->mutex ->unlock ();
}

vector <int> DatabaseGrid::selected_columns () {
	vector <int> out;
	//database ->mutex ->lockForRead ();
	vector <bool> bools ( tab ->columnCount (),false);
	QList<QTableWidgetItem *> indexs = tab ->selectedItems ();
	 for (unsigned int i = 0; i < indexs.size(); i++) {
		bools[indexs.at (i) ->column ()] = true;
	 }
//	database ->mutex ->unlock ();	
	for (unsigned int i = 0; i < bools.size (); i++) {
		if (bools[i]) out.push_back (i);
	}
	return out;
}

void DatabaseGrid::set_mol () {
	ZNMolecule *mol;
	if (database ->get_is_hidden ()) {mol = database ->dummy_mol; }
	else mol = database ->get_molecule (current_number);
    for (unsigned int i=0; i<data -> ddwin->molecules.size (); i++) {
        if (data -> ddwin->molecules[i]->multi) {
            Database_molecule *dm;
            dm = (Database_molecule *) data -> ddwin->molecules[i];
            if (dm->database == database) {
                data -> ddwin->molecules[i] = mol;
                data -> ddwin->set_current_target (i);
                data -> ddwin->gl->draw_molecule (data -> ddwin->target_molecule);
            }
        }
    }
}

void DatabaseGrid::set_row_color (int r, color c) {
	if (r>-1 && r< tab -> rowCount ()) {
		for (unsigned int i = 0; i< tab -> columnCount (); i++) {
			tab -> item (r, i) -> setBackground (QBrush (c));
		}
	}
}



void DatabaseGrid::select_row (int r) {
	color violet = color (190, 20, 255);
	set_row_color (r, violet);
}



void DatabaseGrid::deselect_row (int r) {
	color white = color (255, 255, 255);
	set_row_color (r, white);
}



int DatabaseGrid::real_index_of_line (int i) {
	return i;
}


//////////////////////////////////////////////////////////////////////////////////////////////////


ColorSettingsMenu::ColorSettingsMenu  (QWidget *parent, DDWin *ddw)
:    QWidget (parent)
{
    ddwin = ddw;
    setWindowTitle ("Color ZNMolecule");
    QVBoxLayout *layout = new QVBoxLayout ();
    setLayout (layout);
	MyColorButton *background_color_edit = new MyColorButton (layout, *ddwin->data->background_color);
}
////////////////////////////////////////////////////////////////////////////////////////////////////

BackboneColorMenu::BackboneColorMenu (QWidget *parent, Data *dat)
:    ZNMenu (parent, dat, 0, 0)
{
	
	add_help ();
	setWindowTitle ("Backbone Color");
	MyHideComboBox *_type_hcb = new MyHideComboBox (main_widget () ->layout (), "Color Scheme");
	ZNWidget *_plain_color_w = new ZNWidget ();	
	ZNWidget *_secondary_structure_w = new ZNWidget ();
	main_widget ()->layout () -> addWidget (_plain_color_w);
	main_widget ()->layout () -> addWidget (_secondary_structure_w);
	
	helix_color = color (1.f, 0.f, 0.f, 1.f);
	sheet_color = color (1.f, 1.f, 0.f, 1.f);
	random_color = color (0.3f, 0.f, 1.f, 1.f);

	
	_type_hcb ->insertItem (_plain_color_w, 0, "Single color", "Single color");
	_type_hcb ->insertItem (_secondary_structure_w, 1, "Secondary Structure", "Secondary Structure");
	
	
	MyCompleteColorSettings *solid_color_settings = new MyCompleteColorSettings (_plain_color_w ->layout (),  constant_color);
	QPushButton *color_ok = new QPushButton ("Ok");
	_plain_color_w->layout () -> addWidget (color_ok);
    connect (color_ok, SIGNAL (clicked ()), this, SLOT (color_ok_slot ()));
	
	MyCompleteColorSettings *helix_color_settings = new MyCompleteColorSettings (_secondary_structure_w ->layout (),  helix_color);
	MyCompleteColorSettings *sheet_color_settings = new MyCompleteColorSettings (_secondary_structure_w ->layout (),  sheet_color);
	MyCompleteColorSettings *random_color_settings = new MyCompleteColorSettings (_secondary_structure_w ->layout (),  random_color);
	QPushButton *ss_ok = new QPushButton ("Ok");
	_secondary_structure_w->layout () -> addWidget (ss_ok);
    connect (ss_ok, SIGNAL (clicked ()), this, SLOT (ss_ok_slot ()));
	
}


void BackboneColorMenu::ss_ok_slot () {
	color_backbone_ss (_data ->ddwin ->target_molecule, helix_color, sheet_color, random_color);
}

void BackboneColorMenu::color_ok_slot () {
	color_backbone_color (_data ->ddwin ->target_molecule, constant_color);
}
////////////////////////////////////////////////////////////////////////////////////////////////////

GraphicalObjectsColorMenu::GraphicalObjectsColorMenu (QWidget *parent, Data *dat)
:    ZNMenu (parent, dat, 0, 0), target (0)
{

	add_help ();
	
	hb_acc_color = color (1.f, 0.f, 0.f, 1.f);
	hb_don_color = color (0.f, 0.f, 1.f, 1.f);
	lipo_color = color (0.f, 1.f, 0.f, 1.f);
	
    setWindowTitle ("Color Graphical Object");

	alpha_multiplier = 100.;
	alpha_value = 100.;
    colortype = new QComboBox;
    main_widget ()->layout () -> addWidget (colortype);
	colortype -> insertItem (0,    "Color"  );
    colortype -> insertItem (1,  "Molecule");
    colortype -> insertItem (2,  "Maps");	
	colortype -> insertItem (3,    "Potentials"  );
    colortype -> insertItem (4,   "Alpha Effects" );

	colortype ->setCurrentItem(0);

	
    options = new QStackedWidget;
    main_widget ()->layout () -> addWidget (options);
	
	ZNWidget *color_widget = new ZNWidget ();
	MyCompleteColorSettings *solid_color_settings = new MyCompleteColorSettings (color_widget ->layout (),  constant_color);

	ZNWidget *molecule_widget = new ZNWidget ();
	ZNWidget *maps_widget = new ZNWidget ();
	
	
	maps_cb = new QComboBox;
	maps_widget ->layout() ->addWidget (maps_cb);

	QWidget *square = new QWidget ();
	QVBoxLayout *mapslay = new QVBoxLayout ();
	square -> setLayout(mapslay);
	maps_widget ->layout() ->addWidget (square);
	QHBoxLayout *score_begin_layout = new QHBoxLayout ();
    mapslay -> addLayout (score_begin_layout);
    QHBoxLayout *score_mid_layout = new QHBoxLayout ();
    mapslay -> addLayout (score_mid_layout);
    QHBoxLayout *score_end_layout = new QHBoxLayout ();
    mapslay -> addLayout (score_end_layout);
	
	score_begin_color = _data ->score_begin_color;
	score_mid_color = _data ->score_mid_color;
	score_end_color = _data ->score_end_color;
	score_begin_f = _data ->score_begin;
	score_mid_f = _data ->score_mid;
	score_end_f = _data ->score_end;
	
    MyFloatEditLine *score_begin_line = new MyFloatEditLine (score_begin_layout, "begin", score_begin_f);
    MyFloatEditLine *score_mid_line = new MyFloatEditLine (score_mid_layout, "mid", score_mid_f);
    MyFloatEditLine *score_end_line = new MyFloatEditLine (score_end_layout, "end", score_end_f);
	
	
    MyCompleteColorSettings *score_begin_colorset = new MyCompleteColorSettings (score_begin_layout,  score_begin_color);
	MyCompleteColorSettings *score_mid_colorset = new MyCompleteColorSettings (score_mid_layout, score_mid_color);
    MyCompleteColorSettings *score_end_colorset = new MyCompleteColorSettings (score_end_layout, score_end_color);
	
	
	ZNWidget *potential_widget = new ZNWidget ();
	MyCompleteColorSettings *lipo_color_settings = new MyCompleteColorSettings (potential_widget ->layout (),  lipo_color, "Lipophilic");
	MyCompleteColorSettings *acc_color_settings = new MyCompleteColorSettings (potential_widget ->layout (),  hb_acc_color, "HB acceptor");
	MyCompleteColorSettings *don_color_settings = new MyCompleteColorSettings (potential_widget ->layout (),  hb_don_color, "HB donor");

	
	
	
	ZNWidget *alpha_widget = new ZNWidget ();
	
	molecule_target_cb = new QComboBox;
	molecule_target_cb2 = new QComboBox;
	molecule_target_cb3 = new QComboBox;
	molecule_widget ->layout () ->addWidget(molecule_target_cb);
	alpha_widget ->layout () ->addWidget(molecule_target_cb2);
	potential_widget ->layout () ->addWidget(molecule_target_cb3);
	options ->addWidget (color_widget);
	options ->addWidget (molecule_widget);
	options ->addWidget (maps_widget);
	options ->addWidget (potential_widget);
	options ->addWidget (alpha_widget);
	
	alpha_mol_distance_d = 4.;
	new MyFloatEditLine (alpha_widget ->layout (), "Fade distance", alpha_mol_distance_d);
		QPushButton *fade_but = new QPushButton ("Ok");
	alpha_widget ->layout() ->addWidget (fade_but);
	connect (fade_but, SIGNAL (clicked ()), this, SLOT (fade_slot ()));

	

	

	new MyFloatEditLine (alpha_widget ->layout (), "Multiply Alpha by %", alpha_multiplier, 0, 100);
	QPushButton *perc_mult_but = new QPushButton ("Ok");	
	alpha_widget ->layout() ->addWidget(perc_mult_but);
	connect (perc_mult_but, SIGNAL (clicked ()), this, SLOT (mult_slot ()));
	
	new MyFloatEditLine (alpha_widget ->layout (), "Set Alpha", alpha_value, 0, 100);
	QPushButton *alpha_but = new QPushButton ("Ok");	
	alpha_widget ->layout() ->addWidget(alpha_but);
	connect (alpha_but, SIGNAL (clicked ()), this, SLOT (alpha_slot ()));
	


	
    connect (colortype, SIGNAL(activated(int)), options, SLOT(setCurrentIndex(int)) );
	connect (_data -> ddwin, SIGNAL (targets_updated ()), this, SLOT (update_mols ()));
	connect (_data -> ddwin, SIGNAL (go_updated ()), this, SLOT (update_maps ()));  //temporary... while maps are still considered gos

	
	QPushButton *ok = new QPushButton ("Ok");
	ZNWidget *buttons = new ZNWidget ();
    buttons -> layout () ->addWidget (ok);
	connect (ok, SIGNAL (clicked ()), this, SLOT (ok_slot ()));

	main_widget ()->layout () -> addWidget (buttons);
}



void GraphicalObjectsColorMenu::update_mols () {
	molecule_target_cb ->clear ();
	molecule_target_cb2 ->clear ();
	molecule_target_cb3 ->clear ();
	for (unsigned int i = 0; i < _data -> ddwin -> target ->count (); i++) {
		molecule_target_cb ->insertItem (i, _data -> ddwin -> target ->itemText (i));
		molecule_target_cb2 ->insertItem (i, _data -> ddwin -> target ->itemText (i));
		molecule_target_cb3 ->insertItem (i, _data -> ddwin -> target ->itemText (i));
	}

	
	
}	

void GraphicalObjectsColorMenu::update_maps () {
	maps_cb ->clear ();
	for (unsigned int i = 0; i < _data -> ddwin -> graphical_objects.size (); i++) { //temporary... while maps are still considered gos
		if (_data -> ddwin -> graphical_objects[i] ->is_map ()) {
			maps_cb ->insertItem (i, QString (_data -> ddwin -> graphical_objects[i] ->name.c_str ()));
		}
	}
	
	
	
}	

void GraphicalObjectsColorMenu::mult_slot () {
	if (target) {
		target ->multiply_alpha (alpha_multiplier/100);
		target ->render ();
	}
}

void GraphicalObjectsColorMenu::alpha_slot () {
	if (target) {
		target ->set_alpha (alpha_value/100);
		target ->render ();
	}
}

void GraphicalObjectsColorMenu::fade_slot () {
	if (target) {
		int indx = molecule_target_cb2 ->currentIndex();
		if (indx) {
			ZNMolecule *mol = _data -> ddwin -> molecules [indx];
			target ->alpha_by_mol_distance (mol, alpha_mol_distance_d);
			target ->render ();			
		}
	}
}

void GraphicalObjectsColorMenu::ok_slot () {
	if (target) {
	if ( options -> currentIndex () == 0) {
		target ->color_by_color (constant_color);
		target ->render ();
		
	}
	else if ( options -> currentIndex () == 1) {
		int indx = molecule_target_cb ->currentIndex();
		if (indx) {
			ZNMolecule *mol = _data -> ddwin -> molecules [indx];
			target ->color_by_mol (mol);
			target ->render ();			
		}
	}
		
	else if ( options -> currentIndex () == 2) {
		int indx = maps_cb ->currentIndex();


			if (_data -> ddwin -> graphical_objects [indx] ->is_map ()) {
				Map *map = (Map *) (_data -> ddwin -> graphical_objects [indx]);
				target ->color_by_map (map, score_begin_color, score_mid_color, score_end_color, score_begin_f, score_mid_f, score_end_f);
				target ->render ();			
			}

	}
	else if ( options -> currentIndex () == 3) {
		int indx = molecule_target_cb3 ->currentIndex();
		if (indx) {
			ZNMolecule *mol = _data -> ddwin -> molecules [indx];
			target ->color_by_potential (mol, lipo_color, hb_acc_color, hb_don_color);
			target ->render ();			
		}
	}
		

	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////


ColorMenu::ColorMenu (QWidget *parent, DDWin *ddw)
   :    QWidget (parent)
{
    ddwin = ddw;
    setWindowTitle ("Color Molecule");
    QVBoxLayout *layout = new QVBoxLayout ();
    setLayout (layout);


    colortype = new QComboBox;
    layout -> addWidget (colortype);
    colortype -> insertItem (ELEMENT,  "Element");
    colortype -> insertItem (CHARGE,   "Charge" );
    colortype -> insertItem (SCORE,    "Score"  );
    colortype -> insertItem (COLOR,    "Color"  );

    options = new QStackedWidget;
    layout -> addWidget (options);

    QWidget *element_options = new QWidget ();
    

    QWidget *charge_options = new QWidget ();
    QVBoxLayout *charge_layout = new QVBoxLayout ();
    charge_options -> setLayout (charge_layout);
    QHBoxLayout *charge_begin_layout = new QHBoxLayout ();
    charge_layout -> addLayout (charge_begin_layout);
    QHBoxLayout *charge_end_layout = new QHBoxLayout ();
    charge_layout -> addLayout (charge_end_layout);

    charge_begin_line = new MyFloatEditLine (charge_begin_layout, "begin", ddwin->data->charge_begin);
    charge_end_line = new MyFloatEditLine (charge_end_layout, "end", ddwin->data->charge_end);

/* unused variable
    MyColorButton *charge_begin_color = new MyColorButton (charge_begin_layout, ddwin->data->charge_begin_color);
    MyColorButton *charge_end_color = new MyColorButton (charge_end_layout, ddwin->data->charge_end_color);
*/

    QWidget *score_options = new QWidget ();
    QVBoxLayout *score_layout = new QVBoxLayout ();
    score_options -> setLayout (score_layout);
    QHBoxLayout *score_begin_layout = new QHBoxLayout ();
    score_layout -> addLayout (score_begin_layout);
    QHBoxLayout *score_mid_layout = new QHBoxLayout ();
    score_layout -> addLayout (score_mid_layout);
    QHBoxLayout *score_end_layout = new QHBoxLayout ();
    score_layout -> addLayout (score_end_layout);

    score_begin_line = new MyFloatEditLine (score_begin_layout, "begin", ddwin->data->score_begin);
    score_mid_line = new MyFloatEditLine (score_mid_layout, "mid", ddwin->data->score_mid);
    score_end_line = new MyFloatEditLine (score_end_layout, "end", ddwin->data->score_end);


    MyCompleteColorSettings *score_begin_color = new MyCompleteColorSettings (score_begin_layout,  ddwin->data->score_begin_color);
    MyCompleteColorSettings *score_mid_color = new MyCompleteColorSettings (score_mid_layout, ddwin->data->score_mid_color);
    MyCompleteColorSettings *score_end_color = new MyCompleteColorSettings (score_end_layout, ddwin->data->score_end_color);


    QWidget *color_options = new QWidget ();
    QVBoxLayout *color_layout = new QVBoxLayout ();
    color_options -> setLayout (color_layout);


    MyCompleteColorSettings *select_color_buttons = new MyCompleteColorSettings (color_layout, ddwin->data->constant_color);


    options -> addWidget (element_options);
    options -> addWidget (charge_options);
    options -> addWidget (score_options);
    options -> addWidget (color_options);

    layout -> addSpacing (5);

    QPushButton *multiple = new QPushButton ("Multiple colors");
    layout -> addWidget (multiple);

    QListView *list = new QListView ();
    layout ->addWidget (list);

    QHBoxLayout *buttons = new QHBoxLayout ();
    layout -> addLayout (buttons);

    QPushButton *ok = new QPushButton ("Ok");
    buttons -> addWidget (ok);


    connect (colortype, SIGNAL(activated(int)), options, SLOT(setCurrentIndex(int)) );
    connect (ok, SIGNAL (clicked ()), this, SLOT (ok_slot ()));
}



void ColorMenu::ok_slot () {
    score_begin_line -> set ();    
    score_mid_line -> set ();
    score_end_line -> set ();
    charge_begin_line -> set ();
    charge_end_line -> set ();


    vector <color_mask> masks;
    color_mask mask;
    mask.intensity = 1.0f;
    mask.only_to = 0;
    mask.excluding = 0;
    mask.type = options -> currentIndex ();
    masks.push_back (mask);
	if (is_db_extended (ddwin ->target_molecule)) {
		Database_molecule *dbm = (Database_molecule *) ddwin -> target_molecule;
		Database *db = dbm ->database;
		ddwin -> data -> actions -> apply_color_masks (masks, db);
	}
	else {
		ddwin -> data -> actions -> apply_color_masks (masks, ddwin ->target_molecule);
	}
}



GraphicalObjectsMenu::GraphicalObjectsMenu (QWidget *parent, DDWin *ddw)
   :    QWidget (parent)
{
  //  selected = -1;
    ddwin = ddw;
    setWindowTitle ("Graphical Objects");
    QVBoxLayout *layout = new QVBoxLayout ();
    setLayout (layout);

    list = new QListWidget ();
    layout -> addWidget (list);

    QPushButton *del = new QPushButton ("Delete");
    layout -> addWidget (del);

    connect (del, SIGNAL (clicked ()), this, SLOT (delete_selected_slot ()));

  connect(list,SIGNAL(itemDoubleClicked(QListWidgetItem*)),this,SLOT(show_color_menu(QListWidgetItem*)));
}

void GraphicalObjectsMenu::show_color_menu (QListWidgetItem *item) {
	int row = list ->currentIndex ().row ();
	ddwin ->go_color_menu ->set_target(ddwin ->graphical_objects[row]);
	ddwin ->go_color_menu ->display ();
}


void GraphicalObjectsMenu::update_slot () {
    int row = list -> currentRow ();
    bool selected = list -> isItemSelected (list -> currentItem ());
    list -> clear ();
    for (unsigned int i=0; i < ddwin -> graphical_objects.size (); i++) {
        list -> insertItem (i, tr (ddwin -> graphical_objects [i] -> name .c_str ()));
    }
    list -> setCurrentRow (row);
    list -> setItemSelected (list -> currentItem (), selected);
}



void GraphicalObjectsMenu::delete_selected_slot () {
    int i = list -> currentRow ();
    
    if (i > -1) {
        ddwin -> delete_graphical_object (i);

    }
}



void GraphicalObjectsMenu::paintEvent (QPaintEvent *e) {
    QWidget::paintEvent (e);
    update_slot ();
}



DDSettingsMenu::DDSettingsMenu (QWidget *parent, DDWin *ddw)
   :    QWidget (parent)
{
    ddwin = ddw;
    setWindowTitle ("3D Settings");
    QVBoxLayout *layout = new QVBoxLayout ();
    setLayout (layout);

    focal_d = ddwin ->  gl -> stereo_inter_eye_semi_distance * tan ((90.f -ddwin -> gl -> stereo_toe_in_angle) * PI / 180) ;


    inter_eye_distance = new MyFloatEditLine (layout, "inter-eye semi distance", ddwin->gl->stereo_inter_eye_semi_distance);
    focal_point_distance = new MyFloatEditLine (layout, "focal point distance", focal_d);

    dd_cb = new QComboBox ();
    layout -> addWidget (dd_cb);
    dd_cb ->insertItem(0, "No Stereo" );
    dd_cb ->insertItem(1, "Quad Buffering" );
    dd_cb ->insertItem(2, "Vertical Interlace" );
    dd_cb ->insertItem(3, "Horizontal Interlace" );

    // IJG
    // Need a full set of radio buttons: Use OpenGL, horizontal interlace, vertical interlace, none
    // See DDWin::StereoMode and DDWin::g_stereoMode


    QPushButton *ok = new QPushButton ("Ok");
    layout -> addWidget (ok);


    connect (ok, SIGNAL (clicked ()), this, SLOT (ok_slot ()));
}



void DDSettingsMenu::ok_slot () {
    inter_eye_distance -> set ();
    focal_point_distance -> set ();
    ddwin -> gl -> stereo_toe_in_angle = 90.f - atan (focal_d / ddwin -> gl -> stereo_inter_eye_semi_distance) * 180 / PI;
	switch (dd_cb ->currentIndex ()) {
		case 0:
			ddwin->g_stereoMode = DDWin::StereoMode_None;
			ddwin->g_stencil_mask_needs_redraw = false;
			break;
		case 1:
			ddwin->g_stereoMode = DDWin::StereoMode_ViaOpenGL;
			ddwin->g_stencil_mask_needs_redraw = false;
			break;
		case 2:
//<<<<<<< menu.cc
			ddwin->g_stereoMode = DDWin::StereoMode_HorizontalInterlace;
			ddwin->g_stencil_mask_needs_redraw = true;
//=======
			ddwin->g_stereoMode = DDWin::StereoMode_VerticalInterlace;
//>>>>>>> 1.7
			break;
		case 3:
//<<<<<<< menu.cc
			ddwin->g_stereoMode = DDWin::StereoMode_VerticalInterlace;
			ddwin->g_stencil_mask_needs_redraw = true;
//=======
			ddwin->g_stereoMode = DDWin::StereoMode_HorizontalInterlace;
//>>>>>>> 1.7
			break;
		default:
			ddwin->g_stereoMode = DDWin::StereoMode_None;
			ddwin->g_stencil_mask_needs_redraw = false;
			break;
	}

    // Also need to update ddwin - so it sets up the stencil buffer.
    ddwin -> gl ->refreshStencilBuffer();
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	"My" function
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


MyCheckBox::MyCheckBox (QLayout *parent, bool &v, string str) : QCheckBox ( QString (str.c_str ())) {
    parent -> addWidget (this);
    var = &v;
    get ();
    connect (this, SIGNAL (stateChanged (int)), this, SLOT (set ()) );
} 


void MyCheckBox::get () {
    if (*var) setCheckState (Qt::Checked);
    else setCheckState (Qt::Unchecked);
}


void MyCheckBox::set () {
    if (checkState () == Qt::Checked) *var = TRUE;
    else *var = false;
}



MyColorButton::MyColorButton (QLayout *parent, QColor &col) : QPushButton () {
	parent -> addWidget (this);
	color = &col;
	connect (this, SIGNAL (clicked ()), this, SLOT (my_clicked ()) );
}


void MyColorButton::paintEvent (QPaintEvent *e) {
	QPushButton::paintEvent (e);
    QPixmap pix(48, 48);
    pix.fill(*color);
    this->setIcon(pix);

}


void MyColorButton::my_clicked () {
	QColor new_color = QColorDialog::getColor(*color, this );
	if ( new_color.isValid () ) {
		new_color.setAlphaF (color ->alphaF ());
		*color = new_color;
	}
}



MyComboBox::MyComboBox (QLayout *parent, const char *name) {
	parent -> addWidget (this);
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);
	QLabel *label = new QLabel (name);
	layout -> addWidget (label);

	layout ->setContentsMargins (0, 0, 0, 0);

	_combo_box = new QComboBox ();
	layout -> addWidget (_combo_box);
}



MyCompleteColorSettings::MyCompleteColorSettings (QLayout *parent, color &col, string name) {
	setLayout (new QHBoxLayout ());
	button = new MyColorButton (layout (), col);
	alpha = (int) col.alphaF () * 100;
	slider = new MySlider (layout (), "opacity", alpha, 0, 100);
	if (name!=""){
		QLabel *name_l =  new QLabel (name.c_str ());
		layout () -> addWidget (name_l);
	}
	connect  (slider -> pline -> spinbox, SIGNAL (valueChanged(int)), this, SLOT(setAlpha (int)));
	parent ->addWidget(this);
}


void MyCompleteColorSettings::setAlpha (int i) {
	double v = (double) i;
	v /= 100.;
	button -> color -> setAlphaF (v);
}



MyFloatEditLine::MyFloatEditLine (QLayout *parent, const char *name, double& var, double min, double max)
	: QWidget(){
	parent -> addWidget (this); 
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);
	variable = &var;
	QLabel *label = new QLabel( name);
	layout -> addWidget (label); 
//	label->setMaximumWidth( 200 );
	label->setMinimumWidth( 150 );
	spinbox = new QDoubleSpinBox();
	spinbox ->setMinimum (min);
	spinbox ->setMaximum (max);
	layout -> addWidget (spinbox);
	spinbox->setValue (*variable);
	connect (spinbox, SIGNAL (valueChanged (double)), this, SLOT (set (double)));

	layout ->setContentsMargins (0, 0, 0, 0);
//	spinbox ->setMinimumWidth( 150 );
	spinbox ->setMaximumWidth( 150 );
}


void MyFloatEditLine::set (double d) {
	*variable = d;   
}


void  MyFloatEditLine::set () {
	*variable = spinbox -> value ();   
}


void MyFloatEditLine::set_value (double d) {
	spinbox->setValue (d);    
}



MyGroupBox::MyGroupBox (QLayout *parent, string str, bool checkbox, bool status) : QGroupBox ( QString (str.c_str ()) ) {
	parent -> addWidget (this);
	QVBoxLayout *layout = new QVBoxLayout ();
	this ->setLayout(layout);
	if (checkbox) {
		setCheckable (true);
		if (status) {
			setChecked (true);
		} else {
			setChecked (false);
		}
	}
	layout ->setContentsMargins (5, 5, 5, 5);
 	layout -> setAlignment (Qt::AlignTop);
}

MyGroupBox::MyGroupBox (string str, bool checkbox, bool status) : QGroupBox ( QString (str.c_str ()) ) {
	QVBoxLayout *layout = new QVBoxLayout ();
	this ->setLayout(layout);
	if (checkbox) {
		setCheckable (true);
		if (status) {
			setChecked (true);
		} else {
			setChecked (false);
		}
	}
	layout ->setContentsMargins (5, 5, 5, 5);
 	layout -> setAlignment (Qt::AlignTop);
}


MyHideCheckBox::MyHideCheckBox (QLayout *parent, QWidget *tar, string str) : QCheckBox ( QString (str.c_str ())) {
	_group_box = new MyGroupBox (parent -> layout (), str, true, false);
	_group_box -> hide ();
	parent -> addWidget (this);
	target = tar;
	get ();
	connect (this, SIGNAL (stateChanged (int)), this, SLOT (set ()) );
	connect (_group_box, SIGNAL (clicked (bool)), this, SLOT (set_2 ()) ); 
	_group_box -> layout () -> addWidget (target);
	set ();
} 


void MyHideCheckBox::get () {
    if (target ->isVisible ()) set_checked (); //setCheckState (Qt::Checked);
    else set_unchecked (); //setCheckState (Qt::Unchecked);
}


void MyHideCheckBox::set () {
    if (checkState () == Qt::Checked)  set_checked ();//target -> show ();
    else set_unchecked (); //target -> hide ();
}


void MyHideCheckBox::set_2 () {
	set_unchecked ();
} 


void  MyHideCheckBox::set_checked () {
	hide ();
	_group_box -> show ();
	_group_box ->setChecked (true);
}


void  MyHideCheckBox::set_unchecked () {
	_group_box -> hide ();
	show ();
	_group_box -> hide ();
	setCheckState (Qt::Unchecked);
}



MyHideComboBox::MyHideComboBox (QLayout *parent, const char *name) : MyComboBox (parent, name){
	_widgets.clear ();
	connect (_combo_box, SIGNAL (currentIndexChanged ( int )), this, SLOT (set ( int )) );
}


void MyHideComboBox::insertItem (QWidget *wid, int i, const char *name, const QVariant &data) {

	unsigned int len = _widgets.size ();
	vector<QWidget *>::iterator it;
  	it = _widgets.begin();
	if (i < 1) i = 0;
	else if(i > len) i = len;

	_widgets.insert (it + i, wid);
	_combo_box -> insertItem (i, name, data);
	set (0);
}


void MyHideComboBox::set (int i) {
	if (_combo_box -> count () == _widgets.size ()) {
		for (unsigned int n = 0; n < _widgets.size (); n++) {
			if (n == i) _widgets[n] ->show ();
			else _widgets[n] -> hide ();
		} 
	}
	else cerr << "error in myhidecombobox"<< endl;
}



MyIntegerEditLine::MyIntegerEditLine (QLayout *parent, const char *name, int& var, int min, int max)
	: QWidget(){
	parent -> addWidget (this); 
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);
	variable = &var;
	QLabel *label = new QLabel( name);
	layout -> addWidget (label); 
//	label->setMaximumWidth( 200 );
	label->setMinimumWidth( 150 );
	spinbox = new QSpinBox();
	spinbox ->setMinimum (min);
	spinbox ->setMaximum (max);
	layout -> addWidget (spinbox);
//    stringstream s;
//    s << *variable;
	spinbox->setValue(*variable);
	connect (spinbox, SIGNAL (valueChanged (int)), this, SLOT (set (int)));
	layout ->setContentsMargins (0, 0, 0, 0);
	spinbox ->setMaximumWidth( 150 );
}


void MyIntegerEditLine::set (int d) {
//    istringstream iss (st.toStdString());
//    int i;
//    iss >> i;
	*variable = d;   
}


void MyIntegerEditLine::set ()
{
//	istringstream iss (linedit->text().toStdString());
//	int i;
//	iss >> i;
//	*variable = i;   
	*variable = spinbox -> value (); 
}


void MyIntegerEditLine::set_value (int d) {
//    stringstream s ;
//    s<<v;
//    linedit->setText (QString(s.str().c_str())); 
//    set ();    
	spinbox->setValue (d); 
}



MyLabelf::MyLabelf (QLayout *parent, const char *name, float& var)  
	: QWidget () {
	parent -> addWidget (this);
	setLayout (new QHBoxLayout ()); 
	variable = &var;
	QLabel *name_l = new QLabel ();
	layout () ->addWidget (name_l);
	string nam = name;
	name_l->setText (QString(nam.c_str()));
	label = new QLabel();
	layout ()->addWidget (label); 
	//label->setMaximumWidth( 200 );
	//label->setMinimumWidth( 200 );
	label->setText (QString(double_to_string (*variable).c_str ()));

	update ();
}


void MyLabelf::update () {
//	cerr << variable<<endl;
//		cerr << *variable<<endl;
//	float v = 30;
//	ss << v;
//	cerr << v<<endl;
//	cerr << label<< endl;
//	label ->setNum (3);
	label->setText (QString(double_to_string (*variable).c_str ()));
}


void MyLabelf::set_variable (float *f) {
	variable = f;
}

MyLineEdit::MyLineEdit (QLayout *parent, const char *name) : QWidget() {

	parent -> addWidget (this);
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);

	layout ->setContentsMargins (0, 0, 0, 0);

	tag = name;
	label = new QLabel( name, this ); 
	label->setMaximumWidth( 150 );
	label->setMinimumWidth( 150 );

	linedit = new QLineEdit( this);

	layout -> addWidget (label);
	layout -> addWidget (linedit);

}


MyPushButton::MyPushButton (QLayout *parent, const char *name, int dim, int width) : QWidget() {
	parent ->addWidget (this);
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);

	layout ->setContentsMargins (0, 0, 0, 0);
 	layout -> setAlignment (Qt::AlignTop);

	fbutton = new QPushButton(name, this);
	fbutton ->setMaximumHeight( 24 );
	if (dim != 0 ) fbutton ->setMinimumWidth( dim );
	if (width != 0) fbutton ->setMaximumWidth( width );
	layout -> addWidget (fbutton);
}

MyPushButton::MyPushButton (const char *name, int dim, int width) : QWidget() {
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);

	layout ->setContentsMargins (0, 0, 0, 0);
 	layout -> setAlignment (Qt::AlignTop);

	fbutton = new QPushButton(name, this);
	fbutton ->setMaximumHeight( 24 );
	if (dim != 0 ) fbutton ->setMinimumWidth( dim );
	if (width != 0) fbutton ->setMaximumWidth( width );
	layout -> addWidget (fbutton);
}


MyLineFile::MyLineFile (QLayout *parent, const char *name, int valid) : QWidget() {

	parent -> addWidget (this);
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);

	layout ->setContentsMargins (0, 0, 0, 0);

	tag = name;
	label = new QLabel( name, this ); 
//	label->setMaximumWidth( 150 );
	label->setMinimumWidth( 150 );

	linedit = new QLineEdit( this);
	linedit ->setMinimumWidth( 250 );
//	if (valid == 1) {filetype = "Tripos Mol2 File (*.mol2)";};
//	if (valid == 2) {filetype = "PLANTS Config File (*.pcfg);;All files (*)";};
//	if (valid == 3) {filetype = "List File (*)";};
	fbutton = new QPushButton(  "..." , this);
	fbutton ->setMaximumWidth( 30 );
	fbutton ->setMaximumHeight( 24 );

	control_yes = new QLabel( this ); 
	control_yes ->setPixmap (QPixmap (":icons/V.png") );
	control_yes -> hide();

	control_no = new QLabel( this ); 
	control_no ->setPixmap (QPixmap (":icons/X.png") );
	control_no -> hide();

	layout -> addWidget (label);
	layout -> addWidget (control_no);
	layout -> addWidget (control_yes);
	layout -> addWidget (linedit);
	layout -> addWidget (fbutton);



//	connect (fbutton, SIGNAL( clicked() ), SLOT( set_file () ) );
	if (valid == 1) {
		connect (fbutton, SIGNAL( clicked() ), SLOT( set_file_a () ) );
	}
	if (valid == 2) {
		connect (fbutton, SIGNAL( clicked() ), SLOT( set_file_b () ) );
	}
	if (valid == 3) {
		connect (fbutton, SIGNAL( clicked() ), SLOT( set_file_c () ) );
	}
	if (valid == 4) {
		connect (fbutton, SIGNAL( clicked() ), SLOT( set_file_d () ) );
	}

}

void MyLineFile::set_file_a () {
	QString s = QFileDialog::getOpenFileName(this, tr ("Open file"), "",tr("Tripos Mol2 File (*.mol2)"));
	linedit->clear ();
	linedit->insert (s);
//	set_line ();
}

void MyLineFile::set_file_b () {
	QString s = QFileDialog::getOpenFileName(this, tr ("Open file"), "",tr("PLANTS Config File (*.pcfg);;All files (*)"));
	linedit->clear ();
	linedit->insert (s);
//	set_line ();
}

void MyLineFile::set_file_c () {
	QString s = QFileDialog::getOpenFileName(this, tr ("Open file"), "",tr("List File (*)"));
	linedit->clear ();
	linedit->insert (s);
//	set_line ();
}

void MyLineFile::set_file_d () {
	QString s = QFileDialog::getOpenFileName(this, tr ("Open file"), "",tr("PLANTS (plants*)"));
	linedit->clear ();
	linedit->insert (s);
//	set_line ();
}

/*
QString MyLineFile::ask_file() {
	QString mol_name = QFileDialog::getOpenFileName(this, tr ("Open file"), "",tr("Tripos Mol2 File (*.mol2)"));
	if (valid == 3) {cerr << "www" << endl;};

	return mol_name;
}

void MyLineFile::set_file () {
	QString s = ask_file ();
	linedit->clear ();
	linedit->insert (s);
//	set_line ();
}
*/
string MyLineFile::val () {
	return linedit->text ().toStdString ();
}


MyListView::MyListView (QLayout *parent, int dim, int width, bool button) : QWidget() {

	parent -> addWidget (this);
	QVBoxLayout *layout = new QVBoxLayout ();
	setLayout (layout);

 	layout -> setAlignment (Qt::AlignTop);

	layout ->setContentsMargins (0, 0, 0, 0);

	_lw = new QListWidget ();
	_lw ->sortItems(Qt::AscendingOrder);

	shortcut1 = new QShortcut (this);
	shortcut1 ->setKey(Qt::Key_Delete);
	connect (shortcut1, SIGNAL( activated() ), SLOT( del_list_view_slot() ) );

	shortcut2 = new QShortcut (this);
	shortcut2 ->setKey(Qt::Key_Backspace);
	connect (shortcut2, SIGNAL( activated() ), SLOT( del_list_view_slot() ) );

	layout -> addWidget (_lw);
	if (dim != 0) {
		_lw ->setMaximumHeight( dim );
	}
	else {
		_lw ->setMaximumHeight( 80 );
	}
	if (width != 0) _lw ->setMaximumWidth( width );

	if (button == true) {
		fbutton = new QPushButton("Delete selected", this);
		fbutton ->setMaximumHeight( 24 );
		layout -> addWidget (fbutton);

		connect (fbutton, SIGNAL( clicked() ), SLOT( del_list_view_slot() ) );
	}
}

void MyListView::del_list_view_slot () {
	int _last_row = _lw ->currentRow();
	_lw ->takeItem (_last_row);
	deleting (_last_row);
}


MyListButton::MyListButton (QLayout *parent, const char *name) : QWidget() {
	parent -> addWidget (this);
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);

	layout ->setContentsMargins (0, 0, 0, 0);

	_combo_box = new QComboBox ();
	layout -> addWidget (_combo_box);

	fbutton = new QPushButton( name , this);
//	fbutton ->setMaximumWidth( 30 );
	fbutton ->setMaximumHeight( 24 );
	layout -> addWidget (fbutton);
}


MySlider::MySlider (QLayout *parent, const char *name,  int& var, int vmin, int vmax) : QWidget(){
	parent -> addWidget (this);
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);
	pline = new MyIntegerEditLine (layout, name, var, vmin, vmax);
	layout -> addWidget (pline);
	slider = new QSlider(Qt::Horizontal);
	layout -> addWidget (slider);
	slider->setMinimum ( vmin );
	slider->setMaximum ( vmax );
	slider ->setValue (var);
	connect (slider, SIGNAL(valueChanged(int)), pline, SLOT(set_value(int)) );
	connect (pline ->spinbox, SIGNAL (valueChanged(int)), this, SLOT(setValue (int)));
}


void MySlider::setValue (int i) {
	slider->setValue (i);
}


MyTwoColumn::MyTwoColumn (QLayout *parent) : QWidget () {
	parent -> addWidget (this);
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);

	_left = new ZNWidget ();
	_right = new ZNWidget ();

	layout -> addWidget (_left);
	layout -> addWidget (_right);
 	layout -> setAlignment (Qt::AlignTop);

	layout ->setContentsMargins (0, 0, 0, 0);
	_left ->layout() ->setContentsMargins (0, 0, 5, 0);
 	_left ->layout () -> setAlignment (Qt::AlignTop);
	_right ->layout() ->setContentsMargins (5, 0, 0, 0);
	_right ->layout () -> setAlignment (Qt::AlignTop);
}


My3Column::My3Column (QLayout *parent) : QWidget () {
	parent -> addWidget (this);
	QHBoxLayout *layout = new QHBoxLayout ();
	setLayout (layout);

	_left = new ZNWidget ();
	_center = new ZNWidget ();
	_right = new ZNWidget ();

	layout -> addWidget (_left);
	layout -> addWidget (_center);
	layout -> addWidget (_right);
 	layout -> setAlignment (Qt::AlignTop);

	layout ->setContentsMargins (0, 0, 0, 0);
	_left ->layout() ->setContentsMargins (0, 0, 5, 0);
	_left ->layout () -> setAlignment (Qt::AlignTop);
	_center ->layout() ->setContentsMargins (5, 0, 5, 0);
	_center ->layout () -> setAlignment (Qt::AlignTop);
	_right ->layout() ->setContentsMargins (5, 0, 0, 0);
	_right ->layout () -> setAlignment (Qt::AlignTop);
}

MyGridColumn::MyGridColumn (QLayout *parent, int row, int col) : QWidget () {
	parent -> addWidget (this);
	gridlayout = new QGridLayout ();
	setLayout (gridlayout);

	gridlayout ->setHorizontalSpacing (5);
	gridlayout ->setVerticalSpacing (5);

	gridlayout ->setContentsMargins (0, 0, 0, 0);
}


MyTableWidget::MyTableWidget (QLayout *parent, int row, int height) : QWidget () {
	parent -> addWidget (this);
	QVBoxLayout *layout = new QVBoxLayout ();
	setLayout (layout);
 	layout -> setAlignment (Qt::AlignTop);

	table = new QTableWidget ();
	table ->setColumnCount (row);
	table ->setMaximumHeight( height );
	layout -> addWidget (table);

}

