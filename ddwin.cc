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


//modes
#define NONE 0
#define SELECT 1
#define BUILDER 2
#define MAGIC_PENCIL 3
#define SELECT_EXTEND 4


//color masks
#define ELEMENT 0
#define CHARGE  1
#define SCORE   2
#define COLOR   3



//using namespace OpenBabel;



/*

void DDWin::read_grid_file (string filename, Grid &grid) {
    float min =0;
    float max =0;
    bool line_read = false;
    grid.vec.clear ();
	string line ="";
    ifstream file (filename.c_str ());
	line_read = getline(file, line);

    istringstream iss (line);
    iss >> grid.OX;
    iss >> grid.OY;
    iss >> grid.OZ;
    iss >> grid.NX;
    iss >> grid.NY;
    iss >> grid.NZ;
    iss >> grid.RES;


	for (unsigned int x =0; x<grid.NX; x++) {
    	for (unsigned int y =0; y<grid.NY; y++) {   
        	for (unsigned int z =0; z<grid.NZ; z++) {

	            line_read = getline(file, line);
                istringstream iss (line);
                float val =0;
                iss >> val;   

                if (val) { 
                    if (x==0 && y==0 && z==0) {min =val; max = val;}
                    if (min > val) min =val;
                    if (max < val) max =val;
                    Grid_dot dot;
                    dot.val = val;
                    dot.x() = grid.OX+x*grid.RES;
                    dot.y() = grid.OY+y*grid.RES;
                    dot.z() = grid.OZ+z*grid.RES;
                    grid.vec.push_back (dot);
                }
            }
        }
    }
    grid.min = min; grid.max=max;
    cout<< "min "<<min<<" max "<<max;

}


*/

    GLfloat x = 800.0f;    GLfloat y = 600.0f;
    ArcBallT  ArcBall(x,y);
Point2fT    MousePt;	
const float PI2 = 2.0*3.1415926535f;
GLUquadricObj *quadratic;






DDWin::DDWin (QWidget *parent, Data *dat)
:    QMainWindow (parent), g_stencil_mask(0), g_stencil_mask_needs_redraw(false), g_stencil_mask_frame_counter(0),
g_stereoMode(StereoMode_None)
{
    setAcceptDrops(true) ; 
	last_visited_dir = "";
	edit_lock_counter = 0;
    setMinimumWidth (400);
    dat->set_ddwin (this);
    accel = new Q3Accel (this);
    setup_accel ();
    setWindowTitle( QString(TITLE.c_str()) + QString(VERSION.c_str()) );
    setWindowIcon ((QPixmap (":icons/zeden_ico.png")));
	connect (this, SIGNAL (non_selection_molecules_updated ()), SLOT (emit_targets_updated()));


 //   QGLFormat f;
 //   f.setStereo (true);
//	f.setSampleBuffers (true);

//    QGLFormat::setDefaultFormat (f);
    gl = new MyGl (this);
    assert (gl);

    setCentralWidget (gl);



    undo_view = new QUndoView ();
    undo_view -> setStack (data -> undo_stack);


        ZNMolecule *all = new ZNMolecule;
        assert (all);
	all ->ZNinit();
        all -> SetTitle ("all");
        molecules.push_back (all);
        target_molecule = all;
        current_target = 0;
        backbone_shown = true;
        haptic_mode = false;
        minimizing = false;

		adding_restrains = false;
        rec_movie = false;
        rec_pov_movie = false;
        show_labels = true;
        last_frame = 0;
        mode = NONE;

        resize (800,600);
        builder = new Builder (this);

        set_popups ();

    draw_menu ();
}



void DDWin::closeEvent(QCloseEvent *event)
 {
    QWidget::closeEvent(event);
    cout << "quitting " <<endl;
    assert (data);
    assert (data -> qapp);
	data ->write_preferences ();
    data -> qapp -> quit ();
 }








void DDWin::raytraced_screenshot_slot () {
    string filename;
    filename = "povray.pov";
    write_POV_source (filename);
    string povcommand = "povray";
    float ratio = 1.f;
    int width =  (int) (gl -> width () * ratio);
    int height = (int) (gl -> height () * ratio);
    stringstream command;
    command << povcommand << "  -UV +A -I" << filename << " -W" << width <<" -H" << height <<"&";
    system (command.str ().c_str ());
}


void DDWin::emit_targets_updated () {targets_updated ();}

string DDWin::write_mol2_string (ZNMolecule *mol) {
    return "not implemented";

/*
    stringstream out;
    out << "@<TRIPOS>MOLECULE"<<endl;
    out<<mol->name<<endl;
    out<<mol->atoms.size ()<<" "<<mol->bonds.size ()<<" 0 0"<<endl;
    out<<endl<<endl;
    out<<"@<TRIPOS>ATOM"<<endl;
    for (unsigned int i=0; i<mol->atoms.size (); i++) {
        Atom *at = mol->atoms[i];
        out << at->ID<<" "<<at->mol2Type<<" "<<at-> GetVector ().x()<<" "<<at-> GetVector ().y()<<" "<<at-> GetVector ().z()<<" "<<at->atomType<<" "<<at->residue->number<<" "<<at->residue->name<<" "<<at->charge<<endl;
    }
    out<<"@<TRIPOS>BOND"<<endl;
    for (unsigned int i=0; i<mol->bonds.size (); i++) {
        string type;
        ZNBond *bo = mol->bonds[i];
        if (bo->mol2Type ==1) type = "1";
        else if (bo->mol2Type ==2) type = "2";
        else if (bo->mol2Type ==3) type = "3";
        else if (bo->mol2Type ==4) type = "am";
        else type = "ar";

        out<<bo->number<<" "<<bo->GetBeginAtom ()->ID<<" "<<bo->GetEndAtom ()->ID<<" "<<type<<endl;
    }
    return out.str ();
*/
} 



void DDWin::open_file_slot () {
	QString mol_name = QFileDialog::getOpenFileName(this, tr ("Open file"), last_visited_dir, tr("All Files (*)"));

	if (!mol_name.isEmpty ()) {
		last_visited_dir = QDir (mol_name).path ();
		load_file (mol_name.toStdString());
	}
}

void DDWin::open_session_file_slot () {
	QString mol_name = QFileDialog::getOpenFileName(this, tr ("Open file"), last_visited_dir, tr("Zodiac Session Files (*.zod)"));
		if (!mol_name.isEmpty ())  data ->actions ->load_session (mol_name.toStdString ());
/*
	if (!mol_name.isEmpty ()) {
		last_visited_dir = QDir (mol_name).path ();
		load_file (mol_name.toStdString());
	}
*/
}

void DDWin::set_lists (ZNMolecule *mol) {
    int bl1, bl2, ll, sl;
	int zero = gl ->new_list (); 

    bl1 = gl->new_list ();
	bl2 = gl->new_list ();
    ll = gl->new_list ();
    sl = gl->new_list ();

   // if (mol->multi) {
    //    Database_molecule *dm;
   //     dm = (Database_molecule *) mol;
    //    for (unsigned int i=0; i<dm->database->molecules.size (); i++) {
	//		set_display_lists (dm->database->molecules[i], ll, bl1, bl2, sl);        }
   // }
   // else {
		set_display_lists (mol, ll, bl1, bl2, sl);
    //        }

}


void DDWin::show_grid (Database *dat) {
	for (unsigned int i = 0; i < database_grids.size (); i++) {
		if (database_grids [i] ->database == dat) {
			database_grids[i] -> show ();
			database_grids[i] -> raise ();
		}
	}
}

void DDWin::select (vector <Atom *> atoms) {
	deselect ();
	Selection *sel = new Selection ();
	for (unsigned int i = 0; i < atoms.size (); i++) {
		sel ->select_atom (atoms[i]);
	}
	add_molecule (sel);
}

void DDWin::select (Atom *at) {
    deselect ();
    Selection *sel = new Selection ();
	sel -> select_atom (at);
	add_molecule (sel);
}

void DDWin::load_file (string mol_name) {
	ifstream ifs(mol_name.c_str ());
	OBConversion conv(&ifs);
	OBFormat* inFormat = conv.FormatFromExt(mol_name.c_str ());
	ZNMolecule *mol = new ZNMolecule ();
	if(conv.SetInFormat(inFormat) && conv.Read(mol)) { 
		ZNMolecule *mol2 = new ZNMolecule ();
		if (conv.Read(mol2)) {
			add_database (load_multi_file (mol_name));
		}
		else {
			finalise_molecule (mol);
			CreateZNMoleculeCommand *command = new CreateZNMoleculeCommand (mol, this);
			execute (command); 
			gl -> set_center_of_rotation (get_center (mol));
			gl -> set_center_of_view (get_center (mol));



		}
	}
	else {
cerr << "could not read file "<<mol_name<<endl;
	delete mol;
	}
}


Database *DDWin::load_multi_file (string mol_name) {
    Database *database = new Database ();
    ifstream ifs(mol_name.c_str ());
    OBConversion conv(&ifs);
    OBFormat* inFormat = conv.FormatFromExt(mol_name.c_str ());
    conv.SetInFormat(inFormat) ;
    bool go_on = true;
	int nn = 0;
    while (go_on) {
        Database_molecule *mol = new Database_molecule ();
		go_on = conv.Read (mol);
		if (mol -> NumAtoms ()) {
			nn++;
			database -> add_mol (mol);

			stringstream ss;
			ss  <<nn;
			string name = string ("ZNMolecule ")+ss.str ();
			mol -> SetTitle (name);
			finalise_molecule (mol);
			set_lists (mol);
			if (mol ->NumAtoms () < 100) {
				set_bonds_display_style (mol, 2);
			}
			if (!mol ->NumBonds ()) {
				set_atoms_display_style (mol, 1);
			}
			
			
		}
		else delete mol;
	}
	database -> import_csv (mol_name+".csv");
	return database;

}


void DDWin::deselect () {
	ZNMolecule *selected_mol = 0;
    if (target_molecule->selection) {

		Selection *target_sel = (Selection *) target_molecule;
		set_current_target (0);

		if ( target_sel ->get_molecules ().size () ){

			selected_mol = target_sel ->get_molecules ()[0];
		}
	}

    for (unsigned int i=0; i < molecules.size (); i++) {
        int s = molecules.size ();
        assert ((s-1-(int) i) < molecules.size ());
        assert ((s-1-(int) i) >= 0);
        if (molecules[s-1-i] -> selection) {
            Selection *sel = (Selection *) molecules[s-1-i];
            sel -> deselect ();
            target->removeItem (s-1-i);
            molecules.erase (molecules.begin ()+s-1-i);
     //       delete sel; cannot use this because OBMol distructor would delete all atoms in sel


        }
    }
	if (selected_mol) {
		set_current_target (selected_mol);
	}
	for (unsigned int i=0; i < molecules.size (); i++) {
		gl ->draw_molecule (molecules[i]);
	}
}




Sphere* DDWin::new_sphere (string name) {
    for (unsigned int i=0; i<graphical_objects.size (); i++) {
        if (graphical_objects[i]->name==name) return (Sphere *) graphical_objects[i];
    }
        Sphere *sphere= new Sphere;
        sphere ->set_name (name);
        graphical_objects.push_back (sphere);
        return sphere;


}



void DDWin::add_molecule (ZNMolecule *mol) {
	if (string (mol ->GetTitle ()) == "") {
	    stringstream ssname;
		if (mol ->multi) 		ssname << "Database " << builder ->mcounter++;
		else ssname << "New Molecule " << builder ->mcounter++;
		string s = ssname.str ();
		mol -> SetTitle (s);
	}
    set_lists (mol);
	if (mol ->NumAtoms () < 100) {
		set_bonds_display_style (mol, 2);
	}
	if (!mol ->NumBonds ()) {
		set_atoms_display_style (mol, 1);
	}
    molecules.push_back (mol);
    target->insertItem (1000, QString( mol->GetTitle ()));
	if (!mol ->selection) non_selection_molecules_updated ();
	else targets_updated ();
}


void DDWin::add_database (Database *database) {
	DatabaseGrid *grid = new DatabaseGrid (0, database, data);
	database_grids.push_back (grid);
	ZNMolecule *target = database ->dummy_mol;
	if (database ->get_molecule (0)) {
		target = database ->get_molecule (0);
	}
		CreateZNMoleculeCommand *command = new CreateZNMoleculeCommand (target, this);
		execute (command);
//		gl -> set_center_of_rotation (get_center (target));
//		gl -> set_center_of_view (get_center (target));
	grid ->set_mol();
	

//	grid -> show ();
//	grid -> raise ();

}

void DDWin::new_molecule (string name) {
    ZNMolecule *mol = new ZNMolecule ();
    mol-> SetTitle (name);
    add_molecule (mol);

}



void DDWin::resizeEvent (QResizeEvent *){
//    target->move (width()-400,0);

}

void DDWin::create_menu_actions () {

	openAct = new QAction (tr("&Open File"), this);
	openAct->setShortcut (tr ("Ctrl+O"));
	connect (openAct, SIGNAL (triggered ()), this, SLOT (open_file_slot ()));

	openSessionAct = new QAction (tr("&Open Session File"), this);
	openSessionAct->setShortcut (tr ("Ctrl+E"));
	connect (openSessionAct, SIGNAL (triggered ()), this, SLOT (open_session_file_slot ()));

	saveasAct = new QAction (tr("Save as..."), this);
	saveasAct->setShortcut (tr ("Ctrl+S"));
	connect (saveasAct, SIGNAL (triggered ()), this, SLOT (save_as_slot ()));

	saveSessionAsAct = new QAction (tr("Save Session as..."), this);
	saveSessionAsAct->setShortcut (tr ("Ctrl+A"));
	connect (saveSessionAsAct, SIGNAL (triggered ()), this, SLOT (save_session_as_slot ()));

	screenshotAct = new QAction (tr("&Screenshot"), this);
	screenshotAct->setShortcut (tr ("Ctrl+N"));
	connect (screenshotAct, SIGNAL (triggered ()), this, SLOT (screenshot_slot ()));

    raytracedscreenshotAct = new QAction (tr("&Raytraced Screenshot POVRAY"), this);
    raytracedscreenshotAct->setShortcut (tr ("Ctrl+R"));
    connect (raytracedscreenshotAct, SIGNAL (triggered ()), this, SLOT (raytraced_screenshot_slot ()));

    movieAct = new QAction (tr("Start / save &Movie"), this);
    movieAct->setShortcut (tr ("Ctrl+M"));
    connect (movieAct, SIGNAL (triggered ()), this, SLOT (movie_slot ()));

    raytracedmovieAct = new QAction (tr("Start / save  Raytraced Movie"), this);
  //  screenshotAct->setShortcut (tr ("Ctrl+O"));
    connect (raytracedmovieAct, SIGNAL (triggered ()), this, SLOT (raytraced_movie_slot ()));

    quitAct = new QAction (tr("&Quit"), this);
    quitAct ->setShortcut (tr ("Ctrl+Q"));
    connect (quitAct, SIGNAL (triggered ()), this, SLOT (close ()));
	
	newdatabaseAct = new QAction (tr ("Create new Database"), this);
	connect (newdatabaseAct, SIGNAL (triggered ()), this, SLOT (new_database_slot ()));	

	builderAct = new QAction (tr("&Builder"), this);
	builderAct->setShortcut (tr ("Ctrl+B"));
	connect (builderAct, SIGNAL (triggered ()), this, SLOT (builder_slot ()));
	
	twodwinAct = new QAction (tr("&2D Editor - molsKetch"), this);
	twodwinAct->setShortcut (tr ("Ctrl+K"));
	connect (twodwinAct, SIGNAL (triggered ()), this, SLOT (twodwin_slot ()));

	historyAct = new QAction (tr("History"), this);
	historyAct->setShortcut (tr ("Ctrl+H"));
	connect (historyAct, SIGNAL (triggered ()), this, SLOT (history_slot ()));

	settingsAct = new QAction (tr("&Preferences"), this);
	settingsAct ->setShortcut (tr ("Ctrl+P"));
	connect (settingsAct, SIGNAL (triggered ()), this, SLOT (pref_slot ()));

   wiimoteAct = new QAction (tr("Connect Wiimote for head tracking"), this);
   connect (wiimoteAct, SIGNAL (triggered ()), this, SLOT (wiimote_slot ()));

   wiimote2Act = new QAction (tr("Connect Wiimote for 3D input"), this);
   connect (wiimote2Act, SIGNAL (triggered ()), this, SLOT (wiimote2_slot ()));



    hideHAct = new QAction (tr("Hide Hydrogens"), this);
    connect (hideHAct, SIGNAL (triggered ()), this, SLOT (hide_hydrogens_slot ()));

    hidenpHAct = new QAction (tr("Hide non polar Hydrogens"), this);
    connect (hideHAct, SIGNAL (triggered ()), this, SLOT (hide_nonpolar_hydrogens_slot ()));

    hideallAct = new QAction (tr("Hide All Atoms"), this);
    connect (hideallAct, SIGNAL (triggered ()), this, SLOT (hide_all_atoms_slot ()));

    showallAct = new QAction (tr("Show All Atoms"), this);
    connect (showallAct, SIGNAL (triggered ()), this, SLOT (show_all_atoms_slot ()));

    displaysettingsAct = new QAction (tr("&Display Settings"), this);
    displaysettingsAct->setShortcut (tr ("Ctrl+D"));
    connect (displaysettingsAct, SIGNAL (triggered ()), this, SLOT (display_settings_slot ()));

    sequenceAct = new QAction (tr("&Sequence"), this);
    connect (sequenceAct, SIGNAL (triggered ()), this, SLOT (sequence_slot ()));

    DDsettingsAct = new QAction (tr("3D Settings"), this);
    connect (DDsettingsAct, SIGNAL (triggered ()), this, SLOT (DD_settings_slot ()));
	
	
	enableclippingAct = new QAction (tr("enable clipping"), this);
    connect (enableclippingAct, SIGNAL (triggered ()), this, SLOT (enable_clipping_slot ()));

	disableclippingAct = new QAction (tr("disable clipping"), this);
    connect (disableclippingAct, SIGNAL (triggered ()), this, SLOT (disable_clipping_slot ()));	
	
	DDsettingsAct = new QAction (tr("3D Settings"), this);
    connect (DDsettingsAct, SIGNAL (triggered ()), this, SLOT (DD_settings_slot ()));

    colorAct = new QAction (tr("Color"), this);
    connect (colorAct, SIGNAL (triggered ()), this, SLOT (color_slot ()));
	
	backboneColorAct = new QAction (tr("Backbone Color"), this);
    connect (backboneColorAct, SIGNAL (triggered ()), this, SLOT (backbone_color_slot ()));
	
	
	
	backgroundcolorAct = new QAction (tr("Background Color"), this);
	connect (backgroundcolorAct, SIGNAL (triggered ()), this, SLOT (background_color_slot ()));
	
	
    hapticAct = new QAction (tr("Haptic mode"), this);
	hapticAct -> setIcon (QPixmap (":icons/haptic.png"));
    connect (hapticAct, SIGNAL (triggered ()), this, SLOT (haptic_slot ()));

    dockingAct = new QAction (tr("Docking"), this);
    connect (dockingAct, SIGNAL (triggered ()), this, SLOT (docking_slot ()));
	
    plantsAct = new QAction (tr("Plants"), this);
    connect (plantsAct, SIGNAL (triggered ()), this, SLOT (plants_slot ()));

    gamessAct = new QAction (tr("Create Gamess input"), this);
    connect (gamessAct, SIGNAL (triggered ()), this, SLOT (gamess_slot ()));

    energyAct = new QAction (tr("Compute Energies"), this);
    connect (energyAct, SIGNAL (triggered ()), this, SLOT (compute_energy_slot ()));

    logPAct = new QAction (tr("Predict LogP"), this);
    connect (logPAct, SIGNAL (triggered ()), this, SLOT (logP_slot ()));


    minimiseAct = new QAction (tr("Energy Minimise"), this);
    minimiseAct -> setIcon (QPixmap (":icons/minimise.png"));
    connect (minimiseAct, SIGNAL (triggered ()), this, SLOT (minimise_energy_slot ()));

    partialQAct = new QAction (tr("Gasteiger Partial Charges"), this);
    connect (partialQAct, SIGNAL (triggered ()), this, SLOT (partial_charges_slot ()));

    scoresCharge = new QAction (tr("Scores from Charges"), this);
    connect (scoresCharge, SIGNAL (triggered ()), this, SLOT (scores_from_charges_slot ()));
	
	conformationalSearchAct = new QAction (tr("Conformational Search"), this);
    connect (conformationalSearchAct, SIGNAL (triggered ()), this, SLOT (conformers_slot ()));
	

    addHsAct = new QAction (tr("Add Hydrogens"), this);
    connect (addHsAct, SIGNAL (triggered ()), this, SLOT (add_Hs_slot ()));
	
	duplicateAct = new QAction (tr("Duplicate mol"), this);
    connect (duplicateAct, SIGNAL (triggered ()), this, SLOT (duplicate_slot ()));

    surfaceAct = new QAction (tr("Molecular Surface"), this);
    connect (surfaceAct, SIGNAL (triggered ()), this, SLOT (surface_slot ()));
	
	mapAct = new QAction (tr("Map"), this);
    connect (mapAct, SIGNAL (triggered ()), this, SLOT (map_slot ()));

    sphereAct = new QAction (tr("Sphere"), this);
    connect (sphereAct, SIGNAL (triggered ()), this, SLOT (sphere_slot ()));
	
	backbonetosurfAct = new QAction (tr("Backbone to surface"), this);
    connect (backbonetosurfAct, SIGNAL (triggered ()), this, SLOT (backbone_to_surface_slot ()));

    graphicalobjectsAct = new QAction (tr("Graphical Objects"), this);
    connect (graphicalobjectsAct, SIGNAL (triggered ()), this, SLOT (graphical_objects_slot ()));
	
	
	iodeviceAct = new QAction (tr("I/O device list"), this);
    connect (iodeviceAct, SIGNAL (triggered ()), this, SLOT (iodevice_slot ()));

    aboutAct = new QAction (tr("About"), this);
    connect (aboutAct, SIGNAL (triggered ()), this, SLOT (about_slot ()));



	selectAllAct = new QAction (tr("Select All Atoms"), this);
	selectAllAct -> setIcon (QPixmap (":icons/select_all.png"));
    connect (selectAllAct, SIGNAL (triggered ()), this, SLOT (select_all_slot ()));

	deselectAct = new QAction (tr("Deselect"), this);
	deselectAct -> setIcon (QPixmap (":icons/deselect.png"));
    connect (deselectAct, SIGNAL (triggered ()), this, SLOT (deselect_slot ()));

	
	invertSelectAct = new QAction (tr("Invert Selection"), this);
	invertSelectAct -> setIcon (QPixmap (":icons/invert_selection.png"));
    connect (invertSelectAct, SIGNAL (triggered ()), this, SLOT (invert_selection_slot ()));


	selectCAct = new QAction (tr("Select Carbons"), this);
	selectCAct -> setIcon (QPixmap (":icons/select_C.png"));
    connect (selectCAct, SIGNAL (triggered ()), this, SLOT (select_C_slot ()));
	
	selectHAct = new QAction (tr("Select Hydrogens"), this);
	selectHAct -> setIcon (QPixmap (":icons/select_H.png"));
    connect (selectHAct, SIGNAL (triggered ()), this, SLOT (select_H_slot ()));


	selectsquareAct = new QAction (tr("Select Square"), this);
	selectsquareAct -> setIcon (QPixmap (":icons/select_square.png"));
    connect (selectsquareAct, SIGNAL (triggered ()), this, SLOT (select_square_slot ()));



	buildCAct = new QAction (tr("Carbon"), this);
	buildCAct -> setIcon (QPixmap (":icons/builder_C.png"));
    connect (buildCAct, SIGNAL (triggered ()), this, SLOT (build_C_slot ()));

	buildOAct = new QAction (tr("Oxygen"), this);
	buildOAct -> setIcon (QPixmap (":icons/builder_O.png"));
    connect (buildOAct, SIGNAL (triggered ()), this, SLOT (build_O_slot ()));

	buildNAct = new QAction (tr("Nitrogen"), this);
	buildNAct -> setIcon (QPixmap (":icons/builder_N.png"));
    connect (buildNAct, SIGNAL (triggered ()), this, SLOT (build_N_slot ()));

	buildSAct = new QAction (tr("Sulphur"), this);
	buildSAct -> setIcon (QPixmap (":icons/builder_S.png"));
    connect (buildSAct, SIGNAL (triggered ()), this, SLOT (build_S_slot ()));


	buildsinglebondAct = new QAction (tr("Single bond"), this);
	buildsinglebondAct -> setIcon (QPixmap (":icons/builder_1bond.png"));
    connect (buildsinglebondAct, SIGNAL (triggered ()), this, SLOT (build_single_bond_slot ()));
	
	
	
	builddoublebondAct = new QAction (tr("Double bond"), this);
	builddoublebondAct -> setIcon (QPixmap (":icons/builder_2bond.png"));
    connect (builddoublebondAct, SIGNAL (triggered ()), this, SLOT (build_double_bond_slot ()));
	
	buildtriplebondAct = new QAction (tr("Triple bond"), this);
	buildtriplebondAct -> setIcon (QPixmap (":icons/builder_3bond.png"));
    connect (buildtriplebondAct, SIGNAL (triggered ()), this, SLOT (build_triple_bond_slot ()));
	
	deletebondAct = new QAction (tr("Delete bond"), this);
	deletebondAct -> setIcon (QPixmap (":icons/builder_Xbond.png"));
    connect (deletebondAct, SIGNAL (triggered ()), this, SLOT (delete_bond_slot ()));
	

	deleteatomAct = new QAction (tr("Delete atom"), this);
	deleteatomAct -> setIcon (QPixmap (":icons/builder_del.png"));
//	connect (deleteatomAct, SIGNAL (triggered ()), this, SLOT ( ));

	magicpencilButt = new QPushButton (QIcon (":icons/builder_pencil.png"), "");

	magicpencilButt ->setIconSize (QSize (32, 32));
	magicpencilButt ->resize (32, 32);	
	magicpencilButt ->setCheckable (true);
//	magicpencilAct = new QAction (tr("Magic pencil"), this);
//	magicpencilAct -> setIcon (QPixmap (":icons/builder_pencil.png"));
	connect (magicpencilButt, SIGNAL (toggled (bool)), this, SLOT (magic_pencil_toggle_slot (bool) ));

	buildsmileAct = new QAction (tr("SMILE"), this);
	buildsmileAct -> setIcon (QPixmap (":icons/builder_smile.png"));
//	connect (buildsmileAct, SIGNAL (triggered ()), this, SLOT ( ));


	periodictableAct = new QAction (tr("Periodic Table"), this);
	periodictableAct -> setIcon (QPixmap (":icons/periodic_table.png"));
	connect (periodictableAct, SIGNAL (triggered ()), this, SLOT (periodic_table_slot () ));



    undoAct = data -> undo_stack -> createUndoAction ( this );
    undoAct -> setIcon (QPixmap (":icons/undo.png"));
    undoAct->setShortcut (tr ("Ctrl+Z"));


    redoAct = data -> undo_stack -> createRedoAction ( this );
    redoAct -> setIcon (QPixmap (":icons/redo.png"));
    redoAct->setShortcut (tr ("Shift+Ctrl+Z"));



    colorsAct = new QAction (tr("Color molecule"), this);;
	QPixmap pix(24, 24);
    pix.fill(molecule_color);
    colorsAct->setIcon(pix);
	connect (colorsAct, SIGNAL (triggered ()), this, SLOT (quick_color_slot ()));
	
	colors_chooserAct = new QAction (tr("change color"), this);;
	colors_chooserAct -> setIcon (QPixmap (":icons/arrow.png"));
    connect (colors_chooserAct, SIGNAL (triggered ()), this, SLOT (change_color_slot ()));

	centerAct = new QAction (tr("View"), this);;
	connect (centerAct, SIGNAL (triggered ()), this, SLOT (center_slot ()));


}


void DDWin::periodic_table_slot () {

cerr << "periodic table" <<endl;


}



void DDWin::draw_menu (){


    create_menu_actions ();

    QMenu *file = new QMenu(tr("&File"), this );
    Q_CHECK_PTR( file );
    file -> addAction (openAct);
    file -> addAction (openSessionAct);
    file -> addAction (saveasAct);
    file -> addAction (saveSessionAsAct);
    file -> addAction (screenshotAct);
    file -> addAction (raytracedscreenshotAct);
  //  file -> addAction (movieAct);
  //  file -> addAction (raytracedmovieAct);
    file -> addAction (quitAct);


    QMenu *edit = new QMenu(tr("&Edit"), this );
    Q_CHECK_PTR( edit );
    edit -> addAction (builderAct);
	edit -> addAction (twodwinAct);
	edit -> addAction (newdatabaseAct);
    edit -> addAction (historyAct);

  //  edit -> addAction (wiimoteAct);
 //   edit -> addAction (wiimote2Act);
//	edit -> addAction (iodeviceAct);
	edit -> addAction (settingsAct);

    QMenu *show = new QMenu(tr("&Show"), this );
    Q_CHECK_PTR( show );
    show -> addAction (hideHAct);
    show -> addAction (hidenpHAct);
    show -> addAction (hideallAct);
    show -> addAction (showallAct);


    QMenu *display = new QMenu(tr("&Display"), this );
    Q_CHECK_PTR( display );
    display -> addAction (displaysettingsAct);
    display -> addAction (sequenceAct);
    display -> addAction (DDsettingsAct);
    display -> addAction (enableclippingAct);
	display -> addAction (disableclippingAct);	
    display -> addAction (colorAct);
	display -> addAction (backboneColorAct);
	display -> addAction (backgroundcolorAct);



    QMenu *compute = new QMenu(tr("&Compute"), this );
    Q_CHECK_PTR( compute );
    compute -> addAction (hapticAct);
	compute -> addAction (dockingAct);
    compute -> addAction (plantsAct);
    compute -> addAction (gamessAct);
    compute -> addAction (energyAct);
    compute -> addAction (logPAct);
    compute -> addAction (minimiseAct);
    compute -> addAction (partialQAct);
    compute -> addAction (scoresCharge);
	compute -> addAction (conformationalSearchAct);

    QMenu *utilities = new QMenu(tr("&Utilities"), this );
    Q_CHECK_PTR( utilities );
    utilities -> addAction ( addHsAct );
	utilities -> addAction (duplicateAct);


    QMenu *graph_objects = new QMenu(tr("&Graphical Objects"), this );
    Q_CHECK_PTR( graph_objects );
    graph_objects -> addAction (surfaceAct);
	graph_objects -> addAction (mapAct);
    graph_objects -> addAction (sphereAct);
	graph_objects -> addAction (backbonetosurfAct);
    graph_objects -> addAction (graphicalobjectsAct);


    QMenu *about = new QMenu(tr("&About"), this );
    Q_CHECK_PTR( about );
    about -> addAction (aboutAct);


    QMenuBar *menu = menuBar ();
    menu->addMenu ( file );
    menu->addMenu ( edit );
    menu->addMenu ( show );
    menu->addMenu ( display );
    menu->addMenu ( compute );
    menu->addMenu ( graph_objects );
    menu->addMenu ( utilities );
    menu->addMenu ( about );

    target = new QComboBox( menu );
    target->setMinimumWidth (400);
    target->insertItem(0, "All" );
	targets_updated ();

	current_mode = new QLabel ("");

    toolbar = new QToolBar ("Main");
	target_toolbar = new QToolBar ("Control");
	select_tb = new QToolBar ("Select");
	builder_tb = new QToolBar ("Builder");
	connect (target_toolbar, SIGNAL (orientationChanged (Qt::Orientation)), this, SLOT (resize_target_toolbar (Qt::Orientation)));
    addToolBar (toolbar);
	addToolBar (target_toolbar);
	addToolBar (Qt::RightToolBarArea, select_tb);
	addToolBar (Qt::RightToolBarArea, builder_tb);

	//MyColorButton *colors = new MyColorButton (toolbar ->layout (), molecule_color);
	
	hide_toolbars ();
    toolbar -> addAction (undoAct);
    toolbar -> addAction (redoAct);
    toolbar -> addAction (minimiseAct);
	toolbar -> addAction (hapticAct);
	
	toolbar ->addSeparator ();

	toolbar -> addAction (colorsAct);
	toolbar -> addAction (centerAct);
//	toolbar -> addAction (colors_chooserAct);


   	target_toolbar -> addWidget (target);
	target_toolbar -> addWidget (current_mode);
	
	select_tb -> addAction (selectAllAct);
	select_tb -> addAction (deselectAct);
	select_tb -> addAction (invertSelectAct);
	select_tb -> addAction (selectCAct);
	select_tb -> addAction (selectHAct);
	select_tb -> addAction (selectsquareAct);
 
 
 
	builder_tb -> addAction (buildCAct);
	builder_tb -> addAction (buildOAct);
	builder_tb -> addAction (buildNAct);
	builder_tb -> addAction (buildSAct);
	builder_tb -> addAction (buildsinglebondAct);
	builder_tb -> addAction (builddoublebondAct);
	builder_tb -> addAction (buildtriplebondAct);
	builder_tb -> addAction (deletebondAct);	
	builder_tb -> addAction (deleteatomAct);
	//builder_tb -> addAction (magicpencilAct);
	builder_tb -> addWidget (magicpencilButt);
	builder_tb -> addAction (buildsmileAct);
//	builder_tb -> addAction (periodictableAct);
	
    connect (target, SIGNAL (activated (int)), this, SLOT (set_current_target_slot (int)));

}

void DDWin::quick_color_slot () {
	change_color_slot ();
	if (current_target) {
		ColorAtomCommand *color_atom = new ColorAtomCommand (gl);
		FOR_ATOMS_OF_MOL (at, target_molecule) {
			if (target_molecule ->selection || at ->GetAtomicNum ()==6) {
				color_atom -> add (&*at, molecule_color);
			}
		}
		color_atom->set_name ();
		execute (color_atom);
    }


}


void DDWin::center_slot () {
        if (current_target) {
			gl -> set_center_of_rotation (get_center (target_molecule));
			gl -> set_center_of_view (get_center (target_molecule));
		}
}

void DDWin::change_color_slot () {
	QColor new_color = QColorDialog::getColor(molecule_color, this );
	if ( new_color.isValid () ) {
		molecule_color.setRed (new_color.red ());
		molecule_color.setGreen (new_color.green ());
		molecule_color.setBlue (new_color.blue ());
		molecule_color.setAlpha (new_color.alpha ());
		
	}
	QPixmap pix(24, 24);
    pix.fill(molecule_color);
    colorsAct->setIcon(pix);
}



void DDWin::resize_target_toolbar (Qt::Orientation orient) {
	if (orient == Qt::Horizontal) {
			target -> setMinimumWidth (400);
	}
	else {
	target -> setMinimumWidth (40);
		target -> setMaximumWidth (40);
	} 
}

void DDWin::hide_toolbars () {
//	toolbar -> hide ();
	select_tb -> hide ();
	builder_tb -> hide ();
}


void DDWin::delete_current_molecule () {
	delete_molecule (molecules [current_target]);
	sequence_menu ->del_from_tab (molecules [current_target]);
}

void DDWin::delete_molecule (ZNMolecule *mol) {
    for (unsigned int i=0; i<molecules.size (); i++) {
        if (molecules[i] == mol) {
            target->removeItem (i);
            molecules.erase (molecules.begin ()+i);
            set_current_target (0);
			if (!mol ->selection) non_selection_molecules_updated ();
			else targets_updated ();
            delete mol;
            break;
        }
    }
}


void DDWin::remove_molecule (ZNMolecule *mol) {
    for (unsigned int i=0; i<molecules.size (); i++) {
        if (molecules[i] == mol) {
            target->removeItem (i);
            molecules.erase (molecules.begin ()+i);
            set_current_target (0);

            break;
        }
    }

}



/*

void DDWin::redraw (ZNMolecule *mol) {
    if (mol -> needs_redraw) {
        gl -> draw_molecule (mol);
        mol -> needs_redraw = false;
    }
}

 
 */
void DDWin::recolor_by_score (ZNMolecule *mol) {
    if (mol -> needs_recolor) {
        vector <color_mask> masks;
        color_mask mask;
        mask.intensity = 1.0f;
        mask.only_to = 0;
        mask.excluding = 0;
        mask.type =  2; //score //see menu.cc
        masks.push_back (mask);
        gl -> apply_color_masks (masks, mol, false);
        mol -> needs_recolor = false;
    }
}


void DDWin::end_minimisation () {

    MoveAtomsCommand *command = new MoveAtomsCommand (gl, 0);
  //  FOR_ATOMS_OF_MOL (a, data -> minimize -> minimising_molecule) {
  //      command -> add (&*a, (vect &) a -> GetVector ());
  //  }

  //  execute (command);
    data -> undo_stack -> endMacro ();
    data -> minimize -> deinitialise_minimisation ();
}

void DDWin::delete_graphical_object (int i) {
    assert (i>-1 && i<graphical_objects.size ());
	go_updated ();
    DeleteGraphicalObjectCommand *command = new DeleteGraphicalObjectCommand (graphical_objects[i], this);
    execute (command);  
}






void DDWin::screenshot_slot () {
    QString s = QFileDialog::getSaveFileName(this, 
                    tr ("Save Screenshot"), "",tr("Portable pixel map (*.ppm)"));

    if (!s.isNull()) gl->screenshot (s);

}

void DDWin::movie_slot () {
    if (!rec_movie) {
        rec_movie = TRUE;
        gl->movie_time = QTime::currentTime ();
    }
    else {
        rec_movie = FALSE;        

    QString s = QFileDialog::getSaveFileName(this, 
                    tr ("Save Movie"), "",tr("MPEG (*.mpg)"));

        if (!s.isNull()) {
            ofstream *file = new ofstream("__conf__");
            *file<< "BASE_FILE_FORMAT PPM"<<endl;
             *file<< "INPUT"<<endl;
            *file<< "__movie__*.ppm [0-"<<last_frame-1<<"]"<<endl;
            *file<< "END_INPUT"<<endl;
            *file<< "INPUT_DIR ."<<endl;
            *file<< "IQSCALE 1"<<endl;
            *file<< "PQSCALE 1"<<endl;
            *file<< "BQSCALE 1"<<endl;
            *file<< "RANGE 1"<<endl;
            *file<< "PSEARCH_ALG SUBSAMPLE"<<endl;
            *file<< "GOP_SIZE 1"<<endl;
            *file<< "PATTERN I"<<endl;
            *file<< "PIXEL FULL"<<endl;
            *file<< "INPUT_CONVERT *"<<endl;
            *file<< "SLICES_PER_FRAME 1"<<endl;
            *file<< "BSEARCH_ALG SIMPLE"<<endl;
            *file<< "REFERENCE_FRAME DECODED"<<endl;
            *file<< "FRAME_RATE 25"<<endl;
            *file<< "OUTPUT " << s.toStdString() <<endl;
            system ("ppmtompeg __conf__*");
        }
        
        last_frame = 0;
        system ("rm __*");
    } 
}



void DDWin::raytraced_movie_slot () {
    if (!rec_pov_movie) {
        rec_pov_movie = TRUE;
        gl->movie_time = QTime::currentTime ();
    }
    else {
        rec_pov_movie = FALSE;     
   
        QString s = QFileDialog::getSaveFileName(this, 
                    tr ("Save Movie"), "",tr("MPEG (*.mpg)"));

        if (!s.isNull()) {
            ofstream *file = new ofstream("__conf__");
            *file<< "BASE_FILE_FORMAT PPM"<<endl;
             *file<< "INPUT"<<endl;
            *file<< "__movie__*.ppm [0-"<<last_frame-1<<"]"<<endl;
            *file<< "END_INPUT"<<endl;
            *file<< "INPUT_DIR ."<<endl;
            *file<< "IQSCALE 1"<<endl;
            *file<< "PQSCALE 1"<<endl;
            *file<< "BQSCALE 1"<<endl;
            *file<< "RANGE 1"<<endl;
            *file<< "PSEARCH_ALG SUBSAMPLE"<<endl;
            *file<< "GOP_SIZE 1"<<endl;
            *file<< "PATTERN I"<<endl;
            *file<< "PIXEL FULL"<<endl;
            *file<< "INPUT_CONVERT *"<<endl;
            *file<< "SLICES_PER_FRAME 1"<<endl;
            *file<< "BSEARCH_ALG SIMPLE"<<endl;
            *file<< "REFERENCE_FRAME DECODED"<<endl;
            *file<< "FRAME_RATE 25"<<endl;
            *file<< "OUTPUT " << s.toStdString() <<endl;
            for (unsigned int i=0; i<last_frame; i++) { 
            stringstream ss;
            ss<<i;
            string command = "povray __movie__"+ss.str()+".pov +FP -UV -H480 -W640";
            system (command.c_str());
            }
            system ("ppmtompeg __conf__ ");
        }
        
        last_frame = 0;
        system ("rm __*");
    } 
}


void DDWin::show_atom_properties (Atom* at){
    clicked_atom_menu->set (at);

    clicked_atom_menu->show ();
    clicked_atom_menu->raise ();
}


void DDWin::execute (Command * comm) {
    data -> undo_stack -> push (comm);
}

void DDWin::run_thread (Thread *thread) {
	running_threads.push_back (thread);
	thread ->start ();
}

/*
int DDWin::style_str_to_i (string style){


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
*/


void DDWin::save_as_slot () {
	/*
	QFileDialog save_file(this);

	QStringList filters;
	filters << "Tripos mol2 (*.mol2)"
		<< "InChI (*.inchi)"
		<< "MDL Molfile (*.mol)"
		<< "Protein Data Bank (*.pdb)"
		<< "sdf (*.sdf)"
		<< "smile (*.smi)";
	
		
	save_file.setNameFilters (filters);
	save_file.setAcceptMode (QFileDialog::AcceptSave);
	save_file.exec();

	QString ext = save_file.selectedNameFilter();

	QString filter = save_file.selectedNameFilter();
	filter.truncate (filter.lastIndexOf (')'));
	filter.remove (0, filter.indexOf ('*')+1);
	QString s = save_file.selectedFile();
	if (!s.contains ('.') ) s+= filter;
	 */
	
	QString s = QFileDialog::getSaveFileName(this, 
								 tr ("Save file"), last_visited_dir, tr("All Files (*)"));


	
    if (!s.isNull()) {
        if (current_target) {
			if (!is_db_extended (target_molecule)) data -> actions -> save_as (target_molecule, s.toStdString ());
			else { 
				Database_molecule *dm = (Database_molecule *) target_molecule;
				Database *dat = dm ->database;
				data ->actions ->save_as (dat, s.toStdString ());
			}
		}
		
		
	}
	
}


void DDWin::save_session_as_slot () {
/*
	QFileDialog save_file(this);

	QStringList filters;
	filters << "Zodiac session (*.zod)";
	
	save_file.setNameFilters (filters);
	save_file.setAcceptMode (QFileDialog::AcceptSave);
	save_file.exec();

	QString ext = save_file.selectedNameFilter();

	QString filter = save_file.selectedNameFilter();
	filter.truncate (filter.lastIndexOf (')'));
	filter.remove (0, filter.indexOf ('*')+1);
	QString s = save_file.selectedFile();
		if (!s.contains ('.') ) s+= filter;
	*/
	QString s = QFileDialog::getSaveFileName(this, 
								 tr ("Save file"), last_visited_dir, tr("Zodiac Session files (*.zod)"));

	    if (!s.isNull()) {
 
				data ->actions ->save_session_as (s.toStdString ());
			}


}


void DDWin::set_current_target_slot (int index) {
    current_target = index;
    target_molecule = molecules [index];
    if (target_molecule->multi) {
        Database_molecule *dm;
        dm = (Database_molecule *) target_molecule;
		show_grid (dm ->database);
    }


}

void DDWin::set_current_target (int index) {
    if (index<0) index+=molecules.size ();
    target->setCurrentIndex (index);
    set_current_target_slot (index);

}

void DDWin::set_current_target (ZNMolecule *mol) {
    int index = -1;
    for (unsigned int i=0; i<molecules.size (); i++) {

        if (molecules[i] == mol) {index = i; break;}
    }

    if (index>0) {
        set_current_target (index);
    }

}

void DDWin::scores_from_charges_slot () {
    if (current_target) {
				if (!is_db_extended (target_molecule)) 		data ->actions ->set_scores_from_charges (target_molecule);
			else { 
				Database_molecule *dm = (Database_molecule *) target_molecule;
				Database *dat = dm ->database;
				data ->actions ->set_scores_from_charges (dat);
			}

    }
	else {
        for (unsigned int i=0; i<molecules.size (); i++) {
				data ->actions ->set_scores_from_charges (molecules[i]);
        }
    }

}


void DDWin::new_database_slot () {
	Database *dat = new Database ();
	add_database (dat); 
}




void DDWin::partial_charges_slot () {
    if (current_target) {
        OBGastChrg gastchrg;
        gastchrg.AssignPartialCharges (*target_molecule);
        gl->draw_molecule (target_molecule);
    }
    
    else {
        for (unsigned int i=0; i<molecules.size (); i++) {
            data->mmff->compute_partial_charges (molecules[i]);
            gl->draw_molecule (molecules[i]);
        }
    }
}

void DDWin::atdebug_slot () {
 /*   ofstream *file = new ofstream("/home/nicola/Desktop/debug");
    show_labels = false;
    for (unsigned int i=1; i<346+1; i++) {

        stringstream ss;
        ss<<i;
        load_file ("/home/nicola/MMFF94/hypervalent/hyp"+ss.str()+".mol2"); 
    }
    for (unsigned int m=1; m<molecules.size (); m++) {
        int sum=0;
        for (unsigned int a=0; a<molecules[m]->atoms.size (); a++) {
            sum += molecules[m]->atoms[a]->MMFFtype;
        }
        *file<<sum<<endl;
    }
    file->close ();
    system ("/home/nicola/Desktop/compare.py /home/nicola/Desktop/debug /home/nicola/Desktop/true_atypes");

*/
}

void DDWin::color_slot () {
    color_menu->show ();
    color_menu->raise ();
}

void DDWin::backbone_color_slot () {
    b_color_menu->display ();
}

void DDWin::background_color_slot () {
    color_settings_menu->show ();
    color_settings_menu->raise ();
}


void DDWin::surface_slot () {
    surface_menu -> show ();
    surface_menu -> raise ();
}

void DDWin::map_slot () {
	map_menu  ->display ();
}

void DDWin::sphere_slot () {
    sphere_menu -> show ();
    sphere_menu -> raise ();
}

void DDWin::backbone_to_surface_slot () {
	if (current_target) {
		Surface *surf = new Surface;
		gl ->backbone_to_surface (target_molecule, surf);
		CreateGraphicalObjectCommand *command = new CreateGraphicalObjectCommand (surf, this);
		execute (command);
	}
}

void DDWin::graphical_objects_slot () {
    graphical_objects_menu -> show ();
    graphical_objects_menu -> raise ();
}


void DDWin::haptic_slot () {
    if (current_target) {
            haptic_menu->show ();
            haptic_menu->raise ();
    }
}

void DDWin::compute_energy_slot () {
/*		ZNMolecule *molecule = target_molecule;
		Chemscore *chemscore = new Chemscore;
		MMFF* mmff = new MMFF;
		build_kinematic_chain(molecule);
		Minimize *min = new Minimize (data, mmff, chemscore);
		min ->set_molecule (molecule);
	//	min ->internal_ff ->initialize_internal (molecule, ddwin ->molecules);
		min ->interaction_ff ->initialize_interaction (molecule, molecules, get_center (molecule), 12);


*/

	bool ext = false;
    if (current_target) {
		if (target_molecule -> multi) {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			if (dm -> database -> has_extend_enabled ()) ext = true;
		}
		if (!ext) {
			data->mmff->initialize (target_molecule, molecules);
			float bs, ab, to, sb, vw, op, el;
			stringstream out;
			ab = data->mmff->compute_angle_bendings ();
			to = data->mmff->compute_torsion_interactions ();
			sb = data->mmff->compute_stretch_bend_interactions ();
			vw = data->mmff->compute_van_der_waals_interactions ();
			op = data->mmff->compute_out_of_plane_bendings ();
			el = data->mmff->compute_electrostatic_interactions ();
			bs = data->mmff->compute_bond_stretchings ();
			
			out << "BS = "<<bs<<endl;
			out << "AB = "<<ab<<endl;
			out << "TO = "<<to<<endl;
			out << "SB = "<<sb<<endl;
			out << "OP = "<<op<<endl;
			out << "EL = "<<el<<endl;
			out << "VW = "<<vw<<endl;
			out << "total = "<<data -> actions -> compute_total_energy (target_molecule);
			QMessageBox::information( this, QString("Zodiac"), QString(out.str().c_str())); 
		}
		else {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			Database *dat = dm -> database;
			dat -> add_field ("Energy", data -> actions -> compute_total_energy (dat));
		}
	}
}

void DDWin::logP_slot () {
	bool ext = false;
    if (current_target) {
		if (target_molecule -> multi) {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			if (dm -> database -> has_extend_enabled ()) ext = true;
		}
		if (!ext) {
			stringstream out;
			out << "Predicted logP: ";
			out << data -> actions -> compute_logP (target_molecule);
			QMessageBox::information( this, QString("Zodiac"), QString(out.str().c_str())); 
		}
		else {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			Database *dat = dm -> database;
			dat -> add_field ("cLogP", data -> actions -> compute_logP (dat));
		}
	}
	
}


void DDWin::minimise_energy_slot () {
    if (current_target ) {
        if (!target_molecule -> selection) {
			if (!is_db_extended (target_molecule)) {
				data -> actions -> minimise (target_molecule);
			}
			else {
				Database_molecule *dm = (Database_molecule *) target_molecule;
				Database *dat = dm -> database;
				data -> actions -> minimise (dat);
			}
		}
	}
}	
	

void DDWin::iodevice_slot () {
	iodevice_menu -> display ();
	

}

				
void DDWin::plants_slot () {
	//Plants *plants= new Plants ( 0, "plants conf", data );
	plants_menu -> display ();
}

void DDWin::conformers_slot () {
	conformers_menu -> display ();
}

void DDWin::docking_slot () {
		docking_menu -> display ();
}

void DDWin::gamess_slot () {
	gamess_menu -> display ();
}
				
void DDWin::hide_hydrogens_slot () {
	if (current_target) {
		if (!is_db_extended (target_molecule)) {
			data -> actions -> hide_hydrogens (target_molecule);
		}
		else {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			Database *dat = dm -> database;
			data -> actions -> hide_hydrogens (dat);
		}
	}
				
	else gl -> hide_hydrogens (molecules);
}
				
								
															
void DDWin::hide_nonpolar_hydrogens_slot () {
	if (current_target) {
		if (!is_db_extended (target_molecule)) {
			data -> actions -> hide_nonpolar_hydrogens (target_molecule);
		}
		else {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			Database *dat = dm -> database;
			data -> actions -> hide_nonpolar_hydrogens (dat);
		}
	}
				
	else gl -> hide_hydrogens (molecules);
}
				
					
void DDWin::hide_all_atoms_slot () {
	if (current_target) {
		if (!is_db_extended (target_molecule)) {
			data ->actions -> hide_all_atoms (target_molecule);
		}
		else {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			Database *dat = dm -> database;
			data ->actions -> hide_all_atoms (dat);
		}
	}
				
	else gl -> hide_all_atoms (molecules);
}

void DDWin::show_all_atoms_slot () {
	if (current_target) {
		if (!is_db_extended (target_molecule)) {
			data ->actions -> show_all_atoms (target_molecule);
		}
		else {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			Database *dat = dm -> database;
			data ->actions -> show_all_atoms (dat);
		}
	}
				
	else gl -> show_all_atoms (molecules);
}


void DDWin::add_Hs_slot () {
	if (current_target) {
		bool ok;
		double ph = QInputDialog::getDouble(this, tr("Protonate Molecule"),
                                        tr("pH:"), 7.4, 0, 14, 1, &ok);
		if (ok)
		{
		if (!is_db_extended (target_molecule)) {
			data ->actions ->reprotonate (target_molecule, ph);
		}
		else {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			Database *dat = dm -> database;
			data ->actions ->reprotonate (dat, ph);
			
		}
		}
	}
}


void DDWin::duplicate_slot () {
	if (current_target) {
		ZNMolecule *new_mol = new ZNMolecule (*target_molecule);
		add_molecule (new_mol);
	}
}

/*
void DDWin::color_by_atom_slot () {
    if (current_target) gl->color_by_atom (target_molecule);
    else gl->color_by_atom (molecules);
}



void DDWin::color_by_charge_slot (){
    if (current_target) gl->color_by_charge (target_molecule);
    else gl->color_by_charge (molecules);
}

void DDWin::color_by_score_slot (){
    if (current_target) gl->color_by_score (target_molecule);
    else {
        for (unsigned int i =0; i<molecules.size (); i++) gl->color_by_score (molecules[i]);
    }
}


void DDWin::color_carbons_slot (){
	QColor c = QColorDialog::getColor(Qt::white, this );
    if ( c.isValid() ) {
        float col[3];
        col[0] = (float)(c.red())/255;
        col[1] = (float)(c.green())/255;
        col[2] = (float)(c.blue())/255;
        if (current_target) gl->color_atoms (target_molecule, 6, col);
        else gl->color_atoms (molecules, 6, col);
    }
}
*/
/////////////////////////////////////////////DRAG AND DROP////////////////////////////////////////////

void DDWin::dragEnterEvent( QDragEnterEvent * event) {
  //  cout << "drag enter event"<<endl ;

  //   if (event->mimeData()->hasText())
         event->acceptProposedAction();

}

void DDWin::dragMoveEvent(QDragMoveEvent *event)
{
 //   cout << "drag move event"<<endl ;
    event->acceptProposedAction();
}

void DDWin::dragLeaveEvent(QDragLeaveEvent *event)
{
 //   cout << "drag leave event"<<endl ;
    event->accept();
}


void DDWin::dropEvent(QDropEvent *event)
{
// cout << "drop event" << endl;   
    if (event->mimeData()->hasUrls()) {
        QList<QUrl> urlList = event->mimeData ()->urls();
        for (int i = 0; i < urlList.size()/* && i < 32*/; ++i) {
            QString url = urlList.at(i).toLocalFile ();
			last_visited_dir = QDir (url).path ();
            load_file (url.toStdString ());
        }
	}


    event->acceptProposedAction(); 
}







///////////////////////////////////////////////////////////////////////////////////////////////////
void DDWin::not_impl_slot ()
{
    QMessageBox::about( this, "About Zodiac" ,
			QString(("ZODIAC. Version " + VERSION + "\n"
			"Not Implemented yet.. please be patient.. ").c_str()) );
}


void DDWin::builder_slot () {
    builder_menu->show ();
    builder_menu->raise ();
}

void DDWin::twodwin_slot () {
	data ->twodwin ->show ();
	data ->twodwin ->raise ();
}

void DDWin::history_slot () {
    undo_view -> show ();
    undo_view -> raise ();
}


void DDWin::wiimote_slot () {
    HeadTrackingThread *thread = new HeadTrackingThread (0, this);
	run_thread (thread);
//    thread -> start ();	
}

void DDWin::wiimote2_slot () {
    WiimoteTrackingThread *thread = new WiimoteTrackingThread (0, this);
	run_thread (thread);
//    thread -> start ();	
}


void DDWin::pref_slot () {
	pref_menu ->display ();
}


void DDWin::DD_settings_slot () {
    DDsettings_menu->show ();
    DDsettings_menu->raise ();
}

void DDWin::enable_clipping_slot () {
	if (current_target) {
		set_clippable (target_molecule, true);
	}
}

void DDWin::disable_clipping_slot () {
	if (current_target) {
		set_clippable (target_molecule, false);
	}
}



void DDWin::display_settings_slot () {
	display_menu -> display ();
}

void DDWin::sequence_slot () {
	sequence_menu -> display ();
}

void DDWin::about_slot ()
{
    QMessageBox::about( this, "About Zodiac" ,
			QString(("ZODIAC. I feel good! Version "+VERSION+ "\n"
			"Code by Nicola Zonta. All rights reserved.\n"
			"For bug reports, suggestions or any other feedback please visit our website at www.zeden.org or email me at zontan@cf.ac.uk\n\n").c_str()) );
}




void DDWin::setup_accel () {
    accel->connectItem( accel->insertItem(Qt::Key_Escape), this, SLOT(esc_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_Delete), this, SLOT(del_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_Backspace), this, SLOT(del_pressed()) );
	accel->connectItem( accel->insertItem(Qt::Key_A), this, SLOT(a_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_B), this, SLOT(b_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_C), this, SLOT(c_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_D), this, SLOT(d_pressed()) );	
    accel->connectItem( accel->insertItem(Qt::Key_E), this, SLOT(e_pressed()) );
	accel->connectItem( accel->insertItem(Qt::Key_F), this, SLOT(f_pressed()) );
	accel->connectItem( accel->insertItem(Qt::Key_H), this, SLOT(h_pressed()) );		
	accel->connectItem( accel->insertItem(Qt::Key_I), this, SLOT(i_pressed()) );	
    accel->connectItem( accel->insertItem(Qt::Key_M), this, SLOT(m_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_N), this, SLOT(n_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_O), this, SLOT(o_pressed()) );
	accel->connectItem( accel->insertItem(Qt::Key_R), this, SLOT(r_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_S), this, SLOT(s_pressed()) );

}

void DDWin::deselect_slot () {
	deselect();
}
void DDWin::select_all_slot () {
	if (current_target) {
		ZNMolecule *target = 0;
		if (!target_molecule -> selection)  target = target_molecule;
		else target = ((Selection *) target_molecule) ->get_molecules ()[0];
		deselect ();
		Selection *sel = new Selection;
		FOR_ATOMS_OF_MOL (a, target){
			sel ->select_atom (&*a);
		}
		find_center (sel);
		find_limits (sel);
		add_molecule (sel);
		set_current_target (-1);
	}

}

void DDWin::invert_selection_slot () {
	if (current_target) {
		ZNMolecule *target = 0;
		if (!target_molecule -> selection)  target = target_molecule;
		else target = ((Selection *) target_molecule) ->get_molecules ()[0];


		vector <Atom *> atoms;
		FOR_ATOMS_OF_MOL (a, target){
			if (!get_selected (&*a)) atoms.push_back (&*a);
		}

		deselect ();
		Selection *sel = new Selection;
		sel ->select_atoms (atoms);
		find_center (sel);
		find_limits (sel);
		add_molecule (sel);
		set_current_target (-1);
	}
}

void DDWin::select_C_slot () {
	gl -> select_MW (6);
}


void DDWin::select_H_slot () {
	gl -> select_MW (1);
}


void DDWin::select_square_slot () {
	gl ->select_square_mode = true;
	gl ->setCursor(QCursor (QPixmap (":icons/pointer_select_square.png"), 0, 0));
}



void DDWin::build_C_slot () {
	if (mode == MAGIC_PENCIL) {
	        builder->set_magic_pencil_atomic_number (6);
	}
    else builder->add_atom (6);
}

void DDWin::build_S_slot () {
	if (mode == MAGIC_PENCIL) {
	        builder->set_magic_pencil_atomic_number (16);
	}
    else builder->add_atom (16);
}

void DDWin::build_O_slot () {
	if (mode == MAGIC_PENCIL) {
	        builder->set_magic_pencil_atomic_number (8);
	}
    else builder->add_atom (8);
}

void DDWin::build_N_slot () {
	if (mode == MAGIC_PENCIL) {
	        builder->set_magic_pencil_atomic_number (7);
	}
    else builder->add_atom (7);
}


void DDWin::build_single_bond_slot () {
    builder->add_bond (1);
}

void DDWin::build_double_bond_slot () {
    builder->add_bond (2);
}

void DDWin::build_triple_bond_slot () {
    builder->add_bond (3);
}

void DDWin::delete_bond_slot () {
 //   builder->remove_bond ();
}

//key slots

void DDWin::magic_pencil_toggle_slot (bool checked) {
	if (checked) {
		magic_pencil_start ();
	}
	else {
		magic_pencil_end ();
	}
}
void DDWin::magic_pencil_start () {
		magicpencilButt -> setChecked (true);
        switch_mode (MAGIC_PENCIL);
        builder->set_magic_pencil_atomic_number (6);
}
void DDWin::magic_pencil_end () {
	magicpencilButt -> setChecked (false);
	gl ->unsetCursor();
	gl -> magic_pencil = false;
	builder->last_magic_pencil_atom = NULL;
	switch_mode (BUILDER);
}

void DDWin::del_pressed () {
    if (mode == MAGIC_PENCIL) {
        builder->set_magic_pencil_atomic_number (-2);    

    }
    if (mode == NONE) {
        if (current_target)        delete_current_molecule ();
    }
}

void DDWin::esc_pressed () {
	if (mode == MAGIC_PENCIL) {
		magic_pencil_end ();
	}
	else {
		switch_mode (NONE);
		gl ->unsetCursor();
	}
}



void DDWin::a_pressed () {
	if (mode == SELECT) {
		selectAllAct ->trigger ();
	}
}


void DDWin::b_pressed () {
	if (mode == NONE) switch_mode (BUILDER);

}

void DDWin::c_pressed () {
    if (mode == MAGIC_PENCIL) {
        builder->set_magic_pencil_atomic_number (6);

    }
	if (mode == SELECT) {
		selectCAct ->trigger ();
	}
}

void DDWin::d_pressed () {
	if (mode == SELECT) {
		deselectAct ->trigger ();
	}
}

void DDWin::e_pressed () {
    if (mode == SELECT) {
	switch_mode (SELECT_EXTEND);      
    }
}

void DDWin::f_pressed () {
	if (mode == SELECT_EXTEND) {
		if (target_molecule ->selection) {
			Selection *sel = (Selection *) target_molecule;
			sel ->extend_to_fragment ();
		}

	}
    
}


void DDWin::h_pressed () {
	if (mode == SELECT) {
		selectHAct ->trigger ();
	}
}

void DDWin::m_pressed () {
    if (mode == BUILDER) {
		magic_pencil_start ();
        
    }
}


void DDWin::n_pressed () {
	if (mode == SELECT_EXTEND) {
		if (target_molecule ->selection) {
			Selection *sel = (Selection *) target_molecule;
			sel ->extend_to_neighbours ();
		}

	}
    if (mode == MAGIC_PENCIL) {
        builder->set_magic_pencil_atomic_number (7);    
    }
}

void DDWin::i_pressed () {
	if (mode == SELECT) {
		invertSelectAct ->trigger ();
	}
}

void DDWin::o_pressed () {
    if (mode == MAGIC_PENCIL) {
        builder->set_magic_pencil_atomic_number (8);    
    }
}


void DDWin::r_pressed () {
	if (mode == SELECT_EXTEND) {
		if (target_molecule ->selection) {
			Selection *sel = (Selection *) target_molecule;
			sel ->extend_to_residue ();
		}

	}
    
}


void DDWin::s_pressed () {
	if (mode == NONE) switch_mode (SELECT);
	else  if (mode == MAGIC_PENCIL) {
        builder->set_magic_pencil_atomic_number (16);
	}
	else if (mode == SELECT) {
		selectsquareAct ->trigger ();
	}
    
}

void DDWin::switch_mode (int m) {
	hide_toolbars ();
	if (m == SELECT) {
		mode = SELECT;
		select_tb -> show ();
		mode_string = "Select";
	}
	else if (m == SELECT_EXTEND) {
		mode = SELECT_EXTEND;
		mode_string = "Select - Extend";
	}
	else if (m == BUILDER) {
		mode = BUILDER;
		builder_tb -> show ();
		mode_string = "Builder";
	}
	else if (m == MAGIC_PENCIL) {
		mode = MAGIC_PENCIL;
		builder_tb -> show ();
		mode_string = "Builder";
	}
	else {
		mode = NONE;
		mode_string = "";
	}

	current_mode ->setText (QString(mode_string.c_str ()));

}


void DDWin::lock_editing () {
	if (edit_lock_counter < 0) edit_lock_counter = 0;
	edit_lock_counter ++;
	builderAct -> setEnabled (false);
	minimiseAct -> setEnabled (false);
	builder_menu -> hide ();
	if (mode == MAGIC_PENCIL || mode == BUILDER) esc_pressed ();
}


void DDWin::unlock_editing () {
	edit_lock_counter --;
	if (edit_lock_counter < 1) {
		builderAct -> setEnabled (true);
		minimiseAct -> setEnabled (true);
	}
}








////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

MyGl::MyGl (DDWin *parent)
   :    QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{


    setMouseTracking (true);
    ddwin = parent;
	setAcceptDrops(true) ; 
    init_vars ();
    void    GL ();


}

void MyGl::init_vars (){
   
    needs_GL_update = false;
    next_list = 0;


//    bindingsite_list = new_list ();
//    waters_list = new_list ();

    select_color  = SELECT_COLOR ;



    aromatic_display_style = AROMATIC_DISPLAY_STYLE;


    aromatic_bond_inter_distance =AROMATIC_BOND_INTER_DISTANCE;


    double_bond_stick_radius_scale = DOUBLE_BOND_STICK_RADIUS_SCALE;


    surface_resolution = SURFACE_RESOLUTION;

    fog_begin = FOG_BEGIN;

    stereo_toe_in_angle = 1.14576f;
    stereo_inter_eye_semi_distance = 0.2f;


    head_tracking_x_position = 0.f;
    head_tracking_y_position = 0.f;

    center_of_rotation = vect ();
    center_of_view = vect ();
	view_translations = vect ();

    zbeginy = 0.0f;
    tbeginx = 0.0f;
    tbeginy = 0.0f;


	up_point = vect (0., 1., 0.);
	down_point = vect (0., -1., 0.);
	left_point = vect (-1., 0., 0.);
	right_point = vect (1., 0., 0.);
	
    zoom = false; translate = false; rotate = false; select = false; move = false; spin = false; selection_square = false; magic_pencil = false; select_square_mode = false;
    clicked_atom = new Atom;




    Transform.M[0] = 1.0f; Transform.M[1] = 0.0f; Transform.M[2] = 0.0f; Transform.M[3] = 0.0f;
    Transform.M[4] = 0.0f; Transform.M[5] = 1.0f; Transform.M[6] = 0.0f; Transform.M[7] = 0.0f;
    Transform.M[8] = 0.0f; Transform.M[9] = 0.0f; Transform.M[10] = 1.0f; Transform.M[11] = 0.0f;
    Transform.M[12] = 0.0f; Transform.M[13] = 0.0f; Transform.M[14] = 0.0f; Transform.M[15] = 1.0f;

    Head_Tracking_Transf.M[0] = 1.0f; Head_Tracking_Transf.M[1] = 0.0f; Head_Tracking_Transf.M[2] = 0.0f; Head_Tracking_Transf.M[3] = 0.0f;
    Head_Tracking_Transf.M[4] = 0.0f; Head_Tracking_Transf.M[5] = 1.0f; Head_Tracking_Transf.M[6] = 0.0f; Head_Tracking_Transf.M[7] = 0.0f;
    Head_Tracking_Transf.M[8] = 0.0f; Head_Tracking_Transf.M[9] = 0.0f; Head_Tracking_Transf.M[10] = 1.0f; Head_Tracking_Transf.M[11] = 0.0f;
    Head_Tracking_Transf.M[12] = 0.0f; Head_Tracking_Transf.M[13] = 0.0f; Head_Tracking_Transf.M[14] = 0.0f; Head_Tracking_Transf.M[15] = 1.0f;

    LastRot.M[0] = 1.0f; LastRot.M[1] = 0.0f; LastRot.M[2] = 0.0f;
    LastRot.M[3] = 0.0f; LastRot.M[4] = 1.0f; LastRot.M[5] = 0.0f;
    LastRot.M[6] = 0.0f; LastRot.M[7] = 0.0f; LastRot.M[8] = 1.0f;

 //   Last_Head_Tracking_Rot.M[0] = 1.0f; Last_Head_Tracking_Rot.M[1] = 0.0f; Last_Head_Tracking_Rot.M[2] = 0.0f;
 //   Last_Head_Tracking_Rot.M[3] = 0.0f; Last_Head_Tracking_Rot.M[4] = 1.0f; Last_Head_Tracking_Rot.M[5] = 0.0f;
 //   Last_Head_Tracking_Rot.M[6] = 0.0f; Last_Head_Tracking_Rot.M[7] = 0.0f; Last_Head_Tracking_Rot.M[8] = 1.0f;




    ThisRot.M[0] = 1.0f; ThisRot.M[1] = 0.0f; ThisRot.M[2] = 0.0f;
    ThisRot.M[3] = 0.0f; ThisRot.M[4] = 1.0f; ThisRot.M[5] = 0.0f;
    ThisRot.M[6] = 0.0f; ThisRot.M[7] = 0.0f; ThisRot.M[8] = 1.0f;

	select_pulse = 0;


	startTimer (30); 



}


void MyGl::timerEvent ( QTimerEvent * event ) {
	select_pulse ++; //int that controls the pulsating color for selections. 
	set_clipping_plane (); //rotates the general clipping plane to be parallel to the viewer
	for (unsigned int i = 0; i < ddwin ->molecules.size (); i++) {
		if (get_needs_redraw(ddwin ->molecules[i])) {
			GL_update_molecule(ddwin ->molecules[i]);
		}
		if (get_needs_backbone_redraw (ddwin ->molecules[i])) {
			GL_update_backbone(ddwin ->molecules[i]);
		}
	}
	for (unsigned int i = 0; i < ddwin ->database_grids.size (); i++) {
		if (ddwin ->database_grids[i] ->database ->get_needs_redraw ()) {
			ddwin ->database_grids[i] ->update_graphics ();
		}
	}
	if (!ddwin ->running_threads.size ()) {
		ddwin ->thread_menu ->hide ();
	}
	else {
	ddwin ->thread_menu -> show ();
	for (unsigned int i = 0; i < ddwin ->running_threads.size (); i++) {
		if (ddwin ->running_threads[i] -> isFinished ()) {
			Thread *th = ddwin ->running_threads[i];
			ddwin ->running_threads.erase (ddwin ->running_threads.begin ()+i);
			ddwin ->thread_menu -> clear (i);
			i --;
			delete th;
		}
		else {
			ddwin ->thread_menu -> display_thread (i, ddwin ->running_threads[i]);
		}
	}
	}
	ddwin ->haptic_menu ->update_energy (); //triggers an update of the energy value if haptic simulation is running
	ddwin ->haptic_menu ->maybe_save_result ();
	updateGL (); //redraw the graphic scene
}

void MyGl::set_clipping_plane () {
    GLdouble x1, x2, x3, y1, y2, y3, z1, z2, z3;
    GLint viewport [4];
    GLdouble model [16];
    GLdouble proj [16]; 
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
	double centerz = 0.995;
	
	gluUnProject (viewport[2]/2, viewport[3], centerz, model, proj, viewport, &x1, &y1, &z1);
    gluUnProject (viewport[2], viewport[3]/2, centerz, model, proj, viewport, &x2, &y2, &z2);
    gluUnProject (viewport[2]/2, viewport[3]/2, centerz, model, proj, viewport, &x3, &y3, &z3);

	GLdouble a, b, c, d;
	a = y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2);
	b = z1*(x2-x3) + z2*(x3-x1) + z3*(x1-x2);
	c = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
	d = -(x1* (y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1));
	double eq [4] = {a, b, c, d};
	glClipPlane (GL_CLIP_PLANE0, eq);
}

void MyGl::free_list (int list) {
    if (next_list == list +1) next_list--;
    else empty_lists.push_back (list); 
}



int MyGl::new_list () {
    int out;
    if (empty_lists.size() ) {
        out = empty_lists [empty_lists.size()-1];
        empty_lists.pop_back ();
    
    }
    else {
        out = next_list;
        next_list++;    
    }
    return out;
}





void MyGl::translate_molecule (ZNMolecule *mol, vect v) {
	FOR_ATOMS_OF_MOL(a, mol) {
    sum_to_coordinates(&*a, v);
  }
	vect c = get_center (mol);
		set_center (mol, sum (c, v));
} 




void MyGl::move_molecule (ZNMolecule *mol, float x, float y, float z, bool cut_bool) {
    vect v (x, y, z);
    if (cut_bool) {
        float mod = v.module (); 
        float cut = 0.5;
        if (mod>cut) {
            v.scale_to (cut);
        }
    }
    FOR_ATOMS_OF_MOL(a, mol) {
        sum_to_coordinates(&*a, v);
    }
	vect c = get_center (mol);
	set_center (mol, sum (c, v));
}






void MyGl::initializeGL ()
{
	GLboolean stereo;
	glGetBooleanv(GL_STEREO, &stereo);

    // IJG: Needs to be populated from a tick box in the GUI!!!
    ddwin -> g_stereoMode = DDWin::StereoMode_None;

    // Should we automatically select stereo mode if we discovery stereo
    // is in use?
#if 0
    if (stereo)
    {
        ddwin->g_stereoMode = DDWin::StereoMode_ViaOpenGL;
    }
#endif


  //  glEnable( GL_CULL_FACE );
    glEnable (GL_NORMALIZE);

    glLightfv(GL_LIGHT0,GL_SPECULAR,specular);
	color background_color = *ddwin ->data ->background_color;
	glClearColor (background_color.redF (), background_color.greenF (), background_color.blueF (), background_color.alphaF ());
	GLfloat fogColor[4]= {background_color.redF (), background_color.greenF (), background_color.blueF (), background_color.alphaF ()};	
    glFogfv(GL_FOG_COLOR, fogColor);	
	
	glClearDepth (1.0f);											// Depth Buffer Setup
	glDepthFunc (GL_LEQUAL);										// The Type Of Depth Testing (Less || Equal)

  //  g_openGlStereoAvailable = GL_TRUE;    
    if (ddwin->g_stereoMode == DDWin::StereoMode_ViaOpenGL) {glEnable (GL_STEREO);}

	glEnable (GL_DEPTH_TEST);										// Enable Depth Testing
										// Select Flat Shading (Nice Definition Of Objects)
//	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);				// Set Perspective Calculations To Most Accurate

	quadratic=gluNewQuadric();										// Create A Pointer To The Quadric Object
	gluQuadricNormals(quadratic, GLU_SMOOTH);						// Create Smooth Normals
	gluQuadricTexture(quadratic, GL_TRUE);							// Create Texture Coords

	glEnable(GL_LIGHT0);											// Enable Default Light
	glEnable(GL_LIGHTING);											// Enable Lighting
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,1.0);
 //   glEnable(GL_POINT_SMOOTH);                                       //antialiasing
 //   glEnable(GL_LINE_SMOOTH);
//	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	
	
//    glEnable(GL_POLYGON_SMOOTH);                                      
//	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	
	glShadeModel(GL_SMOOTH);	
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,specular);
    glShadeModel(GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);	
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);			// Enable Alpha Blending (disable alpha testing)
	glEnable(GL_BLEND);	



  //  glFogi(GL_FOG_MODE, GL_LINEAR);		

    glFogf(GL_FOG_DENSITY, 0.005f);				
    glHint(GL_FOG_HINT, GL_DONT_CARE);			
    glFogf(GL_FOG_START, fog_begin);			
    glFogf(GL_FOG_END, 10000.0f);				
    glEnable(GL_FOG);				

}
/*
void MyGl::haptic_timer () {
        QTime t = QTime::currentTime ();
        float diff = time.msecsTo (t);
        if (diff <1) diff = 1;
        float fps = 1000/diff;
        stringstream ss;
        ss<<fps<<" fps"<<endl;
        ddwin->mode_string = ss.str ();
        time = t;
        ZNMolecule * min_mol = ddwin->data->minimize->interaction_ff->target_mol;
        ddwin->data->minimize->haptic_step ();
        draw_molecule (min_mol);
}

*/
/*
void MyGl::timerEvent ( QTimerEvent * ) {
 //   ddwin->mode_string = "";
    float f = 0.001;
 //   if (ddwin->mode == SELECT) {
 //       ddwin->mode_string = "Select";
 //   }
 //   if (ddwin->mode == BUILDER) {
 //       ddwin->mode_string = "Builder";
 //   }


    if (ddwin->haptic_mode) {
        haptic_timer ();
    }

    if (ddwin->minimizing) {
        ddwin->data->minimize->minimize_energy_step ();
    }
    if (ddwin->rec_movie) {
        QTime t = QTime::currentTime ();
        int msecs = movie_time.msecsTo (t);
        if (msecs>1000/30) {
            movie_time = t;
            stringstream ss;
            ss << ddwin->last_frame;
            screenshot (QString("__movie__") + QString(ss.str().c_str()) + QString(".ppm"));
            ddwin->last_frame++;
        }
    }
    if (ddwin->rec_pov_movie) {
        QTime t = QTime::currentTime ();
        int msecs = movie_time.msecsTo (t);
        if (msecs>1000/30) {
            movie_time = t;
            stringstream ss;
            ss << ddwin->last_frame;
            string filename ="__movie__"+ss.str()+".pov";        
            ddwin->write_POV_source (filename);
            ddwin->last_frame++;
        }
    }
}
*/


void MyGl::paintGL () {
	

		color background_color = *ddwin ->data ->background_color;

	GLfloat fogColor[4]= {background_color.redF (), background_color.greenF (), background_color.blueF (), background_color.alphaF ()};
	    glFogfv(GL_FOG_COLOR, fogColor);
		glClearColor (background_color.redF (), background_color.greenF (), background_color.blueF (), background_color.alphaF ());
  //  cout <<"paintGL"<<endl;

	glLineWidth (*ddwin ->data ->line_width);


    GLdouble farx, fary, farz;
    GLdouble nearx, neary, nearz;
    GLdouble frontx, fronty, frontz;
	GLdouble rightx, righty, rightz, upx, upy, upz, downx, downy, downz, leftx, lefty, leftz;
    GLint viewport [4];
    GLdouble model [16];
    GLdouble proj [16]; 
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);

				
    double centerz = 0.995;
	

    gluUnProject (viewport[2]/2, viewport[3], centerz, model, proj, viewport, &upx, &upy, &upz);
    gluUnProject (viewport[2], viewport[3]/2, centerz, model, proj, viewport, &rightx, &righty, &rightz);
    gluUnProject (viewport[2]/2, 0, centerz, model, proj, viewport, &downx, &downy, &downz);
    gluUnProject (0, viewport[3]/2, centerz, model, proj, viewport, &leftx, &lefty, &leftz);

    gluUnProject (viewport[2]/2, viewport[3]/2, 1, model, proj, viewport, &farx, &fary, &farz);
    gluUnProject (viewport[2]/2, viewport[3]/2, 0, model, proj, viewport, &nearx, &neary, &nearz);
    float modfront = sqrt ((nearx-farx)*(nearx-farx)+(neary-fary)*(neary-fary)+(nearz-farz)*(nearz-farz));
    frontx = (nearx-farx)/modfront;
    fronty = (neary-fary)/modfront;
    frontz = (nearz-farz)/modfront;


	up_point = vect (upx, upy, upz);
	down_point = vect (downx, downy, downz);
	left_point = vect (leftx, lefty, leftz);
	right_point = vect (rightx, righty, rightz);
	
	
	
	int numFrameBuffers = 1;

    if (ddwin->g_stereoMode != DDWin::StereoMode_None)
	{
		numFrameBuffers = 2;
	}


	for (int frameBuffer = 0; frameBuffer < numFrameBuffers; frameBuffer++)
	{
		if (ddwin->g_stereoMode == DDWin::StereoMode_ViaOpenGL)
		{
			if (frameBuffer == 0)
			{
				glDrawBuffer(GL_BACK_LEFT);
			}
			else
			{
				glDrawBuffer(GL_BACK_RIGHT);
			}
		}
		else
		{
			glDrawBuffer(GL_BACK);
		}

        // IJG
        if( ddwin->g_stereoMode == DDWin::StereoMode_VerticalInterlace ||
            ddwin->g_stereoMode == DDWin::StereoMode_HorizontalInterlace ) 
        {
            ddwin->g_stencil_mask_frame_counter++;

#if 0
			// Occasionally force a redraw... we never know if (eg) someone else
            // has nuked the stencil buffer...
            if ((ddwin->g_stencil_mask_frame_counter % 20) == 0)
            {
                ddwin->g_stencil_mask_needs_redraw = true;
            }
#endif

            if (ddwin->g_stencil_mask_needs_redraw)
            {
                ddwin->g_stencil_mask_needs_redraw = false;

                glDrawPixels( ddwin->g_stencil_mask_width, ddwin->g_stencil_mask_height, 
                              GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, ddwin->g_stencil_mask );
            }

            // render only every second line
            glEnable(GL_STENCIL_TEST);

            if (frameBuffer == 0)
            {
                // Only screen clear prior to first render - we need to merge both images,
                // so they're rendered to the same buffer.
        		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	

                glStencilFunc(GL_EQUAL, 1, 1);
            }
            else
            {
                glStencilFunc(GL_NOTEQUAL, 1, 1);
            }

            glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP);
        }
        else
        {
            glDisable(GL_STENCIL_TEST);

            // Normal stereo mode or mono mode requires a screen clear before each render
    		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
        }
	

		glDisable (GL_LIGHTING);
		glColor4f (0.5,0.,0.5, 1.);
	//	renderText (20, 20, QString(ddwin->mode_string.c_str()));
		glEnable (GL_LIGHTING);


		glMatrixMode(GL_MODELVIEW);			
		glLoadIdentity();			
	 //   float centerx = translatefx/10; 
	 //   float centery = -translatefy/10; 
	 //   float centerz = -zoomf/10; 
	  //  glTranslatef (centerx, centery, centerz-15.0);	

								
        matrix_transformations (true, frameBuffer);

        // IJG
#if 0
        if( mirror_in_y ) 
        {
            glScalef( 1, -1, 1 );
            glFrontFace( GL_CW );
        } 
        else 
        {
            glFrontFace( GL_CCW );
        }
#endif // 0
#ifdef GL_MULTISAMPLE
	glEnable (GL_MULTISAMPLE);
#endif //GL_MULTISAMPLE
		for (unsigned int i =1; i<ddwin->molecules.size(); i++) {
			if (!ddwin->molecules[i]->selection) {  //selections are not drawn
				if (get_clippable (ddwin->molecules[i])) glEnable (GL_CLIP_PLANE0);
                GLdouble cx, cy, cz;
				vect cent = get_center (ddwin ->molecules[i]);
			    gluProject (cent.x(), cent.y(), cent.z(),model, proj, viewport,&cx, &cy, &cz);
 
            //    glFogf(GL_FOG_START, cz);			
            //    glFogf(GL_FOG_END, cz+0.000002);
         //       GLfloat f;
          //      glGetFloatv (GL_FOG_START, &f);
           //     cerr << f <<endl;
         //       cerr << cz*10000<<endl;
            //    glEnable(GL_FOG);
               // float cx = ddwin -> molecules[i] -> center.x();
             //   float cy = ddwin -> molecules[i] -> center.y();
           //     float cz = ddwin -> molecules[i] -> center.z();
         //       float nz = cx*model[2] + cy*model[6] + cz*model[10] + model[14];
     //           cerr << nz << endl;
       //         glFogf(GL_FOG_START, -nz);			
       //         glFogf(GL_FOG_END, -nz+20);
				glDisable (GL_LIGHTING);
				glCallList (get_line_list (ddwin->molecules[i]));
				glLineWidth (3.f);
				glCallList (get_backbone_list1 (ddwin->molecules[i])); 
				glLineWidth (*ddwin ->data ->line_width);
				glEnable(GL_LIGHTING);
				glCallList (get_stick_list (ddwin->molecules[i])); 
				glCallList (get_backbone_list2 (ddwin->molecules[i])); 
				if (ddwin->show_labels && 0) {
                    FOR_ATOMS_OF_MOL(at, ddwin -> molecules[i]) {
						glColor4f (0.,0.,0.,1.);
						stringstream ss, ss2;

					}
				}
			}
			glDisable (GL_CLIP_PLANE0);
		}

		for (unsigned int i =0; i<ddwin->graphical_objects.size(); i++) {
			if (ddwin->graphical_objects[i]->is_clippable ()) glEnable (GL_CLIP_PLANE0);
			if (ddwin->graphical_objects[i]->mesh) glDisable (GL_LIGHTING);
			glCallList (ddwin->graphical_objects[i]->lst);
			glEnable(GL_LIGHTING);
			glDisable (GL_CLIP_PLANE0);
		}
		glDisable (GL_LIGHTING);
		for (unsigned int i = 0; i<ddwin ->Hbonds.size (); i++) {
			color c (1.f, 0.f, 0.f, ddwin ->Hbonds[i].perc);
			my_line(ddwin ->Hbonds[i].v1, ddwin ->Hbonds[i].v2, c, c);
		}
//		for (unsigned int i = 0; i<ddwin ->vertex_list.size (); i++) {
//			color c (1.f, 0.f, 0.f, 1.f);
//			my_line(ddwin ->vertex_list[i], sum (ddwin ->vertex_list[i], vect(0.1,0.,0.)), c, c);
//		}
		glEnable(GL_LIGHTING);
#ifdef GL_MULTISAMPLE
		glDisable (GL_MULTISAMPLE);
#endif //GL_MULTISAMPLE

		// SELECTION SQUARE
		if (selection_square) {
			GLdouble x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
			gluUnProject (selection_square_x1, viewport[3]-selection_square_y1, 0, model, proj, viewport, &x1, &y1, &z1);
			gluUnProject (selection_square_x1, viewport[3]-selection_square_y2, 0, model, proj, viewport, &x2, &y2, &z2);
			gluUnProject (selection_square_x2, viewport[3]-selection_square_y2, 0, model, proj, viewport, &x3, &y3, &z3);
			gluUnProject (selection_square_x2, viewport[3]-selection_square_y1, 0, model, proj, viewport, &x4, &y4, &z4);
			glDisable (GL_LIGHTING);
            set_color (select_color);
			glBegin (GL_LINE_STRIP);
			glVertex3d (x1, y1, z1);
			glVertex3d (x2, y2, z2);
			glVertex3d (x3, y3, z3);
			glVertex3d (x4, y4, z4);
			glVertex3d (x1, y1, z1);
	    
			glEnd ();
			glEnable(GL_LIGHTING);
		}
	   // MAGIC PENCIL
		if (magic_pencil) {
			glDisable (GL_LIGHTING);
            set_color (select_color);
			glBegin (GL_LINE_STRIP);
			for (unsigned int i=0; i<ddwin->magic_pencil_trail_x.size (); i++) {
				glVertex3d (ddwin->magic_pencil_trail_x[i],ddwin->magic_pencil_trail_y[i],ddwin->magic_pencil_trail_z[i]);
			}
			glEnd ();
			glEnable(GL_LIGHTING);
		}
		if (ddwin->builder->last_magic_pencil_atom) {
			Atom *last = ddwin->builder->last_magic_pencil_atom;
			glDisable (GL_LIGHTING);
			GLdouble lx, ly, lz, lax, lay, laz, lxx, lyy, lzz;
			vect v = get_coordinates(last);
			gluProject (v.x(), v.y(), v.z(),model, proj, viewport,&lax, &lay, &laz);
			gluUnProject (lax, lay, 0, model, proj, viewport, &lxx, &lyy, &lzz);
			gluUnProject (lastx, viewport[3]-lasty, 0, model, proj, viewport, &lx, &ly, &lz);
            set_color (select_color);
			glBegin (GL_LINES);
	//        cout <<lxx <<lyy<<lzz<<endl;
	//        cout <<lx <<ly<<lz<<endl;
			glVertex3d (lxx, lyy, lzz);
			glVertex3d (lx, ly, lz);


			glEnd ();
			glEnable(GL_LIGHTING);
		}
		ddwin ->haptic_menu ->restrain_lock ->lockForRead ();
		for (unsigned int i = 0; i < ddwin ->haptic_menu ->restrains.size (); i++) {
			color c  (0.f, 1.f, 0.f, 1.f);
			glDisable (GL_LIGHTING);
			my_line (get_coordinates (ddwin ->haptic_menu ->restrains[i] ->at1), get_coordinates (ddwin ->haptic_menu ->restrains[i] ->at2), c, c);
			glEnable (GL_LIGHTING);
		}
		ddwin ->haptic_menu ->restrain_lock ->unlock ();
	}
	glFlush ();				

	needs_GL_update = false;
}

void MyGl::refreshStencilBuffer()
{
    // IJG
    if ( (ddwin->g_stereoMode == DDWin::StereoMode_HorizontalInterlace) 
        || (ddwin->g_stereoMode == DDWin::StereoMode_VerticalInterlace) ) //this line read horizontal as well 
    {
        // build up stencil buffer
        int wp2 = width();
        // width needs to be a power of 2 because of how glDrawPixels work.
        // The extra values will be ignored.
        wp2 += 4 - (width() % 4);
        int h = height();

        // Only release the memory if it needs to be reallocated...
        if ( (ddwin->g_stencil_mask_width != wp2)
            || (ddwin->g_stencil_mask_height != h) )
        {
            if ( ddwin -> g_stencil_mask ) 
            {
                free( ddwin -> g_stencil_mask );
            }

            ddwin -> g_stencil_mask = (unsigned char *) malloc( wp2*h );
            ddwin -> g_stencil_mask_width = wp2;
            ddwin -> g_stencil_mask_height = h;
        }

        // Always redraw the stencil contents just in case...
        // we may be the same size window, but could be a different
        // interlace from before.
        if ( ddwin -> g_stereoMode == DDWin::StereoMode_HorizontalInterlace ) 
        {
            for( int i = 0; i < h; i++ )
                for( int j = 0; j < wp2; j++ )
                    ddwin -> g_stencil_mask[i*wp2+j]=(i+1)%2;
        } 
        else if ( ddwin -> g_stereoMode == DDWin::StereoMode_VerticalInterlace ) 
        {
            for( int i = 0; i < h; i++ )
                for( int j = 0; j < wp2; j++ )
                    ddwin -> g_stencil_mask[i*wp2+j]=(j+1)%2;
        }

        ddwin->g_stencil_mask_needs_redraw = true;
    }
    else
    {
        // Don't need a stencil, so release the memory
        if( ddwin -> g_stencil_mask )
		{
			free( ddwin -> g_stencil_mask );
			ddwin -> g_stencil_mask = 0;
		}
    }
}

void MyGl::resizeGL( int width, int height )
{
    GLfloat w = (float) width / (float) height;

  //  cout<<"Viewport resized to "<<width<<" "<<height<<endl;

    glViewport( 0, 0, width, height );
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluPerspective(45.0f, w, 0.1f, 10000.0f);
//    glFrustum( -w, w, -h, h, 1.0, 10000.0 );
    ArcBall.setBounds(width, height);

    refreshStencilBuffer();
}

void MyGl::matrix_transformations (bool stereo, int frameBuffer) {


    if ((ddwin->g_stereoMode != DDWin::StereoMode_None) && stereo)
	{
		if (frameBuffer == 0)
		{
            glRotatef (+stereo_toe_in_angle, 0.f, 1.f, 0.f); 
			glTranslatef(-stereo_inter_eye_semi_distance, 0.0f, 0.0f);
		}
		else
		{
            glRotatef (-stereo_toe_in_angle, 0.f, 1.f, 0.f);
			glTranslatef(+stereo_inter_eye_semi_distance, 0.0f, 0.0f);
		}
	}



	glMultMatrixf(Head_Tracking_Transf.M);
        glTranslatef (head_tracking_x_position, -head_tracking_y_position, 0.f);


        glTranslatef (0, 0, -15.f);
		glTranslatef (view_translations.x(),view_translations.y(),view_translations.z());
		
	//	vect new_center_of_view = rotate_vector (subtract (center_of_view, center_of_rotation));
	//	glTranslatef (new_center_of_view[0],new_center_of_view[1],new_center_of_view[2]);
		

     //   Transform.M[12] = center_of_rotation.x();
     //   Transform.M[13] = center_of_rotation.y();
     //   Transform.M[14] = center_of_rotation.z();

		glMultMatrixf(Transform.M);


		glTranslatef (-center_of_rotation.x(),-center_of_rotation.y(),-center_of_rotation.z());

	//	glTranslatef(-center_of_rotation.x(),-center_of_rotation.y(),-center_of_rotation.z());


}







void MyGl::mousePressEvent ( QMouseEvent * e )
{   
	
    float x = e->x(); float y = e->y();
	Qt::MouseButton button = e->button ();

        if (button == Qt::RightButton)
        {
            if (ddwin->mode == MAGIC_PENCIL) {
                begin_magic_pencil (x, y);
            }
            else {
                select = true;
                begin_zoom (y);
            }
        }
        else if (button == Qt::MidButton)
        {
            if (e->modifiers ()==Qt::ShiftModifier) {
                begin_selection_square (x, y);
            }
            else begin_translate (x, y);
        }
        else if (button == Qt::LeftButton)
        {


            if (e->modifiers () == Qt::ShiftModifier) {
                if (ddwin->current_target!=0)
                    begin_move (x, y);
            }
            else if (e->modifiers () == Qt::ControlModifier) {
                if (ddwin->current_target!=0)
                    begin_spin (x, y);
            }
			else if (select_square_mode) {
				begin_selection_square (x, y);
			}
            else {
                begin_rotate (x, y);
            }
        }

}

void MyGl::mouseReleaseEvent ( QMouseEvent * e )
{
        float x = e->x(); float y = e->y();
		Qt::MouseButton button = e->button ();
        if (button == Qt::RightButton)
        {
            zoom = false;
            if (select){
                 selectf (x,y);
                select = false;
            }
            if (magic_pencil) {
                end_magic_pencil (x, y);
            }    
        }
        else if (button == Qt::MidButton)
        {
            translate = false;
            if (selection_square) end_selection_square ();
        }
        else if (button == Qt::LeftButton)
        {
            rotate = false;
            move = false;
            spin = false;
            if (selection_square) end_selection_square ();

           LastRot = ThisRot;
            

        }
}

void MyGl::mouseMoveEvent ( QMouseEvent * e )
{   

    float x = e->x(); float y = e->y();
    lastx = x; lasty = y;




    if (zoom){
        continue_zoom (y);
        select = false;
    }
    if (translate){
        continue_translate (x, y);
    }
    if (rotate){
        continue_rotate (x, y);
    }
    if (move) {
        continue_move (x, y);
    }
    if (spin) {
        continue_spin (x, y);
    }
    if (selection_square) {
        continue_selection_square (x, y);
    }
    if (magic_pencil)  {
        continue_magic_pencil (x, y);
    }

}

void MyGl::begin_magic_pencil (float x, float y) {
 //   cout <<"begin"<<endl;
    magic_pencil = true;
    Atom *a = select_atom (x, y, 10, 10);
    if (a) {
        ddwin->builder->start_magic_pencil_atom = a;
    }
    else {
        ddwin->builder->start_magic_pencil_atom = NULL;
    }

}

void MyGl::continue_magic_pencil (float x, float y) {
  //  cout <<"continue"<<endl;
 //   cout <<ddwin->magic_pencil_trail_x.size ()<<endl;
    GLint viewport [4];
    GLdouble model [16];
    GLdouble proj [16]; 
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
    GLdouble xx, yy, zz;
    gluUnProject (x, viewport[3]-y, 0, model, proj, viewport, &xx, &yy, &zz);
    ddwin->magic_pencil_trail_x.push_back (xx);
    ddwin->magic_pencil_trail_y.push_back (yy);
    ddwin->magic_pencil_trail_z.push_back (zz);
    ddwin->magic_pencil_trail_wx.push_back (x);
    ddwin->magic_pencil_trail_wy.push_back (y);

}

void MyGl::end_magic_pencil (float x, float y) {


 //   cout <<"end"<<endl;
    GLint viewport [4];
    GLdouble model [16];
    GLdouble proj [16]; 
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
    GLdouble xx, yy, zz;


    if (ddwin->magic_pencil_trail_x.size () < 10) { //click
        if (ddwin->builder->magic_pencil_atomic_number==-2) { //rubber
            Atom *a = select_atom (x, y, 10, 10);
            if (a) { 
                ddwin->builder->delete_atom (a);
            }
        }
        else { //builder
            Atom *a = select_atom (x, y, 10, 10);
            if (a) {  //click on an atom
                if (!a -> IsHydrogen ()) {
                    if (ddwin->builder->last_magic_pencil_atom) { //continuing line
                        if (&*a == ddwin->builder->last_magic_pencil_atom || a -> GetBond (ddwin -> builder -> last_magic_pencil_atom)) { 

                            ddwin->builder->last_magic_pencil_atom = NULL;
                        }
                        
                        else {

                            ddwin->builder->new_bond (a, ddwin->builder->last_magic_pencil_atom);
                        }
                    }
                    ddwin->builder->mutate_atom_to (a, ddwin->builder->magic_pencil_atomic_number);
                    ddwin->builder->last_magic_pencil_atom = a;
                }
            }
            else {  //click on empty spot
                if (ddwin->builder->last_magic_pencil_atom) {
                    GLdouble lx, ly, lz;
                    Atom *last = ddwin->builder->last_magic_pencil_atom;
					vect v = get_coordinates(last);
                    gluProject (v.x(), v.y(), v.z(),model, proj, viewport,&lx, &ly, &lz);
                    gluUnProject (x, viewport[3]-y, lz, model, proj, viewport, &xx, &yy, &zz);
                    Atom *at;
                    vect coord (xx, yy, zz);
                    at = ddwin->builder->add_atom_bonded_to (coord, ddwin->builder->magic_pencil_atomic_number, last);
                    ddwin->builder->last_magic_pencil_atom = at;
                    draw_molecule ((ZNMolecule *) last -> GetParent ());
                }
                else {
                    GLdouble lx, ly, lz;
                    gluProject (0, 0, 0,model, proj, viewport,&lx, &ly, &lz);
                    gluUnProject (x, viewport[3]-y, lz, model, proj, viewport, &xx, &yy, &zz);
                    vect coord (xx, yy, zz);
                    Atom *at;
                    at = ddwin->builder->new_atom (coord, ddwin->builder->magic_pencil_atomic_number);
                    ddwin->builder->last_magic_pencil_atom = at;
                    draw_molecule ((ZNMolecule *) at-> GetParent ());                
                }
            }

        }
    }
    else { //drag


            Atom *a = select_atom (x, y, 10, 10);
            if (a && ddwin->builder->start_magic_pencil_atom) {
                if (a!=ddwin->builder->start_magic_pencil_atom) {
                    ZNBond *bo = a-> GetBond (ddwin->builder->start_magic_pencil_atom);
                    if (bo) { // edit bond
                        if (ddwin->builder->magic_pencil_atomic_number==-2) { //rubber
                            if (bo-> IsTriple ()) ddwin->builder->set_bond (bo, 2);
                            else if (bo-> IsDouble ()) ddwin->builder->set_bond (bo, 1);
                            else if (bo -> IsSingle ()) ddwin -> builder -> delete_bond (bo);
            //                draw_molecule (bo->GetBeginAtom ()-> GetParent ());
                        }
                        else {
                            if (bo -> IsSingle ()) ddwin->builder->set_bond (bo, 2);
                            else if (bo -> IsDouble ()) ddwin->builder->set_bond (bo, 3);
              //              draw_molecule (bo->GetBeginAtom ()-> GetParent ());
                        }
                    }
    
                }
            }
            else {
                float x1 = ddwin->magic_pencil_trail_wx[0];
                float x2 = ddwin->magic_pencil_trail_wx[ddwin->magic_pencil_trail_wx.size ()-1];
                float y1 = ddwin->magic_pencil_trail_wy[0];
                float y2 = ddwin->magic_pencil_trail_wy[ddwin->magic_pencil_trail_wy.size ()-1];
                if ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2) < 100) {
/*
                    vector<OBRing*>::iterator ring;
                    vector<OBRing*> *rlist = (vector<OBRing*>*)ddwin -> target_molecule -> GetData("RingList");
                    for (ring = rlist->begin();ring != rlist->end();++ring) {
                
                        GLdouble cx, cy, cz, ax, ay, az;
                        gluProject (ring -> center.x(), ring- > center.y(), ring->center.z(),model, proj, viewport,&cx, &cy, &cz);

                        gluProject (ring->atoms[0]-> GetVector ().x(), ring->atoms[0]-> GetVector ().y(), ring->atoms[0]-> GetVector ().z(),model, proj, viewport,&ax, &ay, &az);
                        ay = viewport[3]-ay; //mouse coordinates
                        cy = viewport[3]-cy;
                        float r2 = (ax-cx)*(ax-cx)+(ay-cy)*(ay-cy);
                        bool boo = true;
                        for (unsigned int i=0; i<ddwin->magic_pencil_trail_wx.size (); i+=5) {
                            float ccx = ddwin->magic_pencil_trail_wx[i];
                            float ccy = ddwin->magic_pencil_trail_wy[i];
                 //           cout <<(ccx-cx)*(ccx-cx)+(ccy-cy)*(ccy-cy)<<endl;
                            if ((ccx-cx)*(ccx-cx)+(ccy-cy)*(ccy-cy) > r2) {
                                boo = false;
                                break;
                            }
                        }
                    
                        if (boo) {
                            if (ddwin->builder->magic_pencil_atomic_number==-2) { //rubber
                                ddwin->builder -> set_non_aromatic (ddwin->target_molecule->rings[r]);
                            }
                            else  ddwin->builder -> set_aromatic (ddwin->target_molecule->rings[r]);
                        }
                    }
                
    */
                }
                else {
                    ddwin->builder->last_magic_pencil_atom = NULL;
    
                } 
            }
        }

    magic_pencil = false;
    ddwin->magic_pencil_trail_x.clear ();
    ddwin->magic_pencil_trail_y.clear ();
    ddwin->magic_pencil_trail_z.clear ();
    ddwin->magic_pencil_trail_wx.clear ();
    ddwin->magic_pencil_trail_wy.clear ();


}

void MyGl::begin_selection_square (float x, float  y) {
 //   cout<<"begin square"<<endl;
    selection_square = true;
    selection_square_x1 = x;
    selection_square_y1 = y;
    continue_selection_square (x, y);
}

void MyGl::continue_selection_square (float x, float y) {
    selection_square_x2 = x;
    selection_square_y2 = y;

}


void MyGl::end_selection_square () {

    selection_square = false; //drawing the square
	select_square_mode = false; //waiting for a click to start the square
	unsetCursor();
    select_square ();

}





void MyGl::begin_zoom (float y )
{
    zoom = true;
    zbeginy = y;
}

void MyGl::continue_zoom (float y )
{
    view_translations.z() -= (y - zbeginy);
    zbeginy = y;

}

void MyGl::begin_move (float x, float y )
{
    move = true;
    continue_move (x, y);
}

void MyGl::begin_spin (float x, float y )
{
    spin = true;
    sbeginx = x; 
}

void MyGl::begin_translate (float x, float y )
{
    translate = true;
    tbeginx = x; tbeginy = y;
}

void MyGl::continue_translate (float x, float y )
{
//    translatefx += (x - tbeginx);
//    translatefy += (y - tbeginy);
    view_translations.x() +=(x - tbeginx)/10;
    view_translations.y() -=(y - tbeginy)/10;

    tbeginx = x; tbeginy = y;

}

void MyGl::translate_view (vect v) {
	view_translations  = sum (view_translations, v);
}

void MyGl::continue_move (float x, float y )
{
 //   GLdouble camerax, cameray, cameraz;
 //   GLdouble upx, upy, upz;
//    GLdouble rightx, righty, rightz;
    GLdouble centx, centy, centz;
    GLdouble targx, targy, targz;
    GLint viewport [4];
    GLdouble model [16];
    GLdouble proj [16]; 
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);

 //   gluUnProject (viewport[2]/2, viewport[3]/2, 0, model, proj, viewport, &camerax, &cameray, &cameraz);
//   gluUnProject (viewport[2]/2, viewport[3], 0, model, proj, viewport, &upx, &upy, &upz);
//    gluUnProject (viewport[2], viewport[3]/2, 0, model, proj, viewport, &rightx, &righty, &rightz);
	vect center = get_center (ddwin->target_molecule);
    gluProject (center.x(), center.y(), center.z(),model, proj, viewport,&centx, &centy, &centz);
    gluUnProject (x, viewport[3]-y, centz, model, proj, viewport, &targx, &targy, &targz);

    move_molecule (ddwin->target_molecule, targx-center.x(), targy-center.y(), targz-center.z(), false);
//    move_molecule (ddwin->target_molecule, xx, yy, zz);
    draw_molecule (ddwin->target_molecule);

}

void MyGl::continue_spin (float x, float y )
{
    GLdouble axisx, axisy, axisz, farx, fary, farz;
    GLint viewport [4];
    GLdouble model [16];
    GLdouble proj [16]; 
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);

    gluUnProject (viewport[2]/2, viewport[3]/2, 0, model, proj, viewport, &axisx, &axisy, &axisz);
    gluUnProject (viewport[2]/2, viewport[3]/2, 1, model, proj, viewport, &farx, &fary, &farz);

    axisx-=farx; axisy-=fary; axisz-=farz;
    float mod = sqrt(axisx*axisx+axisy*axisy+axisz*axisz);
    axisx/=mod;     axisy/=mod;     axisz/=mod;
    
    float spin = (x - sbeginx)/50;
//    cout << upx<<" "<<upy<<endl;


    vect axis (axisx, axisy, axisz);

    FOR_ATOMS_OF_MOL(a, ddwin -> target_molecule) {
        vect v = get_coordinates (&*a);
        rotate_around_vector (v, axis, get_center (ddwin->target_molecule), spin);
    }
    draw_molecule (ddwin->target_molecule);
    sbeginx = x;
}

void MyGl::begin_rotate (float x, float y )
{
 //   cerr <<"begin rotate";
	rotate = true;										
//	ThisRot = LastRot;	
 //   Point2fT MousePt ;
    MousePt.T[0] =x; MousePt.T[1] = y; 						
	ArcBall.click(&MousePt);
}
void MyGl::continue_rotate (float x, float y )
{
  //  cerr <<"continue rotate";
    Quat4fT     ThisQuat;
  //  Point2fT MousePt ;
    MousePt.T[0] =x; MousePt.T[1] = y; 	
    ArcBall.drag(&MousePt, &ThisQuat);						// Update End Vector && Get Rotation As Quaternion
    Matrix3fSetRotationFromQuat4f(&ThisRot, &ThisQuat);		// Convert Quaternion Into Matrix3fT
    Matrix3fMulMatrix3f(&ThisRot, &LastRot);				// Accumulate Last Rotation Into This One
    Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);	// Set Our Final Transform's Rotation From This One

}

void MyGl::select_MW (unsigned int an) {
    ZNMolecule *selected_mol = NULL;
    if (ddwin -> target_molecule -> selection) {
        Selection *sel = (Selection *) ddwin -> target_molecule; 
        if (sel -> get_molecules ().size () == 1) {
            selected_mol = sel -> get_molecules ()[0];
        }

        if (selected_mol) ddwin -> set_current_target (selected_mol);
        else ddwin -> set_current_target (0);
    }
    ddwin -> deselect ();
    ZNMolecule *mol = ddwin -> target_molecule;

        Selection *sel = new Selection ();
		ddwin ->add_molecule (sel);
     //   ddwin->molecules.push_back (sel);
    //    ddwin->target->insertItem (1000, QString(sel -> GetTitle ())); 
//		ddwin-> emit_targets_updated ();   
        ddwin->set_current_target (-1);
		FOR_ATOMS_OF_MOL (a, mol) {
			if (a-> GetAtomicNum () == an) sel -> select_atom (&*a);
		}

        find_limits (sel);
        find_center (sel);

}



void MyGl::select_square () {
    ZNMolecule *selected_mol = NULL;
    if (ddwin -> target_molecule -> selection) {
        Selection *sel = (Selection *) ddwin -> target_molecule; 
        if (sel -> get_molecules ().size () == 1) {
            selected_mol = sel -> get_molecules ()[0];
        }

        if (selected_mol) ddwin -> set_current_target (selected_mol);
        else ddwin -> set_current_target (0);
    }
    ddwin -> deselect ();
    ZNMolecule *mol = ddwin -> target_molecule;
    GLint viewport[4];
    GLuint *buffer = new GLuint[ddwin->target_molecule -> NumAtoms ()*4];
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix ();
    glLoadIdentity ();
    												
												
    matrix_transformations ();
   
    glSelectBuffer (ddwin->target_molecule -> NumAtoms ()*4, buffer);
    glRenderMode(GL_SELECT);
    glInitNames();                            
    glPushName(0);                              
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();                            
    glLoadIdentity() ;

    glGetIntegerv(GL_VIEWPORT, viewport); 
    int x1, x2, y1, y2;
    if (selection_square_x1<selection_square_x2) {
        x1 = selection_square_x1;
        x2 = selection_square_x2;
    }  
    else {
        x2 = selection_square_x1;
        x1 = selection_square_x2;
    }
    if (selection_square_y1<selection_square_y2) {
        y1 = selection_square_y1;
        y2 = selection_square_y2;
    }  
    else {
        y2 = selection_square_y1;
        y1 = selection_square_y2;
    }
    gluPickMatrix((x1+x2)/2, viewport[3]-(y1+y2)/2, x2-x1, y2-y1, viewport);
    gluPerspective(45.0f, (GLfloat) (viewport[2]-viewport[0])/(GLfloat) (viewport[3]-viewport[1]), 0.1f, 10000.0f);
    glMatrixMode(GL_MODELVIEW);

    draw_atoms_for_selection (ddwin->target_molecule);                                             
    glMatrixMode(GL_PROJECTION);                       
    glPopMatrix();                            
    glMatrixMode(GL_MODELVIEW);      
             

    int hits =  glRenderMode (GL_RENDER);


    if (hits){ 
   //     cerr << hits << endl;
        Selection *sel = new Selection ();

        for (unsigned int i=0; i<hits; i++) {
            Atom * at = mol -> GetAtom (buffer[i*4+3]); 
            assert (at);
            sel -> select_atom (at);
        }

   //     FOR_BONDS_OF_MOL(b, mol) {   
 //            if (get_selected (&*b)) sel -> ZNAddBondToSelection (&*b);
  //      }
		
		find_limits (sel);
        find_center (sel);
		ddwin ->add_molecule (sel);   
        ddwin->set_current_target (-1);


    }

    delete[] buffer;
}


Atom * MyGl::select_atom (float x, float y, int dx, int dy) {
    GLint viewport[4];
    GLuint buffer[512];
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix ();
    glLoadIdentity ();
										
	glLoadIdentity();			
									
    matrix_transformations ();

    glSelectBuffer (510, buffer);
    glRenderMode(GL_SELECT);
    glInitNames();                            
    glPushName(0);                              
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();                            
    glLoadIdentity() ;
      //  GLfloat w = viewport[2] / viewport[3];
      //  GLfloat h = 1.0;
      //  glFrustum( -w, w, -h, h, 1.0, 10000.0 );

              
    glGetIntegerv(GL_VIEWPORT, viewport);   
    gluPickMatrix(x, viewport[3]-y, dx, dy, viewport);
	gluPerspective(45.0f, (GLfloat) (viewport[2]-viewport[0])/(GLfloat) (viewport[3]-viewport[1]), 0.1f, 10000.0f);

    glMatrixMode(GL_MODELVIEW);

    draw_atoms_for_selection (ddwin->target_molecule);                                             
    glMatrixMode(GL_PROJECTION);                       
    glPopMatrix();                            
    glMatrixMode(GL_MODELVIEW);      
             

    int hits =  glRenderMode (GL_RENDER);
    if (hits){ 
        Atom *a = ddwin->target_molecule->GetAtom (buffer[3]);
        return a;

    }
    else return NULL;

 //       if len (hits):
 //           print len (hits)        
 //           print hits[0]
 //           selected = self.master.data.atoms[hits [0][2][0]]
 //           print selected.typ, selected.coords
 //     #      print self.docking_sphere_center_coords 
 //           self.ogl.tkRedraw()



}





void MyGl::selectf (float x, float y) {

    Atom *a = select_atom (x, y, 5, 5);
    if (a) {
		if (ddwin ->adding_restrains) {
			ddwin ->haptic_menu ->add_restrain_atom (a);
		}
		else {
			clicked_atom=a;
			ddwin->show_atom_properties (clicked_atom);
		}
    }

}


void MyGl::draw_ring_line (Ring* ring) {
/*    glPushMatrix ();
    vect center = ring->center;

    rotate_to_plane (ring->bonds[0]->GetBeginAtom ()-> GetVector (),ring->bonds[0]->GetEndAtom ()-> GetVector (), center);


    float slices = 40;
//    float rad1 = stick_rad;
//    float rad2 = stick_rad;
    float angl = 0.0;
    float da = 2*PI/slices;


  
    glBegin (GL_LINE_STRIP);

    ZNBond *this_bond = ring->bonds[0];
    Atom *prev_atom = ring->bonds[0]->GetBeginAtom ();
    Atom *next_atom = ring->bonds[0]->GetEndAtom ();
//    cout <<"next atom "<<ring->bonds[0]->GetEndAtom ()->ID<<endl;
    Atom *start_atom = prev_atom;
  //  float dir_coords [3];
 //   for (unsigned int i=0; i<3; i++) dir_coords[i]=next_atom-> GetVector ()[i];
 //   rotate_around_vector (dir_coords, rot_vec, O, rangle);
 //   float tang;
 //   float test_ang;
 //   dis = sqrt (dir_coords[0]*dir_coords[0]+dir_coords[1]*dir_coords[1]+dir_coords[2]*dir_coords[2]);
//    tang = acos (dir_coords[0]/dis);
 //   if (dir_coords[1]<0) test_ang = tang;
//    else test_ang = PI-tang;

 //   if (test_ang<start_ang) {
 //       Atom *sw;
 //       sw = prev_atom;
 //       prev_atom = next_atom;
 //       next_atom = sw;
 //   }


    float rad = dist (ring->bonds[0]->GetBeginAtom ()-> GetVector (), center) * cos (PI/ring->bonds.size ())-aromatic_bond_inter_distance; //apotema - 
    float curr_ang;
    int maxc = angle (prev_atom-> GetVector (), ring->center, next_atom-> GetVector ())/180*PI/da;
    int c = 0;
    for (unsigned int n=0; n<slices+1; n++){
        angl = n*da+PI/2;
    //    cout <<c<<" "<<maxc<<endl;
        if (c>maxc) {
            c=0;
            prev_atom = next_atom;
        //    next_atom = ring->atoms[0];
            for (unsigned int rb=0; rb<ring->bonds.size (); rb++) {
                if (ring->bonds[rb]!=this_bond && ring->bonds[rb]->GetBeginAtom ()==next_atom) {
                    next_atom = ring->bonds[rb]->GetEndAtom ();
                    this_bond = ring->bonds[rb];
            //        cout<<"next atom "<<ring->bonds[rb]->GetEndAtom ()->ID<<endl;
                    break;
                }
                if (ring->bonds[rb]!=this_bond && ring->bonds[rb]->GetEndAtom ()==next_atom) {
                    next_atom = ring->bonds[rb]->GetBeginAtom ();
                    this_bond = ring->bonds[rb];
          //          cout<<"next atom "<<ring->bonds[rb]->GetBeginAtom ()->ID<<endl;
                    break;
                }
            
            }
            maxc = angle (prev_atom-> GetVector (), ring->center, next_atom-> GetVector ())/180*PI/da;

        }
        float r = (c*next_atom->col.redF () + (maxc-c)*prev_atom->col.redF ())/maxc;
        float g = (c*next_atom->col.greenF ()+(maxc-c)*prev_atom->col.greenF ())/maxc;
        float b = (c*next_atom->col.blueF () +(maxc-c)*prev_atom->col.blueF ())/maxc;
        float a = (c*next_atom->col.alphaF () +(maxc-c)*prev_atom->col.alphaF ())/maxc;

        glColor4f (r, g, b, a);
        glVertex3f (sin (angl)*rad, cos (angl)*rad, 0);
        c++;
    }
    glEnd ();
    glPopMatrix ();

*/
}


void MyGl::draw_ring_stick (Ring* ring) {
/*

    glPushMatrix ();
    vect center = ring->center;
 
    rotate_to_plane (ring->bonds[0]->GetBeginAtom ()-> GetVector (), ring->bonds[0]->GetEndAtom ()-> GetVector (), center);


    float slices = 40;
    float slices2 = 20;
//    float rad1 = stick_rad;
//    float rad2 = stick_rad;
    float angl = 0.0;
    float angl2 = 0.0;
    float da = 2*PI/slices;
    float da2 = 2*PI/slices2;

    float radtwo = stick_rad/2;
    float rad = dist (ring->bonds[0]->GetBeginAtom ()-> GetVector (), center) * cos (PI/ring->bonds.size ())-aromatic_bond_inter_distance; //apotema - /apotema - 
    //equation of plane Ax + By + Cz = 0
    rad-=radtwo;

    ZNBond *this_bond = ring->bonds[0];
    Atom *prev_atom = ring->bonds[0]->GetBeginAtom ();
    Atom *next_atom = ring->bonds[0]->GetEndAtom ();
//    cout <<"next atom "<<ring->bonds[0]->GetEndAtom ()->ID<<endl;
    Atom *start_atom = prev_atom;


    float curr_ang;
    int maxc = angle (prev_atom-> GetVector (), ring->center, next_atom-> GetVector ())/180*PI/da;
    int c = 0;




    vect cent (0.f, 0.f, 0.f);
    vect cent2(0.f, 0.f, 0.f);


    for (unsigned int n=0; n<slices; n++){
        angl = n*da+PI/2;    //+start_ang+PI/2;
        glBegin (GL_QUAD_STRIP);


        if (c>maxc) {
            c=0;
            prev_atom = next_atom;
        //    next_atom = ring->atoms[0];
            for (unsigned int rb=0; rb<ring->bonds.size (); rb++) {
                if (ring->bonds[rb]!=this_bond && ring->bonds[rb]->GetBeginAtom ()==next_atom) {
                    next_atom = ring->bonds[rb]->GetEndAtom ();
                    this_bond = ring->bonds[rb];
            //        cout<<"next atom "<<ring->bonds[rb]->GetEndAtom ()->ID<<endl;
                    break;
                }
                if (ring->bonds[rb]!=this_bond && ring->bonds[rb]->GetEndAtom ()==next_atom) {
                    next_atom = ring->bonds[rb]->GetBeginAtom ();
                    this_bond = ring->bonds[rb];
          //          cout<<"next atom "<<ring->bonds[rb]->GetBeginAtom ()->ID<<endl;
                    break;
                }
            
            }
            maxc = angle (prev_atom-> GetVector (), ring->center, next_atom-> GetVector ())/180*PI/da;

        }
        float r = (c*next_atom->col.redF ()+(maxc-c)*prev_atom->col.redF ())/maxc;
        float g = (c*next_atom-> col.greenF ()+(maxc-c)*prev_atom->col.greenF ())/maxc;
        float b = (c*next_atom-> col.blueF ()+(maxc-c)*prev_atom->col.blueF ())/maxc;
        float a = (c*next_atom-> col.alphaF ()+(maxc-c)*prev_atom-> col.alphaF ())/maxc;

        float r2 = ((c+1)*next_atom-> col.redF ()+(maxc-c-1)*prev_atom-> col.redF ())/maxc;
        float g2 = ((c+1)*next_atom-> col.greenF ()+(maxc-c-1)*prev_atom-> col.greenF ())/maxc;
        float b2 = ((c+1)*next_atom-> col.blueF ()+(maxc-c-1)*prev_atom-> col.blueF ())/maxc;
        float a2 = ((c+1)*next_atom-> col.alphaF ()+(maxc-c-1)*prev_atom-> col.alphaF ())/maxc;








        for (unsigned int n2=0; n2<slices2; n2++) {
            angl2 = n2*da2;
            vect v, v2;
            v.x() = sin (angl)*(rad+radtwo);
            v.y() = cos (angl)*(rad+radtwo);
            v.z() = 0.f;
            v2.x() = sin (angl+da)*(rad+radtwo);
            v2.y() = cos (angl+da)*(rad+radtwo);
            v2.z() = 0.f;
            cent.x() = sin (angl)*(rad);
            cent.y() = cos (angl)*(rad);
            cent.z() = 0.f;
            cent2.x() = sin (angl+da)*(rad);
            cent2.y() = cos (angl+da)*(rad);
            cent2.z() = 0.f;
            vect tang;
            vect tang2;
            tang.x() = (cos (angl));
            tang.y() = (-sin (angl));
            tang.z() = 0.f;
            tang2.x() = (cos (angl+da));
            tang2.y() = (-sin (angl+da));
            tang2.z() = 0.f;


      //      glVertex3fv (v);
       //     glVertex3fv (v2);
            rotate_around_vector (v, tang, cent, angl2); 
            rotate_around_vector (v2, tang2, cent2, angl2);
            glNormal3f (-(v.x() - cent.x()), -(v.y()-cent.y()),-(v.z()-cent.z()));
            glColor4f (r, g, b, a);
            glVertex3f (v.x(), v.y(), v.z());
            glNormal3f (-(v2.x()-cent2.x()), -(v2.y()-cent2.y()),-(v2.z()-cent2.z()));
            glColor4f (r2, g2, b2, a2);
            glVertex3f (v2.x(), v2.y(), v2.z());
            
        }
        glEnd ();  
        c++;      
    }

    glPopMatrix ();

*/
}






void MyGl::draw_bond_line (ZNBond* bond)
{
    glBegin(GL_LINES);
	Atom *at1 = bond->GetBeginAtom ();
	Atom *at2 = bond ->GetEndAtom ();
	vect v1 = get_coordinates (at1);
	vect v2 = get_coordinates (at2);
    setAtomColor(at1);
    glVertex3d(v1.x(), v1.y(), v1.z());
    setAtomColor(at2);
    glVertex3d(v2.x(), v2.y(), v2.z());
    glEnd();
}

void MyGl::draw_aromatic_bond_line (ZNBond* bond) {
    float d=0.08f;
    glBegin(GL_LINES);
	Atom *at1 = bond->GetBeginAtom ();
	Atom *at2 = bond ->GetEndAtom ();
	vect v1 = get_coordinates (at1);
	vect v2 = get_coordinates (at2);
    setAtomColor(at1);
    glVertex3d(v1.x(), v1.y(), v1.z());
    setAtomColor(at2);
    glVertex3d(v2.x(), v2.y(), v2.z());
    glEnd();
	
	
/*    if (bond->hasCoplanar0) {
        Atom *at0 = bond-> GetBeginAtom ();
        Atom *co0 =bond->coplanarAtom0;
         Atom *at1 = bond-> GetEndAtom ();   
        float vecx = co0-> GetVector ().x()-at0-> GetVector ().x();
        float vecy = co0-> GetVector ().y()-at0-> GetVector ().y();
        float vecz = co0-> GetVector ().z()-at0-> GetVector ().z();
        float mod = sqrt (vecx*vecx+vecy*vecy+vecz*vecz);
        vecx/=mod;
        vecy/=mod;
        vecz/=mod;
        glEnable(GL_LINE_STIPPLE);        
        glLineStipple(1, 0xF0F0); 
        glBegin(GL_LINES);
        setAtomColor(at0);
        glVertex3f(at0-> GetVector ().x()+vecx*d,at0-> GetVector ().y()+vecy*d,at0-> GetVector ().z()+vecz*d );
        setAtomColor(at1);
        glVertex3f(at1-> GetVector ().x()+vecx*d,at1-> GetVector ().y()+vecy*d,at1-> GetVector ().z()+vecz*d );
        glEnd();
        glDisable (GL_LINE_STIPPLE);
    }

    else if (bond->hasCoplanar1) {
        Atom *at0 = bond-> GetBeginAtom ();
        Atom *co1 =bond->coplanarAtom1;
        Atom *at1 = bond-> GetEndAtom ();   
        float vecx = co1-> GetVector ().x()-at1-> GetVector ().x(); 
        float vecy = co1-> GetVector ().y()-at1-> GetVector ().y();
        float vecz = co1-> GetVector ().z()-at1-> GetVector ().z();
        float mod = sqrt (vecx*vecx+vecy*vecy+vecz*vecz);
        vecx/=mod;
        vecy/=mod;
        vecz/=mod;
        glEnable(GL_LINE_STIPPLE);         
        glLineStipple(1, 0xF0F0); 
        glBegin(GL_LINES);
        setAtomColor(at0);
        glVertex3f(at0-> GetVector ().x()+vecx*d,at0-> GetVector ().y()+vecy*d,at0-> GetVector ().z()+vecz*d );
        setAtomColor(at1);
        glVertex3f(at1-> GetVector ().x()+vecx*d,at1-> GetVector ().y()+vecy*d,at1-> GetVector ().z()+vecz*d );
        glEnd();
        glDisable (GL_LINE_STIPPLE);
    }*/
}

void MyGl::draw_triple_bond_line (ZNBond * bond) {
    draw_bond_line (bond);
    draw_double_bond_line (bond);
}



void MyGl::draw_double_bond_line (ZNBond* bond) {
    Atom *at0 = bond->GetBeginAtom ();
    Atom *at1 = bond->GetEndAtom ();

    float verts  [4][3];
    compute_double_bond_vertexes (bond, verts);
    glBegin(GL_LINES);
    setAtomColor(at0);
    glVertex3f(verts[0][0], verts[0][1], verts [0][2]);
    setAtomColor(at1);
    glVertex3f(verts[2][0], verts[2][1], verts[2][2]);
    glEnd();

    glBegin(GL_LINES);
    setAtomColor(at0);
    glVertex3f(verts[1][0], verts[1][1], verts[1][2]);
    setAtomColor(at1);
    glVertex3f(verts[3][0], verts[3][1], verts[3][2]);
    glEnd();

}



//    if (found0 && !found1) {
        
//    }



/*
    Atom *at0, *co1, *at1;
    if (bond->hasCoplanar1) {
        at0 = bond-> GetBeginAtom ();
        co1 =bond->coplanarAtom1;
         at1 = bond-> GetEndAtom ();
    }

    else  {
        at1 = bond-> GetBeginAtom ();
        co1 =bond->coplanarAtom0;
        at0 = bond-> GetEndAtom ();  }

  
    float vecx = co1-> GetVector ().x()-at1-> GetVector ().x(); 
    float vecy = co1-> GetVector ().y()-at1-> GetVector ().y();
    float vecz = co1-> GetVector ().z()-at1-> GetVector ().z();
    float mod = sqrt (vecx*vecx+vecy*vecy+vecz*vecz);
    vecx/=mod;
    vecy/=mod;
    vecz/=mod;
    vector <float> vec;
    vector <float> ref;
    vector <float> nor;
    vector <float> par;
    ref.push_back(at1-> GetVector ().x()-at0-> GetVector ().x());
    ref.push_back(at1-> GetVector ().y()-at0-> GetVector ().y());
    ref.push_back(at1-> GetVector ().z()-at0-> GetVector ().z());
    vec.push_back(vecx); vec.push_back(vecy); vec.push_back(vecz);  
    components (vec, ref, par, nor);

    float a1x=at0-> GetVector ().x()+nor[0]*d-par[0]*d;
    float a1y=at0-> GetVector ().y()+nor[1]*d-par[1]*d;
    float a1z=at0-> GetVector ().z()+nor[2]*d-par[2]*d;

    float a2x=at1-> GetVector ().x()+nor[0]*d+par[0]*d;
    float a2y=at1-> GetVector ().y()+nor[1]*d+par[1]*d;
    float a2z=at1-> GetVector ().z()+nor[2]*d+par[2]*d;



        glBegin(GL_LINES);
        setAtomColor(at0);
        glVertex3f(a1x, a1y, a1z);
        setAtomColor(at1);
        glVertex3f(a2x, a2y, a2z);
        glEnd();

    a1x=at0-> GetVector ().x()-nor[0]*d-par[0]*d;
    a1y=at0-> GetVector ().y()-nor[1]*d-par[1]*d;
    a1z=at0-> GetVector ().z()-nor[2]*d-par[2]*d;

    a2x=at1-> GetVector ().x()-nor[0]*d+par[0]*d;
    a2y=at1-> GetVector ().y()-nor[1]*d+par[1]*d;
    a2z=at1-> GetVector ().z()-nor[2]*d+par[2]*d;

        glBegin(GL_LINES);
        setAtomColor(at0);
        glVertex3f(a1x, a1y, a1z);
        setAtomColor(at1);
        glVertex3f(a2x, a2y, a2z);
        glEnd();
    }

*/



void MyGl::draw_bond_stick(ZNBond* bond)
{
 
    float angle, height, mod_of_vector;
    Atom *at1 = bond -> GetBeginAtom ();
    Atom *at2 = bond -> GetEndAtom ();
    vect at1c = get_coordinates(at1);
    vect at2c = get_coordinates (at2);



    vect vector = subtract (at2c, at1c);



    mod_of_vector = vector.module ();

    vector *= 1/mod_of_vector;

    height = mod_of_vector;

    angle = acos (vector.z ()) * 180.0 / PI;
    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  
//    setAtomColor(bond->GetBeginAtom ());
    glPushMatrix ();



    glTranslated (at1c.x (), at1c.y(), at1c.z());
    glRotated (180.0, 0.0, 0.0, 1.0);
    glRotated (angle, vector.y (), -1.0 * vector.x (), 0.0);
  


//    gluCylinder (quadratic,stick_rad,stick_rad,mod_of_vector*0.5,16, 1);

    my_cylinder (*ddwin ->data ->stick_radius, *ddwin ->data ->stick_radius, height, *ddwin ->data ->quality_scale * 5, bond->GetBeginAtom (), bond->GetEndAtom ());

    glPopMatrix ();
    
    Atom *at1p = bond -> GetBeginAtom ();
    Atom *at2p = bond -> GetEndAtom ();


    if (!get_sad (at1p)) {
        set_sad (at1p, true);
        glPushMatrix ();
        setAtomColor(at1);
        glTranslatef (at1c.x(), at1c.y (), at1c.z ());
        gluSphere (quadratic, *ddwin ->data ->stick_radius, *ddwin ->data ->quality_scale * 5, *ddwin ->data ->quality_scale * 5);
        glPopMatrix ();
    }
    if (!get_sad (at2p)) {
        set_sad (at2p, true);
        glPushMatrix ();
        setAtomColor(at2);
        glTranslatef (at2c.x(), at2c.y (), at2c.z ());
        gluSphere (quadratic, *ddwin ->data ->stick_radius, *ddwin ->data ->quality_scale * 5, *ddwin ->data ->quality_scale * 5);
        glPopMatrix ();
    }



  //  setAtomColor(bond->GetEndAtom ());
  //  glPushMatrix ();

//    glTranslatef (bond->GetBeginAtom ()-> GetVector ().x() + 0.5 * (bond->GetEndAtom ()-> GetVector ().x() - bond->GetBeginAtom ()-> GetVector ().x()),
  //				 bond->GetBeginAtom ()-> GetVector ().y()  + 0.5 * (bond->GetEndAtom ()-> GetVector ().y() - bond->GetBeginAtom ()-> GetVector ().y()), 
	//			 bond->GetBeginAtom ()-> GetVector ().z()   + 0.5 * (bond->GetEndAtom ()-> GetVector ().z() - bond->GetBeginAtom ()-> GetVector ().z()));
  //  glRotatef (180.0, 0.0, 0.0, 1.0);
  //  glRotatef (angle, vector[1], -1.0 * vector[0], 0.0);
  

  //  gluCylinder (quadratic,stick_rad,stick_rad,0.5 * mod_of_vector,10, 1);

  //  glPopMatrix ();


}


void MyGl::draw_double_bond_stick (ZNBond* bond) {



    float verts  [4][3];

    float angle, height;
    compute_double_bond_vertexes (bond, verts);
    vect vector (verts[2][0]-verts[0][0], verts[2][1]-verts[0][1], verts[2][2]-verts[0][2]);
    height = vector.module ();
    vector.normalise ();
  

    angle = acos (vector.z()) * 180.0f / PI;
    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  

    glPushMatrix ();


    glTranslatef (verts[0][0], verts[0][1], verts[0][2]);
    glRotatef (180.0f, 0.0f, 0.0f, 1.0f);
    glRotatef (angle, vector.y(), -1.0 * vector.x(), 0.0f);
  
//    gluCylinder (quadratic,stick_rad,stick_rad,mod_of_vector*0.5,16, 1);

    my_cylinder (*ddwin ->data ->stick_radius * *ddwin ->data ->double_bond_stick_scale, *ddwin ->data ->stick_radius * *ddwin ->data ->double_bond_stick_scale, height, *ddwin ->data ->quality_scale * 5, bond->GetBeginAtom (),bond->GetEndAtom ());

    glPopMatrix ();

    vector = vect (verts[3][0]-verts[1][0], verts[3][1]-verts[1][1], verts[3][2]-verts[1][2]);
    height = vector.module ();
    vector.normalise ();


    angle = acos (vector.z()) * 180.0f / PI;
    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  

    glPushMatrix ();


    glTranslatef (verts[1][0], verts[1][1], verts[1][2]);
    glRotatef (180.0, 0.0, 0.0, 1.0);
    glRotatef (angle, vector.y(), -1.0 * vector.x(), 0.0);
  
//    gluCylinder (quadratic,stick_rad,stick_rad,mod_of_vector*0.5,16, 1);

    my_cylinder (*ddwin ->data ->stick_radius * *ddwin ->data ->double_bond_stick_scale, *ddwin ->data ->stick_radius * *ddwin ->data ->double_bond_stick_scale, height, *ddwin ->data ->quality_scale * 5, bond->GetBeginAtom (),bond->GetEndAtom ());

    glPopMatrix ();
  //  ZNMolecule *mol  = (ZNMolecule *) bond -> GetParent ();
    if (CountBonds (bond->GetBeginAtom ())== 1) {
        glPushMatrix ();
        setAtomColor(bond->GetBeginAtom ());
        glTranslatef (verts[0][0], verts[0][1], verts[0][2]);
        gluSphere (quadratic, *ddwin ->data ->stick_radius * *ddwin ->data ->double_bond_stick_scale, *ddwin ->data ->quality_scale * 5, *ddwin ->data ->quality_scale * 5);
        glPopMatrix ();
        glPushMatrix ();
        glTranslatef (verts[1][0], verts[1][1], verts[1][2]);
        gluSphere (quadratic, *ddwin ->data ->stick_radius * *ddwin ->data ->double_bond_stick_scale, *ddwin ->data ->quality_scale * 5, *ddwin ->data ->quality_scale * 5);
        glPopMatrix ();
    }
    if (CountBonds (bond->GetEndAtom ())== 1) {
        glPushMatrix ();
        setAtomColor(bond->GetEndAtom ());
        glTranslatef (verts[2][0], verts[2][1], verts[2][2]);
        gluSphere (quadratic, *ddwin ->data ->stick_radius * *ddwin ->data ->double_bond_stick_scale, *ddwin ->data ->quality_scale * 5, *ddwin ->data ->quality_scale * 5);
        glPopMatrix ();
        glPushMatrix ();
        glTranslatef (verts[3][0], verts[3][1], verts[3][2]);
        gluSphere (quadratic, *ddwin ->data ->stick_radius * *ddwin ->data ->double_bond_stick_scale, *ddwin ->data ->quality_scale * 5, *ddwin ->data ->quality_scale * 5);
        glPopMatrix ();
    }
        
    

}

void MyGl::draw_triple_bond_stick (ZNBond* bond) {
    draw_bond_stick (bond);
    draw_double_bond_stick (bond);

}


void MyGl::draw_aromatic_bond_stick (ZNBond* bond) {
    draw_bond_stick (bond); /*
    float angle, height, mod_of_vector;
    float vecto[3];
    float d=aromatic_bond_inter_distance;
    Atom *at0, *co1, *at1;
    if (bond->hasCoplanar1) {
        at0 = bond-> GetBeginAtom ();
        co1 =bond->coplanarAtom1;
         at1 = bond-> GetEndAtom ();
    }

    else  {
        at1 = bond-> GetBeginAtom ();
        co1 =bond->coplanarAtom0;
        at0 = bond-> GetEndAtom ();  }
  
        float vecx = co1-> GetVector ().x()-at1-> GetVector ().x(); 
        float vecy = co1-> GetVector ().y()-at1-> GetVector ().y();
        float vecz = co1-> GetVector ().z()-at1-> GetVector ().z();
        float mod = sqrt (vecx*vecx+vecy*vecy+vecz*vecz);
        vecx/=mod;
        vecy/=mod;
        vecz/=mod;
        vector <float> vec;
        vector <float> ref;
        vector <float> nor;
        vector <float> par;

        ref.push_back(at0-> GetVector ().x()-at1-> GetVector ().x());
        ref.push_back(at0-> GetVector ().y()-at1-> GetVector ().y());
        ref.push_back(at0-> GetVector ().z()-at1-> GetVector ().z());
        vec.push_back(vecx); vec.push_back(vecy); vec.push_back(vecz);  
        components (vec, ref, par, nor);

    float a1x=at0-> GetVector ().x()+nor[0]*d+par[0]*d;
    float a1y=at0-> GetVector ().y()+nor[1]*d+par[1]*d;
    float a1z=at0-> GetVector ().z()+nor[2]*d+par[2]*d;

    float a2x=at1-> GetVector ().x()+nor[0]*d-par[0]*d;
    float a2y=at1-> GetVector ().y()+nor[1]*d-par[1]*d;
    float a2z=at1-> GetVector ().z()+nor[2]*d-par[2]*d;

    vecto[0] = a2x-a1x;
    vecto[1] = a2y-a1y;
    vecto[2] = a2z-a1z;
    mod_of_vector = sqrt(vecto[0]* vecto[0] + vecto[1]* vecto[1] + vecto[2]* vecto[2]);

    vecto[0] /= mod_of_vector;
    vecto[1] /= mod_of_vector;
    vecto[2] /= mod_of_vector;
  


    height = mod_of_vector;

    angle = acos (vecto[2]) * 180.0 / PI;
    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  


    glPushMatrix ();

    glTranslatef (a1x, a1y, a1z);
    glRotatef (180.0, 0.0, 0.0, 1.0);
    glRotatef (angle, vecto[1], -1.0 * vecto[0], 0.0);


    my_cylinder (stick_rad*double_bond_stick_radius_scale*0.5, stick_rad*double_bond_stick_radius_scale*0.5, mod_of_vector, stick_precision, at0,at1);
    glPopMatrix ();
*/

}



void MyGl::draw_atom_sphere(Atom* atom) {

    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	vect v = get_coordinates(atom);
    setAtomColor(atom);
    glPushMatrix ();
    glTranslatef (v.x(), v.y(), v.z());
  //  my_sphere (sphere_radius,2*sphere_precision,sphere_precision, atom);
    gluSphere (quadratic,*ddwin ->data ->sphere_radius,*ddwin ->data ->quality_scale * 5,*ddwin ->data ->quality_scale * 5);
    glPopMatrix ();

}

void MyGl::draw_atom_sel_sphere(Atom* atom) {
	vect v = get_coordinates(atom);
 //   glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  
 //   setAtomColor(atom);
    glPushMatrix ();
    glTranslatef (v.x(), v.y(), v.z());
    gluSphere (quadratic,*ddwin ->data ->sphere_radius,6,6);
    glPopMatrix ();

}


void MyGl::draw_atom_vdw_sphere(Atom* atom) {
	vect v = get_coordinates(atom);
    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  
    setAtomColor(atom);
    glPushMatrix ();
    glTranslatef (v.x(), v.y(), v.z());
    double vdw = etab.GetVdwRad (atom -> GetAtomicNum ());
    gluSphere (quadratic, vdw, *ddwin ->data ->quality_scale * 7, *ddwin ->data ->quality_scale * 7);
    glPopMatrix ();

}

void MyGl::draw_atom_scaled_vdw_sphere(Atom* atom) {
	vect v = get_coordinates(atom);
    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  
    setAtomColor(atom);
    glPushMatrix ();
    glTranslatef (v.x(), v.y(), v.z());
   // my_sphere (atom->vdw*vdw_scale,sphere_precision,sphere_precision, atom);
    double vdw = etab.GetVdwRad (atom -> GetAtomicNum ());
    gluSphere (quadratic, vdw**ddwin ->data ->vdw_scale,*ddwin ->data ->quality_scale * 5,*ddwin ->data ->quality_scale * 5);
    glPopMatrix ();

}

void MyGl::draw_molecule (ZNMolecule *mol) {
	set_needs_redraw(mol, true);
}

void MyGl::GL_update_molecule (ZNMolecule* mol) {
	lock_geometry_for_read (mol);
	redraw_counter++;
    if (mol->selection) {
        draw_list ( (Selection *) mol);
    }
    else draw_list (mol);
	set_needs_redraw(mol, false);
	unlock_geometry (mol);
}

void MyGl::GL_update_backbone (ZNMolecule* mol) {
	lock_geometry_for_read (mol);
	draw_backbone_list(mol);
	set_needs_backbone_redraw(mol, false);
	unlock_geometry (mol);
}

void MyGl::draw_list (MarchingCubes * cube) {

    float xm = cube->size_x ();
    float ym = cube->size_y ();
    float zm = cube->size_z ();    
    float xmin = cube->xmin ();
    float ymin = cube->ymin ();
    float zmin = cube->zmin ();
    float xmax = cube->xmax ();
    float ymax = cube->ymax ();
    float zmax = cube->zmax ();
//    cout <<xmin<<" "<<ymin<<" "<<zmin<<" "<<xmax<<" "<<" "<<ymax<<" "<<zmax<<endl;

    for (unsigned int i=0; i<cube->ntrigs (); i++) {
        glColor4f (1.,1.,1., 1.);
        glBegin (GL_TRIANGLES);
        Vertex *v1 = cube->vert(cube->trig (i)->v1);
        Vertex *v2 = cube->vert(cube->trig (i)->v2);
        Vertex *v3 = cube->vert(cube->trig (i)->v3);
        glNormal3f (v1->nx, v1->ny, v1->nz);
        glVertex3f (v1->x*(xmax-xmin)/xm+xmin, v1->y*(ymax-ymin)/ym+ymin, v1->z*(zmax-zmin)/zm+zmin);
        glNormal3f (v2->nx, v2->ny, v2->nz);
        glVertex3f (v2->x*(xmax-xmin)/xm+xmin, v2->y*(ymax-ymin)/ym+ymin, v2->z*(zmax-zmin)/zm+zmin);
        glNormal3f (v3->nx, v3->ny, v3->nz);
        glVertex3f (v3->x*(xmax-xmin)/xm+xmin, v3->y*(ymax-ymin)/ym+ymin, v3->z*(zmax-zmin)/zm+zmin);
        glEnd ();
    }
}
/*
void MyGl::draw_list (Grid& grid, int list) {

    glNewList(list,GL_COMPILE);
 //   glBegin (GL_POINTS);
    float min, max;
    min =grid.min;
    max=grid.max;
    glPointSize (20.0);
        glBegin (GL_POINTS);

    for (unsigned i =0; i<grid.vec.size (); i++) {
        float col = (grid.vec[i].val-min)/(max-min);
        glColor4f (2*col-1,0,1-(2*col-1),2*col-1);//col*col*col*0.6);
        if (col >0.5){
        glVertex3f (grid.vec[i].x(), grid.vec[i].y(),grid.vec[i].z());
        }
    }
    glEnd ();

    glEndList ();

    }
*/


void MyGl::draw_list (Selection* sel) {
    find_limits (sel);
    find_center (sel);
    for (unsigned int i=0; i < sel->get_molecules ().size (); i++) {
        draw_list (sel->get_molecules ()[i]);
    }
}


void MyGl::backbone_to_surface (ZNMolecule *mol, Surface *surf) {
	find_backbone_data (mol);
	FOR_RESIDUES_OF_MOL (res, mol) {
		draw_backbone_stick (&*res, surf);
	}
	surf ->render();
}


void MyGl::draw_backbone_list (ZNMolecule *mol) {
	glNewList(get_backbone_list1 (mol),GL_COMPILE);
	FOR_RESIDUES_OF_MOL (r, mol) {
		if (get_backbone_display_style (mol) == 1) {
			find_backbone_data (mol);
			draw_backbone_line (&*r);
		}
	}
	glEndList ();
	glNewList(get_backbone_list2 (mol),GL_COMPILE);
	FOR_RESIDUES_OF_MOL (r, mol) {

		if (get_backbone_display_style (mol) == 2) 		{
			find_backbone_data (mol);
			draw_backbone_stick(&*r);
		}
	}
	glEndList ();
}

void MyGl::draw_list (ZNMolecule* mol) {

    find_limits (mol);

    find_center (mol);


    FOR_ATOMS_OF_MOL(a, mol) {
        set_sad (&*a, false);
    }
    glNewList(get_line_list (mol),GL_COMPILE);

    FOR_BONDS_OF_MOL (b, mol) {
        Atom *at1 = b -> GetBeginAtom ();
        Atom *at2 = b -> GetEndAtom ();

    bool vis = get_visible (&*b);
    int dsi = get_ds (&*b);


    if (vis) {

        int order;
        if (aromatic_display_style == AROMATIC_RINGS && b -> IsAromatic ()) order = 1;
        else order = b -> GetBondOrder ();

  //      cout << order<<" "<<mol->bonds[i]->kekule<<" "<<mol->bonds[i]->mol2Type<< endl;
      	switch (dsi) {
            case LINES:  
                switch (order) {
                case 1:
                    draw_bond_line (&*b);
                    break;
                case 2:
                    draw_double_bond_line (&*b);
                    break;
                case 3:
                    draw_triple_bond_line (&*b);
                    break;
  //              case 4:
   //                 draw_aromatic_bond_line (&*b);
    //                break;
    //            case 5:
    //                if (aromatic_display_style != AROMATIC_RINGS)  draw_aromatic_bond_line (&*b);
    //                else draw_bond_line (&*b);
    //                break;
            }
            break;
        }
    }
    }

  //  vector<OBRing*>::iterator i;
  //  vector<OBRing*> *rlist = (vector<OBRing*>*)mol -> GetData("RingList");
/*
    for (i = rlist->begin();i != rlist->end();++i)
    {

        //// set the center somehow

        if ((*i) -> IsAromatic () && aromatic_display_style == AROMATIC_RING) {
            bool lines = false;
            for (unsigned int ra=0; ra < (*i) -> bonds.size (); ra++) {
                if ((*i) -> bonds[ra] -> displayStyle == LINES) {
                    lines = true;
                    break;
                }
            }
            if (lines) draw_ring_line (*i);
        }
    }
    glEndList ();
*/

    glEndList ();
    glNewList(get_stick_list (mol),GL_COMPILE);


    FOR_BONDS_OF_MOL (b, mol) {        
        Atom *at1 = b -> GetBeginAtom ();
        Atom *at2 = b -> GetEndAtom ();

        bool vis = get_visible (&*b);
        int dsi = get_ds (&*b);
        if (vis) {
          	switch (dsi) {
                case STICKS:
                    int order;
                    if (aromatic_display_style == AROMATIC_RINGS && b -> IsAromatic ()) order = 1;
                    else order = b -> GetBondOrder ();

                    if (order == 1) draw_bond_stick (&*b);                    
                    else if (order == 2) draw_double_bond_stick (&*b);
                    else if (order == 3) draw_triple_bond_stick (&*b);
     //               else if (order == 4) draw_aromatic_bond_stick (&*b);
     //               else if (order == 5) {
     //                   if (aromatic_display_style != AROMATIC_RINGS)  draw_aromatic_bond_stick (&*b);
     //                   else draw_bond_stick (&*b);
            //        }
                break;

            }
        }
    }

    FOR_ATOMS_OF_MOL(a, mol) {
    bool v = get_visible (&*a);
    int dsi = get_ds (&*a);


    if  (v) {  
      	switch (dsi) {
            case  NO_ATOMS:
                break;
            case SPHERES:
                draw_atom_sphere (&*a);
                break;
            case CPK_SPHERES:
                draw_atom_vdw_sphere (&*a);
                break;
            case SCALED_CPK_SPHERES:
                draw_atom_scaled_vdw_sphere (&*a);
                break;
        }
    }

    }

 /*   rlist = (vector<OBRing*>*)mol -> GetData("RingList");
    for (i = rlist->begin();i != rlist->end();++i)
    {
        bool stick = false;

        if ((*i) -> IsAromatic ()) {
            for (unsigned int ra=0; ra< (*i) -> bonds.size (); ra++) {
                if ((*i) -> bonds[ra] -> displayStyle==STICKS) {
                    stick = true;
                    break;
                }
            }
            if (stick && aromatic_display_style==AROMATIC_RINGS) draw_ring_stick (*i);
        }
    }
*/
    glEndList ();

	draw_backbone_list (mol);
}



void MyGl::draw_atoms_for_selection (ZNMolecule *mol){    //only works in rendermode select
  //  glNewList(27,GL_COMPILE);
        int i = 0;
        FOR_ATOMS_OF_MOL(a, mol) {
            i = a -> GetIdx ();    
            glLoadName(i);  
            draw_atom_sel_sphere (&*a); 
        }                 
 //   glEndList ();
}


void MyGl::openGLSetColor (color col) {
    float r, g, b, a;
	r= col.redF ();
	g= col.greenF ();
	b= col.blueF ();
	a= col.alphaF (); 
	glColor4f (r, g, b, a);
}


void MyGl::setAtomColor(Atom* atom)
{

    float r, g, b, a;
    color col = get_color (atom);
    bool sel = get_selected (atom);
    if (sel) {
        float scr, scg, scb, sca;
        scr = select_color.redF ();
        scg = select_color.greenF ();
        scb = select_color.blueF ();
        sca = select_color.alphaF () * (0.5 + 0.5*sin ((float) select_pulse*3/40));

        r = scr * sca + col.redF ()   * ( 1 - sca);
        g = scg * sca + col.greenF () * ( 1 - sca);
        b = scb * sca + col.blueF ()  * ( 1 - sca);
        a = col.alphaF ();
    }
    else {
        r= col.redF ();
        g= col.greenF ();
        b= col.blueF ();
        a= col.alphaF ();    
    }
    glColor4f (r, g, b, a);

}
/*
void MyGl::update_current_color () {
}
*/



void MyGl::hide_hydrogens (ZNMolecule* mol){

    FOR_ATOMS_OF_MOL(a, mol) {
        if (a -> IsHydrogen ()) {
            set_visible (&*a, false);
        }
    }
    draw_molecule (mol);
}


void MyGl::hide_hydrogens (vector <ZNMolecule *> molecules) {
    for (unsigned int n =0; n<molecules.size (); n++ ) hide_hydrogens (molecules[n]);
}

void MyGl::hide_nonpolar_hydrogens (ZNMolecule* mol){
    FOR_ATOMS_OF_MOL(a, mol) {
        if (a -> IsNonPolarHydrogen ()) {
            set_visible (&*a, false);
        }
    }
    draw_molecule (mol);
   
}


void MyGl::hide_nonpolar_hydrogens (vector <ZNMolecule *> molecules) {
    for (unsigned int n =0; n<molecules.size (); n++ ) hide_nonpolar_hydrogens (molecules[n]);
}


void MyGl::show_all_atoms (ZNMolecule* mol){
    FOR_ATOMS_OF_MOL(a, mol) {
            set_visible (&*a, true);
    }
    draw_molecule (mol);
}

void MyGl::show_all_atoms (vector <ZNMolecule *> molecules) {
    for (unsigned int n =0; n<molecules.size (); n++ ) show_all_atoms (molecules[n]);
}

void MyGl::hide_all_atoms (ZNMolecule* mol){
    FOR_ATOMS_OF_MOL(a, mol) {
            set_visible (&*a, false);
    }
    draw_molecule (mol);
}

void MyGl::hide_all_atoms (vector <ZNMolecule *> molecules) {
    for (unsigned int n =0; n<molecules.size (); n++ ) hide_all_atoms (molecules[n]);
}


void MyGl::apply_color_masks (vector <color_mask> masks, ZNMolecule *mol, bool undoable) {
    float c_red, c_green, c_blue, c_alpha,  ce_red, ce_green, ce_blue, ce_alpha, cb_red, cb_green, cb_blue, cb_alpha, sb_red, sb_green, sb_blue, sb_alpha, sm_red, sm_green, sm_blue, sm_alpha, se_red, se_green, se_blue, se_alpha, mid_width;

    mid_width = 0.1f;

    c_red   = ddwin -> data -> constant_color.redF ()  ;
    c_green = ddwin -> data -> constant_color.greenF ();
    c_blue  = ddwin -> data -> constant_color.blueF () ;
    c_alpha = ddwin -> data -> constant_color.alphaF () ;

    cb_red   = ddwin -> data -> charge_begin_color.redF ()  ;
    cb_green = ddwin -> data -> charge_begin_color.greenF ();
    cb_blue  = ddwin -> data -> charge_begin_color.blueF () ;
    cb_alpha  = ddwin -> data -> charge_begin_color.alphaF () ;

    ce_red   = ddwin -> data -> charge_end_color.redF ()  ;
    ce_green = ddwin -> data -> charge_end_color.greenF ();
    ce_blue  = ddwin -> data -> charge_end_color.blueF () ;
    ce_alpha  = ddwin -> data -> charge_end_color.alphaF () ;

    sb_red   = ddwin -> data -> score_begin_color.redF () ;
    sb_green = ddwin -> data -> score_begin_color.greenF ();
    sb_blue  = ddwin -> data -> score_begin_color.blueF () ;
    sb_alpha  = ddwin -> data -> score_begin_color.alphaF () ;


    sm_red   = ddwin -> data -> score_mid_color.redF ()  ;
    sm_green = ddwin -> data -> score_mid_color.greenF ();
    sm_blue  = ddwin -> data -> score_mid_color.blueF () ;
    sm_alpha  = ddwin -> data -> score_mid_color.alphaF () ;  

    se_red   = ddwin -> data -> score_end_color.redF ()  ;
    se_green = ddwin -> data -> score_end_color.greenF ();
    se_blue  = ddwin -> data -> score_end_color.blueF () ;
    se_alpha  = ddwin -> data -> score_end_color.alphaF () ;

    ColorAtomCommand *color_atom = new ColorAtomCommand (this);

    FOR_ATOMS_OF_MOL(at, mol) {
        float newr = 0.f , newg = 0.f , newb = 0.f, newa = 0.f;

        for (unsigned int ms=0; ms < masks.size (); ms++) {
            if (masks[ms].type == ELEMENT) {
                float r=0.f, g=0.f, b=0.f, a = 0.f;
                color col = mol -> get_color_mw (&*at);
                r = col.redF ();
                g = col.greenF ();
                b = col.blueF ();
                a = col.alphaF ();
                r *= masks[ms].intensity;
                g *= masks[ms].intensity;
                b *= masks[ms].intensity;
                a *= masks[ms].intensity;
                assert (!(r > 1.f));
                assert (!(g > 1.f));
                assert (!(b > 1.f));
                assert (!(a > 1.f));

                newr += r;
                newg += g;
                newb += b;
                newa += a;
            }

            else if (masks[ms].type == CHARGE) {
                float r, g, b, a;
                float q = at -> GetPartialCharge ();
                if (q <= ddwin -> data -> charge_begin) {
                    r = cb_red   * masks[ms].intensity;
                    g = cb_green * masks[ms].intensity;
                    b = cb_blue  * masks[ms].intensity;
                    a = cb_alpha  * masks[ms].intensity;
                }
                else if (q >= ddwin -> data -> charge_end) {
                    r = ce_red   * masks[ms].intensity;
                    g = ce_green * masks[ms].intensity;
                    b = ce_blue  * masks[ms].intensity;
                    a = ce_alpha  * masks[ms].intensity;
                }
                else {
                    float perc = (ddwin -> data -> charge_end - q) / (ddwin -> data -> charge_end - ddwin -> data -> charge_begin);
                    r = (cb_red  * perc + ce_red   * (1-perc)) * masks[ms].intensity;
                    g = (cb_green* perc + ce_green * (1-perc)) * masks[ms].intensity;
                    b = (cb_blue * perc + ce_blue  * (1-perc)) * masks[ms].intensity;
                    a = (cb_alpha * perc + ce_alpha  * (1-perc)) * masks[ms].intensity;
                    
                }
                assert (!(r > 1.f));
                assert (!(g > 1.f));
                assert (!(b > 1.f));
                assert (!(a > 1.f));
                newr += r;
                newg += g;
                newb += b;
                newa += a;
    
            }

            else if (masks[ms].type == SCORE) {
                float r, g, b, a;
                float q = get_score (&*at);
				
				
				color col = average_3_colors (q, ddwin ->data ->score_begin_color, ddwin ->data ->score_mid_color, ddwin ->data ->score_end_color,
											   ddwin ->data ->score_begin, ddwin ->data ->score_mid, ddwin ->data ->score_end) ;
				
				r = col.redF () * masks[ms].intensity;
				g = col.greenF () * masks[ms].intensity;
				b = col.blueF () * masks[ms].intensity;
				a = col.alphaF () * masks[ms].intensity;
				
/*                if (q <= ddwin -> data -> score_begin) {
                    r = sb_red   * masks[ms].intensity;
                    g = sb_green * masks[ms].intensity;
                    b = sb_blue  * masks[ms].intensity;
                    a = sb_alpha  * masks[ms].intensity;
                }
                else if (q >= ddwin -> data -> score_end) {
                    r = se_red   * masks[ms].intensity;
                    g = se_green * masks[ms].intensity;
                    b = se_blue  * masks[ms].intensity;
                    a = se_alpha  * masks[ms].intensity;
                }
                else if (q>= ddwin -> data ->score_mid - mid_width && q<=  ddwin -> data ->score_mid + mid_width) {
                    r = sm_red * masks[ms].intensity;
                    g = sm_green *  masks[ms].intensity;
                    b = sm_blue  *  masks[ms].intensity;
                    a = sm_alpha *  masks[ms].intensity;

                }
                else if (q<= ddwin -> data -> score_mid) {
                    float perc = (ddwin -> data -> score_mid - q - mid_width) / (ddwin -> data -> score_mid - ddwin -> data -> score_begin - mid_width);
                    r = (sb_red  * perc + sm_red   * (1-perc)) * masks[ms].intensity;
                    g = (sb_green* perc + sm_green * (1-perc)) * masks[ms].intensity;
                    b = (sb_blue * perc + sm_blue  * (1-perc)) * masks[ms].intensity;
                    a = (sb_alpha * perc + sm_alpha  * (1-perc)) * masks[ms].intensity;
                    
                }
                else {
                    float perc = (ddwin -> data -> score_end - q + mid_width) / (ddwin -> data -> score_end +mid_width- ddwin -> data -> score_mid);
                    r = (sm_red  * perc + se_red   * (1-perc)) * masks[ms].intensity;
                    g = (sm_green* perc + se_green * (1-perc)) * masks[ms].intensity;
                    b = (sm_blue * perc + se_blue  * (1-perc)) * masks[ms].intensity;
                    a = (sm_alpha * perc + se_alpha  * (1-perc)) * masks[ms].intensity;

                }
 
 */
                assert (!(r > 1.f));
                assert (!(g > 1.f));
                assert (!(b > 1.f));
                assert (!(a > 1.f));
                newr += r;
                newg += g;
                newb += b;
                newa += a;
    
            }


            else if (masks[ms].type == COLOR) {
                float r = c_red   * masks[ms].intensity;
                float g = c_green * masks[ms].intensity;
                float b = c_blue  * masks[ms].intensity;
                float a = c_alpha  * masks[ms].intensity;
                assert (!(r > 1.f));
                assert (!(g > 1.f));
                assert (!(b > 1.f));
                newr += r;
                newg += g;
                newb += b;
                newa += a;
            }
        }

        color colo  (newr, newg, newb, newa);
        if (undoable)        color_atom -> add (&*at, colo);

        else {
            set_color (&*at, colo);
        }

    }
    if (undoable) {
        color_atom->set_name ();
        ddwin -> execute (color_atom);
    }
    draw_molecule (mol);
}


void MyGl::set_center_of_rotation (vect v) {
	vect dis = unrotate_vector (subtract (v,center_of_rotation));
	ChangeVectorCommand *command1 = new ChangeVectorCommand (center_of_rotation, v, this,"Set center of rotation");
    ChangeVectorCommand *command2 = new ChangeVectorCommand (view_translations, sum (view_translations,dis), this,"Set center of rotation");
    ddwin -> execute (command1);
	ddwin -> execute (command2);
}
/*
void MyGl::set_center_of_rotation (float x, float y, float z){
    vect v (x, y, z);
    set_center_of_rotation (v);
}
*/

void MyGl::set_center_of_view (vect v) {
	vect dis = subtract (v,center_of_rotation);
	dis = unrotate_vector (dis);
    ChangeVectorCommand *command = new ChangeVectorCommand (view_translations, dis, this, "Set center of view");
    ddwin -> execute (command);
	//view_translations = vect ();

}




void MyGl::screenshot (QString filename){

 GLint viewport [4];
   glGetIntegerv (GL_VIEWPORT, viewport);
   int w = viewport [2]; int h = viewport [3];

   unsigned char* image = new unsigned char[3*w*h];

	FILE* fptr = fopen (filename.toLatin1(), "w" );
   fprintf(fptr, "P6\n");
   fprintf(fptr, "%d %d\n", w, h );
   fprintf(fptr, "255\n" );
   
   glPixelStorei(GL_PACK_ALIGNMENT,1); 
   glReadBuffer(GL_BACK);
   glReadPixels(0,0,w,h,GL_RGB,GL_UNSIGNED_BYTE,image);
   for (int j=h-1;j>=0;j--) {
      for (int i=0;i<w;i++) {
         fputc(image[3*j*w+3*i+0],fptr);
         fputc(image[3*j*w+3*i+1],fptr);
         fputc(image[3*j*w+3*i+2],fptr);
      }
   }
   fclose(fptr);
   delete[] image;





}


/////////////////////////////////////////SLOTS///////////////////////////////////////////////////////////////

void MyGl::head_tracking_update (int x, int y) {
	float scale = 0.1f;
	float xf = x;
	float yf = y;
	xf *= scale;
	yf *= scale;
	float dist = 300.f;
    	head_tracking_x_position = xf;
    	head_tracking_y_position = yf;
	vect new_vec (xf, yf, dist);
	new_vec.normalise ();
	Quat4fT Quat;
	Vector3fT vec1, vec2;
	vec1.s.X = -new_vec.x ();
	vec1.s.Y = new_vec.y ();
	vec1.s.Z = new_vec.z ();

	vec2.s.X = 0.;
	vec2.s.Y = 0.;
	vec2.s.Z = 1.;
	ArcBall.map_vector_on_vector (vec1, vec2, &Quat);

	Matrix3fSetRotationFromQuat4f(&ThisRot, &Quat);		// Convert Quaternion Into Matrix3fT
//	Matrix3fMulMatrix3f(&ThisRot, &Last_Head_Tracking_Rot);				// Accumulate Last Rotation Into This One
	Matrix4fSetRotationFromMatrix3f(&Head_Tracking_Transf, &ThisRot);

	
}

void MyGl::move_camera (float x, float y, float z) {
	vect v (x, y, z);
	view_translations = sum (view_translations, v);

}


void MyGl::move_target (float x, float y, float z) {
	ZNMolecule *mol = ddwin -> target_molecule;
	vect v (x, y, z);
	double rot [9];
	rot [0] = Transform.M [0];
	rot [1] = Transform.M [1];
	rot [2] = Transform.M [2];
	rot [3] = Transform.M [4];
	rot [4] = Transform.M [5];
	rot [5] = Transform.M [6];
	rot [6] = Transform.M [8];
	rot [7] = Transform.M [9];
	rot [8] = Transform.M [10];
	v = rotate_vector_using_matrix_9 (v, rot);
	translate_molecule (mol, v);
	draw_molecule(mol);

}


void MyGl::map_vector_on_vector_world (double x1, double y1, double z1, double x2, double y2, double z2) {
	Quat4fT Quat;
	Vector3fT vec1, vec2;
	vec1.s.X = x1;
	vec1.s.Y = y1;
	vec1.s.Z = z1;

	vec2.s.X = x2;
	vec2.s.Y = y2;
	vec2.s.Z = z2;
	ArcBall.map_vector_on_vector (vec1, vec2, &Quat);


    Matrix3fSetRotationFromQuat4f(&ThisRot, &Quat);		// Convert Quaternion Into Matrix3fT
    Matrix3fMulMatrix3f(&ThisRot, &LastRot);				// Accumulate Last Rotation Into This One
    Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);

	LastRot = ThisRot;

}



void MyGl::map_vector_on_vector_target (double x1, double y1, double z1, double x2, double y2, double z2) {
	quaternion q (1., 0., 0., 0.);
	vect v1 (x1, y1, z1);
	vect v2 (x2, y2, z2);
	double rot [9];
	rot [0] = Transform.M [0];
	rot [1] = Transform.M [1];
	rot [2] = Transform.M [2];
	rot [3] = Transform.M [4];
	rot [4] = Transform.M [5];
	rot [5] = Transform.M [6];
	rot [6] = Transform.M [8];
	rot [7] = Transform.M [9];
	rot [8] = Transform.M [10];
	v1 = rotate_vector_using_matrix_9 (v1, rot);
	v2 = rotate_vector_using_matrix_9 (v2, rot);

	q = map_vector_on_vector_quaternion (v1, v2);
	rotate_molecule (ddwin->target_molecule, q, get_center (ddwin->target_molecule));
	draw_molecule (ddwin->target_molecule);

}

vect MyGl::apply_world_rotation (vect v) {
	double rot [9];
	rot [0] = Transform.M [0];
	rot [1] = Transform.M [1];
	rot [2] = Transform.M [2];
	rot [3] = Transform.M [4];
	rot [4] = Transform.M [5];
	rot [5] = Transform.M [6];
	rot [6] = Transform.M [8];
	rot [7] = Transform.M [9];
	rot [8] = Transform.M [10];
	return rotate_vector_using_matrix_9 (v, rot);
}

vect MyGl::deapply_world_rotation (vect v) {
	double rot [9], inv[9];
	rot [0] = Transform.M [0];
	rot [1] = Transform.M [1];
	rot [2] = Transform.M [2];
	rot [3] = Transform.M [4];
	rot [4] = Transform.M [5];
	rot [5] = Transform.M [6];
	rot [6] = Transform.M [8];
	rot [7] = Transform.M [9];
	rot [8] = Transform.M [10];
	invert_matrix_9 (rot, inv);
	return rotate_vector_using_matrix_9 (v, inv);
}

////////////////////////////////////UTILITIES///////////////////////////////////////////////////////////////


void MyGl::haptic_to_world_coordinates (vect &haptic_p, vect &world_p, float minx, float maxx, float miny, float maxy, float minz, float maxz) {

	vect upv, downv, leftv, rightv;
	get_viewport_points (upv, downv, leftv, rightv);
	float leftx = leftv.x ();
	float lefty = leftv.y ();
	float leftz = leftv.z ();

	float rightx = rightv.x ();
	float righty = rightv.y ();
	float rightz = rightv.z ();

	float upx = upv.x ();
	float upy = upv.y ();
	float upz = upv.z ();

	float downx = downv.x ();
	float downy = downv.y ();
	float downz = downv.z ();

	
    vect up (upx-downx, upy-downy, upz-downz);
    float modup = up.module ();
    up.multiply (1.f/modup);


    vect right (rightx-leftx, righty-lefty, rightz-leftz);
    float modright = right.module ();
    right.multiply (1.f/modright);


    float scale = modup;
    if (modup > modright) scale = modright;
    vect out;
    out = cross_product (right, up);


    vect central_point;



    central_point.x() = leftx + (right.x() * modright)/2;
    central_point.y() = lefty + (right.y() * modright)/2;
    central_point.z() = leftz + (right.z() * modright)/2;

    float hx = (haptic_p.x()-(minx + maxx)/2) / (maxx - minx);
    float hy = (haptic_p.y()-(miny + maxy)/2) / (maxy - miny);
    float hz = (haptic_p.z()-(minz + maxz)/2) / (maxz - minz);

    world_p.x() = central_point.x() + scale * (right.x() * hx + up.x() * hy + out.x() * hz);
    world_p.y() = central_point.y() + scale * (right.y() * hx + up.y() * hy + out.y() * hz);
    world_p.z() = central_point.z() + scale * (right.z() * hx + up.z() * hy + out.z() * hz);

}





void MyGl::world_to_haptic_coordinates (vect &world_p, vect &haptic_p, float minx, float maxx, float miny, float maxy, float minz, float maxz) {

	double rot [9], inv [9];
	rot [0] = Transform.M [0];
	rot [1] = Transform.M [1];
	rot [2] = Transform.M [2];
	rot [3] = Transform.M [4];
	rot [4] = Transform.M [5];
	rot [5] = Transform.M [6];
	rot [6] = Transform.M [8];
	rot [7] = Transform.M [9];
	rot [8] = Transform.M [10];
	invert_matrix_9 (rot, inv);
	world_p = rotate_vector_using_matrix_9 (world_p, inv);

/*
	vect upv, downv, leftv, rightv;
	get_viewport_points (upv, downv, leftv, rightv);

	float leftx = leftv.x ();
	float lefty = leftv.y ();
	float leftz = leftv.z ();

	float rightx = rightv.x ();
	float righty = rightv.y ();
	float rightz = rightv.z ();

	float upx = upv.x ();
	float upy = upv.y ();
	float upz = upv.z ();

	float downx = downv.x ();
	float downy = downv.y ();
	float downz = downv.z ();


    vect up = upv - downv;
    float modup = up.module ();
    up.multiply (1.f/modup);


    vect right (rightx-leftx, righty-lefty, rightz-leftz);
    float modright = right.module ();
    right.multiply (1.f/modright);


    float scale = modup;
    if (modup > modright) scale = modright;
    vect out;
    out = cross_product (right, up);

    vect none;

    vect up_world, right_world, out_world;

    components (world_p, up, up_world, none);
    components (world_p, right, right_world, none);
    components (world_p, out, out_world, none);



    haptic_p.x() = right_world.module ();
    haptic_p.y() = up_world.module ();
    haptic_p.z() = out_world.module ();
*/

}



vect MyGl::rotate_vector (vect v) {
	double rot [9];
	rot [0] = Transform.M [0];
	rot [1] = Transform.M [1];
	rot [2] = Transform.M [2];
	rot [3] = Transform.M [4];
	rot [4] = Transform.M [5];
	rot [5] = Transform.M [6];
	rot [6] = Transform.M [8];
	rot [7] = Transform.M [9];
	rot [8] = Transform.M [10];
	return rotate_vector_using_matrix_9 (v, rot);

}


vect MyGl::unrotate_vector (vect v) { 
	double rot [9], inv [9];
	rot [0] = Transform.M [0];
	rot [1] = Transform.M [1];
	rot [2] = Transform.M [2];
	rot [3] = Transform.M [4];
	rot [4] = Transform.M [5];
	rot [5] = Transform.M [6];
	rot [6] = Transform.M [8];
	rot [7] = Transform.M [9];
	rot [8] = Transform.M [10];
	invert_matrix_9 (rot, inv);
	return rotate_vector_using_matrix_9 (v, inv);
}

void MyGl::draw_backbone_stick (Resid *res, Surface *surf) {
	int n_points = *ddwin ->data ->quality_scale * 9;
	vector <vect> random_points;
	vector <vect> helix_points;
	vector <vect> sheet_points;
	float ha = *ddwin ->data ->backbone_tube_helix_a;
	float hb = *ddwin ->data ->backbone_tube_helix_b;
		float hc = *ddwin ->data ->backbone_tube_helix_c;
		float sa = *ddwin ->data ->backbone_tube_sheet_a;
		float sb = *ddwin ->data ->backbone_tube_sheet_b;
		float sc = *ddwin ->data ->backbone_tube_sheet_c;
	float ra = *ddwin ->data ->backbone_tube_random_a;
	float rb = *ddwin ->data ->backbone_tube_random_b;
	float rc = *ddwin ->data ->backbone_tube_random_c;
	
	for (unsigned int i = 0; i < n_points; i++) {
		float da = 2*M_PI/(n_points-1);
		float angl = i * da;
		random_points.push_back(vect (sin (angl)*rb- (rc*rb *sin (angl)*sin (angl)*sin(angl)), cos (angl)*ra , 0.));
		helix_points.push_back(vect (sin (angl)*hb- (hc*hb *sin (angl)*sin (angl)*sin(angl)), cos (angl)*ha , 0.));
		sheet_points.push_back(vect (sin (angl)*sb- (sc*sb *sin (angl)*sin (angl)*sin(angl)), cos (angl)*sa , 0.));
	}
	vector <vect> *last_shape, *shape, *next_shape;
	color last_col, col, next_col;
	col = get_color (res);
	last_col = col;
	next_col = col;
	last_shape = &random_points;
	next_shape = &random_points;
	shape = &random_points;
	
	if (is_helix (res)) shape = &helix_points;
	else if (is_sheet (res)) shape = &sheet_points;
	
	vect dir = get_backbone_direction (res);
	Resid *prec_res = get_previous_residue(res);
	Resid *follow_res = get_following_residue(res);	
	vect lastdir, nextdir;
	if (prec_res) {
		last_col = get_color (prec_res);
		lastdir = get_backbone_direction(prec_res);
		if (is_helix (prec_res)) last_shape = &helix_points; 
		else if (is_sheet (prec_res)) last_shape = &sheet_points; 
	}
	else lastdir = dir;
	if (follow_res) {
		next_col = get_color (follow_res);
		nextdir = get_backbone_direction(follow_res);
		if (is_helix (follow_res)) next_shape = &helix_points; 
		else if (is_sheet (follow_res)) next_shape = &sheet_points; 
	}
	else nextdir = dir;
	if (dot_product(lastdir, dir) < 0.) lastdir.multiply(-1.);
	if (dot_product(nextdir, dir) < 0.) nextdir.multiply(-1.);	

	lastdir = mean (dir, lastdir);
	nextdir = mean (nextdir, dir);
	lastdir.normalise();
	nextdir.normalise();
	vector <vect> points = get_backbone_points (res);
	add_guide_points_to_backbone (&*res, points);

	color c2 = mean (last_col, col);
	color c1 = mean (next_col, col);
//	if (is_helix (res)) {c1 = c2 = color (1.f, 0.f, 0.f);}
//	else if (is_sheet (res)) {c1 = c2 = color (1.f, 1.f, 0.f);}

	
	
	if (points.size () > 3) {
		float tot = ((float) points.size () -3);

		for (unsigned int i = 3; i < points.size (); i++) {
			
			float dt = ((float) i-3) / tot;
			dt = sin ((dt-0.5f)*M_PI);
			
			dt*= 0.5f;
			dt += 0.5f;

			color cc1, cc2;
			cc1 = color ((float)(c1.redF()*dt+c2.redF()*(1-dt)),((float) c1.greenF()*dt+c2.greenF()*(1-dt)),(float) (c1.blueF()*dt+c2.blueF()*(1-dt)),((float) c1.alphaF()*dt+c2.alphaF()*(1-dt)));
			vect v1, v2;
			v2 = lastdir;
			v1 = nextdir;
			vector <vect> shape1, shape2;
			for (unsigned int n = 0; n < n_points; n++) {
				vect vv2 = mean ((*last_shape)[n], (*shape) [n]);
				vect vv1 = mean ((*next_shape)[n], (*shape) [n]);	
				shape1.push_back (vect (vv1.x()*dt + vv2.x()*(1-dt), vv1.y()*dt + vv2.y()*(1-dt), vv1.z()*dt + vv2.z()*(1-dt) ) );
			}

			vect d (v1.x()*dt + v2.x()*(1-dt), v1.y()*dt + v2.y()*(1-dt), v1.z()*dt + v2.z()*(1-dt)  );

			dt = ((float) i-2) / tot;
			dt = sin ((dt-0.5f) * PI);
			dt*= 0.5f;
			dt += 0.5f;
			
			cc2 = color ((float)(c1.redF()*dt+c2.redF()*(1-dt)),((float) c1.greenF()*dt+c2.greenF()*(1-dt)),(float) (c1.blueF()*dt+c2.blueF()*(1-dt)),((float) c1.alphaF()*dt+c2.alphaF()*(1-dt)));

			for (unsigned int n = 0; n < n_points; n++) {
				vect vv2 = mean ((*last_shape)[n], (*shape) [n]);
				vect vv1 = mean ((*next_shape)[n], (*shape) [n]);	
				shape2.push_back (vect (vv1.x()*dt + vv2.x()*(1-dt), vv1.y()*dt + vv2.y()*(1-dt), vv1.z()*dt + vv2.z()*(1-dt) ) );
			}
			vect d2 (v1.x()*dt + v2.x()*(1-dt), v1.y()*dt + v2.y()*(1-dt), v1.z()*dt + v2.z()*(1-dt)  );
		//	cerr << d<<d2<<lastdir<<nextdir<<endl;
		//	d2.normalise ();

		//	d = d2 = vect (1., 0., 0.);

			my_backbone_ribbon (points[i-3], points [i-2],points [i-1],points [i], d, d2, cc1, cc2, shape1, shape2, surf);

		} 

	}

}

void MyGl::draw_backbone_line (Resid *res) {
	vector <vect> points = get_backbone_points (res);
	color last_col, col, next_col;
	col = get_color (res);
	last_col = col;
	next_col = col;
	Resid *prec_res = get_previous_residue(res);
	Resid *follow_res = get_following_residue(res);	
	if (prec_res) {
		last_col = get_color (prec_res);
	}
	if (follow_res) {
		next_col = get_color (follow_res);
	
	}
	
	
	color c2 = mean (last_col, col);
	color c1 = mean (next_col, col);
	
	
	
	if (points.size () > 1) {
		for (unsigned int i = 1; i < points.size (); i++) {
			float dt = ((float) i-1) / points.size ();
			color cc1 = color ((float)(c1.redF()*dt+c2.redF()*(1-dt)),((float) c1.greenF()*dt+c2.greenF()*(1-dt)),(float) (c1.blueF()*dt+c2.blueF()*(1-dt)),((float) c1.alphaF()*dt+c2.alphaF()*(1-dt)));
			dt = ((float) i) / points.size ();
			color cc2 = color ((float)(c1.redF()*dt+c2.redF()*(1-dt)),((float) c1.greenF()*dt+c2.greenF()*(1-dt)),(float) (c1.blueF()*dt+c2.blueF()*(1-dt)),((float) c1.alphaF()*dt+c2.alphaF()*(1-dt)));

			my_line (points[i-1], points [i], cc1, cc2);
		} 
	}
	
}


void MyGl::my_cylinder (float radone, float radtwo, float lenght,unsigned int slices, Atom* at1, Atom*at2){
	
    glBegin (GL_QUAD_STRIP);
    float angle = 0.0;
    float da = 2*PI/slices;

    for (unsigned int n=0; n<slices+1; n++){
        angle = n*da;
        setAtomColor(at1);
        glNormal3f (sin (angle), cos (angle), 0);
        glVertex3f (sin (angle)*radone, cos (angle)*radone, 0);

        setAtomColor(at2); 
        glNormal3f (sin (angle), cos (angle), 0);       
        glVertex3f (sin (angle)*radtwo, cos (angle)*radtwo, lenght);
    }
    glEnd ();    
}

void MyGl::my_line (vect v1, vect v2, color c1, color c2) {
	glBegin(GL_LINES);
		openGLSetColor(c1);
        glVertex3f (v1.x(), v1.y(), v1.z());
		
		openGLSetColor(c2);
        glVertex3f (v2.x(), v2.y(), v2.z());
	glEnd();
}


void MyGl::my_backbone_ribbon (vect v1, vect v2, vect v3, vect v4, vect dir, vect dir2, color c1, color c2, vector <vect> shape1, vector <vect> shape2, Surface *surf) {
	//Sandri's method
	vect prec_vect = subtract (v2, v1);
	vect cyl_vect = subtract (v3, v2);
	vect post_vect = subtract (v4, v3);
	prec_vect.normalise();
	post_vect.normalise();
	cyl_vect.normalise();
	vect norm1 = mean (prec_vect, cyl_vect);
	vect norm2 = mean (cyl_vect, post_vect);

	
	
	vect par1, pp1, par2, pp2;
	components (dir, norm1, par1, pp1);
	components (dir2, norm2, par2, pp2);	
	
	double m1 [9], m2 [9], m [9];	
	
/*	
	
	vect newz = cyl_vect;
	vect newy = dir;

	vect newx = cross_product(newz, newy);
	newy = cross_product(newz, newx);
	newz.normalise();
	newy.normalise();
	newx.normalise();
	
	
	
	
	m[0] = newx.x ();
	m[3] = newx.y ();
	m[6] = newx.z ();	
	m[1] = newy.x ();
	m[4] = newy.y ();
	m[7] = newy.z ();	
	m[2] = newz.x ();
	m[5] = newz.y ();
	m[8] = newz.z ();
	quaternion q = map_vector_on_vector_quaternion(vect (0., 0., 1.), cyl_vect);
	
	
*/	
	
	
	
	
	
	
	vect newz1 = norm1;
	vect newy1 = pp1;
	newz1.normalise();
	newy1.normalise();
	vect newx1 = cross_product(newy1, newz1);

	m1[0] = newx1.x ();
	m1[3] = newx1.y ();
	m1[6] = newx1.z ();	
	m1[1] = newy1.x ();
	m1[4] = newy1.y ();
	m1[7] = newy1.z ();	
	m1[2] = newz1.x ();
	m1[5] = newz1.y ();
	m1[8] = newz1.z ();	
	
	
	vect newz2 = norm2;
	vect newy2 = pp2;
	newz2.normalise();
	newy2.normalise();
	vect newx2 = cross_product(newy2, newz2);

	m2[0] = newx2.x ();
	m2[3] = newx2.y ();
	m2[6] = newx2.z ();	
	m2[1] = newy2.x ();
	m2[4] = newy2.y ();
	m2[7] = newy2.z ();	
	m2[2] = newz2.x ();
	m2[5] = newz2.y ();
	m2[8] = newz2.z ();	
	
	float angl = 0.f;
	
	glBegin (GL_QUAD_STRIP);
	int slices = shape1.size ();
	vect lastp1, lastp2, lastn1, lastn2;
	SurfVertex *lastv1, *lastv2;

    for (unsigned int n=0; n<slices; n++){

		
		vect p1 = shape1[n];
		vect p2 = shape2[n];
		
		int last_n = n-1;
		if (last_n < 0) last_n = slices-2;
		int next_n = n+1;
		if (next_n >= slices) next_n = 1;
		
		vect tan1 = subtract(shape1[next_n], shape1[last_n]);
		vect tan2 = subtract(shape2[next_n], shape2[last_n]);
		vect n1 = vect (-tan1.y(), tan1.x (), 0);
		vect n2 = vect (-tan2.y(), tan2.x (), 0);		
		
		p1 = rotate_vector_using_matrix_9(p1, m1);
		p2 = rotate_vector_using_matrix_9(p2, m2);	
	
		n1 = rotate_vector_using_matrix_9(n1, m1);
		n2 = rotate_vector_using_matrix_9(n2, m2);
		
	//	p1 = rotate_vector_using_quaternion(p1, q);
//		p2 = rotate_vector_using_quaternion(p2, q);	
		
	//	vect n1 = p1;
	//	vect n2 = p2;
		
		p1 += v2;
		p2 += v3;
		
		if (!surf) {
			openGLSetColor(c1);
			glNormal3f (n1.x (), n1.y (), n1.z ());
			glVertex3f (p1.x (), p1.y (), p1.z ());
		
		
			openGLSetColor(c2); 
			glNormal3f (n2.x (), n2.y (), n2.z ());       
			glVertex3f (p2.x (), p2.y (), p2.z ());
		}
		

		
		
		else {
			int num = surf ->vertices.size ();
			SurfVertex *newv1 = new SurfVertex;
			newv1->n = num;
			newv1->normal = n1;
			newv1->coordinates = p1;
			newv1 ->col = c1;
			
			SurfVertex *newv2 = new SurfVertex;
			newv2->n = num+1;
			newv2->normal = n2;
			newv2->coordinates = p2;
			newv2 ->col = c2;
			if (n>0) {

				SurfFace *face = new SurfFace;	
				face ->v1 = lastv1;
				face ->v2 = lastv2;
				face ->v3 = newv2;				
				SurfFace *face2 = new SurfFace;
				face2 ->v1 = newv2;
				face2 ->v2 = newv1;
				face2 ->v3 = lastv1;

				surf ->faces.push_back (face);
				surf ->faces.push_back (face2);				
			}
			lastp1 = p1;
			lastp2 = p2;
			lastn1 = n1;
			lastn2 = n2;
			lastv1 = newv1;
			lastv2 = newv2;
			surf ->vertices.push_back (lastv1);
			surf ->vertices.push_back (lastv2);
			
		}
    }
    glEnd ();

	
	
	
}



void MyGl::my_cylinder (vect v1, vect v2, vect v3, vect v4, float radone, float radtwo, color c1, color c2, unsigned int slices) {
	vect prec_vect = subtract (v2, v1);
	vect cyl_vect = subtract (v3, v2);
	vect post_vect = subtract (v4, v3);
	prec_vect.normalise();
	post_vect.normalise();
	cyl_vect.normalise();
	vect norm1 = mean (prec_vect, cyl_vect);
	vect norm2 = mean (cyl_vect, post_vect);

	
	glBegin (GL_QUAD_STRIP);
    float angle = 0.0;
    float da = 2*PI/slices;
	quaternion q1 (0., 1., 0., 0.);
	quaternion q2 (0., 1., 0., 0.);
	q1 = map_vector_on_vector_quaternion(vect (0., 0., 1.), norm1);
	q2 = map_vector_on_vector_quaternion(vect (0., 0., 1.), norm2);


	
	quaternion spin_quaternion1 (0., 1., 0., 0.);
	quaternion spin_quaternion2 (0., 1., 0., 0.);		
	vect pp1 = vect (0., radone, 0.);
	pp1 = rotate_vector_using_quaternion(pp1, q1);	
	vect pp2 = vect (0., radtwo, 0.);
	pp2 = rotate_vector_using_quaternion(pp2, q2);	
	
	vect ref (3600., 0., 0.);
	vect origin (0., 0., 0.);
	double spin_angle1 = dihedral(pp1, origin, norm1, sum (norm1, ref))*M_PI/180.;
	spin_quaternion1 = axis_angle_to_quaternion(norm1, spin_angle1);

	double spin_angle2 = dihedral(pp2, origin, norm2, sum (norm2, ref))*M_PI/180.;
	spin_quaternion2 = axis_angle_to_quaternion(norm2, spin_angle2);
	

	
    for (unsigned int n=0; n<slices+1; n++){
        angle = n*da;

		vect p1 = vect ( sin (angle)*radone, cos (angle)*radone, 0.);
		vect p2 = vect ( sin (angle)*radtwo, cos (angle)*radtwo, 0.);
		

		
		p1 = rotate_vector_using_quaternion(p1, q1);
		p2 = rotate_vector_using_quaternion(p2, q2);	
		

//		p1 = rotate_vector_using_quaternion(p1, spin_quaternion1);		
//		p2 = rotate_vector_using_quaternion(p2, spin_quaternion2);	
		
		p1 += v2;
		p2 += v3;



		openGLSetColor(c1);
		vect n1 = subtract (p1,v2);
		vect n2 = subtract (p2,v3);
        glNormal3f (n1.x (), n1.y (), n1.z ());
        glVertex3f (p1.x (), p1.y (), p1.z ());

		
        openGLSetColor(c2); 
        glNormal3f (n2.x (), n2.y (), n2.z ());       
        glVertex3f (p2.x (), p2.y (), p2.z ());
    }
    glEnd ();
    glPopMatrix ();
	
	
	
	
}
	
	
void MyGl::my_cylinder (vect v1, vect v2, float radone, float radtwo, color c1, color c2, unsigned int slices) {
    float angle_dir, height, mod_of_vector;
    vect at1c = v1;
    vect at2c = v2;
	
    vect vector = subtract (at2c, at1c);
    mod_of_vector = vector.module ();
	
    vector *= 1/mod_of_vector;
    height = mod_of_vector;
	
    angle_dir = acos (vector.z ()) * 180.0 / PI;
    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	
	//    setAtomColor(bond->GetBeginAtom ());
	

    glPushMatrix ();
	
	
	
    glTranslated (at1c.x (), at1c.y(), at1c.z());
    glRotated (180.0, 0.0, 0.0, 1.0);
    glRotated (angle_dir, vector.y (), -1.0 * vector.x (), 0.0);
	

	glBegin (GL_QUAD_STRIP);
    float angle = 0.0;
    float da = 2*PI/slices;
	
    for (unsigned int n=0; n<slices+1; n++){
        angle = n*da;
        openGLSetColor(c1);
        glNormal3f (sin (angle), cos (angle), 0);
        glVertex3f (sin (angle)*radone, cos (angle)*radone, 0);
		
        openGLSetColor(c2); 
        glNormal3f (sin (angle), cos (angle), 0);       
        glVertex3f (sin (angle)*radtwo, cos (angle)*radtwo, height);
    }
    glEnd ();
    glPopMatrix ();
	
}



void MyGl::my_sphere (float rad, unsigned int slices, unsigned int stacks, Atom* at){
/*
//    setAtomColor(at);
    float da = 2*PI/slices;
    float x = 0;
    float y = 0;

    float x2 = 0;
    float y2 = 0;

    float vx =0;
    float vy =0;
    float vz =0;

    float vx2 =0;
    float vy2 =0;
    float vz2 =0;

    float px2 =0;
    float py2 =0;
    float pz2 =0;

    float px =0;
    float py =0;
    float pz =0;

    float angle;
    float da2 = PI/stacks;
    float r = 0.5;     float g = 0.5;    float b = 0.5;    float a = 1.;
    float dr = 0;   float dg = 0; float db = 0;

    vector <float> rr;

        for (unsigned int i=0; i<at->bound.size();i++) {
                float bx = at->bound[i]-> GetVector ().x()-at-> GetVector ().x();
                float by = at->bound[i]-> GetVector ().y()-at-> GetVector ().y();
                float bz = at->bound[i]-> GetVector ().z()-at-> GetVector ().z();
                rr.push_back (sqrt (bx*bx+by*by+bz*bz));

        }   

    for (unsigned int n=0; n<slices; n++){


        
        angle = n*da;
        x = sin (angle);
        y = cos (angle);

        x2 = sin (angle+da);
        y2 = cos (angle+da);

        glBegin (GL_QUAD_STRIP);
        for (unsigned int m=0; m<stacks+1; m++){
            r = at -> col.redF ();
            g = at -> col.greenF ();
            b = at -> col.blueF ();
            a = at -> col.alphaF ();

            float angle2 = -PI+m*da2;

            vx = x*sin (angle2); vy = y*sin (angle2); vz= cos(angle2);
            vx2 =x2*sin (angle2); vy2 = y2*sin (angle2); vz2= cos (angle2);
            px = vx*rad; py = vy*rad; pz = vz*rad;
            px2 = vx2*rad; py2 = vy2*rad; pz2 = vz2*rad;
            for (unsigned int i=0; i<at->bound.size(); i++) {
                float bx = at->bound[i]-> GetVector ().x()-at-> GetVector ().x();
                float by = at->bound[i]-> GetVector ().y()-at-> GetVector ().y();
                float bz = at->bound[i]-> GetVector ().z()-at-> GetVector ().z();
                float d=sqrt((bx-px2)*(bx-px2)+(by-py2)*(by-py2)+(bz-pz2)*(bz-pz2));
                float rd = rr[i];
                float dc = (rd-d)/rd;
                if (dc<0) dc=0;

                dr+=(-at-> col.redF ()+at->bound[i]-> col.redF ())*dc;
                dg+=(-at-> col.greenF ()+at->bound[i]-> col.greenF ())*dc;
                db+=(-at-> col.blueF ()+at->bound[i]-> col.blueF ())*dc;
            }

            r+=dr; g+=dg; b+=db;
          
            glColor4f (r, g, b, a);

            glNormal3f ( vx2, vy2, vz2);
            glVertex3f (px2, py2, pz2);

            dr=0; dg=0; db=0; 
            r = at -> col.redF ();
            g = at -> col.greenF ();
            b = at -> col.blueF ();
            a = at -> col.alphaF ();
            for (unsigned int i=0; i<at->bound.size(); i++) {
                if (at->bound[i]->visible==true) {
                    float bx = at->bound[i]-> GetVector ().x()-at-> GetVector ().x();
                    float by = at->bound[i]-> GetVector ().y()-at-> GetVector ().y();
                    float bz = at->bound[i]-> GetVector ().z()-at-> GetVector ().z();
                    float d=sqrt((bx-px)*(bx-px)+(by-py)*(by-py)+(bz-pz)*(bz-pz));
                    float rd = rr[i];
                    float dc = (rd-d)/rd;
                    if (dc<0) dc=0;

                    dr+=(-at-> col.redF ()+at->bound[i]-> col.redF ())*dc;
                    dg+=(-at-> col.greenF ()+at->bound[i]-> col.greenF ())*dc;
                    db+=(-at-> col.blueF ()+at->bound[i]-> col.blueF ())*dc;
                }
            }


            r+=dr; g+=dg; b+=db;

            glColor4f (r, g, b, a);

            glNormal3f (vx, vy, vz); 
            glVertex3f (px, py, pz);

        }

        glEnd ();    
          
    }   
*/
}




void MyGl::compute_double_bond_vertexes (ZNBond *bond, float out [4][3], float d) {

    bool b0a = false;
    bool b0b = false;
    bool b1a = false;
    bool b1b = false;

    if (!d) d = *ddwin ->data ->double_bond_separation/2;
    Atom *at0 = bond -> GetBeginAtom ();
    Atom *at1 = bond -> GetEndAtom ();
    vect c0 = get_coordinates(at0);
    vect c1 = get_coordinates(at1);

    vect v01 = subtract (c1, c0);
    vect v10 = subtract (c0, c1);
    vect v0a = c0;
    vect v0b = c0;
    vect v1a = c1;
    vect v1b = c1;

    int cb0 = CountBonds (at0);
    int cb1 = CountBonds (at1);

    if (0) { //bond -> IsInRing ()) {
    }
    else {
        
        if (cb0 > 3) {
        }
        else if (cb0 == 3) {
            int found = 0;
            Atom *b0 = NULL;
            Atom *b1 = NULL;
            FOR_NBORS_OF_ATOM (n, at0) {
                if (&*n != at1) {
                    if (found == 0) b0 = &*n;
                    else if (found == 1) b1 = &*n;
                    else break;
                    found++;
                }
            }
            assert (b0);
            assert (b1);
            vect veccb0 = get_coordinates (b0);
            vect veccb1 = get_coordinates (b1);

            v0a = subtract (veccb0, c0);  //some form of checking which side of the axis we are is needed to avoid crossing double bonds

            v0b = subtract (veccb1, c0);
            b0a = true;
            b0b = true;

        }
        else if (cb0 == 2) {
        }
        else  {
        }


        if (cb1 > 3) {
        }
        else if (cb1 == 3) {
            int found = 0;
            Atom *b0 = NULL;
            Atom *b1 = NULL;
            FOR_NBORS_OF_ATOM (n1, at1) {
                if (&*n1 != at0) {
                    if (found == 0) b0 = &*n1;
                    else if (found == 1) b1 = &*n1;
                    else break;
                    found++;
                }
            }
            assert (b0);
            assert (b1);
            vect cb0 = get_coordinates(b0);
            vect cb1 = get_coordinates (b1);

            v1a = subtract (cb0, c1);  //some form of checking which side of the axis we are is needed to avoid crossing double bonds
            v1b = subtract (cb1, c1);
            b1a = true;
            b1b = true;

        }
        else if (cb1 == 2) {
        }
        else {
        }



    }
    vect parallel, normal;

    if ((b0a && b0b) && !(b1a && b1b)) {
        components (v0a, v01, parallel, normal);
        v1a = normal;
        normal.multiply (-1);
        v1b = normal;

    }
    else if (!(b0a && b0b) && (b1a && b1b)) {
        components (v1a, v10, parallel, normal);
        v0a = normal;
        normal.multiply (-1.);
        v0b = normal;

    }
    else if ((b0a && b0b) && (b1a && b1b)) {}

    else {
        vect base (1., 1., 1.);
        components (base, v01, parallel, normal);
        v0a = normal;
        v1a = normal;
        normal.multiply (-1.);    
        v0b = normal;
        v1b = normal;
    }

    assert (!isnan (v0a.x()));
    assert (!isnan (v01.x()));


    components (v0a, v01, parallel, normal);
    assert (normal.module ());
    assert (!isnan(normal.module ()));
    v0a.multiply (d/normal.module ());

    components (v0b, v01, parallel, normal);
    v0b.multiply (d/normal.module ());

    components (v1a, v10, parallel, normal);
    v1a.multiply (d/normal.module ());

    components (v1b, v10, parallel, normal);
    v1b.multiply (d/normal.module ());

    assert (!isnan (v0a.x ()));

    float dist1 = dist (v0a, v1a);
    float dist2 = dist (v0a, v1b);
    if (dist2 < dist1) {
	vect swap = v1a;
	v1a = v1b;
	v1b = swap;
    }


    vect o0, o1, o2, o3;
    o0 = sum (c0, v0a);
    o1 = sum (c0, v0b);
    o2 = sum (c1, v1a);
    o3 = sum (c1, v1b);

    assert (!isnan (o0.x ()));
    out[0][0] = o0.x();
    out[0][1] = o0.y();
    out[0][2] = o0.z();

    out[1][0] = o1.x();
    out[1][1] = o1.y();
    out[1][2] = o1.z();


    out[2][0] = o2.x();
    out[2][1] = o2.y();
    out[2][2] = o2.z();

    out[3][0] = o3.x();
    out[3][1] = o3.y();
    out[3][2] = o3.z();



/*
    float c01 [3];
    float c02 [3];
    float c11 [3];
    float c12 [3];
    c01[0]=0.; c01[1]=0.; c01[2]=0.;
    c02[0]=0.; c02[1]=0.; c02[2]=0.;
    c11[0]=0.; c11[1]=0.; c11[2]=0.;
    c12[0]=0.; c12[1]=0.; c12[2]=0.;
    float swapx, swapy, swapz;
    bool found0 = false;
    bool found1 = false;


    if (at0 -> bound.size ()>2) {
        unsigned int oldn = 0;
        for (unsigned int n = 0; n<at0->bound.size(); n++) {
            if (at0->bound[n]->ID!= at1->ID) {
                c01[0] = at0->bound[n]-> GetVector ().x(); 
                c01[1] = at0->bound[n]-> GetVector ().y();
                c01[2] = at0->bound[n]-> GetVector ().z();
                oldn = n;      
                break;
            }   
        }
        for (unsigned int n1 = oldn+1; n1<at0->bound.size(); n1++) {
            if (at0->bound[n1]->ID!= at1->ID) {
                c02[0] = at0->bound[n1]-> GetVector ().x();       
                c02[1] = at0->bound[n1]-> GetVector ().y();
                c02[2] = at0->bound[n1]-> GetVector ().z();
                break;
            }   
        }
        found0 = true;
    }
    
    if (at1->bound.size ()>2) {
        unsigned int oldn = 0;
        for (unsigned int n = 0; n<at1->bound.size(); n++) {
            if (at1->bound[n]->ID!= at0->ID) {
                c11[0] = at1->bound[n]-> GetVector ().x();   
                c11[1] = at1->bound[n]-> GetVector ().y();   
                c11[2] = at1->bound[n]-> GetVector ().z();   
                oldn = n;    
                break;
            }   
        }
        for (unsigned int n1=oldn+1; n1<at1->bound.size(); n1++) {
            if (at1->bound[n1]->ID!= at0->ID) {
                c12[0] = at1->bound[n1]-> GetVector ().x();       
                c12[1] = at1->bound[n1]-> GetVector ().y();       
                c12[2] = at1->bound[n1]-> GetVector ().z();       
                break;
            }   
        }
        found1 = true;
    }

    float refx, refy, refz;

    vect vec, ref, nor, par;



    ref = subtract (at1 -> GetVector (), at0 -> GetVector ());
    ref.multiply (-1);

    if (!found0 && found1) {
    vec.null (); nor.null (); par.null ();
   //     cout <<"not found 0"<<endl;
        vec = vect (c11[0]-at1-> GetVector ().x(), c11[1]-at1-> GetVector ().y(),c11[2]-at1-> GetVector ().z());
        components (vec, ref, par, nor);
        c01[0]=at0-> GetVector ().x()+nor.x();
        c01[1]=at0-> GetVector ().y()+nor.y();
        c01[2]=at0-> GetVector ().z()+nor.z();
        c02[0]=at0-> GetVector ().x()-nor.x();
        c02[1]=at0-> GetVector ().y()-nor.y();
        c02[2]=at0-> GetVector ().z()-nor.z();
    }
  //      ref[0]*=-1;        ref[1]*=-1;        ref[2]*=-1;
    if (found0 && !found1) {
    vec.null (); nor.null (); par.null ();
   //     cout <<"not found 1"<<endl;
        vec = vect (c01[0]-at0-> GetVector ().x(), c01[1]-at0-> GetVector ().y(), c01[2]-at0-> GetVector ().z() );
        components (vec, ref, par, nor);
        c11[0]=at1-> GetVector ().x()+nor.x();
        c11[1]=at1-> GetVector ().y()+nor.y();
        c11[2]=at1-> GetVector ().z()+nor.z();
        c12[0]=at1-> GetVector ().x()-nor.x();
        c12[1]=at1-> GetVector ().y()-nor.y();
        c12[2]=at1-> GetVector ().z()-nor.z();
    }

    


    if ( (c11[0]-c01[0])*(c11[0]-c01[0]) + (c11[1]-c01[1])*(c11[1]-c01[1]) + (c11[2]-c01[2])*(c11[2]-c01[2]) >
(c12[0]-c01[0])*(c12[0]-c01[0]) + (c12[1]-c01[1])*(c12[1]-c01[1]) + (c12[2]-c01[2])*(c12[2]-c01[2])) {
        swapx = c01[0]; swapy = c01[1]; swapz = c01[2]; 
        c01[0] = c02[0]; c01[1] = c02[1]; c01[2] = c02[2];
        c02[0] = swapx; c02[1] = swapy; c02[2] = swapz;
    }
   

    //float v01x, v01y, v01z, v02x, v02y, v02z, v11x, v11y, v11z, v12x, v12y, v12z,





    out[0][0] = c01[0]-at0-> GetVector ().x();
    out[0][1] = c01[1]-at0-> GetVector ().y();
    out[0][2] = c01[2]-at0-> GetVector ().z();

    out[1][0] = c02[0]-at0-> GetVector ().x();
    out[1][1] = c02[1]-at0-> GetVector ().y();
    out[1][2] = c02[2]-at0-> GetVector ().z();

    out[2][0] = c11[0]-at1-> GetVector ().x();
    out[2][1] = c11[1]-at1-> GetVector ().y();
    out[2][2] = c11[2]-at1-> GetVector ().z();

    out[3][0] = c12[0]-at1-> GetVector ().x();
    out[3][1] = c12[1]-at1-> GetVector ().y();
    out[3][2] = c12[2]-at1-> GetVector ().z();


    vec.null (); nor.null (); par.null ();

    vec = vect (out[0][0],out[0][1], out[0][2]);
//    cout << "case 1"<<endl;
    components (vec, ref, par, nor);
    float module = nor.module ();
    out[0][0] *= d / module; 
    out[0][1] *= d / module;
    out[0][2] *= d / module;

    vec.null (); nor.null ();

    vec = vect (out[1][0], out[1][1] ,out[1][2]);
  //  cout << "case 2"<<endl;
    components (vec, ref, par, nor);
    module = nor.module ();
    out[1][0] *= d / module; 
    out[1][1] *= d / module;
    out[1][2] *= d / module;


    vec = vect (out[2][0], out[2][1] ,out[2][2]);

//    cout << "case 3"<<endl;
    components (vec, ref, par, nor);
    module = nor.module ();
    out[2][0] *= d / module; 
    out[2][1] *= d / module;
    out[2][2] *= d / module;

    vec = vect (out[3][0], out[3][1] ,out[3][2]);
//    cout << "case 4"<<endl;
    components (vec, ref, par, nor);
    module = nor.module ();
    out[3][0] *= d / module; 
    out[3][1] *= d / module;
    out[3][2] *= d / module;


    out[0][0] += at0-> GetVector ().x();
    out[0][1] += at0-> GetVector ().y();
    out[0][2] += at0-> GetVector ().z();

    out[1][0] += at0-> GetVector ().x();
    out[1][1] += at0-> GetVector ().y();
    out[1][2] += at0-> GetVector ().z();


    out[2][0] += at1-> GetVector ().x();
    out[2][1] += at1-> GetVector ().y();
    out[2][2] += at1-> GetVector ().z();

    out[3][0] += at1-> GetVector ().x();
    out[3][1] += at1-> GetVector ().y();
    out[3][2] += at1-> GetVector ().z();

*/



}








void set_color (color c) {
    glColor4f (c.redF (), c.greenF (), c.blueF (), c.alphaF ());
}










