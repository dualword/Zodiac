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

#include "constants.h"
#include <qmenubar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include "ddwin.h"
#include <QtOpenGL/qgl.h>
#include <qlabel.h>
#include <Qt3Support/q3grid.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qevent.h>
#include <qdrag.h>
#include <qdir.h>
#include <iostream>
#include <qcolordialog.h>
#include "database.h"

#ifdef WIN32
#include <float.h>
#define isnan _isnan
#endif // WIN32


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


    QGLFormat f;
    f.setStereo (true);

//    QGLFormat::setDefaultFormat (f);
    gl = new MyGl (this);
    assert (gl);

    setCentralWidget (gl);



    undo_view = new QUndoView ();
    undo_view -> setStack (data -> undo_stack);


        Molecule *all = new Molecule;
        assert (all);
        all -> SetTitle ("all");
        molecules.push_back (all);
        target_molecule = all;
        current_target = 0;
        backbone_shown = true;
        haptic_mode = false;
        minimizing = false;

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

string DDWin::write_mol2_string (Molecule *mol) {
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
        Bond *bo = mol->bonds[i];
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




    QString mol_name = QFileDialog::getOpenFileName(this, 
                    tr ("Open file"), last_visited_dir, tr("All Files (*)"));
	if (!mol_name.isEmpty ()) {
		last_visited_dir = QDir (mol_name).path ();
		load_file (mol_name.toStdString());
	}
}

void DDWin::set_lists (Molecule *mol) {
    int bl1, bl2, ll, sl;
	int zero = gl ->new_list (); 

    bl1 = gl->new_list ();
	bl2 = gl->new_list ();
    ll = gl->new_list ();
    sl = gl->new_list ();

    if (mol->multi) {
        Database_molecule *dm;
        dm = (Database_molecule *) mol;
        for (unsigned int i=0; i<dm->database->molecules.size (); i++) {
            dm->database->molecules[i]->backbone_list1 = bl1;
            dm->database->molecules[i]->backbone_list2 = bl2;
            dm->database->molecules[i]->line_list = ll;
            dm->database->molecules[i]->stick_list = sl;
        }
    }
    else {
        mol->backbone_list1 = bl1;
        mol->backbone_list2 = bl2;
        mol->line_list= ll;
        mol->stick_list= sl;
    }

}

void DDWin::select (Atom *at) {
 //   deselect ();
    Selection *sel = new Selection ();
    set_selected (at, true);
    sel -> ZNAddAtom (at);
    sel -> molecules.push_back ((Molecule *) at-> GetParent ());
    molecules.push_back (sel);
    target->insertItem (1000, QString(sel -> GetTitle ()));    
    set_current_target (-1);
    gl->draw_molecule (target_molecule);
}

void DDWin::load_file (string mol_name) {
    ifstream ifs(mol_name.c_str ());
    OBConversion conv(&ifs);
    OBFormat* inFormat = conv.FormatFromExt(mol_name.c_str ());
    Molecule *mol = new Molecule ();
    if(conv.SetInFormat(inFormat) && conv.Read(mol))
    { 
        Molecule *mol2 = new Molecule ();
        if (conv.Read(mol2)) {
            load_multi_file (mol_name);
        }
        else {
            mol -> ZNinit ();  
            CreateMoleculeCommand *command = new CreateMoleculeCommand (mol, this);
            execute (command); 
            gl -> set_center_of_view (mol -> center);
            gl -> set_center_of_rotation (mol -> center);
        }


    }
   else {
    cerr << "could not read file "<<mol_name<<endl;
    delete mol;
    }
}


void DDWin::load_multi_file (string mol_name) {
    Database *database = new Database (data);
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
		ss << nn;
		string name = string ("Molecule ")+ss.str ();
		mol -> SetTitle (name);
        mol -> ZNinit ();  
		}
		else delete mol;
    }
    CreateMoleculeCommand *command = new CreateMoleculeCommand (database -> molecules[0], this);
    execute (command);
    gl -> set_center_of_view (database -> molecules[0] -> center);
    gl -> set_center_of_rotation (database -> molecules[0] -> center);
	database -> set_graphics ();
	database -> import_csv (mol_name+".csv");

}


void DDWin::deselect () {
    cerr << "startdeselect"<<endl;
    if (target_molecule->selection) set_current_target (0);
    for (unsigned int i=0; i < molecules.size (); i++) {
        int s = molecules.size ();
        cerr << "deselect1" << endl;
        assert ((s-1-(int) i) < molecules.size ());
        assert ((s-1-(int) i) >= 0);
        if (molecules[s-1-i] -> selection) {
            Selection *sel = (Selection *) molecules[s-1-i];
            sel -> deselect ();
        cerr << "deselect2" << endl;
            target->removeItem (s-1-i);
            molecules.erase (molecules.begin ()+s-1-i);
     //       delete sel; cannot use this because OBMol distructor would delete all atoms in sel


        }
    }
    cerr << "enddeselect"<<endl;
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



void DDWin::add_molecule (Molecule *mol) {
    set_lists (mol);
    molecules.push_back (mol);
    target->insertItem (1000, QString( mol->GetTitle ()));
	targets_updated ();
}


void DDWin::new_molecule (string name) {
    Molecule *mol = new Molecule ();
    mol-> SetTitle (name);
    add_molecule (mol);
}



void DDWin::resizeEvent (QResizeEvent *){
    target->move (width()-400,0);

}

void DDWin::create_menu_actions () {

    openAct = new QAction (tr("&Open File"), this);
    openAct->setShortcut (tr ("Ctrl+O"));
    connect (openAct, SIGNAL (triggered ()), this, SLOT (open_file_slot ()));

    saveasAct = new QAction (tr("Save as..."), this);
    connect (saveasAct, SIGNAL (triggered ()), this, SLOT (save_as_slot ()));

    screenshotAct = new QAction (tr("&Screenshot"), this);
    screenshotAct->setShortcut (tr ("Ctrl+S"));
    connect (screenshotAct, SIGNAL (triggered ()), this, SLOT (screenshot_slot ()));

    raytracedscreenshotAct = new QAction (tr("&Raytraced Screenshot POVRAY"), this);
    screenshotAct->setShortcut (tr ("Ctrl+R"));
    connect (raytracedscreenshotAct, SIGNAL (triggered ()), this, SLOT (raytraced_screenshot_slot ()));

    movieAct = new QAction (tr("Start / save &Movie"), this);
    screenshotAct->setShortcut (tr ("Ctrl+M"));
    connect (movieAct, SIGNAL (triggered ()), this, SLOT (movie_slot ()));

    raytracedmovieAct = new QAction (tr("Start / save  Raytraced Movie"), this);
  //  screenshotAct->setShortcut (tr ("Ctrl+O"));
    connect (raytracedmovieAct, SIGNAL (triggered ()), this, SLOT (raytraced_movie_slot ()));

    quitAct = new QAction (tr("&Quit"), this);
    screenshotAct->setShortcut (tr ("Ctrl+Q"));
    connect (quitAct, SIGNAL (triggered ()), this, SLOT (close ()));

    builderAct = new QAction (tr("&Builder"), this);
    builderAct->setShortcut (tr ("Ctrl+B"));
    connect (builderAct, SIGNAL (triggered ()), this, SLOT (builder_slot ()));

    historyAct = new QAction (tr("History"), this);
    connect (historyAct, SIGNAL (triggered ()), this, SLOT (history_slot ()));


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

    DDsettingsAct = new QAction (tr("3D Settings"), this);
    connect (DDsettingsAct, SIGNAL (triggered ()), this, SLOT (DD_settings_slot ()));

    colorAct = new QAction (tr("Color"), this);
    connect (colorAct, SIGNAL (triggered ()), this, SLOT (color_slot ()));

    hapticAct = new QAction (tr("Haptic mode"), this);
	hapticAct -> setIcon (QPixmap (":icons/haptic.png"));
    connect (hapticAct, SIGNAL (triggered ()), this, SLOT (haptic_slot ()));

    dockingAct = new QAction (tr("Docking"), this);
    connect (dockingAct, SIGNAL (triggered ()), this, SLOT (plants_slot ()));

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

    addHsAct = new QAction (tr("Add Hydrogens"), this);
    connect (addHsAct, SIGNAL (triggered ()), this, SLOT (add_Hs_slot ()));

    surfaceAct = new QAction (tr("Molecular Surface"), this);
    connect (surfaceAct, SIGNAL (triggered ()), this, SLOT (surface_slot ()));

    sphereAct = new QAction (tr("Sphere"), this);
    connect (sphereAct, SIGNAL (triggered ()), this, SLOT (sphere_slot ()));

    graphicalobjectsAct = new QAction (tr("Graphical Objects"), this);
    connect (graphicalobjectsAct, SIGNAL (triggered ()), this, SLOT (graphical_objects_slot ()));

    aboutAct = new QAction (tr("About"), this);
    connect (aboutAct, SIGNAL (triggered ()), this, SLOT (about_slot ()));


//toolbars



    undoAct = data -> undo_stack -> createUndoAction ( this );
    undoAct -> setIcon (QPixmap (":icons/undo.png"));
    undoAct->setShortcut (tr ("Ctrl+Z"));


    redoAct = data -> undo_stack -> createRedoAction ( this );
    redoAct -> setIcon (QPixmap (":icons/redo.png"));
    redoAct->setShortcut (tr ("Shift+Ctrl+Z"));




}



void DDWin::draw_menu (){


    create_menu_actions ();

    QMenu *file = new QMenu(tr("&File"), this );
    Q_CHECK_PTR( file );
    file -> addAction (openAct);
    file -> addAction (saveasAct);
    file -> addAction (screenshotAct);
    file -> addAction (raytracedscreenshotAct);
  //  file -> addAction (movieAct);
  //  file -> addAction (raytracedmovieAct);
    file -> addAction (quitAct);


    QMenu *edit = new QMenu(tr("&Edit"), this );
    Q_CHECK_PTR( edit );
    edit -> addAction (builderAct);

    edit -> addAction (historyAct);

  //  edit -> addAction (wiimoteAct);
  //  edit -> addAction (wiimote2Act);

    QMenu *show = new QMenu(tr("&Show"), this );
    Q_CHECK_PTR( show );
    show -> addAction (hideHAct);
    show -> addAction (hidenpHAct);
    show -> addAction (hideallAct);
    show -> addAction (showallAct);


    QMenu *display = new QMenu(tr("&Display"), this );
    Q_CHECK_PTR( display );
    display -> addAction (displaysettingsAct);
    display -> addAction (DDsettingsAct);
    display -> addAction (colorAct);



    QMenu *compute = new QMenu(tr("&Compute"), this );
    Q_CHECK_PTR( compute );
    compute -> addAction (hapticAct);
    compute -> addAction (dockingAct);
    compute -> addAction (energyAct);
    compute -> addAction (logPAct);
    compute -> addAction (minimiseAct);
    compute -> addAction (partialQAct);
    compute -> addAction (scoresCharge);


    QMenu *utilities = new QMenu(tr("&Utilities"), this );
    Q_CHECK_PTR( utilities );
    utilities -> addAction ( addHsAct );


    QMenu *graph_objects = new QMenu(tr("&Graphical Objects"), this );
    Q_CHECK_PTR( graph_objects );
    graph_objects -> addAction (surfaceAct);
    graph_objects -> addAction (sphereAct);
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
    target->move (400,0);
    target->setMinimumWidth (400);
    target->insertItem(0, "All" );
	targets_updated ();

    QToolBar *toolbar = new QToolBar ();
    addToolBar (toolbar);
    toolbar -> addAction (undoAct);
    toolbar -> addAction (redoAct);
    toolbar -> addAction (minimiseAct);
	toolbar -> addAction (hapticAct);
   	toolbar -> addWidget (target);
	
 
    connect (target, SIGNAL (activated (int)), this, SLOT (set_current_target_slot (int)));

}

void DDWin::delete_current_molecule () {
    delete_molecule (molecules [current_target]);
}

void DDWin::delete_molecule (Molecule *mol) {
    for (unsigned int i=0; i<molecules.size (); i++) {
        if (molecules[i] == mol) {
            target->removeItem (i);
            molecules.erase (molecules.begin ()+i);
            set_current_target (0);
			targets_updated ();
            delete mol;
            gl->updateGL ();
            break;
        }
    }
}


void DDWin::remove_molecule (Molecule *mol) {
    for (unsigned int i=0; i<molecules.size (); i++) {
        if (molecules[i] == mol) {
            target->removeItem (i);
            molecules.erase (molecules.begin ()+i);
            set_current_target (0);
            gl->updateGL ();
            break;
        }
    }

}


void DDWin::redraw (Molecule *mol) {
    if (mol -> needs_redraw) {
        gl -> draw_molecule (mol);
        mol -> needs_redraw = false;
    }
}

void DDWin::recolor_by_score (Molecule *mol) {
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
        gl -> draw_molecule (mol);
        mol -> needs_redraw = false;
    }
}


void DDWin::end_minimisation () {
	cerr << "end" << endl;
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

void DDWin::save_as_slot () {
    QString s = QFileDialog::getSaveFileName(this, 
											 tr ("Save As"), "",tr("TRIPOS mol2 (*.mol2)"));
	
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



void DDWin::set_current_target_slot (int index) {
    current_target = index;
    target_molecule = molecules [index];
    if (target_molecule->multi) {
        Database_molecule *dm;
        dm = (Database_molecule *) target_molecule;
		dm -> database -> set_graphics ();
    //    browser_menu->target = dm->database;
   //     browser_menu->current_number = dm->number-1;

    }


}

void DDWin::set_current_target (int index) {
    if (index<0) index+=molecules.size ();
    target->setCurrentIndex (index);
    set_current_target_slot (index);

}

void DDWin::set_current_target (Molecule *mol) {
    int index = -1;
    for (unsigned int i=0; i<molecules.size (); i++) {
        if (molecules[i] == mol) index = i;
    }
    
    if (index>0) {
        target->setCurrentIndex (index);
        set_current_target_slot (index);
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


void DDWin::surface_slot () {
    surface_menu -> show ();
    surface_menu -> raise ();
}

void DDWin::sphere_slot () {
    sphere_menu -> show ();
    sphere_menu -> raise ();
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
				
				
void DDWin::plants_slot () {
	Plants *plants=new Plants ( 0, "plants conf", data );
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
		if (!is_db_extended (target_molecule)) {
			data ->actions ->reprotonate (target_molecule);
		}
		else {
			Database_molecule *dm = (Database_molecule *) target_molecule;
			Database *dat = dm -> database;
			data ->actions ->reprotonate (dat);
			
		}
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



void DDWin::history_slot () {
    undo_view -> show ();
    undo_view -> raise ();
}


void DDWin::wiimote_slot () {
    HeadTrackingThread *thread = new HeadTrackingThread (0, this);
    thread -> start ();	
}

void DDWin::wiimote2_slot () {
    WiimoteTrackingThread *thread = new WiimoteTrackingThread (0, this);
    thread -> start ();	
}


void DDWin::DD_settings_slot () {
    DDsettings_menu->show ();
    DDsettings_menu->raise ();
}


void DDWin::display_settings_slot () {
    dsetpopup->show ();
    dsetpopup->raise ();
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
    accel->connectItem( accel->insertItem(Qt::Key_M), this, SLOT(m_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_N), this, SLOT(n_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_O), this, SLOT(o_pressed()) );
    accel->connectItem( accel->insertItem(Qt::Key_S), this, SLOT(s_pressed()) );

}



//key slots

void DDWin::del_pressed () {
    if (mode == MAGIC_PENCIL) {
        builder->set_magic_pencil_atomic_number (-2);    
        gl->updateGL ();
    }
    if (mode == NONE) {
        if (current_target)        delete_current_molecule ();
    }
}

void DDWin::esc_pressed () {
    mode = NONE; 
    mode_string = "";
    gl -> magic_pencil = false;
    builder->last_magic_pencil_atom = NULL;
    gl->updateGL ();

}



void DDWin::a_pressed () {
    cout <<"a pressed"<<endl;
}


void DDWin::b_pressed () {
    mode = BUILDER;
    mode_string = "Builder";
    gl->updateGL ();
}

void DDWin::c_pressed () {
    if (mode == MAGIC_PENCIL) {
        builder->set_magic_pencil_atomic_number (6);    
        gl->updateGL ();
    }
}


void DDWin::m_pressed () {
    if (mode == BUILDER) {
        mode = MAGIC_PENCIL;
        builder->set_magic_pencil_atomic_number (6); 
        gl->updateGL ();
    }
}


void DDWin::n_pressed () {
    if (mode == MAGIC_PENCIL) {
        builder->set_magic_pencil_atomic_number (7);    
        gl->updateGL ();
    }
}



void DDWin::o_pressed () {
    if (mode == MAGIC_PENCIL) {
        builder->set_magic_pencil_atomic_number (8);    
        gl->updateGL ();
    }
}


void DDWin::s_pressed () {
    mode = SELECT;
    mode_string = "Select";
    gl->updateGL ();
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
   :    QGLWidget(parent)
{


    setMouseTracking (true);
    ddwin = parent;
	setAcceptDrops(true) ; 
    init_vars ();
 //   time = QTime::currentTime ();
 //   startTimer( TIMER_INTERVAL );
    void    GL ();


}

void MyGl::init_vars (){
   
    needs_GL_update = false;
    next_list = 0;


//    bindingsite_list = new_list ();
//    waters_list = new_list ();

    select_color  = SELECT_COLOR ;
    background_color  = BACKGROUND_COLOR ;


    bindingsite_color  = BINDINGSITE_COLOR ;


    water_color  = WATER_COLOR ;

    current_color = CURRENT_COLOR;


    sphere_radius = SPHERE_RADIUS;
    stick_rad = STICK_RAD;
    aromatic_display_style = AROMATIC_DISPLAY_STYLE;

    double_bond_inter_distance =DOUBLE_BOND_INTER_DISTANCE;
    aromatic_bond_inter_distance =AROMATIC_BOND_INTER_DISTANCE;

    
    vdw_scale =VDW_SCALE;
    vdw_precision = VDW_PRECISION;
    stick_precision = STICK_PRECISION;
    sphere_precision = SPHERE_PRECISION;
    double_bond_stick_radius_scale = DOUBLE_BOND_STICK_RADIUS_SCALE;


    surface_resolution = SURFACE_RESOLUTION;

    fog_begin = FOG_BEGIN;

    stereo_toe_in_angle = 1.14576f;
    stereo_inter_eye_semi_distance = 0.2f;


    head_tracking_x_position = 0.f;
    head_tracking_y_position = 0.f;

    center_of_rotation = vect ();
    center_of_view = vect ();

    zbeginy = 0.0f;
    tbeginx = 0.0f;
    tbeginy = 0.0f;


	up_point = vect (0., 1., 0.);
	down_point = vect (0., -1., 0.);
	left_point = vect (-1., 0., 0.);
	right_point = vect (1., 0., 0.);
	
    zoom = false; translate = false; rotate = false; select = false; move = false; spin = false; selection_square = false; magic_pencil = false;
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










void MyGl::move_molecule (Molecule *mol, float x, float y, float z, bool cut_bool) {
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

	glClearColor (background_color.redF (), background_color.greenF (), background_color.blueF (), background_color.alphaF ());							// Black Background
	glClearDepth (1.0f);											// Depth Buffer Setup
	glDepthFunc (GL_LEQUAL);										// The Type Of Depth Testing (Less || Equal)

  //  g_openGlStereoAvailable = GL_TRUE;    
    if (ddwin->g_stereoMode == DDWin::StereoMode_ViaOpenGL) {glEnable (GL_STEREO); cerr <<"stereo enabled"<<endl;}

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
    glEnable(GL_POINT_SMOOTH);                                       //antialiasing
    glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
//    glEnable(GL_POLYGON_SMOOTH);                                      
//	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	
	glShadeModel(GL_SMOOTH);	
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,specular);
    glShadeModel(GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);	
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);			// Enable Alpha Blending (disable alpha testing)
	glEnable(GL_BLEND);	


    GLfloat fogColor[4]= {1.0f, 1.0f, 1.0f, 1.0f};
  //  glFogi(GL_FOG_MODE, GL_LINEAR);		
    glFogfv(GL_FOG_COLOR, fogColor);			
    glFogf(GL_FOG_DENSITY, 0.005f);				
    glHint(GL_FOG_HINT, GL_DONT_CARE);			
 //   glFogf(GL_FOG_START, fog_begin);			
 //   glFogf(GL_FOG_END, 10000.0f);				
   // glEnable(GL_FOG);				

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
        Molecule * min_mol = ddwin->data->minimize->interaction_ff->target_mol;
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
  //  cout <<"paintGL"<<endl;

//    GLdouble camerax, cameray, cameraz;
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
		renderText (20, 20, QString(ddwin->mode_string.c_str()));
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

		for (unsigned int i =1; i<ddwin->molecules.size(); i++) {
			if (!ddwin->molecules[i]->selection) {  //selections are not drawn
                GLdouble cx, cy, cz;
			    gluProject (ddwin -> molecules[i] -> center.x(), ddwin -> molecules[i] -> center.y(), ddwin -> molecules[i] -> center.z(),model, proj, viewport,&cx, &cy, &cz);
 
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
				glCallList (ddwin->molecules[i]->line_list);
				glCallList (ddwin->molecules[i]->backbone_list1); 
				glEnable(GL_LIGHTING);
				glCallList (ddwin->molecules[i]->stick_list);
				glCallList (ddwin->molecules[i]->backbone_list2);
				if (ddwin->show_labels && 0) {
                    FOR_ATOMS_OF_MOL(at, ddwin -> molecules[i]) {
						glColor4f (0.,0.,0.,1.);
						stringstream ss, ss2;

					}
				}
			}
		}

		for (unsigned int i =0; i<ddwin->graphical_objects.size(); i++) {
			if (ddwin->graphical_objects[i]->mesh) glDisable (GL_LIGHTING);
			glCallList (ddwin->graphical_objects[i]->lst);
			glEnable(GL_LIGHTING);
		}
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




        Transform.M[12] = center_of_view.x() -center_of_rotation.x();
        Transform.M[13] = center_of_view.y()-center_of_rotation.y();
        Transform.M[14] = center_of_view.z()-15-center_of_rotation.z();

	glMultMatrixf(Transform.M);


    GLint viewport [4];
    GLdouble model [16];
    GLdouble proj [16]; 
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
    GLdouble xxu, yyu, zzu, xxc, yyc, zzc, xxr, yyr, zzr;
    gluUnProject (viewport[2]/2, viewport[3], 0, model, proj, viewport, &xxu, &yyu, &zzu);
    gluUnProject (viewport[2]/2, viewport[3]/2, 0, model, proj, viewport, &xxc, &yyc, &zzc);
    gluUnProject (viewport[2], viewport[3]/2, 0, model, proj, viewport, &xxr, &yyr, &zzr);


    float newyx = xxu - xxc; 
    float newyy = yyu - yyc;
    float newyz = zzu - zzc;

    float newxx = xxr - xxc; 
    float newxy = yyr - yyc;
    float newxz = zzr - zzc;

		glTranslatef(-center_of_rotation.x(),-center_of_rotation.y(),-center_of_rotation.z());


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


           LastRot = ThisRot;
            

        }
}

void MyGl::mouseMoveEvent ( QMouseEvent * e )
{   

    float x = e->x(); float y = e->y();
    lastx = x; lasty = y;

    if (ddwin->builder->last_magic_pencil_atom) {
        updateGL ();
    }


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
    updateGL ();
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
                    draw_molecule ((Molecule *) last -> GetParent ());
                }
                else {
                    GLdouble lx, ly, lz;
                    gluProject (0, 0, 0,model, proj, viewport,&lx, &ly, &lz);
                    gluUnProject (x, viewport[3]-y, lz, model, proj, viewport, &xx, &yy, &zz);
                    vect coord (xx, yy, zz);
                    Atom *at;
                    at = ddwin->builder->new_atom (coord, ddwin->builder->magic_pencil_atomic_number);
                    ddwin->builder->last_magic_pencil_atom = at;
                    draw_molecule ((Molecule *) at-> GetParent ());                
                }
            }
            updateGL ();
        }
    }
    else { //drag


            Atom *a = select_atom (x, y, 10, 10);
            if (a && ddwin->builder->start_magic_pencil_atom) {
                if (a!=ddwin->builder->start_magic_pencil_atom) {
                    Bond *bo = a-> GetBond (ddwin->builder->start_magic_pencil_atom);
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
    updateGL ();

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
    updateGL ();
}


void MyGl::end_selection_square () {

    selection_square = false;
    select_square ();
    updateGL ();
}





void MyGl::begin_zoom (float y )
{
    zoom = true;
    zbeginy = y;
}

void MyGl::continue_zoom (float y )
{
    center_of_view.z() -= (y - zbeginy);
    zbeginy = y;
    updateGL ();
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
    center_of_view.x() +=(x - tbeginx)/10;
    center_of_view.y() -=(y - tbeginy)/10;

    tbeginx = x; tbeginy = y;
    updateGL ();
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
    gluProject (ddwin->target_molecule->center.x(), ddwin->target_molecule->center.y(),ddwin->target_molecule->center.z(),model, proj, viewport,&centx, &centy, &centz);
    gluUnProject (x, viewport[3]-y, centz, model, proj, viewport, &targx, &targy, &targz);

   // upx-=camerax; upy-=cameray; upz-=cameraz;
  //  rightx-=camerax; righty-=cameray; rightz-=cameraz;
  //  float modup = sqrt(upx*upx+upy*upy+upz*upz);
 //   float modright = sqrt(rightx*rightx+righty*righty+rightz*rightz);
 //   upx/=modup;     upy/=modup;     upz/=modup;
 //   rightx/=modright;    righty/=modright;    rightz/=modright;
    
  //  float xmove = (x - mbeginx)/50;
 //   float ymove = -(y - mbeginy)/50;
//    cout << upx<<" "<<upy<<endl;

//    float xx = xmove*rightx+ymove*upx;
//    float yy = xmove*righty+ymove*upy;
//    float zz = xmove*rightz+ymove*upz;

    if (ddwin->haptic_mode && ddwin->data->minimize->haptic_dof_mode == 0){
        if (ddwin->data->minimize->fragments.size ()) {
            ddwin->data->minimize->fragments[0].translation.x() = targx;
            ddwin->data->minimize->fragments[0].translation.y() = targy;
            ddwin->data->minimize->fragments[0].translation.z() = targz;
        }
    }

    else {
    move_molecule (ddwin->target_molecule, targx-ddwin->target_molecule->center.x(), targy-ddwin->target_molecule->center.y(), targz-ddwin->target_molecule->center.z(), false);
//    move_molecule (ddwin->target_molecule, xx, yy, zz);
    draw_molecule (ddwin->target_molecule);
    }
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
        rotate_around_vector (v, axis, ddwin->target_molecule->center, spin);
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


    updateGL ();
}


void MyGl::select_square () {
    Molecule *selected_mol = NULL;
    if (ddwin -> target_molecule -> selection) {
        Selection *sel = (Selection *) ddwin -> target_molecule; 
        if (sel -> molecules.size () == 1) {
            selected_mol = sel -> molecules[0];
        }

        if (selected_mol) ddwin -> set_current_target (selected_mol);
        else ddwin -> set_current_target (0);
    }
    ddwin -> deselect ();
    Molecule *mol = ddwin -> target_molecule;
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

    cerr << "select1"<<endl;
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
             
    cerr << "select2"<<endl;
    int hits =  glRenderMode (GL_RENDER);


    if (hits){ 
        cerr << hits << endl;
        Selection *sel = new Selection ();
        ddwin->molecules.push_back (sel);
        ddwin->target->insertItem (1000, QString(sel -> GetTitle ())); 
		ddwin-> emit_targets_updated ();   
        ddwin->set_current_target (-1);
        for (unsigned int i=0; i<hits; i++) {
            Atom * at = mol -> GetAtom (buffer[i*4+3]); 
            assert (at);
            set_selected (at, true);
            sel -> ZNAddAtomToSelection (at);
        }

        FOR_BONDS_OF_MOL(b, mol) {   
 //            if (get_selected (&*b)) sel -> ZNAddBondToSelection (&*b);
        }
        sel->add_mol (mol);
        sel->find_limits ();
        sel->find_center ();

    }
    cerr << "draw mol "<< mol -> GetTitle ()<<endl;
    draw_molecule (mol);
    cerr << "finished draw mol "<< mol -> GetTitle ()<<endl;
    cerr << "select 3"<<endl;
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
        clicked_atom=a;
        ddwin->show_atom_properties (clicked_atom);
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

    Bond *this_bond = ring->bonds[0];
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

    Bond *this_bond = ring->bonds[0];
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






void MyGl::draw_bond_line (Bond* bond)
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

void MyGl::draw_aromatic_bond_line (Bond* bond) {
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

void MyGl::draw_triple_bond_line (Bond * bond) {
    draw_bond_line (bond);
    draw_double_bond_line (bond);
}



void MyGl::draw_double_bond_line (Bond* bond) {
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



void MyGl::draw_bond_stick(Bond* bond)
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

    my_cylinder (stick_rad, stick_rad, height, stick_precision, bond->GetBeginAtom (), bond->GetEndAtom ());

    glPopMatrix ();
    
    Atom *at1p = bond -> GetBeginAtom ();
    Atom *at2p = bond -> GetEndAtom ();


    if (!get_sad (at1p)) {
        set_sad (at1p, true);
        glPushMatrix ();
        setAtomColor(at1);
        glTranslatef (at1c.x(), at1c.y (), at1c.z ());
        gluSphere (quadratic, stick_rad, stick_precision, stick_precision);
        glPopMatrix ();
    }
    if (!get_sad (at2p)) {
        set_sad (at2p, true);
        glPushMatrix ();
        setAtomColor(at2);
        glTranslatef (at2c.x(), at2c.y (), at2c.z ());
        gluSphere (quadratic, stick_rad, stick_precision, stick_precision);
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


void MyGl::draw_double_bond_stick (Bond* bond) {



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

    my_cylinder (stick_rad/2, stick_rad/2, height, stick_precision, bond->GetBeginAtom (),bond->GetEndAtom ());

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

    my_cylinder (stick_rad/2, stick_rad/2, height, stick_precision, bond->GetBeginAtom (),bond->GetEndAtom ());

    glPopMatrix ();
  //  Molecule *mol  = (Molecule *) bond -> GetParent ();
    if (CountBonds (bond->GetBeginAtom ())== 1) {
        glPushMatrix ();
        setAtomColor(bond->GetBeginAtom ());
        glTranslatef (verts[0][0], verts[0][1], verts[0][2]);
        gluSphere (quadratic, stick_rad/2, stick_precision, stick_precision);
        glPopMatrix ();
        glPushMatrix ();
        glTranslatef (verts[1][0], verts[1][1], verts[1][2]);
        gluSphere (quadratic, stick_rad/2, stick_precision, stick_precision);
        glPopMatrix ();
    }
    if (CountBonds (bond->GetEndAtom ())== 1) {
        glPushMatrix ();
        setAtomColor(bond->GetEndAtom ());
        glTranslatef (verts[2][0], verts[2][1], verts[2][2]);
        gluSphere (quadratic, stick_rad/2, stick_precision, stick_precision);
        glPopMatrix ();
        glPushMatrix ();
        glTranslatef (verts[3][0], verts[3][1], verts[3][2]);
        gluSphere (quadratic, stick_rad/2, stick_precision, stick_precision);
        glPopMatrix ();
    }
        
    

}

void MyGl::draw_triple_bond_stick (Bond* bond) {
    draw_bond_stick (bond);
    draw_double_bond_stick (bond);

}


void MyGl::draw_aromatic_bond_stick (Bond* bond) {
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
    gluSphere (quadratic,sphere_radius,sphere_precision,sphere_precision);
    glPopMatrix ();

}

void MyGl::draw_atom_sel_sphere(Atom* atom) {
	vect v = get_coordinates(atom);
 //   glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  
 //   setAtomColor(atom);
    glPushMatrix ();
    glTranslatef (v.x(), v.y(), v.z());
    gluSphere (quadratic,sphere_radius,6,6);
    glPopMatrix ();

}


void MyGl::draw_atom_vdw_sphere(Atom* atom) {
	vect v = get_coordinates(atom);
    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  
    setAtomColor(atom);
    glPushMatrix ();
    glTranslatef (v.x(), v.y(), v.z());
    double vdw = etab.GetVdwRad (atom -> GetAtomicNum ());
    gluSphere (quadratic, vdw, vdw_precision, vdw_precision);
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
    gluSphere (quadratic, vdw*vdw_scale,sphere_precision,sphere_precision);
    glPopMatrix ();

}



void MyGl::draw_molecule (Molecule* mol) {
    if (mol->selection) {
        draw_list ( (Selection *) mol);
    }
    else draw_list (mol);
    updateGL ();
}

void MyGl::draw_backbone (Molecule* mol) {
  /*  glNewList(mol->backbone_list1,GL_COMPILE);
    for (unsigned int i=0; i<mol->residues.size();i++) {
        if (mol->residues[i]->backbone_visible == true) {
            if (mol->residues[i]->backbone_style == BACKBONE_LINE) {
                set_color (mol -> residues [i] -> backbone_color);
                glBegin(GL_LINE_STRIP);
                for (unsigned int j=0; j<mol->residues[i]->backbone_list.size(); j++) {
                    glVertex3f(mol->residues[i]->backbone_list[j].x(), mol->residues[i]->backbone_list[j].y(), mol->residues[i]->backbone_list[j].z());
                }
                glEnd();
            }
        }
    }
    glEndList();
    glNewList(mol->backbone_list2,GL_COMPILE);
    for (unsigned int i=0; i<mol->residues.size();i++) {
        if (mol->residues[i]->backbone_visible == true) {
            if (mol->residues[i]->backbone_style == BACKBONE_STICK) {
                backbone_cylinder (mol->residues[i]);
            }
        }


    }
    glEndList();
*/
}

//void MyGl::draw_surface (Molecule* mol) {
 //   SurfaceMenu *sur_men = new SurfaceMenu (0, "surf", ddwin->data, mol);

//}


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
    sel->find_limits ();
    sel->find_center ();
    for (unsigned int i=0; i < sel->molecules.size (); i++) {
        draw_list (sel->molecules[i]);
    }
}

void MyGl::draw_backbone_list (Molecule *mol) {
	glNewList(mol->backbone_list1,GL_COMPILE);
	FOR_RESIDUES_OF_MOL (r, mol) {
		draw_backbone_line (&*r);
	}
	glEndList ();
	glNewList(mol->backbone_list2,GL_COMPILE);
	FOR_RESIDUES_OF_MOL (r, mol) {
		if (0) draw_backbone_stick(&*r);
	}
	glEndList ();
}

void MyGl::draw_list (Molecule* mol) {
 //   cerr << 0<< endl;
    mol-> find_limits ();
 //   cerr << 0<< endl;
    mol-> find_center ();

  //  cerr << 0<< endl;
    FOR_ATOMS_OF_MOL(a, mol) {
        set_sad (&*a, false);
    }
    glNewList(mol->line_list,GL_COMPILE);

    FOR_BONDS_OF_MOL (b, mol) {
        Atom *at1 = b -> GetBeginAtom ();
        Atom *at2 = b -> GetEndAtom ();
 //   cerr << b -> HasData ("visible");
 //   cerr << "here" << endl;
    bool vis = get_visible (&*b);
    int dsi = get_ds (&*b);

   // cerr << dsi;
    if (vis) {
    //    cerr << "yess"<<endl;
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
 //   cerr <<"3 ";
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
 //   cerr <<"4 ";
    glEndList ();
    glNewList(mol->stick_list,GL_COMPILE);


    FOR_BONDS_OF_MOL (b, mol) {        
        Atom *at1 = b -> GetBeginAtom ();
        Atom *at2 = b -> GetEndAtom ();
 //   cerr << b -> HasData ("visible");
 //   cerr << "here" << endl;
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
  //  cerr <<"5 ";
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
  // cerr <<"6 ";
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



void MyGl::draw_atoms_for_selection (Molecule *mol){    //only works in rendermode select
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
        sca = select_color.alphaF ();

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
void MyGl::set_conf (int conf){

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


    ddwin->at_disp->setCurrentIndex (atom_mode);
    ddwin->bo_disp->setCurrentIndex (bond_mode);

}







void MyGl::hide_hydrogens (Molecule* mol){

    FOR_ATOMS_OF_MOL(a, mol) {
        if (a -> IsHydrogen ()) {
            set_visible (&*a, false);
        }
    }
    draw_molecule (mol);
}


void MyGl::hide_hydrogens (vector <Molecule *> molecules) {
    for (unsigned int n =0; n<molecules.size (); n++ ) hide_hydrogens (molecules[n]);
}

void MyGl::hide_nonpolar_hydrogens (Molecule* mol){
    FOR_ATOMS_OF_MOL(a, mol) {
        if (a -> IsNonPolarHydrogen ()) {
            set_visible (&*a, false);
        }
    }
    draw_molecule (mol);
   
}


void MyGl::hide_nonpolar_hydrogens (vector <Molecule *> molecules) {
    for (unsigned int n =0; n<molecules.size (); n++ ) hide_nonpolar_hydrogens (molecules[n]);
}


void MyGl::show_all_atoms (Molecule* mol){
    FOR_ATOMS_OF_MOL(a, mol) {
            set_visible (&*a, true);
    }
    draw_molecule (mol);
}

void MyGl::show_all_atoms (vector <Molecule *> molecules) {
    for (unsigned int n =0; n<molecules.size (); n++ ) show_all_atoms (molecules[n]);
}

void MyGl::hide_all_atoms (Molecule* mol){
    FOR_ATOMS_OF_MOL(a, mol) {
            set_visible (&*a, false);
    }
    draw_molecule (mol);
}

void MyGl::hide_all_atoms (vector <Molecule *> molecules) {
    for (unsigned int n =0; n<molecules.size (); n++ ) hide_all_atoms (molecules[n]);
}


void MyGl::apply_color_masks (vector <color_mask> masks, Molecule *mol, bool undoable) {
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
                if (q <= ddwin -> data -> score_begin) {
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
    ChangeVectorCommand *command = new ChangeVectorCommand (center_of_rotation, v, this,"Set center of rotation");
    ddwin -> execute (command);
    
}
/*
void MyGl::set_center_of_rotation (float x, float y, float z){
    vect v (x, y, z);
    set_center_of_rotation (v);
}
*/

void MyGl::set_center_of_view (vect v) {
    ChangeVectorCommand *command = new ChangeVectorCommand (center_of_view, v, this, "Set center of view");
    ddwin -> execute (command);

}
/*
void MyGl::set_center_of_view (float x, float y, float z){
    vect v (x, y, z);
    set_center_of_view (v);
 //   cout << "center of view set at: "<<x<<" "<<y<<" "<<z<<endl;
}
*/




void MyGl::screenshot (QString filename){

 GLint viewport [4];
   glGetIntegerv (GL_VIEWPORT, viewport);
   int w = viewport [2]; int h = viewport [3];
//	cerr << w << " " << h << endl;
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

void MyGl::set_wireframe (){
    set_conf (WIREFRAME);
}

void MyGl::set_stick (){
    set_conf (STICK);
}

void MyGl::set_cpk (){
    set_conf (CPK);
}

void MyGl::set_ballandline (){
    set_conf (BALLANDLINE);
}

void MyGl::set_ballandstick (){
    set_conf (BALLANDSTICK);
}


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

	updateGL ();
}

void MyGl::move_camera (float x, float y, float z) {
    center_of_view.x() += x;
    center_of_view.y() += y;
    center_of_view.z() += z;
    updateGL ();
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

	float rot [9], inv [9];
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







vector <vect> MyGl::get_backbone_points (Resid *res) {
	bool prec = false, fol = false;
	vect prec_vect = vect (0, 0,0);
	vect fol_vect = vect (0, 0,0);
	vector <vect> out;
	Molecule *mol = NULL;
	FOR_ATOMS_OF_RESIDUE (a, res) {
		mol = (Molecule *) &*a -> GetParent ();
	}
	int indx =  res -> GetIdx ();
	if (indx > 0) {
		Resid *prec_res = mol -> GetResidue (indx - 1);
		if (res ->GetChainNum () == prec_res ->GetChainNum ()) {
			
			FOR_ATOMS_OF_RESIDUE(a, prec_res) {
				QString atomID = QString(prec_res->GetAtomID(&*a).c_str());
				atomID.trimmed();
				if (atomID == "C") {
					prec = true;
					prec_vect = get_coordinates(&*a);
				}
			}		}		
	}
	if (indx < mol -> NumResidues()-1) {
		Resid *fol_res = mol -> GetResidue (indx +1);
		if (res ->GetChainNum () == fol_res ->GetChainNum ()) {
			
			FOR_ATOMS_OF_RESIDUE(a, fol_res) {
				QString atomID = QString(fol_res->GetAtomID(&*a).c_str());
				atomID.trimmed();
				if (atomID == "N") {
					fol = true;
					fol_vect = get_coordinates(&*a);
				}
			}	
		}		
	}
	
	FOR_ATOMS_OF_RESIDUE(a, res) {
		QString atomID = QString(res->GetAtomID(&*a).c_str());
		atomID.trimmed();
		if (atomID == "N") {
			if (prec) {
				out.push_back(mean (prec_vect, get_coordinates(&*a)));
			}
			out.push_back(get_coordinates(&*a));
		}
	}
	FOR_ATOMS_OF_RESIDUE(a, res) {
		QString atomID = QString(res->GetAtomID(&*a).c_str());
		atomID.trimmed();
		if (atomID == "CA") {
			out.push_back(get_coordinates(&*a));
		}
	}
	FOR_ATOMS_OF_RESIDUE(a, res) {
		QString atomID = QString(res->GetAtomID(&*a).c_str());
		atomID.trimmed();
		if (atomID == "C") {
			out.push_back(get_coordinates(&*a));
			if (fol) {
				out.push_back(mean (fol_vect, get_coordinates(&*a)));
			}
		}
	}
	
	out = smooth_list (out);
	out = smooth_list (out);
	out = smooth_list (out);
	out = smooth_list (out);
	out = smooth_list (out);
	
	
	return out;
}




vector <vect> MyGl::smooth_list (vector <vect> lis){
	vector <vect> out;
	if (lis.size () > 1) {
		vect beg = lis[0];
		vect end = lis[lis.size ()-1];
		out.push_back(beg);
		for (unsigned int i = 1; i<lis.size (); i++) {
			out.push_back (mean (lis [i], lis [i-1]));
		}
		out.push_back (end);
	}
	return out;
}
						   

void MyGl::draw_backbone_stick (Resid *res) {
	vector <vect> points = get_backbone_points (res);
	color c1 = color (0.f, 1.f, 0.f);
	color c2 = color (0.f, 1.f, 0.f);
	float rad1 = 0.2;
	float rad2 = 0.2;
	if (points.size () > 1) {
		for (unsigned int i = 1; i < points.size (); i++) {
			my_cylinder (points[i-1], points [i], rad1, rad2, c1, c2, stick_precision);
		} 
	}

}

void MyGl::draw_backbone_line (Resid *res) {
	vector <vect> points = get_backbone_points (res);
	color c1 = color (0.f, 1.f, 0.f);
	color c2 = color (0.f, 1.f, 0.f);
	float rad1 = 0.2;
	float rad2 = 0.2;
	if (points.size () > 1) {
		for (unsigned int i = 1; i < points.size (); i++) {
			my_line (points[i-1], points [i], c1, c2);
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




void MyGl::compute_double_bond_vertexes (Bond *bond, float out [4][3], float d) {

    bool b0a = false;
    bool b0b = false;
    bool b1a = false;
    bool b1b = false;

    if (!d) d = double_bond_inter_distance/2;
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
        //    cerr << "cb0 "<<cb0.x()<<" "<<cb0.y()<<" "<<cb0.z()<<endl;
            v0a = subtract (veccb0, c0);  //some form of checking which side of the axis we are is needed to avoid crossing double bonds
       //     cerr << "v0a "<<v0a.x()<<" "<<v0a.y()<<" "<<v0a.z()<<endl;
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
        //    cerr << cb0.x()<<" "<< cb0.y()<<" "<< cb0.z()<<" "<< cb1.x()<<" "<< cb1.y()<<" "<< cb1.z()<<endl;
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










