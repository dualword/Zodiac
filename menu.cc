
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

#define ELEMENT 0
#define CHARGE  1
#define SCORE   2
#define COLOR   3

#define FLOAT  32

void DDWin::disp_ok () {
    for (unsigned int i=0;i<at_opts.size();i++) at_opts[i]->set ();
	
    gl->aromatic_display_style = style_str_to_i (ar_disp->currentText().toStdString());
	int at_st = style_str_to_i (at_disp->currentText().toStdString());
	int bo_st = style_str_to_i (bo_disp->currentText().toStdString());
	
    if (current_target) {
        if (!target_molecule -> selection) {
			bool ext = false;
			if (target_molecule -> multi) {
				Database_molecule *dm = (Database_molecule *) target_molecule;
				if (dm -> database -> has_extend_enabled ()) ext = true;
			}
			if (!ext) {
				data -> actions -> change_display_style (target_molecule, at_st, bo_st);
			}
			else {
				Database_molecule *dm = (Database_molecule *) target_molecule;
				Database *dat = dm -> database;
				data -> actions -> change_display_style (dat, at_st, bo_st);
			}
			
			
		}
		else {
			for (unsigned int m=0;m<molecules.size();m++){
				ChangeDisplayStyleCommand *command = new ChangeDisplayStyleCommand (gl);
				Molecule *mol = molecules[m];
				FOR_ATOMS_OF_MOL(a, mol) {
					int new_style = style_str_to_i (at_disp->currentText().toStdString());
					command -> add (&*a, new_style);
				}
				
				FOR_BONDS_OF_MOL (b, mol) {
					int new_style = style_str_to_i (bo_disp->currentText().toStdString());
					command -> add (&*b, new_style);
					
				}
				command -> set_name ();
				data -> ddwin -> execute (command);
				
				gl -> draw_molecule (mol);
				
			}
		}
		
	}
}

void DDWin::set_popups (){

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
    at_opts.push_back (new MyFloatEditLine (bondsvb, "Double Bond inter Distance", gl->double_bond_inter_distance));
    at_opts.push_back (new MyFloatEditLine (bondsvb, "Aromatic Bond inter Distance", gl->aromatic_bond_inter_distance));
    at_opts.push_back (new MyFloatEditLine (atomsvb, "VdW scale", gl->vdw_scale));
//    at_opts.push_back (new MyFloatEditLine (atomsvb, "VdW Precision", gl->vdw_precision));
    at_opts.push_back (new MyFloatEditLine (bondsvb, "Double bond scale", gl->double_bond_stick_radius_scale));
    
    QPushButton *ok = new QPushButton("Ok");
    butt2 -> addWidget (ok);
    connect (ok, SIGNAL( clicked() ), SLOT( disp_ok() ) );





  //  browser_menu = new BrowserMenu (0, this);
    color_menu = new ColorMenu (0, this);
    builder_menu = new BuilderMenu (0, builder);
    haptic_menu = new HapticMenu (0, data->minimize);
    clicked_atom_menu = new Clicked_atomMenu (0, this);
    surface_menu = new SurfaceMenu (0, data);
    sphere_menu = new SphereMenu (0, data);
 //   grid_menu = new GridMenu (0, data);
    graphical_objects_menu = new GraphicalObjectsMenu (0, this);
    DDsettings_menu = new DDSettingsMenu (0, this);
}













void SurfaceMenu::draw_surface () {

    if (data -> ddwin -> current_target) {
        Molecule *mol = data -> ddwin -> target_molecule;
        mesh = (gtype->currentIndex ()==1);
        float a = alpha/100;
        surface = new Surface ();
        surface -> set_molecule (mol);
        surface -> list = data -> ddwin -> gl -> new_list ();
        surface -> name = string ("Surface ") + mol -> GetTitle ();
        surface -> mesh = mesh;


        SurfaceThread *thread = new SurfaceThread (0, surface, data -> ddwin);

        thread -> a = a;
        thread -> res = res;
        thread -> start ();


        connect (thread, SIGNAL (finished ()), this, SLOT (add_surface ()));

    }
}


void SurfaceMenu::add_surface () {
    surface -> render ();
    CreateGraphicalObjectCommand *command = new CreateGraphicalObjectCommand (surface, data -> ddwin);
    data -> ddwin -> execute (command);
}

SurfaceMenu::SurfaceMenu (QWidget *parent, Data *dat )
   :    QWidget(parent)
{
    data = dat;
  //  molecule = mol;
    res = data->ddwin->gl->surface_resolution;


    QVBoxLayout *vbox = new QVBoxLayout (this);
    this->setMinimumSize (500, 250);   
    this->setMaximumSize (500, 250);
    this->setMinimumSize (500, 250);   
    this->setMaximumSize (500, 250);
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

//    QSlider *slider = new QSlider( Horizontal, this, "slider" );
 //   alpha_p = new MyFloatEditLine (vbox, "Opacity", alpha);
//    connect (slider, SIGNAL(valueChanged(int)), alpha_p, SLOT(set_value(int)) );
//    connect ();
    alpha_s = new MySlider (vbox, "Opacity",alpha, 0, 100);
    alpha_s->slider->setValue (100);
    QPushButton *ok = new QPushButton ("Ok");
    vbox -> addWidget (ok);
    connect (ok, SIGNAL (clicked ()) ,SLOT (draw_surface ()));


}

///////////////////////////////////////////////////////////////////////


SphereMenu::SphereMenu (QWidget *parent, Data *dat )
   :    QWidget(parent)
{
    data = dat;
  //  molecule = mol;

    center.null ();
    QVBoxLayout *vbox = new QVBoxLayout (this);

    setWindowTitle("Spheres");

    cent_x = new MyFloatEditLine (vbox, "Center x", center.x());
    cent_y = new MyFloatEditLine (vbox, "Center y", center.y());
    cent_z = new MyFloatEditLine (vbox, "Center z", center.z());
    rad = new MyFloatEditLine (vbox, "Radius", radius);
    alpha_s = new MySlider (vbox, "Opacity",alpha, 0, 100);
    alpha_s->slider->setValue (100);
    MyColorButton *colorbutt = new MyColorButton (vbox, col);


    QPushButton *ok = new QPushButton ("Ok");
    vbox -> addWidget (ok);
    connect (ok, SIGNAL (clicked ()) ,SLOT (draw_sphere ()));
}

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









///////////////////////////////////////////////////////////////////////77
/*
GridMenu::GridMenu (QWidget *parent, Data *dat )
   :    QWidget(parent)
{
    data = dat;
  //  molecule = mol;

    setWindowTitle("Grids");
    QVBoxLayout *vbox = new QVBoxLayout (this);
    res_le = new MyFloatEditLine (vbox, "resolution", resolution);

    QPushButton *ok = new QPushButton ("Ok");
    vbox -> addWidget (ok);
    connect (ok, SIGNAL (clicked ()) ,SLOT (draw_sphere ()));
}

*/



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

BuilderMenu::BuilderMenu (QWidget *parent, Builder *build )
   :    QWidget (parent)
{
    builder = build;


    Q3VBox *vbox = new Q3VBox (this);
    this->setMinimumSize (500, 250);   
    this->setMaximumSize (500, 250);
    this->setMinimumSize (500, 250);   
    this->setMaximumSize (500, 250);
    this->setWindowTitle("Builder");


    QPushButton *C = new QPushButton ("C", vbox);
    connect (C, SIGNAL (clicked ()) ,SLOT (add_C ()));

    QPushButton *N = new QPushButton ("N", vbox);
    connect (N, SIGNAL (clicked ()) ,SLOT (add_N ()));

    QPushButton *O = new QPushButton ("O", vbox);
    connect (O, SIGNAL (clicked ()) ,SLOT (add_O ()));

    QPushButton *single_b = new QPushButton ("-", vbox);
    connect (single_b, SIGNAL (clicked ()) ,SLOT (single_bond ()));

    QPushButton *double_b = new QPushButton ("=", vbox);
    connect (double_b, SIGNAL (clicked ()) ,SLOT (double_bond ()));

    QPushButton *triple_b = new QPushButton ("#", vbox);
    connect (triple_b, SIGNAL (clicked ()) ,SLOT (triple_bond ()));

    QPushButton *no_b = new QPushButton ("|", vbox);
    connect (no_b, SIGNAL (clicked ()) ,SLOT (no_bond ()));


    smiles = new QLineEdit (vbox);
    QPushButton *smilesb = new QPushButton ("Add smiles", vbox);
    connect (smilesb, SIGNAL (clicked ()) ,SLOT (add_smiles ()));
    

 //   QPushButton *ok = new QPushButton ("Ok", vbox);
 //   connect (ok, SIGNAL (clicked ()) ,SLOT (draw_surface ()));

}

void BuilderMenu::add_smiles () {
    string str = smiles -> text ().toStdString ();
    OBConversion conv;
    Molecule *mol = new Molecule ();
    conv.SetInFormat ("SMI");
    conv.ReadString (mol, str);
    mol -> AddHydrogens ();
    CreateMoleculeCommand *command = new CreateMoleculeCommand (mol, builder -> ddwin);
    builder -> ddwin -> execute (command);
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


void BuilderMenu::add_C () {
    builder->add_atom (6);
}


void BuilderMenu::add_N () {
    builder->add_atom (7);
}


void BuilderMenu::add_O () {
    builder->add_atom (8);
}


/*

SphereMenu::SphereMenu (QWidget *parent, const char *name, Data *dat)
   :    Q3VBox(parent, name)
{

    data = dat;
    Q3VBox *vbox = new Q3VBox (this);
    this->setMinimumSize (500, 250);   
    this->setMaximumSize (500, 250);
    this->setMinimumSize (500, 250);   
    this->setMaximumSize (500, 250);
    this->setWindowTitle("Surfaces");



}

*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HapticMenu::HapticMenu (QWidget *parent, Minimize *min )
   :    QTabWidget (parent)
{
    minimize = min;
    this->setWindowTitle("Haptic");

    QWidget *settings_widget = new QWidget (this);
    QVBoxLayout *setmainlayout = new QVBoxLayout;
    settings_widget -> setLayout (setmainlayout);

    dofmode = new QComboBox;
 //   dofmode->insertItem(1,  "6" );
    dofmode->insertItem(2,  "3 N" );
    setmainlayout -> addWidget (dofmode);

    interff = new QComboBox();
    interff->insertItem(1,  "MMFF" );
//    interff->insertItem(2,  "PLP" );
 //   interff->insertItem(3,  "Chemscore" );
    setmainlayout -> addWidget (interff);
    MyCheckBox *automove = new MyCheckBox (setmainlayout, minimize -> haptic_thread -> automove, "automove");
    MyCheckBox *color_by_scor = new MyCheckBox (setmainlayout, minimize -> haptic_thread -> color_by_score, "color by score");

    QHBoxLayout *button_lay = new QHBoxLayout;
    setmainlayout -> addLayout (button_lay);
    QPushButton *Ok_b = new QPushButton ("Start");
    connect (Ok_b, SIGNAL (clicked ()) ,SLOT (Ok ()));
    button_lay -> addWidget (Ok_b);

    QPushButton *end_b = new QPushButton ("Stop");
    connect (end_b, SIGNAL (clicked ()) ,SLOT (end ()));


    button_lay -> addWidget (end_b);

    addTab( settings_widget, "Settings" );
    

    QWidget *energy_widget = new QWidget (this);
    addTab( energy_widget, "Energy" );

  //  total_E = new MyLabelf(energy_widget, "Total Energy", &minimize->total_E);
  //  internal_E = new MyLabelf(energy_widget, "Internal Energy", &minimize->total_internal_E);
  //  interaction_E = new MyLabelf(energy_widget, "Interaction Energy", &minimize->total_interaction_E);

}





void HapticMenu::Ok () {
	if (!minimize -> haptic_thread -> isRunning ()) {
		delete minimize->interaction_ff;
		string iff = interff->currentText ().toStdString ();
		if (iff == "MMFF") minimize->interaction_ff = new MMFF ();
		else if (iff == "PLP") minimize->interaction_ff = new PLP ();
		else if (iff == "Chemscore") minimize->interaction_ff = new Chemscore ();
		else minimize->interaction_ff = new MMFF ();

		string dofm =dofmode->currentText ().toStdString ();
		if (dofm == "6") minimize -> haptic_dof_mode = 0;
		else if (dofm == "3 N") minimize -> haptic_dof_mode = 2;
		else minimize -> haptic_dof_mode = 2;
		minimize -> start_haptic_mode ();
	}

}

void HapticMenu::end () {
    minimize->stop_haptic_mode ();

}

void HapticMenu::update () {
    total_E->update ();
    internal_E->update ();
    interaction_E->update ();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Clicked_atomMenu::Clicked_atomMenu (QWidget *parent, DDWin* ddw) 
   :    QTabWidget (parent)
    {
    ddwin = ddw;

    Q3VBox *atomselpopupv = new Q3VBox (this);
    setMinimumSize (350, 255);   
    setMaximumSize (350, 255);
    setWindowTitle ("Clicked Atom Properties");
   // atomselpopupv->setMinimumSize (350, 255);   
  //  atomselpopupv->setMaximumSize (350, 255);    
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



    QPushButton *setascenterb = new QPushButton ("Set As Center of view", atomselpopupv);
    connect (setascenterb, SIGNAL (clicked ()),SLOT (set_clicked_atom_as_center_of_view ()));
    QPushButton *setasrotcenterb = new QPushButton ("Set As Center of Rotation", atomselpopupv);
    connect (setasrotcenterb, SIGNAL (clicked ()),SLOT (set_clicked_atom_as_center_of_rotation ()));
    addTab( atomselpopupv, "Atom" );


    Q3VBox *resv = new Q3VBox (this);
    Q3Grid *resgrd = new Q3Grid (5,resv);
    (void) new QLabel("Residue Name", resgrd );
    resna = new QLabel("", resgrd );
    (void) new QLabel("      ", resgrd );
    (void) new QLabel("Residue Number", resgrd );
    resnu = new QLabel("", resgrd );
    resna->setFrameStyle( Q3Frame::Panel | Q3Frame::Sunken );
    resnu->setFrameStyle( Q3Frame::Panel | Q3Frame::Sunken );
    addTab( resv, "Residue" );



    Q3VBox *molv = new Q3VBox (this);
    QPushButton *addh = new QPushButton ("Add hydrogens", molv);
    connect (addh, SIGNAL (clicked ()),SLOT (add_Hs ()));
    addTab( molv, "Molecule" );




}


void Clicked_atomMenu::update (){
    Atom *at = clicked_atom;


    set_value (aplid, at -> GetIdx ());
    set_value (aplat, string (etab.GetSymbol (at -> GetAtomicNum ())));
    set_value (aplq, at -> GetPartialCharge ());
    set_value (aplfc, at-> GetFormalCharge ());
    set_value (aplx, at-> GetVector ().x());
    set_value (aply, at-> GetVector ().y());
    set_value (aplz, at-> GetVector ().z());
  //  set_value (resna, at->residue->name);
  //  set_value (resnu, at->residue->number);
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

    ddwin->gl->set_center_of_view ((vect&) clicked_atom-> GetVector ());
}

void Clicked_atomMenu::set_clicked_atom_as_center_of_rotation () {

    ddwin->gl->set_center_of_rotation ((vect&) clicked_atom-> GetVector ());
 //   set_center_of_view (clicked_atom-> GetVector ()[0],clicked_atom-> GetVector ()[1],clicked_atom-> GetVector ()[2]);
}


void Clicked_atomMenu::add_Hs () {
    Molecule *mol = (Molecule *) clicked_atom -> GetParent ();
    mol -> DeleteHydrogens ();
    mol -> AddHydrogens ();
    mol -> ZNinit ();
    ddwin -> gl -> draw_molecule (mol);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////


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
    if (current_number+1 < target->molecules.size ()) {
    current_number++;
    set_mol ();
    }
}

void BrowserMenu::last_slot () {
    current_number = target->molecules.size ()-1;
    set_mol ();
}

void BrowserMenu::set_mol () {
	cerr << "set_mol"<< endl;
    for (unsigned int i=0; i<ddwin->molecules.size (); i++) {
        if (ddwin->molecules[i]->multi) {
            Database_molecule *dm;
            dm = (Database_molecule *) ddwin->molecules[i];
            if (dm->database == target) {
                ddwin->molecules[i] = target->molecules[current_number];
                ddwin->set_current_target (i);
                ddwin->gl->draw_molecule (ddwin->target_molecule);
            }
        }
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////


DatabaseGrid::DatabaseGrid  (QWidget *parent, Database *datb, Data *dat) 
:QWidget (parent)

{
    QVBoxLayout *layout = new QVBoxLayout ();
    setLayout (layout);
	ext_bool = false;
		MyCheckBox *extend = new MyCheckBox (layout, ext_bool, "Extend actions to database");
		tab = new QTableWidget(this);
			tab -> setSortingEnabled (false);
    layout -> addWidget (tab);

	current_number = 0;
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
	

	layout -> addLayout (hbox);

	data = dat;
	db = datb;
	
	add (db);
	//browser -> target = db;
    //browser ->current_number = 0;
	show ();
	connect (tab, SIGNAL (cellDoubleClicked (int, int)), this, SLOT (manage_double_click (int, int)));

	tab -> setSortingEnabled (true);
}

void DatabaseGrid::add_field (DatabaseField *df) {
	tab -> setSortingEnabled (false);
	tab -> insertColumn (tab -> columnCount ());
int 	column = tab -> columnCount () -1;
	QTableWidgetItem *header = new QTableWidgetItem (QString (df -> name.c_str ()));
	tab -> setHorizontalHeaderItem (column, header);
	for (unsigned int i = 0; i < tab -> rowCount (); i++) {
		stringstream ss;
		ss << df -> data [i];
		QTableWidgetItem *item = new QTableWidgetItem (/*QString (ss.str ().c_str ())*/);
		QVariant string_var (tr(ss.str ().c_str ()));
		item ->setData (FLOAT, QVariant (string_var.toDouble ()));
	//	item ->setData (Qt::TextAlignmentRole, Qt::AlignRight); 
		item ->setData (Qt::DisplayRole, string_var.toDouble ());
		tab ->setItem (i, column, item);
	}
	tab -> setSortingEnabled (true);
}

void DatabaseGrid::manage_number_changed (const QString str) {
	istringstream ss (str.toStdString ());
	int i;
	ss >> i;
	set_current_molecule (i-1);
}

void DatabaseGrid::manage_double_click (int r, int c) {
	set_current_molecule (r);
}


void DatabaseGrid::set_current_molecule (int i) {
	if (i>-1 && i< db -> molecules.size ()) {
		deselect_row (current_number);
		stringstream ss;
		ss << i+1;
		le -> setText (QString (ss.str ().c_str ()));
		current_number = i;
		select_row (current_number);
		set_mol ();
	}
}
void DatabaseGrid::add (Database *db) {
	tab -> setSortingEnabled (false);
	tab -> setColumnCount (1);
	tab ->setRowCount (db -> molecules.size ());
	QTableWidgetItem *name = new QTableWidgetItem (QString ("Name"));
	tab ->setHorizontalHeaderItem (0, name);
	for (unsigned int i=0; i< db -> molecules.size (); i++) {
		Database_molecule *mol = db -> molecules [i];
	//	stringstream ss;
	//	ss << mol -> number;
	//	string str ("molecule "+ss.str ());
	//	QTableWidgetItem *header = new QTableWidgetItem (QString (str.c_str ()));
		QTableWidgetItem *num_widget = new QTableWidgetItem ();
		num_widget ->setData (Qt::DisplayRole, mol -> number);
	//	tab -> setVerticalHeaderItem (i, header);
		tab -> setItem (i, 0, num_widget);
	}
	tab -> setSortingEnabled (true);
}

void DatabaseGrid::first_slot () {
    set_current_molecule (0);
}

void DatabaseGrid::prev_slot () {
	set_current_molecule (current_number - 1);

}



void DatabaseGrid::next_slot () {
    set_current_molecule (current_number + 1);
}

void DatabaseGrid::last_slot () {
    set_current_molecule ( db->molecules.size ()-1);

}


void DatabaseGrid::set_mol () {
    for (unsigned int i=0; i<data -> ddwin->molecules.size (); i++) {
        if (data -> ddwin->molecules[i]->multi) {
            Database_molecule *dm;
            dm = (Database_molecule *) data -> ddwin->molecules[i];
            if (dm->database == db) {
                data -> ddwin->molecules[i] = db->molecules[real_index_of_line (current_number)];
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
	QTableWidgetItem *num_it = tab ->item (i, 0);
	string text_st = num_it -> text ().toStdString (); 
	istringstream ss (text_st);
	int num;
	ss >> num;
	return num;
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
    MyColorButton *charge_begin_color = new MyColorButton (charge_begin_layout, ddwin->data->charge_begin_color);
    MyColorButton *charge_end_color = new MyColorButton (charge_end_layout, ddwin->data->charge_end_color);


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
    MyColorButton *select_color_buttons = new MyColorButton (color_layout, ddwin->data->constant_color);



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
    MyCheckBox *stereo = new MyCheckBox (layout, ddwin ->g_stereoAvailable, "enable stereo");

    QPushButton *ok = new QPushButton ("Ok");
    layout -> addWidget (ok);


    connect (ok, SIGNAL (clicked ()), this, SLOT (ok_slot ()));
}

void DDSettingsMenu::ok_slot () {
    inter_eye_distance -> set ();
    focal_point_distance -> set ();
    ddwin -> gl -> stereo_toe_in_angle = 90.f - atan (focal_d / ddwin -> gl -> stereo_inter_eye_semi_distance) * 180 / PI;
    cerr << ddwin -> gl -> stereo_toe_in_angle << " "<<ddwin -> gl -> stereo_inter_eye_semi_distance<< endl;

    ddwin -> gl -> updateGL ();
}









///////////////////////////////POPUPMENUS/////////////////////////////////////////////////////////////

MyFloatEditLine::MyFloatEditLine (QLayout *parent, const char *name, float& var)
    : QWidget(){
    parent -> addWidget (this); 
    QHBoxLayout *layout = new QHBoxLayout ();
    setLayout (layout);
    variable = &var;
    QLabel *label = new QLabel( name);
    layout -> addWidget (label); 
    label->setMaximumWidth( 200 );
    label->setMinimumWidth( 200 );
    linedit = new QLineEdit();
    layout -> addWidget (linedit);
    stringstream s;
    s << *variable;
    linedit->setText (QString(s.str().c_str()));
    connect (linedit, SIGNAL (textChanged ( const QString & )), this, SLOT (set (const QString &)));
}

MyFloatEditLine::MyFloatEditLine (QLayout *parent, const char *name, double& var)
    : QWidget(){
    parent -> addWidget (this); 
    QHBoxLayout *layout = new QHBoxLayout ();
    setLayout (layout);
    variable = (float *) &var;
    QLabel *label = new QLabel( name);
    layout -> addWidget (label); 
    label->setMaximumWidth( 200 );
    label->setMinimumWidth( 200 );
    linedit = new QLineEdit();
    layout -> addWidget (linedit);
    stringstream s;
    s << *variable;
    linedit->setText (QString(s.str().c_str()));
    connect (linedit, SIGNAL (textChanged ( const QString & )), this, SLOT (set (const QString &)));
}






void MyFloatEditLine::set (const QString &st) {
    istringstream iss (st.toStdString());
    float f;
    iss >> f;
    *variable = f;   
}

void  MyFloatEditLine::set ()
{
    istringstream iss (linedit->text().toStdString());
    float f;
    iss >> f;
    *variable = f;   
}

void MyFloatEditLine::set_value (int v) {
    stringstream s ;
    s<<v;
    linedit->setText (QString(s.str().c_str())); 
    set ();    
}


MyLabelf::MyLabelf (QWidget *parent, const char *name, float* var)  
   : Q3HBox( parent, name ) {
    variable = var;
    QLabel *name_l = new QLabel ("name", this);
    string nam = name;
    name_l->setText (QString(nam.c_str()));
    label = new QLabel( name, this ); 
    label->setMaximumWidth( 200 );
    label->setMinimumWidth( 200 );


    update ();


}

void MyLabelf::update () {

    stringstream ss;
    ss << *variable;
    label->setText (QString(ss.str().c_str()));

}

void MyLabelf::set_variable (float *f) {
    variable = f;
}

MySlider::MySlider (QLayout *parent, const char *name,  float& var, int vmin, int vmax) : QWidget(){

    parent -> addWidget (this);
    QHBoxLayout *layout = new QHBoxLayout ();
    setLayout (layout);
    pline = new MyFloatEditLine (layout, name, var);
    layout -> addWidget (pline);
	slider = new QSlider(Qt::Horizontal);
    layout -> addWidget (slider);
    slider->setMinimum ( vmin );
    slider->setMaximum ( vmax );
    connect (slider, SIGNAL(valueChanged(int)), pline, SLOT(set_value(int)) );
    connect (pline->linedit, SIGNAL (textChanged(const QString&)), this, SLOT(setValue (const QString&)));
    

}

void MySlider::setValue (const QString& s) {
    int i;
    stringstream ss (s.toStdString());
    ss>>i;
    slider->setValue (i);
}







MyColorButton::MyColorButton (QLayout *parent, QColor &col) : QPushButton () {
    parent -> addWidget (this);
    color = &col;
    connect (this, SIGNAL (clicked ()), this, SLOT (my_clicked ()) );

}


MyCompleteColorSettings::MyCompleteColorSettings (QLayout *parent, color &col) {
    button = new MyColorButton (parent, col);
    alpha = (int) col.alphaF () * 100;
    slider = new MySlider (parent, "opacity", alpha, 0, 100);
    connect  (slider -> pline -> linedit, SIGNAL (textChanged(const QString&)), this, SLOT(setAlpha (const QString&)));
}


void MyCompleteColorSettings::setAlpha (const QString &s) {
    stringstream ss (s.toStdString ());
    float v = 0.f;
    ss >> v;
    v /= 100;
    button -> color -> setAlphaF (v);
  //  cerr <<button ->color -> alphaF ()<<endl;
}



void MyColorButton::paintEvent (QPaintEvent *e) {
    QPushButton::paintEvent (e);

    QPainter painter(this);
    int border = 5;
    painter.fillRect (border, border, width ()-2*border, height ()-2*border, *color);

}

void MyColorButton::my_clicked () {
	QColor new_color = QColorDialog::getColor(*color, this );
    if ( new_color.isValid () ) {
    *color = new_color;
    }
}


//////////////////////////////////////////////////////

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

