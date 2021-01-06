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

#include "myline.h"
#include "plants_mainwin.h"
#include <qvalidator.h>
#include <Qt3Support/q3vbox.h>
#include <Qt3Support/q3hbox.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qdatetime.h>
#include <Qt3Support/q3buttongroup.h>
#include <qcheckbox.h>
#include <Qt3Support/q3groupbox.h>

#include <qapplication.h>
#include <qcombobox.h>
#include <qpushbutton.h>
#include <Qt3Support/q3grid.h>
#include <Qt3Support/q3frame.h>
#include <qlayout.h>
#include <qfiledialog.h>
#include <iostream>
#include <fstream>
#include "constants.h"
#include <Qt3Support/q3accel.h>
#include <qpainter.h>
#include <qimage.h>
#include <qstyle.h>
#include "graphical_object.h"





MyListBoxItem::MyListBoxItem(string s, Q3ListBox *par)
: Q3ListBoxText(QString(s.c_str()))
{
	valid=check_mol2 (s);
	tex = s;
	parent = par;
	//        setCustomHighlighting( TRUE );
}


bool MyListBoxItem::check_mol2 (string line){
	istringstream iss (line);
	string filename;
	iss >> filename;
	iss >> filename; 
    return true;
//	Molecule *mol = new Molecule ();
	//ifstream ifs (filename.c_str ());
//	bool b = false;
//	mol->readMultiMOL2 (filename, ifs, b);
//	unsigned int asize =mol->atoms.size ();
//	if (asize>1) {
//		return true;
	}    
	//else return false;}



void MyListBoxItem::paint( QPainter *painter )
{
	// evil trick: find out whether we are painted onto our listbox
	//   bool in_list_box = listBox() && listBox()->viewport() == painter->device();

	//    QRect r ( 0, 0, -10+width( listBox() ), height( listBox() ) );
	// if ( in_list_box && isSelected() ) {
	//      painter->fillRect( r, sel_col ); 

	//   }
	if (!valid) {
        QBrush br = parent->palette ().brush (ERR_COL);
		painter->fillRect( 2, 2, -10+width( listBox() ), -4+height( listBox() ), br);
		//       valid=false;
	}

	painter->drawText(5,-4+height( listBox()), QString(tex.c_str()));
	//   if ( in_list_box && isCurrent() )
	//       listBox()->style().drawPrimitive( QStyle::PE_FocusRect, painter, r, listBox()->colorGroup() );
}











PlantsMainWin::PlantsMainWin( Plants *parent)
:    QTabWidget (parent)
{
	window = parent;
	setupwin();




}

void PlantsMainWin::setupwin()
{

	/////////////////////////////////////INPUT OUTPUT////////////////////////////////////////


	Q3VBox *iotab = new Q3VBox( this );
	iotab->setMargin( 5 );

	Q3ButtonGroup *inputgr = new Q3ButtonGroup( 1, Qt::Horizontal, "Input Options", iotab );

	pfile = new MyLineF (inputgr, this, "Protein File", 1);

	connect(pfile->linedit, SIGNAL (textChanged(const QString &)), SLOT (load_protein () ));
	lfile = new MyLineF (inputgr, this, "Ligand File", 1);
	connect(lfile->linedit, SIGNAL (textChanged(const QString &)), SLOT (load_ligand () ));


	llist = new MyLineF (inputgr, this, "Ligand List", 3);
	connect(lfile->linedit, SIGNAL( returnPressed() ), SLOT( add_lf_slot() )) ;
	connect(llist->linedit, SIGNAL( returnPressed() ), SLOT( add_ll_slot() )) ;

	//    connect(lfile->adbutt, SIGNAL( clicked() ), SLOT( add_lf_slot() )) ;
	//   connect(llist->adbutt, SIGNAL( clicked() ), SLOT( add_ll_slot() )) ;
	iflb = new Q3ListBox (inputgr);

	check_ligands ();
	Q3Accel *ifa = new Q3Accel( iflb );        
	ifa->connectItem( ifa->insertItem(Qt::Key_Delete), 
		this,                  
		SLOT(del_ifile()) );

	connect (iflb, SIGNAL(selected (int)), SLOT (ligand_selected(int))); 

	Q3ButtonGroup *bsgr = new Q3ButtonGroup( 1, Qt::Horizontal, "Binding Site options", inputgr );
	bsfl = new QPushButton( "Bindingsite From ligand",bsgr);
	bsfl->setCheckable( FALSE );
	connect (bsfl, SIGNAL (clicked ()), SLOT (load_bs_from_ligand_slot () ));  

	bsx = new MyLine (bsgr, "Binding Site Center X", 2);
	bsy = new MyLine (bsgr, "Binding Site Center Y", 2);
	bsz = new MyLine (bsgr, "Binding Site Center Z", 2);
	bsr = new MyLine (bsgr, "Binding Site Radius", 2);

	connect (bsx->linedit, SIGNAL (textChanged(const QString &)), SLOT (check_bindingsite () ));
	connect (bsy->linedit, SIGNAL (textChanged(const QString &)), SLOT (check_bindingsite () ));
	connect (bsz->linedit, SIGNAL (textChanged(const QString &)), SLOT (check_bindingsite () ));
	connect (bsr->linedit, SIGNAL (textChanged(const QString &)), SLOT (check_bindingsite () ));


	Q3ButtonGroup *outputgr = new Q3ButtonGroup( 1, Qt::Horizontal, "Output Options", iotab );

	odir = new MyLine (outputgr, "Output Dir", 0);
	connect (odir->linedit, SIGNAL (textChanged(const QString &)), SLOT (check_outdir () ));
	Q3Grid *grid = new Q3Grid( 2, outputgr);

	wpc = new QPushButton( "Write Protein Conformation", grid);
	wpc->setCheckable( TRUE );


	wrs = new QPushButton( "Write Rescored Structures", grid);
	wrs->setCheckable( TRUE );


	wmm2 = new QPushButton( "Write Multi Mol2", grid);
	wmm2->setCheckable( TRUE );


	wrl = new QPushButton( "Write Ranking links", grid);
	wrl->setCheckable( TRUE );


	wrmm2 = new QPushButton( "Write Ranking multi Mol2", grid);
	wrmm2->setCheckable( TRUE );


	wpb = new QPushButton( "Write Protein Bindingsite", grid);
	wpb->setCheckable( TRUE );

	wps = new QPushButton( "Write Protein Splitted", grid);
	wps->setCheckable( TRUE );

	wpas = new QPushButton( "Write Per Atom Scores", grid);
	wpas->setCheckable( TRUE );


	addTab( iotab, "Input && Output" );



	////////////////////////////////SEARCH ALGORITHM////////////////////////////////////////////
	Q3VBox *algtab = new Q3VBox( this );
	algtab->setMargin( 5 );


	Q3ButtonGroup *sa = new Q3ButtonGroup( 1, Qt::Horizontal, "Search Algorithm", algtab );

	aants = new MyLine (sa, "Number of Ants ", 1);
	aevap = new MyLine (sa, "Evaporation Factor", 2);
	asigma = new MyLine (sa, "Sigma Scaling Factor", 2);

	Q3Grid *sag = new Q3Grid (3, sa);
	flipab = new QPushButton( "Flip amide bonds", sag);
	flipab->setCheckable( TRUE );


	flipn = new QPushButton( "Flip planar N", sag);
	flipn->setCheckable( TRUE );

	ffbp = new QPushButton( "Force Flipped Bond Planarity", sag);
	ffbp->setCheckable( TRUE );




	Q3ButtonGroup *sf = new Q3ButtonGroup( 1, Qt::Horizontal, "Scoring Function", algtab );

	Q3Grid *sfg = new Q3Grid (2, sf);
	sfchoice = new QComboBox( sfg );
	sfchoice->insertItem(-1,  "chemplp" );
	sfchoice->insertItem(-1,  "plp" );




	cho = new QPushButton( "Weak CHO", sfg);
	cho->setCheckable( TRUE );

	(void) new QLabel("Ligand Internal Clash", sfg );

	lintra = new QComboBox(sfg );
	lintra->insertItem(-1,  "lj" );
	lintra->insertItem(-1,  "clash" );

	rigidl = new QPushButton( "Rigid Ligand", sfg);
	rigidl->setCheckable( TRUE );
	rigidl->setChecked( FALSE );

	rigida = new QPushButton( "Rigid All", sfg);
	rigida->setCheckable( TRUE );
	rigida->setChecked( FALSE );

	Q3ButtonGroup *cpweights = new Q3ButtonGroup( 1, Qt::Horizontal, "ChemPLP weights", sf );


	chbw = new MyLine (cpweights, "Charged HB", 2);
	cmw = new MyLine (cpweights, "Charged Metal", 2);
	hbw = new MyLine (cpweights, "HB", 2);
	hbchow = new MyLine (cpweights, "HB CHO", 2);
	mw = new MyLine (cpweights, "Metal", 2);
	plpw = new MyLine (cpweights, "PLP", 2);
	iw = new MyLine (cpweights, "Intercept", 2);



	Q3ButtonGroup *ca = new Q3ButtonGroup( 1, Qt::Horizontal, "Cluster Algorithm", algtab );

	crmsd = new MyLine (ca, "Cluster RMSD", 2);
	cs = new MyLine (ca, "Cluster Structure", 1); 



	addTab( algtab, "Algorithm" );


	//////////////////////////////////CONSTRAINTS//////////////////////////////////////////




	Q3VBox *contab = new Q3VBox( this );
	contab->setMargin( 5 );
	contab->setSpacing( 5 );

	Q3ButtonGroup *constr = new Q3ButtonGroup( 1, Qt::Horizontal, "Constraints", contab );
	Q3Grid *cgr = new Q3Grid (2,constr);

	addphbcb = new QPushButton("Add Protein HB Constraint", cgr);
	connect (addphbcb, SIGNAL( clicked() ), SLOT( addphbc() ) );
	//  frame // // // // // // // // protein hb // // // // // //
	phbcpopup = new QWidget;

	QVBoxLayout *phbcpopupv = new QVBoxLayout;
    phbcpopup -> setLayout (phbcpopupv);
//	phbcpopup->setMinimumSize (350, 250);   
//	phbcpopup->setMaximumSize (350, 250);
//	phbcpopupv->setMinimumSize (350, 250);   
//	phbcpopupv->setMaximumSize (350, 250);
	QHBoxLayout *anb = new QHBoxLayout;
    phbcpopupv -> addLayout (anb);
	QLabel *lab = new QLabel("Atom number");
    anb -> addWidget (lab);
	anlined = new QComboBox;
    anb ->addWidget (anlined);
	pat.push_back (anlined);
	anlined->setFocus();
	//  connect( anlined, SIGNAL( returnPressed() ), phbcpopup, SLOT( hide() ) );
	//  tmpE2->setGeometry(10,10, 260, 30);
	QHBoxLayout *phbwei = new QHBoxLayout;
    phbcpopupv -> addLayout (phbwei);
	QLabel *lab2 =  new QLabel("Weight");
    phbwei -> addWidget (lab2);
	phbweil = new QLineEdit;
    phbwei -> addWidget (phbweil);
	phbweil->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,phbweil ) );
	QPushButton *setphbcb = new QPushButton("Set constraint");
    phbcpopupv -> addWidget (setphbcb);
	connect ( setphbcb, SIGNAL( clicked() ), SLOT( add_phbc_slot() ) );





	addshcb = new QPushButton("Add Shape Constraint", cgr);
	connect (addshcb, SIGNAL( clicked() ), SLOT( addshc() ) );
	//  frame // // // // // // // // shape // // // // // //
	shcpopup = new QWidget;

	Q3VBox *shcpopupv = new Q3VBox (shcpopup);
	shcpopup->setMinimumSize (350, 250);   
	shcpopup->setMaximumSize (350, 250);
	shcpopupv->setMinimumSize (350, 250);   
	shcpopupv->setMaximumSize (350, 250);
	Q3HBox *shfilehb = new Q3HBox (shcpopupv);
	(void) new QLabel("Shape File", shfilehb );
	shfl = new QLineEdit( shfilehb);
	QPushButton *shfb = new QPushButton("", shfilehb);
	shfb->setMaximumWidth (20);
	shfb->setMinimumWidth (20);
	connect ( shfb, SIGNAL( clicked() ), SLOT( set_sf_slot() ) );
	shfl->setFocus();
	Q3HBox *shwei = new Q3HBox (shcpopupv);
	(void) new QLabel("Weight", shwei );
	shweil = new QLineEdit( shwei);
	shweil->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,shweil ) );
	QPushButton *setshcb = new QPushButton("Set constraint", shcpopupv);
	connect ( setshcb, SIGNAL( clicked() ), SLOT( add_shc_slot() ) );




	addsdcb = new QPushButton("Add Surface Distance Constraint", cgr);
	connect (addsdcb, SIGNAL( clicked() ), SLOT( addsdc() ) );

	//  frame // // // // // // // // surface distance // // // // // //
	sdcpopup = new QWidget;

	Q3VBox *sdcpopupv = new Q3VBox (sdcpopup);
	sdcpopup->setMinimumSize (350, 250);   
	sdcpopup->setMaximumSize (350, 250);
	sdcpopupv->setMinimumSize (350, 250);   
	sdcpopupv->setMaximumSize (350, 250);
	Q3HBox *sdmmb = new Q3HBox (sdcpopupv);

	(void) new QLabel("From", sdmmb );    sdbbfl = new QLineEdit( sdmmb);
	(void) new QLabel("To", sdmmb );    sdbbtl = new QLineEdit( sdmmb);
	sdbbfl->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,sdbbfl ) ); 
	sdbbtl->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,sdbbtl ) );
	Q3HBox *sdanb = new Q3HBox (sdcpopupv);
	(void) new QLabel("Atom Number", sdanb );
	sdanl = new QComboBox(sdanb);
	lat.push_back (sdanl);
	sdanl->setFocus();
	Q3HBox *sdwei = new Q3HBox (sdcpopupv);
	(void) new QLabel("Weight", sdwei );
	sdweil = new QLineEdit( sdwei);
	sdweil->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,sdweil ) );
	QPushButton *setsdb = new QPushButton("Set constraint", sdcpopupv);
	connect ( setsdb, SIGNAL( clicked() ), SLOT( add_sdc_slot() ) );



	addidcb = new QPushButton("Add Intra Distance Constraint", cgr);
	connect (addidcb, SIGNAL( clicked() ), SLOT( addidc() ) );
	//  frame // // // // // // // // intra distance // // // // // //
	idcpopup = new QWidget;

	Q3VBox *idcpopupv = new Q3VBox (idcpopup);
	idcpopup->setMinimumSize (350, 250);   
	idcpopup->setMaximumSize (350, 250);
	idcpopupv->setMinimumSize (350, 250);   
	idcpopupv->setMaximumSize (350, 250);
	Q3HBox *idmmb = new Q3HBox (idcpopupv);

	(void) new QLabel("From", idmmb );    idbbfl = new QLineEdit( idmmb);
	(void) new QLabel("To", idmmb );    idbbtl = new QLineEdit( idmmb);
	idbbfl->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,idbbfl ) ); 
	idbbtl->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,idbbtl ) );
	Q3HBox *idanb = new Q3HBox (idcpopupv);
	(void) new QLabel("Ligand Atom", idanb );
	idanl = new QComboBox( idanb);
	idanl->setFocus();
	lat.push_back (idanl);
	Q3HBox *idanb2 = new Q3HBox (idcpopupv);
	(void) new QLabel("Ligand Atom", idanb2 );
	idanl2 = new QComboBox(idanb2);
	lat.push_back (idanl2);  

	Q3HBox *idwei = new Q3HBox (idcpopupv);
	(void) new QLabel("Weight", idwei );
	idweil = new QLineEdit( idwei);
	idweil->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,idweil ) );
	QPushButton *setidb = new QPushButton("Set constraint", idcpopupv);
	connect ( setidb, SIGNAL( clicked() ), SLOT( add_idc_slot() ) );


	addldcb = new QPushButton("Add Ligand Protein Distance Constraint", cgr);
	connect (addldcb, SIGNAL( clicked() ), SLOT( addldc() ) );
	//  frame // // // // // // // // ligand protein distance // // // // // //
	ldcpopup = new QWidget;
	Q3VBox *ldcpopupv = new Q3VBox (ldcpopup);
	ldcpopup->setMinimumSize (350, 250);   
	ldcpopup->setMaximumSize (350, 250);
	ldcpopupv->setMinimumSize (350, 250);   
	ldcpopupv->setMaximumSize (350, 250);
	Q3HBox *ldmmb = new Q3HBox (ldcpopupv);

	(void) new QLabel("From", ldmmb );   ldbbfl = new QLineEdit( ldmmb);
	(void) new QLabel("To", ldmmb );    ldbbtl = new QLineEdit( ldmmb);
	ldbbfl->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,ldbbfl ) ); 
	ldbbtl->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,ldbbtl ) );
	Q3HBox *ldanb = new Q3HBox (ldcpopupv);
	(void) new QLabel("Ligand Atom", ldanb );
	ldanl = new QComboBox(ldanb);
	lat.push_back (ldanl);

	ldanl->setFocus();
	Q3HBox *ldanb2 = new Q3HBox (ldcpopupv);
	(void) new QLabel("Protein Atom", ldanb2 );
	ldanl2 = new QComboBox( ldanb2);
	pat.push_back (ldanl2);
	Q3HBox *ldwei = new Q3HBox (ldcpopupv);
	(void) new QLabel("Weight", ldwei );
	ldweil = new QLineEdit( ldwei);
	ldweil->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,ldweil ) );
	QPushButton *setldb = new QPushButton("Set constraint", ldcpopupv);
	connect ( setldb, SIGNAL( clicked() ), SLOT( add_ldc_slot() ) );


	constrlb = new Q3ListBox (constr);

	Q3Accel *cna = new Q3Accel( constrlb );        
	cna->connectItem( cna->insertItem(Qt::Key_Delete), 
		this,                  
		SLOT(del_constr()) );

	addTab( contab, "Constraints" );

	Q3VBox *flextab = new Q3VBox( this );
	flextab->setMargin( 5 );
	flextab->setSpacing( 5 );

	Q3ButtonGroup *flex = new Q3ButtonGroup( 1, Qt::Horizontal, "Flexibility", flextab );
	Q3HBox *flexhb = new Q3HBox (flex);     Q3HBox *flexhb2 = new Q3HBox (flex);
	flscl = new QComboBox(  flexhb );
	pres.push_back (flscl);
	addflscb = new QPushButton("Add Flexible side chain", flexhb);
	connect (addflscb, SIGNAL( clicked() ), SLOT( add_flsc_slot() ) );


	fixpl = new QComboBox( flexhb2 );
	addfixpb = new QPushButton("Add Fixed Protein Bond", flexhb2);
	flexlb = new Q3ListBox (flex);

	Q3Accel *fla = new Q3Accel( flexlb );        
	fla->connectItem( fla->insertItem(Qt::Key_Delete), 
		this,                  
		SLOT(del_flex()) );


	connect (addfixpb, SIGNAL( clicked() ), SLOT( add_fixp_slot() ) );

	ipsw = new MyLine (flex, "Intra Protein Score Weight", 2);

	addTab( flextab, "Flexibility" );


	//}
	//////////////////////////////////WATER//////////////////////////////////////////
	//void PlantsMainWin::setupWater()


	//{   


	Q3VBox *watertab = new Q3VBox( this );
	watertab->setMargin( 5 );
	watertab->setSpacing( 5 );  
	Q3HBox *watercoor = new Q3HBox (watertab);
	Q3HBox *waterop = new Q3HBox (watertab);


	(void) new QLabel("X", watercoor );    xwc = new QLineEdit( watercoor);
	(void) new QLabel("Y", watercoor );    ywc = new QLineEdit( watercoor);
	(void) new QLabel("Z", watercoor );    zwc = new QLineEdit( watercoor);
	(void) new QLabel("Radius", waterop );    wr = new QLineEdit( waterop);

	xwc->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,xwc ) ); 
	ywc->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,ywc ) );
	zwc->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,zwc ) );
	wr->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,wr ) );
	addwaterb = new QPushButton("Add Water Molecule", waterop);
	connect ( addwaterb, SIGNAL( clicked() ), SLOT( add_water_slot() ) );

	wref = new MyLineF (watertab, this, "Water Reference", 1);



	waterlb = new Q3ListBox (watertab);
	Q3Accel *wta = new Q3Accel( waterlb );        
	wta->connectItem( wta->insertItem(Qt::Key_Delete), 
		this,                  
		SLOT(del_water()) );


	Q3ButtonGroup *wweights = new Q3ButtonGroup( 1, Qt::Horizontal, "Water molecule weights", watertab );

	wphb = new MyLine (wweights, "Water-Protein HB", 2);
	wlhb = new MyLine (wweights, "Water-Ligand HB", 2);
	wwhb = new MyLine (wweights, "Water-Water HB", 2);
	nwhbp = new MyLine (wweights, "No HB Penalty", 2);


	addTab( watertab, "Water" );



}



////////////////////////////////POP UP MENU////////////////////////////////////////////////




void PlantsMainWin::addphbc(){
	phbcpopup->show();
}

void PlantsMainWin::addshc(){
	shcpopup->show();
}

void PlantsMainWin::addsdc(){
	sdcpopup->show();
}

void PlantsMainWin::addidc(){
	idcpopup->show();
}

void PlantsMainWin::addldc(){
	ldcpopup->show();
}


////////////////////////////////SLOTS////////////////////////////////////////////////


void PlantsMainWin::add_phbc_slot (){
	string an = anlined->currentText ().toStdString ();
	string phbw = phbweil->displayText ().toStdString ();
	if ((an!="") && (phbw!="")){
		istringstream ss (an);
		string anone;
		ss >> anone;
		QString st = "chemplp_protein_hb_constraint ";

		constrlb->insertItem(st + QString(" ") + QString(anone.c_str()) + QString(" ")
			+ QString(phbw.c_str()));
	}
}

void PlantsMainWin::add_shc_slot (){
	string fil = shfl->displayText ().toStdString ();
	string we = shweil->displayText ().toStdString ();
	if (fil!="" && we !=""){
		QString st = "shape_constraint ";
		constrlb->insertItem(st + QString(" ") + QString(fil.c_str()) + QString(" ")
			+ QString(we.c_str()));
	}
}

void PlantsMainWin::add_sdc_slot (){
	string from = sdbbfl->displayText ().toStdString ();
	string to = sdbbtl->displayText ().toStdString ();
	string at = sdanl->currentText ().toStdString ();
	string we = sdweil->displayText ().toStdString ();
	if (from!="" && to!="" && at!="" && we!="") {
		QString st = "surface_distance_constraint ";
		istringstream ss (at);
		string atone;
		ss >> atone;
		constrlb->insertItem(st + QString(" ")
			+ QString(from.c_str()) + QString(" ")
			+ QString(to.c_str()) + QString(" ")
			+ QString(atone.c_str()) + QString(" ")
			+ QString(we.c_str()));
	}
}

void PlantsMainWin::add_idc_slot (){
	string from = idbbfl->displayText ().toStdString ();
	string to = idbbtl->displayText ().toStdString ();
	string la = idanl->currentText ().toStdString ();
	string pa = idanl2->currentText ().toStdString ();
	string we = idweil->displayText ().toStdString ();

	if (from!="" && to !="" && la!="" && pa!="" && we!=""){
		QString st = "ligand_intra_distance_constraint ";
		istringstream ss (la);
		string laone;
		ss >> laone;
		istringstream ss2 (pa);
		string paone;
		ss2 >> paone;

//		constrlb->insertItem(st+" "+from+" "+to+" "+laone+" "+paone+" "+we);
	}
}

void PlantsMainWin::add_ldc_slot (){
	string from = ldbbfl->displayText ().toStdString ();
	string to = ldbbtl->displayText ().toStdString ();
	string la = ldanl->currentText ().toStdString ();
	string pa = ldanl2->currentText ().toStdString ();
	string we = ldweil->displayText ().toStdString ();

	if (from!="" && to !="" && la!="" && pa!="" && we!=""){
		QString st = "protein_ligand_distance_constraint ";
		istringstream ss (la);
		string laone;
		ss >> laone;
		istringstream ss2 (pa);
		string paone;
		ss2 >> paone;

//		constrlb->insertItem(st+" "+from+" "+to+" "+laone+" "+paone+" "+we);
	}
}

void PlantsMainWin::add_flsc_slot (){
	string res = flscl->currentText ().toStdString ();
	if (res!=""){
		istringstream ss (res);
		string resone;
		ss >> resone;
		QString st = "flexible_protein_side_chain_number ";
		//flexlb->insertItem(st+resone);
	}
}

void PlantsMainWin::add_fixp_slot (){
	QString st = "fixed_protein_bond ";
	st.append (fixpl->currentText());
	//flexlb->insertItem(st);
}

void PlantsMainWin::add_water_slot (){
	string x = xwc->displayText ().toStdString ();
	string y = ywc->displayText ().toStdString ();
	string z = zwc->displayText ().toStdString ();
	string r = wr->displayText ().toStdString ();
	if (x!="" && y!="" && z!="" && r!=""){

		QString st = "water_molecule ";
		//waterlb->insertItem(st+" "+x+" "+y+" "+z+" "+r);
	}
}

void PlantsMainWin::set_sf_slot (){
	QString s = QFileDialog::getOpenFileName(this,	tr("open file dialog"), "",
		tr("Tripos Mol2 File (*.mol2)"));
	shfl->insert(s);

}

void PlantsMainWin::add_lf_slot (){
	QString val = lfile->linedit->displayText ();
	if (val!=""){
		QString st = "ligand_file ";
		//    MyListBoxItem *a = new MyListBoxItem (st+val);
//		iflb->insertItem(new MyListBoxItem (st+val, iflb));
		check_ligands ();
	}
}

void PlantsMainWin::add_ll_slot (){
	QString val = llist->linedit->displayText ();
	if (val!=""){
		QString st = "ligand_list ";
//		iflb->insertItem(new MyListBoxItem (st+val, iflb));
		check_ligands ();
	}
}

void PlantsMainWin::load_bs_from_ligand_slot (){
	QString s = QFileDialog::getOpenFileName(this, tr ("open file"), " ",
		tr ("Tripos Mol2 File (*.mol2)"));
	Molecule *mol = window->data->check_mol2 (s.toStdString());
//	if (mol->valid) {
	//	vect bs = find_mass_center (mol->atoms); ///readd radius calculation
	//	set_bindingsite (bs.x(), bs.y(), bs.z(), 12.f);
//	}
}

void PlantsMainWin::load_protein () {
/*	string pname = pfile->val ();
	vector <string> patlist, preslist;
	if (pname.find (".mol2")!=string::npos) {

		Molecule *prot;
		prot =  window->data->check_mol2 (pname);
		if (prot->valid) {
			//       window->data->protein = prot;
			//      window->data->ddwin->load_protein (prot);

			for (unsigned int i=0;i<prot->atoms.size ();i++){
				int n = prot->atoms[i]->ID;
				string name = prot->atoms[i]->name;
				string ssname =prot->atoms[i]->substructureName;
				int ssnum =prot->atoms[i]->substructureNumber;
				stringstream s, ssn;
				s << n; ssn << ssnum;
				patlist.push_back ( s.str()+" "+name+"  "+ssname+ssn.str());
				pfile->linedit->setBackgroundRole (base);
			}
			preslist = prot->res_strings;

		}        
		else {
			pfile->linedit->setBackgroundRole (ERR_COL);
			//    window->data->ddwin->hide_protein ();
		}

	}
	else {
		pfile->linedit->setBackgroundRole (ERR_COL);
		//       window->data->ddwin->hide_protein ();
	}
	update_boxes (pat, patlist);  
	update_boxes (pres, preslist);
*/    
}

void PlantsMainWin::load_ligand () {
/*	string lname = lfile->val ();
	vector <string> latlist;
	if (lname.find (".mol2")!=string::npos) {

		Molecule *lig;
		lig =  window->data->check_mol2 (lname);
//		if (lig->valid) {
			//     window->data->ligand=lig;

			for (unsigned int i=0;i<lig->atoms.size ();i++){
				int n = lig->atoms[i]->ID;
				string name = lig->atoms[i]->name;
				string ssname =lig->atoms[i]->substructureName;
				int ssnum =lig->atoms[i]->substructureNumber;
				stringstream s, ssn;
				s << n; ssn << ssnum;
				latlist.push_back ( s.str()+" "+name+"  "+ssname+ssn.str());
			}
			add_lf_slot ();
			window->data->ligands.push_back (lig);
			//          window->data->ddwin->load_ligand (lig);
	//	}

	}
	update_boxes (lat, latlist);
*/
}

void PlantsMainWin::ligand_selected (int index){
	//   cout <<index <<" "<< window->data->ligands.size ()<<endl;
	//   if ((unsigned int) index < window->data->ligands.size ()){
	//       window->data->ddwin->load_ligand (window->data->ligands[index]);
	//   }
}

void PlantsMainWin::set_bindingsite (float x, float y, float z, float r) {

/*
    stringstream ss, ss2, ss3, ss4;
	ss << x;
	bsx->ins (ss.str ());
	ss2 << y;
	bsy->ins (ss2.str ());
	ss3 << z;
	bsz->ins (ss3.str ());
	ss4 << r;
	bsr->ins (ss4.str ());
*/
}

void PlantsMainWin::set_bindingsite (float x, float y, float z) {
/*
	stringstream ss, ss2, ss3;
	ss << x;
	bsx->ins (ss.str ());
	ss2 << y;
	bsy->ins (ss2.str ());
	ss3 << z;
	bsz->ins (ss3.str ());
*/
}



void PlantsMainWin::check_bindingsite () {
	string x = bsx->val ();
	string y = bsy->val ();
	string z = bsz->val ();
	string r = bsr->val ();
	if (x=="") { bsx->linedit->setBackgroundRole (ERR_COL);}
	else {bsx->linedit->setBackgroundRole (base);}

	if (y=="") { bsy->linedit->setBackgroundRole (ERR_COL);}
	else {bsy->linedit->setBackgroundRole (base);}

	if (z=="") { bsz->linedit->setBackgroundRole (ERR_COL);}
	else {bsz->linedit->setBackgroundRole (base);}

	if (r=="") { bsr->linedit->setBackgroundRole (ERR_COL);}
	else {bsr->linedit->setBackgroundRole (base);}
	if (x!="" && y!="" && z!="" && r!=""){


		//       window->data->ddwin->load_bindingsite (getfloat (x), getfloat (y), getfloat (z), getfloat (r));
		Sphere *sph=window->data->ddwin->new_sphere ("Plants bindingsite");      
		sph->set_center (getfloat (x), getfloat (y), getfloat (z));
		sph->set_radius (getfloat (r));
		sph->set_color (1.0f, 1.0f, 1.0f, 0.4f);
		sph->render_as_surface (); 
		//      window->data->ddwin->gl->updateGL ();  
	}
	else {
		//       window->data->ddwin->hide_bindingsite ();    
	}


}
void PlantsMainWin::check_ligands () {

	int n = iflb->numRows ();

	if (n  <1) {
		iflb->setBackgroundRole (ERR_COL);
//		unsetPalette ();

	} 
	else {

		iflb->setBackgroundRole (QPalette::Base);
		//   for (unsigned int i=0;i<iflb->count();i++) {
		//      cout << iflb->item(iflb->count()-1)->valid;

		//     if (iflb->count() == 2) {it->valid =true;}



	}
}

void PlantsMainWin::check_outdir () {

	string outdir = odir->val ();

	if (outdir=="") {
		odir->linedit->setBackgroundRole (ERR_COL);

	} 
	else {

		odir->linedit->setBackgroundRole (QPalette::Base);

	}
}

void PlantsMainWin::del_ifile () {
	iflb->removeItem (iflb->currentItem ());
	if (iflb->currentItem ()<window->data->ligands.size ()) 
		window->data->ligands.erase (window->data->ligands.begin()+iflb->currentItem ());
	check_ligands ();
}

void PlantsMainWin::del_constr () {
	constrlb->removeItem (constrlb->currentItem ());
}

void PlantsMainWin::del_flex () {
	flexlb->removeItem (flexlb->currentItem ());
}

void PlantsMainWin::del_water () {
	waterlb->removeItem (waterlb->currentItem ());
}


///////////////////////////////////////////////////////////////////////////////////////////
void PlantsMainWin::update_boxes (vector<QComboBox*> &vec, vector<string>& values){
	for (unsigned int i=0;i<vec.size (); i++){
		vec[i]->clear ();
		for (unsigned int j=0; j<values.size (); j++){
			vec[i]->insertItem (-1, QString(values[j].c_str()));


		}

	}
}


float PlantsMainWin::getfloat (string s) {
	istringstream iss (s);
	float f;
	iss>>f;
	return f;
}
