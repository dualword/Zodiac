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


#ifndef MENU_H
#define MENU_H


#include "constants.h"
#include "ZNdata.h"
#include <qmenubar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qlabel.h>
#include <Qt3Support/q3grid.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qevent.h>
#include <Qt3Support/q3dragobject.h>
#include <qdir.h>
#include <qcolordialog.h>
#include <qtabwidget.h>
#include <qslider.h>
#include "ZNmolecule.h"
#include "qstackedwidget.h"
#include <qlistwidget.h>
#include "command.h"
#include <QRadioButton>
#include "thread.h"
#include <QTableWidget>


#define BIG 999999999

class MySlider; class MyFloatEditLine; class MyIntegerEditLine; class Builder; class Minimize; class DDWin; class MyLabelf; class MyCompleteColorSettings; class MyColorButton; class DatabaseField; class MyComboBox; class MyHideComboBox; class Data; class MyHideCheckBox; class MyGroupBox; class MyTwoColumn; class MyCheckBox; class My3Column; class MyGridColumn; class MyLineFile; class MyLineEdit; class MyListView; class MyPushButton; class MyListButton; class MyTableWidget;
class Thread; class HapticThread;


class ZNWidget : public QWidget {
	Q_OBJECT

public:
	ZNWidget (QWidget *parent = 0);
};



class ZNAdvancedWidget : public ZNWidget {
	Q_OBJECT

public:
	ZNAdvancedWidget (const char *advance="Advanced options", bool show=true);
	QLayout *basic_layout () {return _public ->layout ();};
	QLayout *advanced_layout () {return _private ->layout ();};
//	QWidget *basic_widget () {return(QWidget *) _public;};
//	QWidget *advanced_widget () {return (QWidget *) _private;};

private:

protected:
//	MyHideCheckBox *_hide_check_box;
	MyGroupBox *_group_box;
	ZNWidget *_private, *_public;

	protected slots:
		void hide_advanced_slot();
};



class ZNMenu : public QMainWindow {
      Q_OBJECT

public:
	ZNMenu (QWidget *parent, Data *dat, bool tabs=false, bool advanced = false);
//	~ZNMenu () {~QWidget ();};
	QWidget *main_widget () {return _main_widget;};
	void addMenu (QMenu *menu) {menuBar () ->addMenu (menu);};
	void display () {if (display_requirements_met ()) {show (); raise ();}};

private:

protected:
	Data *_data;
	//QMenuBar *_menu_bar;
	virtual bool display_requirements_met () {return true;};
	bool check_for_a_set_target ();
	QAction *_about_action;
	QWidget *_main_widget;
	QTabWidget *_tabs;
	virtual void add_menu () {};
	void add_help ();

	protected slots:
		virtual void _show_about () {cerr<<"about"<<endl;};
		void _load_bot_fel (string& buffer, const char *reference, MyFloatEditLine *target, const char *sep=" ");
		void _load_bot_iel (string& buffer, const char *reference, MyIntegerEditLine *target, const char *sep=" ");
		void _load_bot_lf (string& buffer, const char *reference, MyLineFile *target);
		void _load_bot_chb (string& buffer, const char *reference, MyCheckBox *target);
		void _load_bot_box (string& buffer, const char *reference, MyGroupBox *target);
		void _load_bot_cb (string& buffer, const char *reference, MyComboBox *target);
		void _load_bot_le (string& buffer, const char *reference, MyLineEdit *target);
};



class SurfaceMenu : public QWidget {
   Q_OBJECT

public:

    SurfaceMenu ( QWidget *parent, Data *dat);
    Data *data;
    Surface *surface;


private:
//    ZNMolecule *molecule;
    QComboBox *stype;
    QComboBox *gtype;
	
	QComboBox *near_to;
    MyFloatEditLine *resolution;
    MyFloatEditLine *alpha_p, *near_to_dle;
    MySlider *alpha_s;
    double res;
    float color [4];
    int alpha;
    int type;
    bool mesh;
	double near_to_f;

private slots:
    void add_surface ();
    void draw_surface ();
	void update_near_to ();

};



class MapMenu : public ZNMenu {
	Q_OBJECT
	public:
	Map *map;
	MapMenu (QWidget *parent, Data *dat);
	MyGroupBox *_site_box;
	MyFloatEditLine *_site_x_el;
	MyFloatEditLine *_site_y_el;
	MyFloatEditLine *_site_z_el;
	MyFloatEditLine *_site_r_el;
	MyFloatEditLine *_threshold_el;
	MyFloatEditLine *_resolution_el;
	MyLineEdit *_name_el;
	MyComboBox *_type_cb;	
	double _x, _y, _z, _r, _threshold, _resolution;
	MyPushButton *_compute_map_bt;
	color _solid_color;
	private:
	private slots:
	bool display_requirements_met () ;
	void compute_map_slot ();
	void add_map ();
};


class SphereMenu : public QWidget {
   Q_OBJECT

public:

    SphereMenu ( QWidget *parent, Data* dat);
    Data *data;

private:
    MyFloatEditLine *cent_x;
    MyFloatEditLine *cent_y;
    MyFloatEditLine *cent_z;
    MyFloatEditLine *rad;
    MySlider *alpha_s;
    color col;

    double x, y, z;
    int alpha;
    double radius;

	private slots:
	    void draw_sphere ();
};


class SequenceMenu : public ZNMenu {
	Q_OBJECT

public:
	SequenceMenu (QWidget *parent, Data *dat);

	public slots:
		void update_tab ();
		void del_from_tab (ZNMolecule *mol);

private:
	QTableWidget *tab;
	int row, column;
	QAction *_residue_indice, *_residue_uid, *_save_sequence;
	QActionGroup *aminoacidGroup;

	private slots:
		void update_graphics (int row, int r, QString residue);
		void clicked_cell (int x, int y);
		void selected_cells ();
		void pressed_cell (QTableWidgetItem *item);
		void add_menu ();
		void add_help ();
		void _residue_indice_slot ();
		void _residue_uid_slot ();
		void _save_sequence_slot ();

};


class DisplayMenu : public ZNMenu {
	Q_OBJECT

public:
	DisplayMenu (QWidget *parent, Data *dat);
	MyTwoColumn *values_tc;
	MyGridColumn *disp_gc, *_dsetbutts_gc;
	QComboBox *at_disp, *bo_disp, *ar_disp;
	MyComboBox *backbone_disp;
	QPushButton *dset_wireframe, *dset_stick, *dset_cpk, *dset_ball_and_stick, *dset_ball_and_line;
	MyPushButton *ok_p;
	MyGroupBox *backbone_box, *display_mode_box, *atoms_box, *bonds_box;

private:
        vector <MyFloatEditLine*> at_opts, bo_opts;
        int style_str_to_i (string style);

	private slots:
	        void disp_ok ();
	        void set_wireframe ();
	        void set_stick ();
	        void set_cpk ();
	        void set_ballandline ();
	        void set_ballandstick ();
	        void set_conf (int conf);
};

class BuilderMenu : public ZNMenu {
	Q_OBJECT

public:
	BuilderMenu (QWidget *parent, Data *dat, Builder *build);
	MyPushButton *_C, *_N, *_O, *_S, *_F, *_H, *_P, *_Cl, *_I, *_Br, *_H_b, *_He_b, *_Li_b, *_Be_b, *_B_b, *_C_b, *_N_b, *_O_b, *_F_b, *_Ne_b, *_Na_b, *_Mg_b, *_Al_b, *_Si_b, *_P_b, *_S_b, *_Cl_b, *_Ar_b, *_K_b, *_Ca_b, *_Sc_b, *_Ti_b, *_V_b, *_Cr_b, *_Mn_b, *_Fe_b, *_Co_b, *_Ni_b, *_Cu_b, *_Zn_b, *_Ga_b, *_Ge_b, *_As_b, *_Se_b, *_Br_b, *_Kr_b, *_Rb_b, *_Sr_b, *_Y_b, *_Zr_b, *_Nb_b, *_Mo_b, *_Tc_b, *_Ru_b, *_Rh_b, *_Pd_b, *_Ag_b, *_Cd_b, *_In_b, *_Sn_b, *_Sb_b, *_Te_b, *_I_b, *_Xe_b, *_Cs_b, *_Ba_b, *_Lu_b, *_Hf_b, *_Ta_b, *_W_b, *_Re_b, *_Os_b, *_Ir_b, *_Pt_b, *_Au_b, *_Hg_b, *_Tl_b, *_Pb_b, *_Bi_b, *_Po_b, *_At_b, *_Rn_b, *_Fr_b, *_Ra_b, *_Lr_b, *_Rf_b, *_Db_b, *_Sg_b, *_Bh_b, *_Hs_b, *_Mt_b, *_Ds_b, *_Rg_b, *_Uub_b, *_Uut_b, *_Uuq_b, *_Uup_b, *_Uuh_b, *_Uus_b, *_Uuo_b, *_La_b, *_Ce_b, *_Pr_b, *_Nd_b, *_Pm_b, *_Sm_b, *_Eu_b, *_Gd_b, *_Tb_b, *_Dy_b, *_Ho_b, *_Er_b, *_Tm_b, *_Yb_b, *_Ac_b, *_Th_b, *_Pa_b, *_U_b, *_Np_b, *_Pu_b, *_Am_b, *_Cm_b, *_Bk_b, *_Cf_b, *_Es_b, *_Fm_b, *_Md_b, *_No_b, *_no_bond_b, *_single_bond_b, *_double_bond_b, *_triple_bond_b, *_smiles_b, *_ok_b, *_benzene_b, *_ring3_b, *_ring4_b, *_ring5_b, *_ring6_b, *_ring7_b, *_ring8_b, *_furan_b, *_furanO_b, *_Ala_b, *_Arg_b, *_Asn_b, *_Asp_b, *_Cys_b, *_Glu_b, *_Cln_b, *_Gly_b, *_His_b, *_Ile_b, *_Leu_b, *_Lys_b, *_Met_b, *_Phe_b, *_Pro_b, *_Ser_b, *_Thr_b, *_Trp_b, *_Tyr_b, *_Val_b, *_CO_b, *_NCO_b, *_COOH_b, *_PO3_b, *_CCd_b, *_CCt_b, *_CN_b, *_SO2_b, *_NO2_b;
	Builder *builder;
	QLineEdit *smiles;
	MyTwoColumn *_smiles_tc;
	MyGridColumn *_periodictable_gc, *_atoms_gc, *_bonds_gc, *_fragments_gc, *_rings_gc, *_basic_gc, *_aminoacids_gc;
	MyGroupBox *_atoms_box, *_rings_box, *_fragments_box, *_bonds_box, *_aminoacids_box, *_nucleotides_box, *_heterocycles_box, *_smiles_box, *_del_box;

	private slots:
		void add_smiles ();
		void add_benzene ();
		void add_mol (string str);
		void add_fragment (string str);
		void add_ring3 ();
		void add_ring4 ();
		void add_ring5 ();
		void add_ring6 ();
		void add_ring7 ();
		void add_ring8 ();
		void add_furan ();
		void add_furanO ();
		void add_CO ();
		void add_NCO ();
		void add_COOH ();
		void add_CCd ();
		void add_CCt ();
		void add_CN ();
		void add_PO3 ();
		void add_NO2 ();
		void add_SO2 ();
		void add_H ();
		void add_He ();
		void add_Li ();
		void add_Be ();
		void add_B ();
		void add_C ();
		void add_N ();
		void add_O ();
		void add_F ();
		void add_Ne ();
		void add_Na ();
		void add_Mg ();
		void add_Al ();
		void add_Si ();
		void add_P ();
		void add_S ();
		void add_Cl ();
		void add_Ar ();
		void add_K ();
		void add_Ca ();
		void add_Sc ();
		void add_Ti ();
		void add_V ();
		void add_Cr ();
		void add_Mn ();
		void add_Fe ();
		void add_Co ();
		void add_Ni ();
		void add_Cu ();
		void add_Zn ();
		void add_Ga ();
		void add_Ge ();
		void add_As ();
		void add_Se ();
		void add_Br ();
		void add_Kr ();
		void add_Rb ();
		void add_Sr ();
		void add_Y ();
		void add_Zr ();
		void add_Nb ();
		void add_Mo ();
		void add_Tc ();
		void add_Ru ();
		void add_Rh ();
		void add_Pd ();
		void add_Ag ();
		void add_Cd ();
		void add_In ();
		void add_Sn ();
		void add_Sb ();
		void add_Te ();
		void add_I ();
		void add_Xe ();
		void add_Cs ();
		void add_Ba ();
		void add_Lu ();
		void add_Hf ();
		void add_Ta ();
		void add_W ();
		void add_Re ();
		void add_Os ();
		void add_Ir ();
		void add_Pt ();
		void add_Au ();
		void add_Hg ();
		void add_Tl ();
		void add_Pb ();
		void add_Bi ();
		void add_Po ();
		void add_At ();
		void add_Rn ();
		void add_Fr ();
		void add_Ra ();
		void add_Lr ();
		void add_Rf ();
		void add_Db ();
		void add_Sg ();
		void add_Bh ();
		void add_Hs ();
		void add_Mt ();
		void add_Ds ();
		void add_Rg ();
		void add_Uub ();
		void add_Uut ();
		void add_Uuq ();
		void add_Uup ();
		void add_Uuh ();
		void add_Uus ();
		void add_Uuo ();
		void add_La ();
		void add_Ce ();
		void add_Pr ();
		void add_Nd ();
		void add_Pm ();
		void add_Sm ();
		void add_Eu ();
		void add_Gd ();
		void add_Tb ();
		void add_Dy ();
		void add_Ho ();
		void add_Er ();
		void add_Tm ();
		void add_Yb ();
		void add_Ac ();
		void add_Th ();
		void add_Pa ();
		void add_U ();
		void add_Np ();
		void add_Pu ();
		void add_Am ();
		void add_Cm ();
		void add_Bk ();
		void add_Cf ();
		void add_Es ();
		void add_Fm ();
		void add_Md ();
		void add_No ();

		void single_bond ();
		void double_bond ();
		void triple_bond ();
		void no_bond ();

	        void add_Ala ();
        	void add_Arg ();
        	void add_Asn ();
        	void add_Asp ();
        	void add_Cys ();
        	void add_Glu ();
        	void add_Cln ();
        	void add_Gly ();
        	void add_His ();
        	void add_Ile ();
        	void add_Leu ();
        	void add_Lys ();
        	void add_Met ();
        	void add_Phe ();
        	void add_Pro ();
        	void add_Ser ();
        	void add_Thr ();
        	void add_Trp ();
        	void add_Tyr ();
        	void add_Val ();

};



class HapticMenu : public ZNMenu {
   Q_OBJECT

public:
	QReadWriteLock *restrain_lock;
	QComboBox *interff;
	QComboBox *dofmode;
	HapticMenu ( QWidget *parent, Data* dat);
	Minimize *minimize;
	HapticThread *haptic_thread;
	MyLabelf *interaction_E;
	MyListView *restrain_list;
	vector <ForceFieldInteraction *> restrains;
	void update_energy ();
	void maybe_save_result ();
	bool automove_b, color_by_score_b, saving;
	double cluster_RMSD, E_tolerance, last_k, last_dist, mult;
	void add_restrain_atom (Atom *at);
	MyFloatEditLine *restrain_dist;
	MyFloatEditLine *restrain_k;
	QLabel *label;
	Atom *last_atom;
	private slots:
	void Ok ();
	void end ();
	void add_restrain ();
	void clear_restrains ();
	void delete_restrain (int);
	void update_k (double);
	void update_dist (double);
	void update_spinboxes (int);

	void user_save_current_pose ();
};

/*

class BrowserMenu : public QWidget {
	Q_OBJECT

public:
	BrowserMenu ( QWidget *parent, DDWin* ddwin);
	int current_number;

	Database *target;
	DDWin *ddwin;

	void set_mol ();
	void set_target (Database *db); 

	private slots:
		void first_slot ();
		void prev_slot ();
		void next_slot ();
	void last_slot ();
};

*/

class Clicked_atomMenu : public ZNMenu {
	Q_OBJECT

public:
	Clicked_atomMenu ( QWidget *parent, Data* dat);

        //Q3Frame *atomselpopup;
        QLabel *aplid, *aplat, *aplq;
        QLineEdit *aplfc, *aplx, *aply, *aplz;
        QLabel *resna, *resnu, *aptype;
		QLabel *idl;
		MyIntegerEditLine *formal_charge_le;

        Atom *clicked_atom;

        void set (Atom *at);
        void update ();

		int formal_charge;
	private slots:
		void set_value (QLabel *lab, float val) ;
		void set_value (QLabel *lab, string val) ;
		void set_value (QLineEdit *lab, float val) ;

	void add_Hs ();
	void set_clicked_atom_as_center_of_view ();
	void set_clicked_atom_as_center_of_rotation ();
};



class ColorSettingsMenu : public QWidget {

public:
	DDWin *ddwin;
	ColorSettingsMenu ( QWidget *parent, DDWin* ddwin);
};



class ColorMenu : public QWidget {
	Q_OBJECT

public:
	DDWin *ddwin;
	ColorMenu ( QWidget *parent, DDWin* ddwin);

	QComboBox *colortype;
	QStackedWidget *options;

private:
	QColor constant_color;
	MyFloatEditLine *score_begin_line, *score_mid_line, *score_end_line, *charge_begin_line, *charge_end_line;

	private slots:
		void ok_slot ();
};

class BackboneColorMenu : public ZNMenu {
	Q_OBJECT
	
public:
	DDWin *ddwin;
	BackboneColorMenu ( QWidget *parent, Data *dat);
	
	QComboBox *colortype;
	QStackedWidget *options;
	
private:
	color constant_color, helix_color, sheet_color, random_color;
	MyFloatEditLine *score_begin_line, *score_mid_line, *score_end_line, *charge_begin_line, *charge_end_line;
	
	private slots:
	void color_ok_slot ();
	void ss_ok_slot ();
};


class GraphicalObjectsColorMenu : public ZNMenu {
	Q_OBJECT
	
public:
	DDWin *ddwin;
	GraphicalObjectsColorMenu ( QWidget *parent, Data *dat);
	
	QComboBox *colortype, *molecule_target_cb, *molecule_target_cb2, *molecule_target_cb3, *maps_cb;
	QStackedWidget *options;
	void set_target (GraphicalObject *go) {target = go;};
private:
	color constant_color;
	GraphicalObject *target;
		double score_begin_f, score_end_f, score_mid_f;
	color score_begin_color, score_end_color, score_mid_color, lipo_color, hb_acc_color, hb_don_color;
	double alpha_mol_distance_d, alpha_value, alpha_multiplier ;
	private slots:
	void ok_slot ();
	void update_mols ();
	void update_maps ();
	void mult_slot ();
	void alpha_slot ();
	void fade_slot ();

	
};



class DDSettingsMenu : public QWidget {
	Q_OBJECT

public:
	DDWin *ddwin;
	DDSettingsMenu ( QWidget *parent, DDWin* ddwin);
	double focal_d;  

private:
	MyFloatEditLine *inter_eye_distance, *focal_point_distance;
	QComboBox *dd_cb;

	private slots:
		void ok_slot ();
};



class GraphicalObjectsMenu : public QWidget {
	Q_OBJECT

public:
	DDWin *ddwin;
	GraphicalObjectsMenu ( QWidget *parent, DDWin* ddwin);

private:
	QListWidget *list;
	void paintEvent ( QPaintEvent * ) ;

	public slots:
		void update_slot ();
		void delete_selected_slot ();
	void show_color_menu (QListWidgetItem *item);
};




class DatabaseGrid : public ZNMenu {
	Q_OBJECT

public:
	DatabaseGrid ( QWidget *parent, Database *db, Data *dat);
	Database *database;
	Data *data;
	QTableWidget *tab;
	//BrowserMenu *browser;

	void add_menu ();
	inline bool has_extend_enabled () {return ext_bool;};
	int real_index_of_line (int);
	void set_database (Database *db);
	void update_graphics ();
	vector <int> selected_columns ();
	
	void set_mol ();

private:
	MyCheckBox *hide_cb;
	void setup_actions ();
	bool ext_bool;
	QLineEdit *le;
	int current_number;

//	void set_no_mol ();

	QAction *sortupAct, *sortdownAct, *newFieldAct, *addTargetMolAct, *loadCsvAct, *mergedatabaseAct, *univocalnamesAct;
	QAction *FiTAct, *calcAct;
	QToolBar *toolbar;
	private slots:
		void manage_hide (bool b);
		void select_row (int r);
		void deselect_row (int r);
		void set_row_color (int r, color c);
		void set_current_molecule (int i);
		void manage_number_changed (const QString str);
		void manage_double_click (int r); 
		void add_target_molecule_slot ();
		void merge_database_slot ();
		void univocal_names_slot ();
		void new_field_slot ();
	void first_slot ();
	void prev_slot ();
	void next_slot ();
	void last_slot ();
	
	void sort_up_slot ();
	void sort_down_slot ();
	
	void FiT_slot ();
	void calc_slot ();
	void load_csv_slot ();
	void update_cell (int row, int column);
};



class IODeviceMenu : public ZNMenu {
	Q_OBJECT

public:
	IODeviceMenu (QWidget *parent, Data *dat);
//	~IODeviceMenu () {~ZNMenu ();};

private:
	QListView *_input_listview, *_output_listview;
	
protected:

	protected slots:
};


class DockingMenu : public ZNMenu {
	Q_OBJECT
	public:
	DockingMenu (QWidget *parent, Data *dat);
	MyGroupBox *_binding_site_box;
				MyFloatEditLine *_site_x_el;
			MyFloatEditLine *_site_y_el;
			MyFloatEditLine *_site_z_el;
			MyFloatEditLine *_site_r_el;
			double _x, _y, _z, _r;
	MyPushButton *_start_docking_bt;
	private:
	bool display_requirements_met ();
	private slots:
		void start_docking_slot ();
};

class ThreadWidget : public ZNWidget {
	Q_OBJECT
	public:
	ThreadWidget (QWidget *parent = 0);
	QLabel *_name, *_action;
	QProgressBar *_progress_bar;
};

class ThreadMenu : public ZNMenu {
	Q_OBJECT
	public:
	ThreadMenu (QWidget *parent, Data *dat);
	void display_thread (int n, Thread *thread);
	void clear (int i = 0);
	private:
	vector <ThreadWidget *> widgets;	
};



///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Preferences menu
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


class PrefMenu : public ZNMenu {
	Q_OBJECT

public:
	PrefMenu (QWidget *parent, Data *dat);

	public slots:


private:

protected:
	MyLineFile *plants_file;
	MyLineFile *gamess_file;
	MyPushButton *_save_preferences_b;
	QLabel *label;

	protected slots:
		void _save_preferences ();



};


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	PLANTS menu FE
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


const string VERSION_PLANTS = "1.1";



class PLANTSMenu : public ZNMenu {
	Q_OBJECT

public:
	PLANTSMenu (QWidget *parent, Data *dat);
//	vector <QComboBox*> lat, pat, pres;
//	OBMol *prot;
	vect d;
	float f;
	string plants_exe;

private:


protected:
	bool runPlants;
	void add_menu ();
	bool check_mol2 (string line);
	bool mol;
// main tab
//	input options
		MyGroupBox *_input_box;
		MyTwoColumn *_input_tc;
		MyLineFile *load_protein_file;
		MyLineFile *load_ligand_file;
		MyLineFile *load_ligand_list_file;
		MyListView *input_file_list;
		MyListView *input_file;
		double _zero_value;
		double _min_radius;
//		binding site
			MyGroupBox *_binding_site_box;
			MyPushButton *_binding_site_from_ligand_b;
			MyTwoColumn *_binding_site_tc;
			MyFloatEditLine *_site_x_el;
			MyFloatEditLine *_site_y_el;
			MyFloatEditLine *_site_z_el;
			MyFloatEditLine *_site_radius_el;
			vector3 coordinates;
			double site_x;
			double site_y;
			double site_z;
//	output options
		MyGroupBox *_output_box;
		MyLineEdit *_output_dir;
		My3Column *_output_tc;
		MyGroupBox *_output_adv_box;
		bool _write_protein_conformations_value;
		MyCheckBox *_write_protein_conformations_chb;
		bool _write_protein_bindingsite_value;
		MyCheckBox *_write_protein_bindingsite_chb;
		bool _write_protein_splitted_value;
		MyCheckBox *_write_protein_splitted_chb;
		bool _write_rescored_structures_value;
		MyCheckBox *_write_rescored_structures_chb;
		bool _write_multi_mol2_value;
		MyCheckBox *_write_multi_mol2_chb;
		bool _write_ranking_links_value;
		MyCheckBox *_write_ranking_links_chb;
		bool _write_ranking_multi_mol2_value;
		MyCheckBox *_write_ranking_multi_mol2_chb;
		bool _write_per_atom_scores_value;
		MyCheckBox *_write_per_atom_scores_chb;
		bool _write_merged_ligand_value;
		MyCheckBox *_write_merged_ligand_chb;
		bool _write_merged_protein_value;
		MyCheckBox *_write_merged_protein_chb;
		bool _write_merged_water_value;
		MyCheckBox *_write_merged_water_chb;
		bool _keep_original_mol2_description_value;
		MyCheckBox *_keep_original_mol2_description_chb;
//		bool _merge_multi_conf_output_value;
		MyGroupBox *_merge_multi_conf_output_box;
		const char *_merge_multi_conf_character_value;
		MyLineEdit *_merge_multi_conf_character;
		int _merge_multi_conf_after_characters_value;
		MyIntegerEditLine *_merge_multi_conf_after_characters_el;
		MyTwoColumn *_merge_multi_conf_output_tc;
// Algorithm tab
//	search options
		MyGroupBox *_search_box;
		MyComboBox *_search_speed_cb;
		My3Column *_search_tc;
		int _aco_ants_value;
		MyIntegerEditLine *_aco_ants_el;
		double _aco_evap_value;
		MyFloatEditLine *_aco_evap_el;
		double _aco_sigma_value;
		MyFloatEditLine *_aco_sigma_el;
		MyTwoColumn *_flip_tc;
		bool _flip_amide_bonds_value;
		MyCheckBox *_flip_amide_bonds_chb;
		bool _flip_planar_n_value;
		MyCheckBox *_flip_planar_n_chb;
		bool _force_flipped_bonds_planarity_value;
		MyCheckBox *_force_flipped_bonds_planarity_chb;
		bool _force_planar_bond_rotation_value;
		MyCheckBox *_force_planar_bond_rotation_chb;
		MyComboBox *_rescore_mode_cb;
//	scoring options
		MyGroupBox *_scoring_box;
		MyComboBox *_algorithm_type_cb;
		double _outside_binding_site_penalty_value;
		MyFloatEditLine *_outside_binding_site_penalty_el;
		MyTwoColumn *_score_inter_tc;
		My3Column *_score_intra_tc;
		MyGroupBox *_search_adv_box;
		MyTwoColumn *_search_adv_tc;
		bool _enable_sulphur_acceptors_value;
		MyCheckBox *_enable_sulphur_acceptors_chb;
		bool _chemplp_clash_include_HH_value;
		MyCheckBox *_chemplp_clash_include_HH_chb;
		double _chemplp_clash_include_14_value;
		MyFloatEditLine *_chemplp_clash_include_14_el;
		MyComboBox *_ligand_intra_score_cb;
		MyGroupBox *_scoring_inter_box;
		MyGroupBox *_scoring_intra_box;
//		plp/plp5
			MyGroupBox *_plp_weights_box;
			My3Column *_plp_weights_tc;
			double _plp_steric_e_value;
			MyFloatEditLine *_plp_steric_e_el;
			double _plp_burpolar_e_value;
			MyFloatEditLine *_plp_burpolar_e_el;
			double _plp_hbond_e_value;
			MyFloatEditLine *_plp_hbond_e_el;
			double _plp_metal_e_value;
			MyFloatEditLine *_plp_metal_e_el;
			double _plp_repulsive_weight_value;
			MyFloatEditLine *_plp_repulsive_weight_el;
			double _plp_tors_weight_value;
			MyFloatEditLine *_plp_tors_weight_el;
//		chemplp
			MyGroupBox *_chemplp_weights_box;
			My3Column *_chemplp_weights_tc;
			bool _chemplp_weak_cho_value;
			MyCheckBox *_chemplp_weak_cho_chb;
			double _chemplp_charged_hb_weight_value;
			MyFloatEditLine *_chemplp_charged_hb_weight_el;
			double _chemplp_charged_metal_weight_value;
			MyFloatEditLine *_chemplp_charged_metal_weight_el;
			double _chemplp_hbond_weight_value;
			MyFloatEditLine *_chemplp_hbond_weight_el;
			double _chemplp_hbond_cho_weight_value;
			MyFloatEditLine *_chemplp_hbond_cho_weight_el;
			double _chemplp_metal_weight_value;
			MyFloatEditLine *_chemplp_metal_weight_el;
			double _chemplp_plp_weight_value;
			MyFloatEditLine *_chemplp_plp_weight_el;
			double _chemplp_plp_steric_e_value;
			MyFloatEditLine *_chemplp_plp_steric_e_el;
			double _chemplp_plp_burpolar_e_value;
			MyFloatEditLine *_chemplp_plp_burpolar_e_el;
			double _chemplp_plp_hbond_e_value;
			MyFloatEditLine *_chemplp_plp_hbond_e_el;
			double _chemplp_plp_metal_e_value;
			MyFloatEditLine *_chemplp_plp_metal_e_el;
			double _chemplp_plp_repulsive_weight_value;
			MyFloatEditLine *_chemplp_plp_repulsive_weight_el;
			double _chemplp_tors_weight_value;
			MyFloatEditLine *_chemplp_tors_weight_el;
			double _chemplp_lipo_weight_value;
			MyFloatEditLine *_chemplp_lipo_weight_el;
			double _chemplp_intercept_weight_value;
			MyFloatEditLine *_chemplp_intercept_weight_el;
	MyTwoColumn *_cluster_docking_tc;
//	cluster options
		MyGroupBox *_cluster_box;
		double _cluster_rmsd_value;
		MyFloatEditLine *_cluster_rmsd_el;
		int _cluster_structures_value;
		MyIntegerEditLine *_cluster_structures_el;
//
	MyGroupBox *_docking_box;
	bool _rigid_ligand_value;
	MyCheckBox *_rigid_ligand_chb;
	bool _rigid_all_value;
	MyCheckBox *_rigid_all_chb;
// Costraints tab
	MyGroupBox *_costraints_box;
	MyGroupBox *_protein_hb_costraint_box;
	MyGroupBox *_shape_costraint_box;
	MyGroupBox *_surface_distance_costraint_box;
	MyGroupBox *_ligand_intra_distance_costraint_box;
	MyGroupBox *_protein_ligand_distance_costraint_box;
	MyTwoColumn *_protein_hb_costraint_tc;
	MyTwoColumn *_shape_costraint_tc;
	MyListView *costraints_file_list;
	MyLineFile *load_shape_file;
	MyTwoColumn *_costraints_tc;
	MyFloatEditLine *_shape_weight;
	MyPushButton *_shape_weight_b;
	const char *title_shape;
	MyComboBox *_protein_hb_costraint_cb;
	MyFloatEditLine *_protein_hb_costraint_fel;
	MyPushButton *_protein_hb_costraint_b;
	MyFloatEditLine *_to_surface_distance;
	MyFloatEditLine *_from_surface_distance;
	MyPushButton *_surface_distance_b;
	MyComboBox *_ligand_surface_distance_cb;
	MyFloatEditLine *_weight_surface_distance_fel;
	MyPushButton *_ligand_intra_distance_b;
	MyPushButton *_protein_ligand_distance_b;
	MyFloatEditLine *_from_ligand_intra_distance;
	MyFloatEditLine *_to_ligand_intra_distance;
	MyComboBox *_number_a_ligand_intra_distance_cb;
	MyComboBox *_number_b_ligand_intra_distance_cb;
	MyFloatEditLine *_weight_ligand_intra_distance;
	MyFloatEditLine *_from_protein_ligand_distance;
	MyFloatEditLine *_to_protein_ligand_distance;
	MyComboBox *_number_a_protein_ligand_distance_cb;
	MyComboBox *_number_b_protein_ligand_distance_cb;
	MyFloatEditLine *_weight_protein_ligand_distance;
	MyTwoColumn *_ligand_intra_distance_costraint_tc;
	MyTwoColumn *_protein_ligand_distance_costraint_tc;
	MyTwoColumn *_surface_distance_costraint_tc;
// Flexibility tab
	MyGroupBox *_flexibility_box;
	MyTwoColumn *_flexibility_tc;
	MyListView *flexibility_file_list;
	double _intra_protein_score_weight_value;
	MyFloatEditLine *_intra_protein_score_weight_el;
	MyListButton *_flexible_protein_side_chain_string;
	MyListButton *_fix_protein_bond;
	MyListButton *_flexible_protein_side_chain_number;
// Water tab
	MyGroupBox *_water_box;
	MyGroupBox *_water_weights_box;
	MyTwoColumn *_water_weights_tc;
	double _water_protein_hb_weight_value;
	MyFloatEditLine *_water_protein_hb_weight_el;
	double _water_ligand_hb_weight_value;
	MyFloatEditLine *_water_ligand_hb_weight_el;
	double _water_water_hb_weight_value;
	MyFloatEditLine *_water_water_hb_weight_el;
	double _no_water_ligand_hb_penalty_value;
	MyFloatEditLine *_no_water_ligand_hb_penalty_el;
	double _water_enable_penalty_value;
	MyPushButton *_water_site_b;
	MyFloatEditLine *_water_enable_penalty_el;
	MyLineFile *load_water_file;
	MyTwoColumn *_water_site_tc;
	MyFloatEditLine *_water_site_x_el;
	MyFloatEditLine *_water_site_y_el;
	MyFloatEditLine *_water_site_z_el;
	MyFloatEditLine *_water_site_radius_el;
	vector3 water_coordinates;
	double _water_min_radius;
	double water_site_x;
	double water_site_y;
	double water_site_z;
//
	QAction *_load_action;
	QAction *_run_action;
	QAction *_save_action;
	QAction *_restore_action;

	string saved_file;

	private slots:
		void plp_weight (int i);
		void set_binding_site ();
		void set_water_site ();
//		void set_binding_site_values (float x, float y, float z, float r);
		void update_boxes (MyListView *parent, const char *mol_name, int type);

	protected slots:
		void hide_group_box2 ();
		void set_protein_hb_constraint_slot();
		void set_ligand_intra_distance_slot ();
		void set_protein_ligand_distance_slot ();
		void set_shape_constraint_slot ();
		void set_surface_distance_slot ();
		void hide_group_box ();
		void _load_slot ();
		void _save_slot ();
		void _run_slot ();
		void _restore_slot ();
		void _show_about ();
		void add_fscs_slot();
		void add_fscn_slot();
		void add_fpb_slot();
		void load_protein ();
		void load_ligand ();
		void load_shape ();
		void load_ligand_list ();
		void _err_mol ();
		void _err_multi ();
		void _load_bot_list_view (string& buffer, const char *reference, MyListView *target);
		void _load_bot_water (string& buffer, const char *reference);
		void _load_bot_bindingsite (string& buffer, const char *reference);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////


class ConformersMenu : public ZNMenu {
	Q_OBJECT

public:
	ConformersMenu (QWidget *parent, Data *dat);
//	~ConformersMenu () {~ZNMenu ();};

private:

	ZNWidget *_stochastic_widget, *_systematic_widget;
	MyHideComboBox *_type_hcb;
protected:
	bool display_requirements_met () {return check_for_a_set_target ();};
	int _results, _timlim;
	MyIntegerEditLine *_results_iel, *_timlim_iel;
protected slots:
	void start_thread ();

};



///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	GAMESS menu FE
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


const string VERSION_GAMESS = "0.1";

const string GAMESS;


class GamessMenu : public ZNMenu {
	Q_OBJECT

public:
	GamessMenu (QWidget *parent, Data *dat);
//	~GamessMenu () {~ZNMenu ();};

private:

protected:
// main tab
	MyGroupBox *_base_box;
	MyComboBox *_run_type_cb;
	QComboBox *_semi_emp_cb;
	QComboBox *_H_F_cb;
	MyComboBox *_mult_cb, *_guess_cb;
	MyComboBox *_scftyp_cb;
	MyHideComboBox *_gbasis_hcb;
	MyGroupBox *_solv_box;
	QRadioButton *radio1;
//	MyComboBox *_prova_cb;
// solvent group
	MyHideCheckBox *_enable_solv;
	MyComboBox *_solvent_cb;
	MyFloatEditLine *_tabs_el;
	QWidget *_hide;
	double _temp;
	MyHideCheckBox *_hide_check_box;
	ZNWidget *_private_solv, *_public_solv;
	int _solv;
// ielpot tab
	MyGroupBox *_elpot_box, *_pdc_box, *_grid_box, *_points_box;
	MyComboBox *_where_cb, *_output_cb, *_ptsel_cb, *_constr_cb;
	MyIntegerEditLine *_layer_iel, *_maxpdc_iel;
	MyFloatEditLine *_vdwscl_fel, *_vdwinc_fel;
	int _layer, _maxpdc;
	double _vdwscl, _vdwinc;
// system tab
	MyGroupBox *_run_box;
	MyComboBox *_exetyp_cb;
	MyIntegerEditLine *_timlim_iel;
	MyIntegerEditLine *_memory_iel;
	int _timlim;
	int _memory;
	int _memory2;
// tabs
//	ZNAdvancedWidget *_elpot_tab;
	QWidget *_system_tab;
	void add_menu (); 
	QAction *_load_action;
	QAction *_test_action;
	QAction *_save_action;
	QAction *_restore_action;
	int _gbasis;
	int _gbasis_load;

	bool display_requirements_met () {return check_for_a_set_target ();};

	private slots:
		void hide_temp (int i);
		void hide_solv (int i);
		void hide_elpot (int i);
		void _load_bot_iel (string& buffer, const char *reference, MyIntegerEditLine *target);
		void _load_bot_fel (string& buffer, const char *reference, MyFloatEditLine *target);
		void _load_bot_box (string& buffer, const char *reference, MyGroupBox *target);
		void _load_bot_cb (string& buffer, const char *reference, MyComboBox *target);
		void _load_bot_gbasis (string& buffer);
		void _load_bot_ngauss (string& buffer);
		void _load_bot_ndfunc (string& buffer);
		void _load_bot_npfunc (string& buffer);
		void _load_bot_diffsp (string& buffer);
		void _load_bot_diffs (string& buffer);

	protected slots:
		void _load_slot ();
		void _save_slot ();
		void _restore_slot ();
		void _show_about ();
};


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	"My" classes
//
///////////////////////////////////////////////////////////////////////////////////////////////////////


class MyCheckBox : public QCheckBox {
	Q_OBJECT

public:
	MyCheckBox (QLayout *parent, bool &var, string string);

private:
	bool *var;

	private slots:
		void get ();
		void set ();
};



class MyColorButton : public QPushButton {
	Q_OBJECT

public:
	MyColorButton (QLayout *parent, QColor &color);

	void paintEvent ( QPaintEvent * ) ;
	QColor *color;

	private slots:
		void my_clicked ();
};



class MyComboBox : public QWidget {
	Q_OBJECT

public:
	MyComboBox (QLayout *parent, const char *name);
	QComboBox *combo_box () {return _combo_box;};
	void insertItem (int index, const char *name);
	QVariant currentData () {return _combo_box ->itemData (_combo_box ->currentIndex ()); };
	QComboBox *_combo_box;

protected:

};



class MyCompleteColorSettings : public QWidget {
    Q_OBJECT

public:
    MyCompleteColorSettings (QLayout *parent, color &col, string name = "");
    MyColorButton *button;
    MySlider *slider;
    int alpha;

	private slots:
	    void setAlpha (int i);
};


class MyFloatEditLine : public QWidget {
	Q_OBJECT

public:
	MyFloatEditLine (QLayout *parent, const char *name, double& var, double min = -BIG, double max = BIG);
	
	void set ();
	QDoubleSpinBox *spinbox;
	double* variable;

	double currentData () {return spinbox ->value ( ); };

	public slots:
		void set_value (double d);
		void set (double d);
};



class MyGroupBox : public QGroupBox {
	Q_OBJECT

public:
	MyGroupBox (QLayout *parent, string str, bool checkbox = false, bool status = false);
	MyGroupBox (string str, bool checkbox = false, bool status = false);
};



class MyHideCheckBox : public QCheckBox {
	Q_OBJECT

public:
	MyHideCheckBox (QLayout *parent, QWidget *tar, string string);

private:
	QWidget *target;
	MyGroupBox *_group_box;
	void set_checked ();
	void set_unchecked ();

	private slots:
	void get ();
	void set ();
	void set_2 ();
};



class MyHideComboBox : public MyComboBox {
	Q_OBJECT

public:
	MyHideComboBox (QLayout *parent, const char *name);
	void insertItem ( QWidget *wid, int index, const char *name, const QVariant &data = QVariant ());

protected:
	vector <QWidget *> _widgets;

	protected slots:
		void set (int i);
};



class MyIntegerEditLine : public QWidget {
	Q_OBJECT

public:
	MyIntegerEditLine (QLayout *parent, const char *name, int& var, int min = -BIG, int max = BIG);

	void set ();
	QSpinBox *spinbox;
	int* variable;

	public slots:
    		void set_value (int val);
//		void set (const QString &);
	void set (int d);
};




class MyLabelf : public QWidget {
   Q_OBJECT

public:
    MyLabelf (QLayout *parent, const char *name, float& var);
//    ~MySlider ();
//private:
    QLabel *label;
    float* variable;
    void update ();
    void set_variable (float *f);
};


class MyLineEdit : public QWidget {
	Q_OBJECT

public:
	MyLineEdit (QLayout *parent, const char *name );
	QLineEdit *linedit;
	QLabel *label;
protected:
	const char *tag;
};



class MyLineFile : public QWidget {
	Q_OBJECT

public:
	MyLineFile (QLayout *parent, const char *name, int valid);
//	const char* get_value ();
//	QPushButton *adbutt; 
	QLineEdit *linedit;
	QLabel *label;
	QLabel *control_yes;
	QLabel *control_no;
//	void ins (QString s);
	string val ();

protected:
	const char *tag;
	QPushButton *fbutton; 
	char *filetype;
	int valid;

	protected slots:
//		QString ask_file ();
		void set_file_a (); 
		void set_file_b ();
		void set_file_c ();
		void set_file_d ();
};


class MyListButton : public QWidget {
	Q_OBJECT

public:
	MyListButton (QLayout *parent, const char *name);
	QPushButton *fbutton;
	QComboBox *_combo_box;

protected:


};

class MyPushButton : public QWidget {
	Q_OBJECT

public:
	MyPushButton (QLayout *parent, const char *name, int dim=0, int width=0);
	MyPushButton (const char *name, int dim=0, int width=0);
	QPushButton *fbutton;

protected:

};


class MyListView : public QWidget {
	Q_OBJECT

public:
	MyListView (QLayout *parent, int dim=0, int width=0, bool button=true);
	QListWidget *_lw;
	QPushButton *fbutton;
	QShortcut *shortcut1, *shortcut2;


	protected slots:
		void del_list_view_slot ();
	signals:
	void deleting (int r);

};


class MySlider : public QWidget {
   Q_OBJECT

public:
    MySlider (QLayout *parent, const char *name, int& var, int vmin, int vmax);
//    ~MySlider ();
//private:
    MyIntegerEditLine *pline;
    QSlider *slider;

    public slots:
    void setValue (int i);
};


class MyTwoColumn : public QWidget {
	Q_OBJECT

public:
	MyTwoColumn (QLayout *parent);
	QLayout *left_layout () {return _left ->layout ();};
	QLayout *right_layout () {return _right ->layout ();};
	ZNWidget *_left, *_right;

private:

protected:

};

class My3Column : public QWidget {
	Q_OBJECT

public:
	My3Column (QLayout *parent);
	QLayout *left_layout () {return _left ->layout ();};
	QLayout *center_layout () {return _center ->layout ();};
	QLayout *right_layout () {return _right ->layout ();};

private:

protected:
	ZNWidget *_left, *_center, *_right;
};

class MyGridColumn : public QWidget {
	Q_OBJECT

public:
	MyGridColumn (QLayout *parent, int row, int col);
	QGridLayout *gridlayout;

private:

protected:

};


class MyTableWidget : public QWidget {
	Q_OBJECT

public:
	MyTableWidget (QLayout *parent, int row, int height);
	QTableWidget *table;

private:

protected:

};

#endif
