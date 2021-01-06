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

#ifndef DDWIN_H
#define DDWIN_H


#include "constants.h"
#include "ZNdata.h"
#include <QMainWindow>
#include "arcball.h"
#include "ZNmolecule.h"
#include <qlabel.h>
#include <qcombobox.h>
#include "plants.h"
#include "MarchingCubes.h"
#include "graphical_object.h"
#include "minimize.h"
#include "builder.h"
#include <Qt3Support/q3accel.h>
#include <QtGui>
#include <qdatetime.h>
#include "menu.h"
#include <QUndoView>
#include "command.h"
#include "obabel_includes.h"
#include "thread.h"




#include <qmenubar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>

#include <qlabel.h>
#include <Qt3Support/q3grid.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qevent.h>
#include <qdrag.h>
#include <qdir.h>
#include <qcolordialog.h>
#include "database.h"



#include "function.h"
#include "ils.h"

#ifdef WIN32
#include <float.h>
#define isnan _isnan
#endif // WIN32






static bool g_stereoAvailable = FALSE;

//class vect;
class MyGl; class Data; class Builder; class BuilderMenu; class HapticMenu;  class Clicked_atomMenu; class ColorMenu; class ColorSettingsMenu; class GamessMenu; class PLANTSMenu; class PrefMenu;
class BackboneColorMenu; class SurfaceMenu; class SphereMenu; class GraphicalObjectsMenu; class GraphicalObjectsColorMenu; class DDSettingsMenu; class GridMenu; class ConformersMenu; class Command; class DatabaseGrid; class IODeviceMenu; class DockingMenu; class DisplayMenu; class SequenceMenu;
class ThreadMenu; class MapMenu;
class MyFloatEditLine;
class Grid;
class Thread;

typedef struct {
	vect coordinates;
} myVertex;



class DDWin : public QMainWindow {
	Q_OBJECT

public:


		vector <HBond> Hbonds;
		vector <vect> vertex_list;
		int mode;
		QComboBox *target;
		QLabel *current_mode;
		QUndoView *undo_view;
		DDWin (QWidget *parent, Data *dat);
		
		string mode_string;
		
		enum StereoMode { StereoMode_None, StereoMode_VerticalInterlace, StereoMode_HorizontalInterlace, StereoMode_ViaOpenGL};
		
		StereoMode g_stereoMode;
		unsigned char *g_stencil_mask;
		bool g_stencil_mask_needs_redraw;
		unsigned int g_stencil_mask_frame_counter;
		unsigned int g_stencil_mask_height;
		unsigned int g_stencil_mask_height_tile;
		unsigned int g_stencil_mask_width;
		unsigned int g_stencil_mask_width_tile;
		
		vector <float > magic_pencil_trail_x;
		vector <float > magic_pencil_trail_y;
		vector <float > magic_pencil_trail_z;
		vector <float > magic_pencil_trail_wx;
		vector <float > magic_pencil_trail_wy;
		
		int current_target;
		ZNMolecule *target_molecule;
		
		vector <ZNMolecule *> molecules;
		vector <DatabaseGrid *> database_grids;
		vector <GraphicalObject *> graphical_objects;
		vector <Thread *> running_threads;
		
        void run_thread (Thread *);
		bool show_labels;
		bool  backbone_shown, haptic_mode, minimizing, rec_movie, rec_pov_movie, adding_restrains;
		Data* data;
		Builder* builder;
		
        void load_file (string mol_name);
        Database *load_multi_file (string mol_name);
        Sphere* new_sphere (string name);
		
        void new_molecule (string name);
        void add_molecule (ZNMolecule *mol);
		void add_database (Database *dat);
		
        void show_atom_properties (Atom *at);
		
/*
        QComboBox *at_disp;
        QComboBox *bo_disp;
        QComboBox *ar_disp;
	   QComboBox *backbone_disp;
*/		
        int last_frame;
        void write_POV_source (string filename);
        void deselect ();
		
        MyGl *gl ;  
        void set_current_target (int index); 
        void set_current_target (ZNMolecule * mol); 
        void set_lists (ZNMolecule *mol);
		void show_grid (Database *dat); 
		
        void select (Atom *at);
		void select (vector <Atom *> atoms);
		
        void delete_current_molecule ();
        void delete_molecule (ZNMolecule *mol);
		
        void remove_molecule (ZNMolecule *mol);
		
        void delete_graphical_object (int i);
		
        string write_mol2_string (ZNMolecule *mol);
	BackboneColorMenu *b_color_menu;
		ConformersMenu *conformers_menu;
        HapticMenu *haptic_menu;
		GamessMenu *gamess_menu;
		PLANTSMenu *plants_menu;
		DisplayMenu *display_menu;
		SequenceMenu *sequence_menu;
		DockingMenu *docking_menu;
		ThreadMenu *thread_menu;
		IODeviceMenu *iodevice_menu;
        SurfaceMenu *surface_menu;
		GraphicalObjectsColorMenu *go_color_menu;
        GridMenu *grid_menu;
		SphereMenu *sphere_menu;
		MapMenu *map_menu;
        GraphicalObjectsMenu *graphical_objects_menu;
		
        void execute (Command *comm);
		
		
		void lock_editing ();
		void unlock_editing ();
		
		
	private:
		int edit_lock_counter;
		void closeEvent(QCloseEvent *event);
		color molecule_color;
		
		QAction *openAct;
		QAction *openSessionAct;
		QAction *saveasAct;
		QAction *saveSessionAsAct;
		QAction *screenshotAct;
		QAction *raytracedscreenshotAct;
		QAction *movieAct;
		QAction *raytracedmovieAct;
		QAction *quitAct;
		QAction *newdatabaseAct;
		
		QAction *builderAct;
		QAction *twodwinAct;
		QAction *historyAct;
		QAction *wiimoteAct;
		QAction *wiimote2Act;
		QAction *settingsAct;
		
		QAction *hideHAct;
		QAction *hidenpHAct;
		QAction *showallAct;
		QAction *hideallAct;
		
		QAction *displaysettingsAct;
		QAction *sequenceAct;
		QAction *colorAct;
		QAction *backboneColorAct;
		QAction *backgroundcolorAct;
		QAction *DDsettingsAct;
		QAction *enableclippingAct;
		QAction *disableclippingAct;
		
		QAction *hapticAct;    
		QAction *dockingAct;
		QAction *plantsAct;
		QAction *gamessAct;
		QAction *energyAct;
		QAction *logPAct;
		QAction *minimiseAct;
		QAction *partialQAct;
		QAction *scoresCharge;
		QAction *conformationalSearchAct;
		
		QAction *surfaceAct;
		QAction *mapAct;
		QAction *sphereAct;
		QAction *backbonetosurfAct;
		QAction *graphicalobjectsAct;
		
		QAction *addHsAct;
		QAction *duplicateAct;
		
		QAction *iodeviceAct;
		
		QAction *aboutAct;
		
		
		QAction *undoAct;
		QAction *redoAct;        
		
		QAction *colorsAct, *colors_chooserAct;
		QAction *centerAct;
		
		QAction *selectAllAct;
		QAction *deselectAct;
		QAction *invertSelectAct;
		QAction *selectCAct;
		QAction *selectHAct;
		QAction *selectsquareAct;
		
		
		
		QAction *buildCAct;
		QAction *buildSAct;
		QAction *buildNAct;
		QAction *buildOAct;
		QAction *buildsinglebondAct;
		QAction *builddoublebondAct;
		QAction *buildtriplebondAct;
		QAction *deletebondAct;

		QAction *deleteatomAct;
	//	QAction *magicpencilAct;
	
		QPushButton *magicpencilButt;
		QAction *buildsmileAct;
		QAction *periodictableAct;

		QToolBar *toolbar, *target_toolbar, *select_tb, *builder_tb;

		
        Q3Accel *accel;
        void setup_accel ();
		
		
        int style_str_to_i (string style);
        QWidget *dsetpopup;
        BuilderMenu *builder_menu;
        PrefMenu *pref_menu;
        ColorMenu *color_menu;
		ColorSettingsMenu *color_settings_menu;
        DDSettingsMenu *DDsettings_menu;
        Clicked_atomMenu * clicked_atom_menu;
		
		
        void draw_menu ();
		void hide_toolbars ();
        void create_menu_actions ();
        void set_popups ();
		
        
        void dragEnterEvent( QDragEnterEvent * );
        void dragMoveEvent( QDragMoveEvent * );
        void dragLeaveEvent( QDragLeaveEvent * );
        void dropEvent( QDropEvent * );
		
        string POV_vector (vect v);
        string POV_color (color c);
        string POV_bond (ZNBond *bo, bool clippable = false);
        string POV_atom (Atom *at,  bool clippable = false);
        string POV_ring (Ring *ring);
        string POV_surface (Surface *surf, bool clippable = false);
		string POV_map (Map *map, bool clippable = false);
        string POV_sphere (Sphere *sph, bool clippable = false);
		string POV_backbone (ZNMolecule *mol, bool clippable = false);
		
		string POV_stick (vect v1, vect v2, color c1, color c2, float rad, float a1size, float a2size, bool clippable = false);
		
		public:
		QString last_visited_dir;
		
		private slots:
		
		void resize_target_toolbar (Qt::Orientation);
		void resizeEvent (QResizeEvent *);
		
		void open_file_slot ();
		void open_session_file_slot ();
		void save_as_slot (); 
		void save_session_as_slot (); 
		void screenshot_slot ();
		void movie_slot ();
		void raytraced_movie_slot ();
		
		void new_database_slot ();
		void builder_slot ();
		void twodwin_slot ();
		
		void history_slot ();
		void wiimote_slot ();
		void wiimote2_slot ();
		void pref_slot ();
		
		void hide_hydrogens_slot ();
		void hide_nonpolar_hydrogens_slot ();
		void show_all_atoms_slot ();
		void hide_all_atoms_slot ();
		
		void display_settings_slot ();
		void sequence_slot ();
		void DD_settings_slot ();
		void disable_clipping_slot ();
		void enable_clipping_slot ();
		void color_slot ();
	void backbone_color_slot ();
		void background_color_slot ();
		
		void surface_slot ();
		void map_slot ();
		void sphere_slot ();
		void backbone_to_surface_slot ();
		void graphical_objects_slot ();
		
		void add_Hs_slot ();
		void duplicate_slot ();
		

        void atdebug_slot ();    
        void about_slot ();
        void partial_charges_slot ();
        void scores_from_charges_slot ();
		void conformers_slot ();
        void compute_energy_slot ();
        void logP_slot ();
        void minimise_energy_slot ();
        void not_impl_slot ();
		//     void color_by_atom_slot ();
		//     void color_by_charge_slot ();
		//     void color_by_score_slot ();
		//     void color_carbons_slot ();
		
        void haptic_slot ();
        void plants_slot ();
		void docking_slot ();
		void gamess_slot ();
		void iodevice_slot ();
		
        void set_current_target_slot (int index);  
		void invert_selection_slot ();
		void deselect_slot ();
		void select_all_slot ();
        void select_C_slot ();
		void select_H_slot ();
		
		
		
		void select_square_slot ();
		
		void build_C_slot ();
		void build_S_slot ();
		void build_O_slot ();
		void build_N_slot ();
		
	
		void build_single_bond_slot ();
		void build_double_bond_slot ();
		void build_triple_bond_slot ();
		void delete_bond_slot ();
		
		void magic_pencil_toggle_slot (bool checked);
		void magic_pencil_start ();
		void magic_pencil_end ();

		void periodic_table_slot ();
		
	void switch_mode (int m);
		//key slots
        void del_pressed ();
        void esc_pressed ();
        void a_pressed ();
        void b_pressed ();
        void c_pressed ();
        void d_pressed ();		
        void e_pressed ();
        void f_pressed ();
        void h_pressed ();		
		void i_pressed ();
        void m_pressed ();
        void n_pressed ();
        void o_pressed ();
		void r_pressed ();
        void s_pressed ();
		
		
		//toolbar slots
		
		
		
		public slots:
        void raytraced_screenshot_slot ();
		//  void redraw (ZNMolecule *mol);
        void recolor_by_score (ZNMolecule *mol);
        void end_minimisation ();
		void change_color_slot ();
		void quick_color_slot ();
		void center_slot ();
	void emit_go_updated () {go_updated();};
		void emit_targets_updated ();
		
	signals:
		void non_selection_molecules_updated ();
		void targets_updated ();
		void go_updated ();
	};





























class MyGl : public QGLWidget
	{
		Q_OBJECT
    public:
        MyGl (DDWin *parent);
		void timerEvent ( QTimerEvent * event );
		void set_clipping_plane ();
		bool needs_GL_update;
		void matrix_transformations (bool stereo = false, int frameBuffer = 0);
		
		void draw_molecule (ZNMolecule* mol);
		void GL_update_molecule (ZNMolecule *mol);
		void GL_update_backbone (ZNMolecule *mol);
		void draw_bindingsite (float x,float  y,float  z, float r);
		void screenshot (QString filename);
		void set_center_of_rotation (vect v);
		void set_center_of_view (vect v);
		
		
		void refreshStencilBuffer();
		
		color select_color;

		void translate_view (vect v);
		
		int redraw_counter;
		int aromatic_display_style;

		double double_bond_inter_distance;
		double aromatic_bond_inter_distance;


		double double_bond_stick_radius_scale;
		int fog_begin;
		double surface_resolution;
		double stereo_toe_in_angle;
		double stereo_inter_eye_semi_distance;
		bool select_square_mode;
		
		float head_tracking_x_position, head_tracking_y_position;
		
		vect apply_world_rotation (vect v);
		vect deapply_world_rotation (vect v);
				
		void hide_hydrogens (ZNMolecule* mol);
		void hide_nonpolar_hydrogens (ZNMolecule* mol);
		void show_all_atoms (vector <ZNMolecule*> mols);
		void hide_all_atoms (vector <ZNMolecule*> mols);
		void hide_hydrogens (vector <ZNMolecule*> mols);
		void hide_nonpolar_hydrogens (vector <ZNMolecule*> mols);
		void show_all_atoms (ZNMolecule* mol);
		void hide_all_atoms (ZNMolecule* mol);
		
		void draw_list (ZNMolecule* mol);
		void draw_backbone_list (ZNMolecule *mol);
		void draw_list (Selection* sel);
		void draw_list (MarchingCubes * cube);
		//    void draw_list (Grid& grid, int list);
		
		
        void haptic_to_world_coordinates (vect &haptic_p, vect &world_p, float minx, float maxx, float miny, float maxy, float minz, float maxz);
        void world_to_haptic_coordinates (vect &world_p, vect &haptic_p, float minx, float maxx, float miny, float maxy, float minz, float maxz);
		vect rotate_vector (vect v);
		vect unrotate_vector (vect v);
		
		
		void move_molecule (ZNMolecule *mol, float x, float y, float z, bool cut_bool = true);
		
		
		void translate_molecule (ZNMolecule *mol, vect v);

		
        unsigned int next_list;
        vector <int> empty_lists;
		
		void free_list (int list);
		int new_list ();
		int lastx, lasty;
		
		
		void compute_double_bond_vertexes (ZNBond *bond, float out[4][3], float d=0);
		
        Matrix4fT Transform, Head_Tracking_Transf;
		
        void update_current_color ();
        QTime time, movie_time;
        bool  magic_pencil;
		void apply_color_masks (vector <color_mask> masks, ZNMolecule *mol, bool undoable = TRUE);
		
		inline void get_viewport_points (vect &up, vect &down, vect &left, vect &right) {up = up_point; down = down_point; left = left_point; right = right_point;}
		public:
				void select_MW (unsigned int an = 6);
				
						
        vect center_of_rotation;
        vect center_of_view, view_translations;
		        Matrix3fT LastRot, ThisRot; //Last_Head_Tracking_Rot;
				void backbone_to_surface (ZNMolecule *mol, Surface *surf);
    private:
		int select_pulse;
		vect up_point, down_point, left_point, right_point; // to define the current view. updated by every paintGL call
		
        vector <float> vw_haptic_force;
        vector <float> el_haptic_force;
        vector <float> vw_haptic_torque;
        vector <float> el_haptic_torque;
		
        float el_haptic_e, vw_haptic_e;
		
        Atom *clicked_atom;
		
        float lasta2, lasta2increment;
        bool    isClicked, isDragging;																		

		
        Point2fT    MousePt;
        bool zoom, translate, rotate, select, move, spin, selection_square;
        int selection_square_x1, selection_square_y1, selection_square_x2, selection_square_y2;

        float zbeginy, tbeginx, tbeginy, mbeginx, mbeginy, sbeginx;
		
        DDWin *ddwin ;
        void initializeGL (); 
        void init_vars ();
        void paintGL();
        void resizeGL( int width, int height );
        void mousePressEvent ( QMouseEvent * e );
        void mouseReleaseEvent ( QMouseEvent * e );
        void mouseMoveEvent ( QMouseEvent * e );
        void begin_zoom (float y);
        void continue_zoom (float y);
        void begin_magic_pencil (float x, float y);
        void continue_magic_pencil (float x, float y);
        void end_magic_pencil (float x, float y);
        void begin_translate (float x, float y);
        void begin_move (float x, float y);
        void continue_translate (float x, float y);
        void continue_move (float x, float y);
        void begin_spin (float x, float y);
        void continue_spin (float x, float y);
        void begin_rotate (float x, float y);
        void continue_rotate (float x, float y);
        void selectf (float x, float y);
        Atom *select_atom (float x, float y, int dx, int dy);
        void select_square ();

		
        void begin_selection_square (float x, float y);
        void continue_selection_square (float x, float y);
        void end_selection_square ();
		
		
		
		
        void draw_atoms_for_selection (ZNMolecule* mol);
        void draw_bonds_for_selection (ZNMolecule* mol);
        void draw_bond_stick(ZNBond* bond);
        void draw_double_bond_stick(ZNBond* bond);
        void draw_bond_line(ZNBond* bond);
        void draw_double_bond_line(ZNBond* bond);
        void draw_triple_bond_line (ZNBond* bond);
        void draw_aromatic_bond_line(ZNBond* bond);
        void draw_aromatic_bond_stick(ZNBond* bond);
        void draw_triple_bond_stick (ZNBond* bond); 
        void draw_ring_line (Ring* ring);
        void draw_ring_stick (Ring* ring);
		
        void draw_atom_sphere(Atom* atom);
        void draw_atom_sel_sphere(Atom* atom);
        void draw_atom_vdw_sphere(Atom* atom);
        void draw_atom_scaled_vdw_sphere(Atom* atom);
		
		
        void draw_backbone_line (Resid *res);
        void draw_backbone_stick (Resid *res, Surface *surf = NULL);


		
		
        void setAtomColor(Atom* atom);
		void openGLSetColor (color col);
		
	
		
		void my_line (vect v1, vect v2, color c1, color c2);
		void my_cylinder (float rada, float radb, float lenght, unsigned int slices, Atom* at1, Atom *at2);
		void my_backbone_ribbon (vect v1, vect v2, vect v3, vect v4, vect dir, vect dir2, color c1, color c2, vector <vect> shape1, vector <vect> shape2, Surface *surf) ;
		void my_cylinder (vect v1, vect v2, float radone, float radtwo, color c1, color c2, unsigned int slices) ;
		void my_cylinder (vect v1, vect v2, vect v3, vect v4, float radone, float radtwo, color c1, color c2, unsigned int slices) ;
        void my_sphere (float rad, unsigned int slices, unsigned int stacks, Atom* at);
		
		
		
		//      void dragEnterEvent( QDragEnterEvent * );
		//     void dragMoveEvent( QDragMoveEvent * );
		//     void dragLeaveEvent( QDragLeaveEvent * );
		//    void dropEvent( QDropEvent * );
		
		
		private slots:
	
		void head_tracking_update (int x, int y);
		void move_camera (float x, float y, float z);
		void move_target (float x, float y, float z);
		void map_vector_on_vector_world (double x1, double y1, double z1, double x2, double y2, double z2);
		void map_vector_on_vector_target (double x1, double y1, double z1, double x2, double y2, double z2);
		
		//   void hide_hydrogens ();
		//   void hide_nonpolar_hydrogens ();
		//   void show_all ();
		
		
		
		
		
	};

void set_color (color c);


#endif
