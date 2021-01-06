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

#include "ZNdata.h"
#include <QMainWindow>
#include <QtOpenGL/qgl.h>
#include "arcball.h"
#include "molecule.h"
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
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/molchrg.h>
#include <openbabel/groupcontrib.h>

static bool g_stereoAvailable = FALSE;

//class vect;
class MyGl; class Data; class Builder; class BuilderMenu; class HapticMenu; class BrowserMenu; class Clicked_atomMenu; class ColorMenu;
class SurfaceMenu; class SphereMenu; class GraphicalObjectsMenu; class DDSettingsMenu; class GridMenu; class Command;
class MyFloatEditLine;
class Grid;


typedef struct {
 vect coordinates;
} myVertex;

class DDWin : public QMainWindow 
{
   Q_OBJECT


public:
    int mode;
    QComboBox *target;
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
    Molecule *target_molecule;

    vector <Molecule *> molecules;
    vector <GraphicalObject *> graphical_objects;

        
    bool show_labels;
    bool  backbone_shown, haptic_mode, minimizing, rec_movie, rec_pov_movie;
    Data* data;
    Builder* builder;

        void load_file (string mol_name);
        void load_multi_file (string mol_name);
        Sphere* new_sphere (string name);

        void new_molecule (string name);
        void add_molecule (Molecule *mol);


        void show_atom_properties (Atom *at);


        QComboBox *at_disp;
        QComboBox *bo_disp;
        QComboBox *ar_disp;

        int last_frame;
        void write_POV_source (string filename);
        void deselect ();

        MyGl *gl ;  
        void set_current_target (int index); 
        void set_current_target (Molecule * mol); 
        void set_lists (Molecule *mol);

        void select (Atom *at);

        void delete_current_molecule ();
        void delete_molecule (Molecule *mol);

        void remove_molecule (Molecule *mol);

        void delete_graphical_object (int i);

        string write_mol2_string (Molecule *mol);

        HapticMenu * haptic_menu;
        SurfaceMenu * surface_menu;
        GridMenu *grid_menu;
	SphereMenu *sphere_menu;
        GraphicalObjectsMenu *graphical_objects_menu;

        void execute (Command *comm);
		
		
	void lock_editing ();
	void unlock_editing ();
	

private:
	int edit_lock_counter;
    void closeEvent(QCloseEvent *event);

    QAction *openAct;
    QAction *saveasAct;
    QAction *screenshotAct;
    QAction *raytracedscreenshotAct;
    QAction *movieAct;
    QAction *raytracedmovieAct;
    QAction *quitAct;

    QAction *builderAct;
    QAction *historyAct;
    QAction *wiimoteAct;
    QAction *wiimote2Act;

    QAction *hideHAct;
    QAction *hidenpHAct;
    QAction *showallAct;
    QAction *hideallAct;

    QAction *displaysettingsAct;
    QAction *colorAct;
    QAction *DDsettingsAct;

    QAction *hapticAct;    
    QAction *dockingAct;
    QAction *energyAct;
    QAction *logPAct;
    QAction *minimiseAct;
    QAction *partialQAct;
    QAction *scoresCharge;

    QAction *surfaceAct;
    QAction *sphereAct;
    QAction *graphicalobjectsAct;

    QAction *addHsAct;

    QAction *aboutAct;

//toolbar

    QAction *undoAct;
    QAction *redoAct;        


        Q3Accel *accel;
        void setup_accel ();


        int style_str_to_i (string style);
        QWidget *dsetpopup;
      //  BrowserMenu * browser_menu;
        BuilderMenu * builder_menu;
        ColorMenu * color_menu;
        DDSettingsMenu *DDsettings_menu;
        Clicked_atomMenu * clicked_atom_menu;


        void draw_menu ();
        void create_menu_actions ();
        void set_popups ();
 

        vector <MyFloatEditLine*> at_opts, bo_opts;


        
        void dragEnterEvent( QDragEnterEvent * );
        void dragMoveEvent( QDragMoveEvent * );
        void dragLeaveEvent( QDragLeaveEvent * );
        void dropEvent( QDropEvent * );

        string POV_vector (vect v);
        string POV_color (color c);
        string POV_bond (Bond *bo);
        string POV_atom (Atom *at);
        string POV_ring (Ring *ring);
        string POV_surface (Surface *surf);
        string POV_sphere (Sphere *sph);
	
	string POV_stick (vect v1, vect v2, color c1, color c2, float rad, float a1size, float a2size);


		QString last_visited_dir;

private slots:


    void resizeEvent (QResizeEvent *);

    void open_file_slot ();
    void save_as_slot (); 
    void screenshot_slot ();
    void movie_slot ();
    void raytraced_movie_slot ();

    void builder_slot ();

    void history_slot ();
    void wiimote_slot ();
   void wiimote2_slot ();
   
    void hide_hydrogens_slot ();
    void hide_nonpolar_hydrogens_slot ();
    void show_all_atoms_slot ();
    void hide_all_atoms_slot ();

    void display_settings_slot ();
    void DD_settings_slot ();
    void color_slot ();

    void surface_slot ();
    void sphere_slot ();
    void graphical_objects_slot ();

    void add_Hs_slot ();


        void disp_ok ();

        void atdebug_slot ();    
        void about_slot ();
        void partial_charges_slot ();
        void scores_from_charges_slot ();
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

        void set_current_target_slot (int index);  
        
//key slots
        void del_pressed ();
        void esc_pressed ();
        void a_pressed ();
        void b_pressed ();
        void c_pressed ();
        void m_pressed ();
        void n_pressed ();
        void o_pressed ();
        void s_pressed ();


//toolbar slots



public slots:
        void raytraced_screenshot_slot ();
        void redraw (Molecule *mol);
        void recolor_by_score (Molecule *mol);
        void end_minimisation ();

		void emit_targets_updated ();

signals:
	void targets_updated ();
};





























class MyGl : public QGLWidget
{
   Q_OBJECT
    public:
        MyGl (DDWin *parent);
	bool needs_GL_update;
    void matrix_transformations (bool stereo = false, int frameBuffer = 0);

    void draw_molecule (Molecule* mol);
    void draw_bindingsite (float x,float  y,float  z, float r);
    void screenshot (QString filename);
    void set_center_of_rotation (vect v);
    void set_center_of_view (vect v);

    void refreshStencilBuffer();

    color select_color;
    color background_color;
    color bindingsite_color;
    color water_color;
    color current_color;
    int aromatic_display_style;

    float sphere_radius;
    float  stick_rad;
    float double_bond_inter_distance;
    float aromatic_bond_inter_distance;
    float vdw_scale;
    int vdw_precision, stick_precision, sphere_precision;
    float double_bond_stick_radius_scale;
    int fog_begin;
    float surface_resolution;
    float stereo_toe_in_angle;
    float stereo_inter_eye_semi_distance;


    float head_tracking_x_position, head_tracking_y_position;





    void draw_backbone (Molecule* mol);
    void draw_surface (Molecule* mol);
    void hide_hydrogens (Molecule* mol);
    void hide_nonpolar_hydrogens (Molecule* mol);
    void show_all_atoms (vector <Molecule*> mols);
    void hide_all_atoms (vector <Molecule*> mols);
    void hide_hydrogens (vector <Molecule*> mols);
    void hide_nonpolar_hydrogens (vector <Molecule*> mols);
    void show_all_atoms (Molecule* mol);
    void hide_all_atoms (Molecule* mol);

    void draw_list (Molecule* mol);
	void draw_backbone_list (Molecule *mol);
    void draw_list (Selection* sel);
    void draw_list (MarchingCubes * cube);
//    void draw_list (Grid& grid, int list);


        void haptic_to_world_coordinates (vect &haptic_p, vect &world_p, float minx, float maxx, float miny, float maxy, float minz, float maxz);
        void world_to_haptic_coordinates (vect &world_p, vect &haptic_p, float minx, float maxx, float miny, float maxy, float minz, float maxz);

    void move_molecule (Molecule *mol, float x, float y, float z, bool cut_bool = true);

        unsigned int next_list;
        vector <int> empty_lists;

    void free_list (int list);
    int new_list ();
    int lastx, lasty;


    void compute_double_bond_vertexes (Bond *bond, float out[4][3], float d=0);

        Matrix4fT Transform, Head_Tracking_Transf;

        void update_current_color ();
        QTime time, movie_time;
        bool  magic_pencil;
    void apply_color_masks (vector <color_mask> masks, Molecule *mol, bool undoable = TRUE);
	
	inline void get_viewport_points (vect &up, vect &down, vect &left, vect &right) {up = up_point; down = down_point; left = left_point; right = right_point;}
public slots:


    private:

	vect up_point, down_point, left_point, right_point; // to define the current view. updated by every paintGL call

        vector <float> vw_haptic_force;
        vector <float> el_haptic_force;
        vector <float> vw_haptic_torque;
        vector <float> el_haptic_torque;

        float el_haptic_e, vw_haptic_e;

        Atom *clicked_atom;
 
        float lasta2, lasta2increment;
        bool    isClicked, isDragging;																		
        Matrix3fT LastRot, ThisRot; //Last_Head_Tracking_Rot;

        Point2fT    MousePt;
        bool zoom, translate, rotate, select, move, spin, selection_square;
        int selection_square_x1, selection_square_y1, selection_square_x2, selection_square_y2;
        vect center_of_rotation;
        vect center_of_view;
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




        void draw_atoms_for_selection (Molecule* mol);
        void draw_bonds_for_selection (Molecule* mol);       
        void draw_bond_stick(Bond* bond);
        void draw_double_bond_stick(Bond* bond);
        void draw_bond_line(Bond* bond);
        void draw_double_bond_line(Bond* bond);
        void draw_triple_bond_line (Bond* bond);
        void draw_aromatic_bond_line(Bond* bond);
        void draw_aromatic_bond_stick(Bond* bond);
        void draw_triple_bond_stick (Bond* bond); 
        void draw_ring_line (Ring* ring);
        void draw_ring_stick (Ring* ring);

        void draw_atom_sphere(Atom* atom);
        void draw_atom_sel_sphere(Atom* atom);
        void draw_atom_vdw_sphere(Atom* atom);
        void draw_atom_scaled_vdw_sphere(Atom* atom);


        void draw_backbone_line (Resid *res);
        void draw_backbone_stick (Resid *res);
	vector <vect> get_backbone_points (Resid *res) ;
	vector <vect> smooth_list (vector <vect> lis);
        void setAtomColor(Atom* atom);
		void openGLSetColor (color col);


        void set_conf (int conf);



		void my_line (vect v1, vect v2, color c1, color c2);
		void my_cylinder (float rada, float radb, float lenght, unsigned int slices, Atom* at1, Atom *at2);
		void my_cylinder (vect v1, vect v2, float radone, float radtwo, color c1, color c2, unsigned int slices) ;
        void my_sphere (float rad, unsigned int slices, unsigned int stacks, Atom* at);



  //      void dragEnterEvent( QDragEnterEvent * );
   //     void dragMoveEvent( QDragMoveEvent * );
   //     void dragLeaveEvent( QDragLeaveEvent * );
    //    void dropEvent( QDropEvent * );


    private slots:
        void set_wireframe ();
        void set_stick ();
        void set_cpk ();
        void set_ballandline ();
        void set_ballandstick ();

	void head_tracking_update (int x, int y);
	void move_camera (float x, float y, float z);

     //   void hide_hydrogens ();
     //   void hide_nonpolar_hydrogens ();
     //   void show_all ();





};

        void set_color (color c);


#endif
