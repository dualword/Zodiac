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
#include <qmenubar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <QtOpenGL/qgl.h>
#include <qlabel.h>
#include <Qt3Support/q3grid.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qevent.h>
#include <Qt3Support/q3dragobject.h>
#include <qdir.h>
#include <iostream>
#include <qcolordialog.h>
#include <qtabwidget.h>
#include <qslider.h>
#include "builder.h"
#include "ddwin.h"
#include "molecule.h"
#include "qstackedwidget.h"
#include <qlistwidget.h>
#include "command.h"
#include <QRadioButton>
#include "thread.h"
#include <QTableWidget>
#include "database.h"



class MySlider; class MyFloatEditLine; class Builder; class Minimize; class DDWin; class MyLabelf; class MyCompleteColorSettings; class MyColorButton;
class DatabaseField;





class SurfaceMenu : public QWidget 
{
   Q_OBJECT

public:

    SurfaceMenu ( QWidget *parent, Data* dat);
    Data *data;
    Surface *surface;

private slots:
    void add_surface ();

private:
//    Molecule *molecule;
    QComboBox *stype;
    QComboBox *gtype;
    MyFloatEditLine *resolution;
    MyFloatEditLine *alpha_p;
    MySlider *alpha_s;
    float res;
    float color [4];
    float alpha;
    int type;
    bool mesh;


private slots:
    void draw_surface ();

};




class SphereMenu : public QWidget 
{
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

    vect center;
    float alpha;
    float radius;



private slots:
    void draw_sphere ();

};






class BuilderMenu : public QWidget
{
   Q_OBJECT

public:

    BuilderMenu ( QWidget *parent, Builder* builder);
    Builder *builder;
    QLineEdit *smiles;

private slots:
    void add_smiles ();
    void add_C ();
    void add_N ();
    void add_O ();

    void single_bond ();
    void double_bond ();
    void triple_bond ();
    void no_bond ();




};



class HapticMenu : public QTabWidget 
{
   Q_OBJECT

public:
    QComboBox *interff;
    QComboBox *dofmode;
    HapticMenu ( QWidget *parent, Minimize* min);
    Minimize *minimize;

    MyLabelf *total_E, *internal_E, *interaction_E;
    void update ();

private slots:
    void Ok ();
    void end ();




};



class BrowserMenu : public QWidget
{
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

class Clicked_atomMenu : public QTabWidget 
{
   Q_OBJECT

public:
    Clicked_atomMenu ( QWidget *parent, DDWin* ddwin);


    DDWin *ddwin;
        Q3Frame *atomselpopup;
        QLabel *aplid, *aplat, *aplq;
        QLineEdit *aplfc, *aplx, *aply, *aplz;
        QLabel *resna, *resnu;

        Atom *clicked_atom;

        void set (Atom *at);
        void update ();

private slots:
    void set_value (QLabel *lab, float val) ;
    void set_value (QLabel *lab, string val) ;
    void set_value (QLineEdit *lab, float val) ;


    void add_Hs ();
    void set_clicked_atom_as_center_of_view ();
    void set_clicked_atom_as_center_of_rotation ();


};


class ColorMenu : public QWidget 
{
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




class DDSettingsMenu : public QWidget 
{
   Q_OBJECT

public:
    DDWin *ddwin;
    DDSettingsMenu ( QWidget *parent, DDWin* ddwin);
    float focal_d;  

private:
    MyFloatEditLine *inter_eye_distance, *focal_point_distance;

private slots:
    void ok_slot ();

};









class GraphicalObjectsMenu : public QWidget 
{
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

};




/*

class GridMenu : public QWidget 
{
   Q_OBJECT

public:

    GridMenu ( QWidget *parent, Data* dat);
    Data *data;


private:
    MyFloatEditLine *res_le;

    vect center;
    float resolution;




private slots:


};

*/
class DatabaseGrid : public QWidget
{
   Q_OBJECT

public:

    DatabaseGrid ( QWidget *parent, Database* db, Data *dat);
    Database *db;
	Data *data;
	QTableWidget *tab;
	BrowserMenu *browser;
	inline bool has_extend_enabled () {return ext_bool;};
	void add_field (DatabaseField *);
	int real_index_of_line (int);

private:
	bool ext_bool;
	QLineEdit *le;
	int current_number;
	void add (Database *db);
	void set_mol ();

private slots:
	void select_row (int r);
	void deselect_row (int r);
	void set_row_color (int r, color c);
	void set_current_molecule (int i);
	void manage_number_changed (const QString str);
	void manage_double_click (int r, int c); 
    void first_slot ();
    void prev_slot ();
    void next_slot ();
    void last_slot ();



private slots:


};















class MySlider : QWidget {
   Q_OBJECT
public:
    MySlider (QLayout *parent, const char *name, float& var, int vmin, int vmax);
//    ~MySlider ();
//private:
    MyFloatEditLine *pline;
    QSlider *slider;
    public slots:
    void setValue (const QString& s);
};


class MyFloatEditLine : public QWidget
{
    Q_OBJECT

public:
    MyFloatEditLine (QLayout *parent, const char *name, float& var);
    MyFloatEditLine (QLayout *parent, const char *name, double& var);

 //   const char* get_value ();
 //   void ins (QString s);
    void set ();
    QLineEdit *linedit;
    float* variable;
public slots:
    void set_value (int val);
    void set (const QString &);

};


class MyCompleteColorSettings : public QWidget {
    Q_OBJECT
public:
    MyCompleteColorSettings (QLayout *parent, color &col);
    MyColorButton *button;
    MySlider *slider;
    float alpha;
private slots:
    void setAlpha (const QString &);
};


class MyLabelf : Q3HBox {
   Q_OBJECT
public:
    MyLabelf (QWidget *parent, const char *name, float* var=NULL);
//    ~MySlider ();
//private:
    QLabel *label;
    float* variable;
    void update ();
    void set_variable (float *f);
};



class MyColorButton : QPushButton {
    Q_OBJECT
public:
    MyColorButton (QLayout *parent, QColor &color);

  void paintEvent ( QPaintEvent * ) ;
    QColor *color;

private slots:
    void my_clicked ();
};


class MyCheckBox : QCheckBox {
    Q_OBJECT
public:
     MyCheckBox (QLayout *parent, bool &var, string string);

private:
    bool *var;
private slots:
    void get ();
    void set ();
};

#endif
