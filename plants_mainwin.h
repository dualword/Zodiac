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


#ifndef PLANTS_MAINWIN_H
#define PLANTS_MAINWIN_H

#include "constants.h"
#include <qtabwidget.h>
#include <qstring.h>
#include <qfileinfo.h>
#include <Qt3Support/q3frame.h>
#include <qpushbutton.h>
#include <Qt3Support/q3listbox.h>
#include <qcombobox.h>
#include "myline.h"
#include "plants.h"
#include "constants.h"


using namespace std;

class MyListBoxItem : public Q3ListBoxText
{

public:
    bool valid;
    MyListBoxItem(string s, Q3ListBox *par);


protected:
bool check_mol2 (string line);
    Q3ListBox *parent;
    string tex;
    virtual void paint( QPainter * );
    virtual int width( const Q3ListBox* ) const { return (parent->width ()-5); }
    virtual int height( const Q3ListBox* ) const { return 20; }
};



class Plants;
class PlantsMainWin : public QTabWidget 
{
    Q_OBJECT

public:
    PlantsMainWin( Plants *parent);
    Plants *window;
    vector <QComboBox*> lat, pat, pres; 
    

protected:
    QString filename;
    QFileInfo fileinfo;


    void setupwin();

 private slots:
    void add_water_slot();
    void addphbc();
    void addshc();
    void addsdc();
    void addidc();
    void addldc();
    void add_flsc_slot();
 //   void askfil(const char *dir, const char *type, QWidget *parent, const char *title, const char *text);
    void add_phbc_slot ();
    void add_shc_slot ();
    void add_sdc_slot ();
    void add_idc_slot ();
    void add_ldc_slot ();
    void add_lf_slot ();
    void add_ll_slot ();
    void set_sf_slot ();
    void add_fixp_slot ();
    void load_bs_from_ligand_slot ();
    void ligand_selected (int index);

    void load_protein ();
    void load_ligand ();
    void check_outdir ();
    void check_bindingsite ();


    void del_ifile ();
    void del_constr ();
    void del_flex ();
    void del_water ();

 private:

    void update_boxes (vector<QComboBox*> &vec, vector<string>& values);


    QWidget* phbcpopup; 
    QWidget* shcpopup;  
    QWidget* sdcpopup; 
    QWidget* idcpopup;
    QWidget* ldcpopup;
    QWidget* flscpopup;
    QPushButton *addwaterb;
    QPushButton *addphbcb;
    QPushButton *addshcb;
    QPushButton *addsdcb;
    QPushButton *addidcb;
    QPushButton *addldcb;
    QPushButton *addflscb;
    QPushButton *addfixpb;

    QComboBox *anlined;
    QLineEdit *phbweil;
    QLineEdit *shfl;
    QLineEdit *shweil;


 
    QComboBox *sdanl;
    QLineEdit *sdbbfl;
    QLineEdit *sdbbtl;
    QLineEdit *sdweil;
    QComboBox *idanl;
    QComboBox *idanl2;
    QLineEdit *idbbfl;
    QLineEdit *idbbtl;
    QLineEdit *idweil;
    QComboBox *ldanl;
    QComboBox *ldanl2;
    QLineEdit *ldbbfl;
    QLineEdit *ldbbtl;
    QLineEdit *ldweil;

    QComboBox *flscl, *fixpl;

    QPushButton *bsfl;


    QLineEdit *xwc, *ywc, *zwc, *wr;


    MyLineF *lfile, *llist;


public:

    MyLineF *pfile, *wref;
    MyLine *aants, *aevap, *asigma, *chbw, *cmw, *hbw, *hbchow, *mw, *plpw, *iw, *crmsd, *cs, *odir, *ipsw, *wphb, *wlhb, *wwhb, *nwhbp, *bsx, *bsy, *bsz, *bsr;
    QPushButton *flipab, *flipn, *ffbp;
    QPushButton *wpc, *wrs, *wmm2, *wrl, *wrmm2, *cho, *rigidl, *rigida, *wpb, *wps, *wpas;
    QComboBox *sfchoice, *lintra;
    Q3ListBox *iflb, *constrlb, *flexlb, *waterlb;
    void check_ligands ();
    void set_bindingsite (float x, float y, float z);
    void set_bindingsite (float x, float y, float z, float r);
    float getfloat (string s);
};

#endif
