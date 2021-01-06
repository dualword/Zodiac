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
#include <qlineedit.h>
#include <qvalidator.h>
#include <qlabel.h>
#include <qfiledialog.h>
#include <iostream>


using namespace std;

MyLine::MyLine( QWidget *parent, const char *name, int valid)
    : Q3HBox( parent, name )
{


    QLabel *label = new QLabel( name, this ); 
    label->setMaximumWidth( 250 );
    label->setMinimumWidth( 250 );
    linedit = new QLineEdit( this);
    if (valid == 1)
        linedit->setValidator( new QIntValidator( linedit ) );
    if (valid == 2)
        linedit->setValidator( new QDoubleValidator(  -999.0, 999.0, 2,linedit ) );

}

 void MyLine::ins (QString s)
{   
    linedit-> clear ();
    linedit->insert (s);
 
}

 string MyLine::val ()
{
   return linedit->text ().toStdString ();
   
}

////////////////////////////////////////////////////////////////////////////////////////////////////


MyLineF::MyLineF( QWidget *parent, QWidget *master, const char *name, int valid, int butt)
    : Q3HBox( parent, name )
{
static const char * p1_xpm[] = {
"16 16 3 1",
" 	c None",
".	c #000000000000",
"X	c #FFFFFFFF0000",
"                ",
"                ",
"         ....   ",
"        .XXXX.  ",
" .............. ",
" .XXXXXXXXXXXX. ",
" .XXXXXXXXXXXX. ",
" .XXXXXXXXXXXX. ",
" .XXXXXXXXXXXX. ",
" .XXXXXXXXXXXX. ",
" .XXXXXXXXXXXX. ",
" .XXXXXXXXXXXX. ",
" .XXXXXXXXXXXX. ",
" .XXXXXXXXXXXX. ",
" .............. ",
"                "};
    QPixmap p1( p1_xpm );

    tag = name;
    QLabel *label = new QLabel( name, this ); 
    label->setMaximumWidth( 150 );
    label->setMinimumWidth( 150 );
    linedit = new QLineEdit( this);
    if (valid == 1) {filetype = "Tripos Mol2 File (*.mol2)";};
    if (valid == 2) {filetype = "PLANTS Config File (*.pcfg);;All files (*)";};
    if (valid == 3) {filetype = "List File (*)";};
    QPushButton *fbutt = new QPushButton(  "" , this);
 //   fbutt->setPixmap( p1);
    connect (fbutt, SIGNAL( clicked() ), SLOT( set_file() ) );

    if (butt){
            adbutt = new QPushButton( "", this );

    adbutt->setMaximumWidth( 20 );
        }
    
}

const char* MyLineF::get_value()
{

 //   QString a = linedit->displayText();
   // return a.latin1 ();
    return "a";

}

QString MyLineF::ask_file()
{

    QString mol_name = QFileDialog::getOpenFileName(this, 
                    tr ("Open file"), "",tr("Tripos Mol2 File (*.mol2)"));

	return mol_name;
 }
 void MyLineF::set_file ()
{
    QString s = ask_file ();
    linedit->clear ();
    linedit->insert (s);
 //   set_line ();
}

 void MyLineF::set_line ()
{   QString st = tag; 
    st.append ("   | ");
    st.append (linedit->displayText ());
 //   master->iflb->insertItem(st);

}

 void MyLineF::ins (QString s)
{   
    linedit-> clear ();
    linedit->insert (s);
 
}

 string MyLineF::val ()
{
   return linedit->text ().toStdString ();
   
}
