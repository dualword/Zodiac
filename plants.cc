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
#include "plants.h"
#include <qmenubar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include "constants.h"
#include "plantsconfig.h"






Plants::Plants (QWidget *parent, const char *name, Data *dat )
   :    Q3VBox(parent, name)//, filename( _filename ), fileinfo( filename )
{


    create_menu ();


    setWindowTitle(QString("PLANTS FE v-") + QString(VERSION.c_str()));
    data = dat;
    mainwin = new PlantsMainWin ( this);
    config = new Configuration (this);
    config->Default();
    config->Apply();
    show ();
}



void Plants::create_menu () {
    create_menu_actions ();

//    QPixmap p1 ( );
//    QPixmap p2  ( );
//    QPixmap p3  ( );


    QMenu *file = new QMenu(tr("&File"), this );
    Q_CHECK_PTR( file );

    file -> addAction (openAct);
    file -> addAction (saveasAct);

    QMenu *settings = new QMenu(tr("&Settings"), this );
    Q_CHECK_PTR( settings );

    settings -> addAction (resetAct);

    QMenu *about = new QMenu(tr("?"), this );
    Q_CHECK_PTR( about );
    about -> addAction (aboutAct);

    QMenuBar *menu = new QMenuBar( this );
    Q_CHECK_PTR( menu );
    menu-> addMenu( file );
    menu-> addMenu( settings );
    menu-> addMenu( about );


}



void Plants::create_menu_actions () {
    openAct = new QAction (tr("&Open File"), this);
    openAct->setShortcut (tr ("Control+O"));
    connect (openAct, SIGNAL (triggered ()), this, SLOT (load ()));

    saveasAct = new QAction (tr("&Open File"), this);
    saveasAct->setShortcut (tr ("Control+O"));
    connect (saveasAct, SIGNAL (triggered ()), this, SLOT (save ()));

    resetAct = new QAction (tr("Reset to &Default"), this);
    resetAct->setShortcut (tr ("Control+D"));
    connect (resetAct, SIGNAL (triggered ()), this, SLOT (reset ()));


    aboutAct = new QAction (tr("&Open File"), this);
    aboutAct->setShortcut (tr ("Control+O"));
    connect (aboutAct, SIGNAL (triggered ()), this, SLOT (about ()));

}



void Plants::load (){
     QString s = QFileDialog::getOpenFileName(this, 
                    tr ("Open file"), "",tr("PLANTS Config File (*.pcfg);;All files (*)"));
    if (!s.isNull()) {
        config->Default();
        config->Read(s.toStdString ().c_str ());
        config->Apply();
    }
    
}

void Plants::reset (){
        config->Default();
        config->Apply();

}

void Plants::save (){
     QString s = QFileDialog::getSaveFileName(this, 
                    tr ("Save configuration"), "",tr("PLANTS Config File (*.pcfg)"));

    config->Get();
    config->Write(s.toStdString ().c_str ());
}

void Plants::about ()
{
    QMessageBox::about( this, QString("About PLANTS Front End") ,
			QString("PLANTS FE. Version ") + QString(VERSION.c_str()) + QString("\n")
			+ QString("Code by Nicola Zonta. All rights reserved.\n")
			+ QString("For bug reports, suggestions || any other feedback please email me at zontan@cf.ac.uk\n\nMAY THE ANT BE WITH YOU!") );
}


