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


#ifndef MYLINE_H
#define MYLINE_H

#include <Qt3Support/q3hbox.h>
#include <qlineedit.h>
#include <iostream>
#include <qpushbutton.h>
#include <Qt3Support/q3tabdialog.h>
#include <string>


using namespace std ;


class QLineEdit;

class MyLine : public Q3HBox
{
    Q_OBJECT

public:
    MyLine (QWidget *parent, const char *name, int valid);

 //   const char* get_value ();
    void ins (QString s);
    string val ();


    QLineEdit *linedit;

};


class MyLineF : public Q3HBox
{
    Q_OBJECT

public:
    MyLineF( QWidget *parent, QWidget *master, const char *name, int valid, int butt=0);
    const char* get_value ();
    QPushButton *adbutt; 
    QLineEdit *linedit;
    void ins (QString s);
    string val ();



protected:
    const char *tag;
    QPushButton *fbutton; 
    char *filetype;
//    MainWin *master;

protected slots:
   QString ask_file ();
    void set_file ();  
    void set_line ();


    
};




#endif
