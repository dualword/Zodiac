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
#ifndef PLANTS_H
#define PLANTS_H


#include <qfileinfo.h>
#include "plants_mainwin.h"
#include <Qt3Support/q3vbox.h>
#include "plantsconfig.h"
#include "ZNdata.h"
class Data;
class Configuration;


class Plants : public Q3VBox 
{
   Q_OBJECT

public:

    Plants ( QWidget *parent, const char *name, Data* dat );
    PlantsMainWin *mainwin;
    Data *data;
    Configuration *config;
private :

    void create_menu ();
    void create_menu_actions ();

    QAction *openAct, *saveasAct, *resetAct, *aboutAct;

private slots:

    void reset ();
    void save ();

public slots:
    void about ();
    void load ();

};
#endif
