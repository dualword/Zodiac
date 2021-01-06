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

#ifndef DATA_H
#define DATA_H

#include "MMFF.h"
#include "PLP.h"
#include "chemscore.h"
#include "molecule.h"
#include "ddwin.h"
#include <qapplication.h>
#include "minimize.h"
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <QUndoStack>
#include "actions.h"



class Myline;
class MylineF;

//using namespace std;

class Actions;
class Data;

class DDWin;
class Minimize;




class Data {

public:
    Data (QApplication *master);//, DDWin *ddw);
 
    void set_ddwin (DDWin *ddw);
 //   void set_FF (MMFF *mmff, TriposFF *triposff, PLP *plp);
    QApplication *qapp;
    DDWin *ddwin;
	Actions *actions;
    MMFF *mmff;
    Minimize *minimize ;
 //   PLP *plp;
    Molecule *protein;
    vector <Molecule*> ligands;
    vector <vector <float> > waters;

    QUndoStack * undo_stack;



    color constant_color;
    color charge_begin_color, charge_end_color;
    color score_begin_color, score_mid_color, score_end_color;


    float charge_begin, charge_end;
    float score_begin, score_mid, score_end;


    Molecule* check_mol2 (string filename);



};

bool is_db_extended (Molecule *mol);


#endif
