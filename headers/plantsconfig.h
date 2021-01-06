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
#ifndef PLANTSCONFIG_H
#define PLANTSCONFIG_H


#include "plantsconfig.h"
#include "plants.h"
#include "plants_mainwin.h"

class PlantsMainWin;
class Plants;

class Confline
{
public:
    Confline (char type,  string filestr,string defval, MyLine* targetl);
    Confline (char type,  string filestr, string defval,QPushButton* targetb);
    Confline (char type,  string filestr, string defval, QLineEdit* targetle);
    Confline (char type,  string filestr, string defval, MyLineF* targetlf);
    Confline (char type,  string filestr, string defval, QComboBox* targetcb);
    Confline (char type,  string filestr, string defval, Q3ListBox* targetlb);
    Confline (char type,  string filestr, string defval, PlantsMainWin* targetlb);

    void Default ();
    void Apply ();
    void Get ();
    int getInt (string mystr);


   char type;
    string defval, filestr, current;
    MyLine* targetl;
    QPushButton* targetb;
    QLineEdit* targetle;
    MyLineF* targetlf;
    QComboBox* targetcb;
    Q3ListBox* targetlb;
    PlantsMainWin* targetmw;
    string boolstr (bool boo);

    vector<string> vect;


 //   void Read ();
 //   void Write ();

};














class Configuration 
{

public:
    Configuration (Plants *plant);
//    void setdata (Data *datas);
    void Apply ();
    void Default ();
    void Write (const char * filename);
    void Read (const char * filename);
    void Get ();
//void Load (File* file);


private:
    Plants *plants; 




//    vector <string> ifiles, vflex, vwater, vconstr;
    vector <Confline*> conf;
    



};


#endif
