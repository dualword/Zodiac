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

#ifndef WIIMOTE_H
#define WIIMOTE_H
#include <vector>
#include "maths.h"

using namespace std;

#ifdef CWIID
#include "/usr/local/include/cwiid.h"

typedef cwiid_wiimote_t Wiimote;

#elif WIIUSE
#include <wiiuse.h>
typedef wiimote_t Wiimote;


#else
typedef int Wiimote;



#endif //CWIID


int wii ();

Wiimote *init_wiimote (Wiimote *wiimote);
bool is_wiimote_closed (Wiimote *wiimote);

vect get_Acc_data (Wiimote *wiimote);


void switch_led (Wiimote *wiimote, unsigned int n, bool on);
void set_rumble (Wiimote *wiimote, bool b);
void get_IR_data (Wiimote *wiimote, vector <int> &x, vector <int> &y);

bool is_A_pressed (Wiimote *wiimote);
bool is_B_pressed (Wiimote *wiimote);
bool is_1_pressed (Wiimote *wiimote);
bool is_2_pressed (Wiimote *wiimote);

void poll (Wiimote *wiimote);



#endif //WIIMOTE_H
