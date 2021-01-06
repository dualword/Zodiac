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

#ifdef LINUX
#include "/usr/local/include/cwiid.h"




int wii ();

cwiid_wiimote_t *init_wiimote (cwiid_wiimote_t *wiimote);
void get_IR_data (cwiid_wiimote_t *wiimote, vector <int> &x, vector <int> &y);

void set_led_state (cwiid_wiimote_t *wiimote, unsigned int n);
bool is_A_pressed (cwiid_wiimote_t *wiimote);



#else
class cwiid_wiimote_t;
inline int wii () {return 0;};

inline cwiid_wiimote_t *init_wiimote (cwiid_wiimote_t *wiimote) {return 0;};
inline void get_IR_data (cwiid_wiimote_t *wiimote, std::vector <int> &x, std::vector <int> &y) {};

inline void set_led_state (cwiid_wiimote_t *wiimote, unsigned int n) {};
inline bool is_A_pressed (cwiid_wiimote_t *wiimote) {return false;};
inline bool cwiid_close (cwiid_wiimote_t *wiimote) {return true;};
#endif //LINUX

#endif //WIIMOTE_H
