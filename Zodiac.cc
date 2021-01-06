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

#include "qwidget.h"
#include <qapplication.h>
#include "ZNdata.h"
#include <string>
#include "constants.h"
#include "ddwin.h"
#include "MMFF.h"

#ifdef HAPTICS
#include "haptics.h"
#endif //HAPTICS

int main( int argc, char **argv )
{
    QApplication a( argc, argv );

#ifdef HAPTICS
	// - haptics - //
	// Call the haptics initialization function
	HapticsClass::StaticInit();
	// - haptics - //
#endif //HAPTICS

    Data dat (&a);
    DDWin ddwin (0, &dat);




    ddwin.show ();

    if (argc >1) {
        for (unsigned int i=1; i<argc; i++) {
            ddwin.load_file (argv[i]);
        }
    }
    int exitCode = a.exec();

#ifdef HAPTICS
	// - haptics - //
	gHaptics.uninit();
	// - haptics - //
#endif //HAPTICS

	return exitCode;

}
