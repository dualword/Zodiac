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


#ifdef HAPTICS
#include <HD/hd.h>
#endif //HAPTICS

#include "qwidget.h"
#include <qapplication.h>
#include "ZNdata.h"
#include <string>
#include "constants.h"
#include "ddwin.h"
#include "MMFF.h"

#include <QTranslator>
#include <QLocale>
#ifdef HAPTICS
//#include "haptics.h"



static HHD ghHD = HD_INVALID_HANDLE;
static HDSchedulerHandle gSchedulerCallback = HD_INVALID_HANDLE;
HDSchedulerHandle gCallbackHandle = 0;
HDCallbackCode HDCALLBACK Haptics(void *);
#endif //HAPTICS

int main( int argc, char **argv )
{

// wiimote_kickstart() ;


    QApplication a( argc, argv );

  //QString locale = QLocale::system().name();
 // QTranslator translator;
 // translator.load(QString("zodiac_") + locale);
 // a.installTranslator(&translator);




    Data dat (&a);

    DDWin ddwin (0, &dat);
	dat.load_haptic_device ();

#ifdef HAPTICS
	// - haptics - //
	// Call the haptics initialization function
//	HapticsClass::StaticInit();
	// - haptics - //


    std::cout << "haptics callback" << "\n";
    gSchedulerCallback = hdScheduleAsynchronous(
        Haptics, &dat, HD_MAX_SCHEDULER_PRIORITY);

    HDErrorInfo error;
    ghHD = hdInitDevice(HD_DEFAULT_DEVICE);
    if (HD_DEVICE_ERROR(error = hdGetError()))
    {

        fprintf(stderr, "\nPress any key to quit.\n");
        getchar();
        //exit(-1);
    }
    printf("Found device %s\n",hdGetString(HD_DEVICE_MODEL_TYPE));
    hdEnable(HD_FORCE_OUTPUT);
    hdEnable(HD_MAX_FORCE_CLAMPING);
    hdStartScheduler ();
  //  current_force_x = new double (0.);
   // current_force_y = new double (0.);
   // current_force_z = new double (0.);
    if (HD_DEVICE_ERROR(error = hdGetError()))
    {

        fprintf(stderr, "\nPress any key to quit.\n");
        return -1;
    }
#endif //HAPTICS

    ddwin.show ();

    if (argc >1) {
        for (unsigned int i=1; i<argc; i++) {
            ddwin.load_file (argv[i]);
        }
    }
    int exitCode = a.exec();

#ifdef HAPTICS
	// - haptics - //
//	gHaptics.uninit();
    hdStopScheduler();
    hdUnschedule(gSchedulerCallback);

    if (ghHD != HD_INVALID_HANDLE)
    {
        hdDisableDevice(ghHD);
        ghHD = HD_INVALID_HANDLE;
    }
	// - haptics - //
#endif //HAPTICS

	return exitCode;

}



#ifdef HAPTICS
HDCallbackCode HDCALLBACK Haptics(void *data)
{

    hdBeginFrame(ghHD);
     Data *dat = (Data *) data; 

    static HDboolean bRenderForce = 0;
    HDErrorInfo error;

    HDint nCurrentButtons, nLastButtons;
	double position [3];
	double angles [3];


    hdGetDoublev(HD_CURRENT_GIMBAL_ANGLES,angles);
    hdGetDoublev(HD_CURRENT_POSITION, position);
	dat -> current_position_x = position[0];
	dat -> current_position_y = position[1];
	dat -> current_position_z = position[2];




	dat -> current_pitch = angles[2];
	dat -> current_roll = angles[0];
	dat -> current_yaw = angles[1];


//	cerr << dat ->current_pitch << endl;
double force [3];
force [0] = dat -> current_force_x;
force [1] = dat -> current_force_y;
force [2] = dat -> current_force_z;

	hdSetDoublev(HD_CURRENT_FORCE, force);
    hdEndFrame(ghHD);
	return HD_CALLBACK_CONTINUE;

}
#endif // HAPTICS


