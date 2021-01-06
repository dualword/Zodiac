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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "wiimote.h"


#ifdef CWIID

//cwiid_mesg_callback_t cwiid_callback;

#define toggle_bit(bf,b)	\
	(bf) = ((bf) & b)		\
	       ? ((bf) & ~(b))	\
	       : ((bf) | (b))


void set_rumble (Wiimote *wiimote, bool b) {}

void set_led_state(Wiimote *wiimote, unsigned char led_state);
void set_rpt_mode(Wiimote *wiimote, unsigned char rpt_mode);
void print_state(struct cwiid_state *state);

cwiid_err_t err;
void err(Wiimote *wiimote, const char *s, va_list ap)
{
	if (wiimote) printf("%d:", cwiid_get_id(wiimote));
	else printf("-1:");
	vprintf(s, ap);
	printf("\n");
}


Wiimote *init_wiimote (Wiimote *wiimote) {
	bdaddr_t bdaddr;	/* bluetooth device address */
	unsigned char mesg = 0;
	unsigned char led_state = 0;
	unsigned char rpt_mode = 0;
	unsigned char rumble = 0;
	int exit = 0;

	cwiid_set_err(err);

	/* Connect to address given on command-line, if present */

		bdaddr = *BDADDR_ANY;


	/* Connect to the wiimote */
	printf("Put Wiimote in discoverable mode now (press 1+2)...\n");
	if (!(wiimote = cwiid_open(&bdaddr, 0))) {
		fprintf(stderr, "Unable to connect to wiimote\n");
		return NULL;
	}
//	cerr << wiimote << endl;


	//turns IR on
	toggle_bit(rpt_mode, CWIID_RPT_IR);
	set_rpt_mode(wiimote, rpt_mode);

	//turns ACC on	
	toggle_bit(rpt_mode, CWIID_RPT_ACC);
	set_rpt_mode(wiimote, rpt_mode);


	//enable button reporting
	toggle_bit(rpt_mode, CWIID_RPT_BTN);
	set_rpt_mode(wiimote, rpt_mode);


	toggle_bit(led_state, CWIID_LED1_ON);
	set_led_state(wiimote, led_state);
	return wiimote;
}


void get_IR_data (Wiimote *wiimote, vector <int> &x, vector <int> &y) {
    	struct cwiid_state state;	/* wiimote state */
	vector <int> xx;
	vector <int> yy;
	if (cwiid_get_state(wiimote, &state)) {
	   fprintf(stderr, "Error getting state\n");
	}
	for (unsigned int i = 0; i < CWIID_IR_SRC_COUNT; i++) {
		if (state.ir_src[i].valid) {
			xx.push_back (state.ir_src[i].pos[CWIID_X]);
			yy.push_back (state.ir_src[i].pos[CWIID_Y]);
		}
	}

	x = xx;
	y = yy;
}


vect get_Acc_data (Wiimote *wiimote) {
    	struct cwiid_state state;	/* wiimote state */
	vector <int> xx;
	vector <int> yy;
	if (cwiid_get_state(wiimote, &state)) {
	   fprintf(stderr, "Error getting state\n");
	}
	 vect out (state.acc[CWIID_X]-125, state.acc[CWIID_Y]-125, state.acc[CWIID_Z]-125);
	return out;
}

bool is_A_pressed (Wiimote *wiimote) {
    	struct cwiid_state state;	/* wiimote state */
	if (cwiid_get_state(wiimote, &state)) {
	   fprintf(stderr, "Error getting state\n");
	}

//	cerr << state.buttons << endl;
	return (state.buttons == CWIID_BTN_A);
}


bool is_1_pressed (Wiimote *wiimote) {
    	struct cwiid_state state;	/* wiimote state */
	if (cwiid_get_state(wiimote, &state)) {
	   fprintf(stderr, "Error getting state\n");
	}

//	cerr << state.buttons << endl;
	return (state.buttons == CWIID_BTN_1);
}


bool is_2_pressed (Wiimote *wiimote) {
    	struct cwiid_state state;	/* wiimote state */
	if (cwiid_get_state(wiimote, &state)) {
	   fprintf(stderr, "Error getting state\n");
	}

//	cerr << state.buttons << endl;
	return (state.buttons == CWIID_BTN_2);
}



bool is_B_pressed (Wiimote *wiimote) {
    	struct cwiid_state state;	/* wiimote state */
	if (cwiid_get_state(wiimote, &state)) {
	   fprintf(stderr, "Error getting state\n");
	}

//	cerr << state.buttons << endl;
	return (state.buttons == CWIID_BTN_B);
}



void set_led_state (Wiimote *wiimote, unsigned int n) {
	unsigned char state;
	switch (n) {
	case 1:
		state = CWIID_LED1_ON;
		break;
	case 2:
		state = CWIID_LED2_ON;
		break;
	case 3:
		state = CWIID_LED3_ON;
		break;
	case 4:
		state = CWIID_LED4_ON;
		break;
	default:
		state = CWIID_LED1_ON;
		break;
	}
    unsigned char led_state;
    toggle_bit(led_state, state);
    set_led_state(wiimote, led_state);

}

void set_led_state(Wiimote *wiimote, unsigned char led_state)
{
	if (cwiid_set_led(wiimote, led_state)) {
		fprintf(stderr, "Error setting LEDs \n");
	}
}

void switch_led (Wiimote *wiimote, unsigned int n, bool on) {}
	
void set_rpt_mode(Wiimote *wiimote, unsigned char rpt_mode)
{
	if (cwiid_set_rpt_mode(wiimote, rpt_mode)) {
		fprintf(stderr, "Error setting report mode\n");
	}
}

bool is_wiimote_closed (Wiimote *wiimote) {
	return cwiid_close (wiimote);
}

void poll (Wiimote *wiimote) {};

#elif WIIUSE

void set_rumble (Wiimote *wiimote, bool b) {
	wiiuse_rumble (wiimote, b);
}

void switch_led (Wiimote *wiimote, unsigned int n, bool on) {
	switch (n) {
	case 1:
		wiiuse_set_leds(wiimote, WIIMOTE_LED_1);
		break;
	case 2:
		wiiuse_set_leds(wiimote, WIIMOTE_LED_2);
		break;
	case 3:
		wiiuse_set_leds(wiimote, WIIMOTE_LED_3);
		break;
	case 4:
		wiiuse_set_leds(wiimote, WIIMOTE_LED_4);
		break;
	
	}
}
void get_IR_data (Wiimote *wiimote, vector <int> &x, vector <int> &y) {
	vector <int> xx, yy;
		int i = 0;



		/* go through each of the 4 possible IR sources */

		for (; i < 4; ++i) {

			/* check if the source is visible */

			if (wiimote->ir.dot[i].visible) {
				xx.push_back (wiimote->ir.dot[i].x);
				yy.push_back (wiimote->ir.dot[i].y);



			}
		}
	x = xx;
	y = yy;

}

vect get_Acc_data (Wiimote *wiimote){
	float roll = wiimote ->orient.roll;
	float pitch = wiimote ->orient.pitch;
	float yaw = wiimote ->orient.yaw;
	return vect (roll, pitch, yaw);
}
bool is_A_pressed (Wiimote *wiimote) {return IS_HELD(wiimote, WIIMOTE_BUTTON_A);}
bool is_B_pressed (Wiimote *wiimote) {return IS_HELD(wiimote, WIIMOTE_BUTTON_B);}
bool is_1_pressed (Wiimote *wiimote) {return IS_HELD(wiimote, WIIMOTE_BUTTON_ONE);}
bool is_2_pressed (Wiimote *wiimote) {return IS_HELD(wiimote, WIIMOTE_BUTTON_TWO);}
bool is_wiimote_closed (Wiimote *wiimote) {return true;};

Wiimote *init_wiimote (Wiimote *wiimote) {
	Wiimote **wiimotes =  wiiuse_init(1);
	int found = wiiuse_find(wiimotes, 1, 5);
	if (found) {
		int connected = wiiuse_connect(wiimotes, 1);

		if (connected) {
			Wiimote *wiimote = wiimotes[0];
			wiiuse_motion_sensing (wiimote, 1);
			return wiimote;
		}
		else { return 0;}
	}

	else {

		printf("Failed to connect to any wiimote.\n");

		return 0;
	}

};


void poll (Wiimote *wiimote) {
	Wiimote *wiimotes [1];
	wiimotes [0] = wiimote;
	wiiuse_poll (wiimotes, 1); 
};

#else //CWIID

void switch_led (Wiimote *wiimote, unsigned int n, bool on) {}
void get_IR_data (Wiimote *wiimote, vector <int> &x, vector <int> &y) {}
vect get_Acc_data (Wiimote *wiimote){return vect ();};
bool is_A_pressed (Wiimote *wiimote) {return false;}
bool is_B_pressed (Wiimote *wiimote) {return false;}
bool is_1_pressed (Wiimote *wiimote) {return false;}
bool is_2_pressed (Wiimote *wiimote) {return false;}
bool is_wiimote_closed (Wiimote *wiimote) {return true;};
Wiimote *init_wiimote (Wiimote *wiimote) {return 0;};
void poll (Wiimote *wiimote) {};
void set_rumble (Wiimote *wiimote, bool b) {};

#endif //CWIID
