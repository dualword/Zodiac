
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

#ifndef OPTIMISER_H
#define OPTIMISER_H

#include "function.h"
#include <algorithm>
//#include <QTime>



typedef struct {
    vector  <float> state;
    float score;
} conformation;


class Optimiser {
	public:
	Optimiser () : time_limit (0) { iteration_limit = -1; optFunc = NULL; iterations_n = 1;};
	Optimiser (Function *f) : time_limit (0) {optFunc = f;  iteration_limit = -1; iterations_n = 1;};	
	~Optimiser () {};
	virtual void init (vector<Variable*>& vars) {};
	virtual float run () = 0;
	virtual void set_iteration_limit (int limit) {iteration_limit = limit; iterations_n = 0;};
	vector <conformation *> &get_results () {return results;};
	void set_time_limit (int t) {time_limit = t;};
	

protected:
		Function* optFunc;
		vector <conformation *> results;
		QTime start_time;
		int time_limit, iteration_limit;
		int iterations_n;
		bool reached_limits () {
	//		cerr << "iteration_limit "<<iteration_limit<<" iteration_n "<< iterations_n<<" time_limit "<<time_limit<<" elapsed_time "<<start_time.secsTo ( QTime::currentTime () )<<endl;
			bool t = true; bool it = true;
			if (!time_limit) {
				t = false;
			} 
			else {
				int elapsed_time = start_time.secsTo ( QTime::currentTime () );
				if (elapsed_time > time_limit) t = true;
				else t = false;
			}
			if (iteration_limit < 1) {
				it = false;
			} 
			else {
				if (iterations_n > iteration_limit) it = true;
				else it = false;
			}
		//	if (t ) cerr << "reachced time limits"<<endl;
		//	if (it ) cerr << "reachced iterations limits "<<iterations_n << " "<<iteration_limit<<endl;			
				return (t || it);
		}
	void start_timer () {
		start_time = QTime::currentTime ();
	}
};



struct conformationScoreCompare
{
    bool operator()(const conformation* a, const conformation* b)
    {
           return a->score < b->score;
    }
};

#endif //OPTIMISER_H


