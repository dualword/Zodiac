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


#ifndef ILS_H
#define ILS_H

#include "optimiser.h"
#include "nms.h"

class ILS : public Optimiser {
public:
	ILS () : Optimiser () {};	
	ILS (Function* f) : Optimiser (f) 
	{
		optFunc = f;
				set_iteration_limit (0);
		
	};	
	~ILS () {};
	
	float run () 
	{
		bool optimising = true;
		minimiser = new NMS (optFunc);

	//	minimiser -> set_iteration_limit (20);
		float variation_range = (float)M_PI;
		unsigned int i = 0, j = 0;
		vector<float> bestSolution;
		float bestSolutionValue;
		float tmpSolutionValue;
		vector<Variable*>& vars = optFunc->access_variables();
		for (j = 0; j < vars.size(); j++) {
			bestSolution.push_back(*vars[j]->value);
		}
		bestSolutionValue = optFunc->evaluate();
		start_timer ();
		while (!reached_limits ()) {
			for (j = 0; j < vars.size(); j++) {
				do {
					float variation_range =  (vars[j]->max_val - vars[j]->min_val)/2;

					*vars[j]->value = bestSolution[j] -variation_range + 2 * variation_range*(float)rand()/(float)RAND_MAX;
										//cerr << *vars[j]->value << " "<<vars[j]->min_val <<" "<< vars[j]->max_val<<endl;
				}	
				while ((*vars[j]->value < vars[j]->min_val) || (*vars[j]->value > vars[j]->max_val));
			}
			if (optimising) {
				minimiser -> init (vars);
				minimiser -> run ();
				
			}
			bool accept = true;
			for (j = 0; j < vars.size(); j++) {
				if ((*vars[j]->value < vars[j]->min_val) || (*vars[j]->value > vars[j]->max_val)) {
					accept = false;
					break;
				}
			}
			
			if (accept) {
			tmpSolutionValue = optFunc->evaluate();
			conformation *solution = new conformation;
			solution ->score = tmpSolutionValue;
			for 	(j = 0; j < vars.size(); j++) {
				solution ->state.push_back (*vars[j]->value);
			}	
			results.push_back (solution);
			if ( tmpSolutionValue < bestSolutionValue) {
				for (j = 0; j < vars.size(); j++) {
					bestSolution[j] = (*vars[j]->value);
				}
				bestSolutionValue = tmpSolutionValue;
				cout << "BEST: " << tmpSolutionValue << endl;
			}
			iterations_n ++;
			}
		}
		for (j = 0; j < vars.size(); j++) {
			*vars[j]->value = bestSolution[j];
		}
		tmpSolutionValue = optFunc->evaluate();	
		
		sort (results.begin (), results.end (), conformationScoreCompare ());
	/*	for (j = 0; j < results.size(); j++) {
			cerr << results[j] ->score<<"    ";
			for (i=0; i< results[j] ->state.size (); i++) {
				cerr << results[j]->state[i];
			}
			cerr << endl;
		}*/
		return tmpSolutionValue;
	};
	
private:
	Optimiser *minimiser;
	
	
};

#endif //ILS_H