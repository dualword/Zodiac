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


#ifndef PSO_H
#define PSO_H

#include "optimiser.h"
#include "nms.h"


class PSO : public Optimiser {
public:
	PSO () : Optimiser () {};	
	PSO (Function* f) : Optimiser (f) 
	{
		optFunc = f;
		
	};	
	~PSO () {};
	
	float run () 
	{
		
		
		unsigned int i;
		unsigned int j;
		vector<Variable*>& vars = optFunc->access_variables();
		
		float bestSolutionValue = numeric_limits<float>::max();
		vector<float> maxV(vars.size());
		unsigned int nofFunc = 0;
		float r1;
		float r2;
		
		vector<vector<float> > x;
		
		bool optimising = true;
		minimiser = new NMS (optFunc);
		minimiser ->set_iteration_limit (20);
		
		vector<float> fx;
		vector<float> bestState;
		vector<float> currentState;
		vector<vector<float> > v;
		
		vector<vector<float> > p;
		vector<float> fp;
		
		unsigned int nofParticles = 20;
		
		x.resize(nofParticles);
		v.resize(nofParticles);
		p.resize(nofParticles);
		fx.resize(nofParticles);
		fp.resize(nofParticles);
		
		bestState.resize(vars.size());
		
		
		for (i = 0; i < vars.size(); i++) {
			if (i < 3) {
				maxV[i] = 2.0;
			}
			else {
				maxV[i] = 50.0;
			}
		}
		
		for (i = 0; i < nofParticles; i++) {
			x[i].resize(vars.size());
			v[i].resize(vars.size());
			p[i].resize(vars.size());
			
			for (j = 0; j < vars.size(); j++) {
				x[i][j] = (float)rand()/(float)RAND_MAX * (vars[j]->max_val - vars[j]->min_val) + vars[j]->min_val;
			}
			p[i] = x[i];
			
			for (j = 0; j < vars.size(); j++) {
				*vars[j]->value = x[i][j];
			}
			
			fx[i] = optFunc->evaluate();	
			fp[i] = fx[i];
			
			nofFunc++;
			
			
			if (fx[i] < bestSolutionValue) {
				bestSolutionValue = fx[i];
				bestState = x[i];
			}
			
			for (j = 0; j < vars.size(); j++) {
				v[i][j] = -maxV[j] + (float)rand()/(float)RAND_MAX * (2.0 * maxV[j]);
			}
		}
		
		//unsigned int iterations = 0;
		start_timer ();

		while (!reached_limits ()) {
			iterations_n ++;
			currentState = bestState;
			for (i = 0; i < nofParticles; i++) {
				
				if (fx[i] < fp[i]) {
					fp[i] = fx[i];
					p[i] = x[i];
				}
				
				for (j = 0; j < vars.size(); j++) {
					
					r1 = (float)rand()/(float)RAND_MAX;
					r2 = (float)rand()/(float)RAND_MAX;
					
					
					v[i][j] = iw*v[i][j] + ci*r1*(p[i][j]-x[i][j]) + cg*r2*(currentState[j]-x[i][j]);
					
					if (v[i][j] > maxV[j]) {
						v[i][j] = maxV[j];
					}
					if (v[i][j] < -maxV[j]) {
						v[i][j] = -maxV[j];
					}
					
					x[i][j] += v[i][j];
				}
				for (j = 0; j < vars.size(); j++) {
					*vars[j]->value = x[i][j];
				}
				
				if (optimising) {
					minimiser -> init (vars);
					minimiser -> run ();
					
				}
				fx[i] = optFunc->evaluate();	
				conformation *solution = new conformation;
				solution ->score = fx[i];
				for 	(j = 0; j < vars.size(); j++) {
					solution ->state.push_back (*vars[j]->value);
				}	
				results.push_back (solution);
				for (j = 0; j < vars.size(); j++) {
					x[i][j] = *vars[j]->value;
				}
				
				nofFunc++;
				
				if (fx[i] < bestSolutionValue) {
					bestState = x[i];
					bestSolutionValue = fx[i];
				}
				
			}
			
			for (j = 0; j < vars.size(); j++) {
				*vars[j]->value = bestState[j];
			}
			
		}
		sort (results.begin (), results.end (), conformationScoreCompare ());
		return bestSolutionValue;
		
	};
	
private:
	float ci;
	float cg;
	float iw;
//	PSOlocalWeight = 1.0;
//	PSOglobalWeight = 1.0;
//	PSOinertiaWeight = 0.6;
	Optimiser *minimiser;
};


#endif //PSO_H