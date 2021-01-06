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

#ifndef FUNCTION_H
#define FUNCTION_H

//#include "constants.h"
#include <vector>

using namespace std;

typedef struct {
	float* value;
	float min_val, max_val;
} Variable;

class Function {
	public:
	Function () {evaluations = 0;};
	~Function () {};
	virtual float evaluate () = 0;
	vector <Variable*>& access_variables() {return _variables;};

protected:
	vector <Variable *> _variables;
	unsigned int evaluations;
	
};


class Square : public Function {
	public:
	Square () 
	{
	//	values.resize(1);
		Variable* v = new Variable;
	//	v->value = &values[0];
		*v->value = 2.0;
		_variables.push_back(v);
		
	};
	
	Square (vector<float>& vv) 
	{
		unsigned int i;
	//	values.resize(vv.size());
		for (i = 0; i < vv.size(); i++) {
			Variable* v = new Variable;
			v->value = &vv[i];
			*v->value = vv[i];
			_variables.push_back(v);
		}
		
	};
	
	~Square () {
		unsigned int i;
		
		for (i = 0; i < _variables.size(); i++) {
			if (_variables[i]) {
				delete _variables[i];
			}
		}
	};
	
	float evaluate ()
	{ 
		evaluations++;
		float v = 0.0;
		unsigned int i;
		for (i = 0; i < _variables.size(); i++) {
			v += *_variables[i]->value * *_variables[i]->value;
		}
		return v;
	};

private:
//	vector<float> values;
	
};




#endif //FUNCTION_H
