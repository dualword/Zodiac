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

#include "constants.h"

class DataGrid {
private:
	float resolution;
    float *_data;
	int _x_size, _y_size, _z_size;
	float x_length, y_length, z_length;
	float _min_x, _min_y, _min_z;
public:
	DataGrid (vect min_corn, vect max_corner, float res) {
		resolution = res;
		_min_x = min_corn.x ();
		_min_y = min_corn.y ();
		_min_z = min_corn.z ();
		x_length = max_corner.x() - _min_x;
		y_length = max_corner.y() - _min_y;
		z_length = max_corner.z() - _min_z;
		if (x_length < 0.) x_length = - x_length;
		if (y_length < 0.) y_length = - y_length;
		if (z_length < 0.) z_length = - z_length;
		_x_size = (int) (x_length / resolution) + 1 + 1; 		
		_y_size = (int) (y_length / resolution) + 1 + 1;
		_z_size = (int) (z_length / resolution) + 1 + 1;
		_data = new float [_x_size *_y_size * _z_size];
		for (unsigned int i = 0; i < _x_size  * _y_size  * _z_size  ; i++) _data[i] = 0.;
	};
	~DataGrid () {delete[] _data;};
	void set_value (int i, int j, int k, float val) {if (i*j*k < (_x_size -1) * (_y_size-1) * (_z_size-1)) {_data[i*(_z_size*_y_size) + j * (_z_size) + k  ] = val;}};
	int x_size () {return _x_size;};
	int y_size () {return _y_size;};
	int z_size () {return _z_size;};
	
	
	int get_i (float x) {
		int i = 0;
		float newx = x - _min_x;
		if (0 <= newx  && newx<= ( x_length)) i = (int) (newx / resolution);
	//			cerr << "get i " << x << " " << newx << " " << i << endl;
		return i;
	};

	int get_j (float y) {
		int j = 0;
		float newy = y - _min_y;
		if (0 <= newy && newy<= (y_length)) j = (int) (newy / resolution);
		//cerr << _min_y << "  " << resolution << " "<<newy << " " << y<<" "<<y_length<<" "<<j<<endl;
		return j;
	};
	int get_k (float z) {
		int k = 0;
		float newz = z - _min_z;

		if (0 <= newz  && newz<= ( z_length)) k = (int) (newz / resolution);
				//		cerr << _min_z << "  " << resolution << " "<<newz << " " << k<<" "<<z_length<<endl;
		return k;
	};
	float get_x (int i) {
		float x = _min_x + i*resolution;
	//				cerr << "get x " << i << " " << _min_x + x_length << " " << x << endl;
		if (_min_x <= x  && x<= (_min_x + x_length)) return x;

		else return _min_x;
	};
	float get_y (int j) {
		float y = _min_y + j*resolution;
		if (_min_y <= y && y<= _min_y + y_length) return y;
		else return _min_y;
	};
	float get_z (int k) {
		float z = _min_z + k*resolution;
		if (_min_z <= z && z<= _min_z + z_length) return z;
		else return _min_z;
	};
	float get_value (float x, float y, float z) {
	
		if (x < _min_x) x = _min_x;
		if (y < _min_y) y = _min_y;
		if (z < _min_z) z = _min_z;

		if (x > _min_x + x_length) x = _min_x + x_length;
		if (y > _min_y + y_length) y = _min_y + y_length;
		if (z > _min_z + z_length) z = _min_z + z_length;
							
		int i = get_i (x);
		int j = get_j (y);
		int k = get_k (z);
		float xl = get_x (i);
		float yl = get_y (j);
		float zl = get_z (k);

		int ii, jj, kk;
		float dx, dy, dz;
		if (i < _x_size -1) 
			ii = i + 1;
		
		else ii = i;
		if (j < _y_size -1) jj = j + 1;
		else ii = i;
		if (k < _z_size -1) kk = i + 1;
		else kk = k;	
		
		dx = (x -xl) / resolution;
		dy = (y - yl) / resolution;
		dz = (z - zl) / resolution;
		
	//	cerr << x << " " << i << " " << xl << endl;
	//	cerr << y << " " << j << " " << yl << endl;
	//	cerr << z << " " << k << " " << zl << endl;
										
		float V000 = _data [i*(_z_size*_y_size) + j * (_z_size) + k ];
		float V001 = _data [i*(_z_size*_y_size) + j * (_z_size) + kk ];
		float V010 = _data [i*(_z_size*_y_size) + jj * (_z_size) + k ];
		float V011 = _data [i*(_z_size*_y_size) + jj * (_z_size) + kk ];
		float V100 = _data [ii*(_z_size*_y_size) + j * (_z_size) + k ];
		float V101 = _data [ii*(_z_size*_y_size) + j * (_z_size) + kk ];
		float V110 = _data [ii*(_z_size*_y_size) + jj * (_z_size) + k ];
		float V111 = _data [ii*(_z_size*_y_size) + jj * (_z_size) + kk ];
		
	//	cerr << "values "<< dx<<" " << dy<<" " <<dz << " "<< V000 << " "<< V001 << " "<< V010 << " "<< V011 << " "<< V100 << " "<< V101 << " "<< V110 << " "<< V111 << " "<<endl;
		
		return	V000 * (1 - dx) * (1 - dy) * (1 - dz) +
		V100 * dx       * (1 - dy) * (1 - dz) +
		V010 * (1 - dx) * dy       * (1 - dz) +
		V001 * (1 - dx) * (1 - dy) * dz +
		V101 * dx       * (1 - dy) * dz +
		V011 * (1 - dx) * dy       * dz +
		V110 * dx       * dy       * (1 - dz) +
		V111 * dx       * dy       * dz; 
	};
	

	
};

