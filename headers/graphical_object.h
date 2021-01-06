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

#ifndef GRAPHICAL_OBJECT_H
#define GRAPHICAL_OBJECT_H



#include "constants.h"
#include "MarchingCubes.h"
#include "FF.h"

using namespace std;
class MarchingCubes;


typedef struct {
    SurfVertex* v1;
    SurfVertex* v2;
    SurfVertex* v3;
} SurfFace;

class Map;

class GraphicalObject {
public:
    GraphicalObject ();
    string name;
    int lst;
    inline void set_name (string nam) {name = nam;};
    bool mesh;
	bool _clippable;
    virtual bool is_surface ();
    virtual bool is_sphere ();
    virtual bool is_grid ();
	virtual bool is_map ();
	virtual void render () {};
	virtual 	void color_by_color (color c) {};
	virtual void color_by_mol (ZNMolecule *mol, float a = 1.) {};
	virtual void color_by_potential (ZNMolecule *mol, color lipo, color hb_acc, color hb_don, float threshold = -4.f) {};
	virtual void color_by_map (Map *map, color c1, color c2, color c3, float f1, float f2, float f3) {}	
	virtual void alpha_by_mol_distance (ZNMolecule *mol, float distance = 4.f) {};
	virtual void multiply_alpha (double d) {};
	virtual void set_alpha (double d) {};
	bool is_clippable () {return _clippable;}
	void set_clippable (bool b) {_clippable = b;}
};

class Surface : public GraphicalObject {

public:
    Surface ();
    void color_by_atom (float a = 1.f);
	void color_by_color (color c);
	void color_by_mol (ZNMolecule *mol, float a = 1.f);
	void color_by_map (Map *map, color c1, color c2, color c3, float f1, float f2, float f3);
	void color_by_potential (ZNMolecule *mol, color lipo, color hb_acc, color hb_don, float threshold = -4.f);
	void alpha_by_mol_distance (ZNMolecule *mol, float distance = 4.f);
	virtual void multiply_alpha (double d) ;
	virtual void set_alpha (double d) ;
    void render ();
    bool is_surface ();

	float near_to_dist;
    vector <SurfFace *> faces;
    vector <SurfVertex *> vertices;
    void set_molecule (ZNMolecule *mo, ZNMolecule *near_to);
	ZNMolecule *near_to;
    ZNMolecule *molecule;
    cutoffGrid<Atom*> *grid;
    cutoffGrid<Atom*> *near_to_grid;
	
protected:
	void fade_vertex_alpha_by_atom (SurfVertex *vert, float d);
    void color_vertex_by_atom (SurfVertex *vert, float a); 
    void color_vertex_by_atom_smooth (SurfVertex *vert, float a);
	void color_vertex_by_color (SurfVertex *vert, color c);

	void set_vertex_alpha (SurfVertex *vert, float a);
	void multiply_vertex_alpha (SurfVertex *vert, float a);
	void color_vertex_by_map (SurfVertex *vert, Map *map, color c1, color c2, color c3, float f1, float f2, float f3);

    void render_as_mesh ();
    void render_as_surface ();

};


class Map : public Surface {
public:
    Map ();
	void clean ();
	bool is_map ();
	void render ();
	int type;
//	vector <SurfFace *> faces;
//    vector <SurfVertex *> vertices;
	vect site_center;
	float site_radius;
//	MarchingCubes *_mc;
	ZNMolecule *molecule;
	float threshold, resolution;
	MarchingCubes *cube;
	color solid_color;
	float get_value (float x, float y, float z) {if (cube) return cube -> get_interpolated_data(x, y, z); else cerr <<"no no no"<<endl; return 0.f;};

};


class Sphere : public GraphicalObject {

public:
    bool is_sphere ();
    inline void set_color (float r, float g, float b, float a) {col = color (r, g, b, a);};
    inline void set_color (color c) {col = c;};
//    inline void set_center (float x, float y, float z) {center = vect (x, y, z);};
    void set_center (vect cent) {center = cent;};
    void set_radius (float r) {radius = r;};

    void render ();
    void render_as_surface ();
    

    vect center;
    color col;
    float radius;

private:



};



class Grid : public GraphicalObject {
    public:
    Grid ();
    bool is_grid ();
    
    vect origin;
    vect end;
    float res;
    float threshold;
    void load ();
    ZNMolecule *probe;
    vector <ZNMolecule *> env;

    void init_ff ();
    inline void set_env (vector <ZNMolecule *> mols) {env = mols;};
    float (*function)(float, float, float);
    float null (float x, float y, float z);
    MarchingCubes *cube;
    ForceField *ff;
    
};

#endif
