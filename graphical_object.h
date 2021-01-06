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

#include <string>
#include <vector>
#include "molecule.h"
#include "cutoffGrid.h"
#include "constants.h"
#include "MarchingCubes.h"
#include "FF.h"

using namespace std;
class MarchingCubes;

class SurfVertex
 {
    public:
  //  SurfVertex ();
    inline vect& GetVector () {return coordinates;};
    vect coordinates;
    color col;
    vect normal;
    int n;
}; 

typedef struct {
    SurfVertex* v1;
    SurfVertex* v2;
    SurfVertex* v3;
} SurfFace;



class GraphicalObject {
public:
    GraphicalObject ();
    string name;
    int list;
    inline void set_name (string nam) {name = nam;};
    bool mesh;
    virtual bool is_surface ();
    virtual bool is_sphere ();
    virtual bool is_grid ();
};

class Surface : public GraphicalObject {

public:
    Surface ();
    void color_by_atom (float a = 1.f);
    void render ();
    bool is_surface ();


    vector <SurfFace *> faces;
    vector <SurfVertex *> vertices;
    void set_molecule (Molecule *mo);
    Molecule *molecule;
    cutoffGrid<Atom*> *grid;
private:

    void color_vertex_by_atom (SurfVertex *vert, float a); 
    void color_vertex_by_atom_smooth (SurfVertex *vert, float a); 


    void render_as_mesh ();
    void render_as_surface ();

};

class Sphere : public GraphicalObject {

public:
    bool is_sphere ();
    inline void set_color (float r, float g, float b, float a) {col = color (r, g, b, a);};
    inline void set_color (color c) {col = c;};
    inline void set_center (float x, float y, float z) {center = vect (x, y, z);};
    inline void set_center (vect cent) {center = cent;};
    inline void set_radius (float r) {radius = r;};

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
    float treshold;
    void load ();
    Molecule *probe;
    vector <Molecule *> env;

    void init_ff ();
    inline void set_env (vector <Molecule *> mols) {env = mols;};
    float (*function)(float, float, float);
    float null (float x, float y, float z);
    MarchingCubes *cube;
    ForceField *ff;
    
};

#endif
