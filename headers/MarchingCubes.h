/**
 * @file    MarchingCubes.h
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.2
 * @date    12/08/2002
 *
 * @brief   MarchingCubes Algorithm
 */
//________________________________________________


#ifndef MARCHINGCUBES_H
#define MARCHINGCUBES_H

#if !defined(WIN32) || defined(__CYGWIN__)
#pragma interface
#endif // WIN32


#include <iostream>;
using namespace std;
//_____________________________________________________________________________
// types
/** unsigned char alias */
typedef unsigned char uchar ;
/** signed char alias */
typedef   signed char schar ;
/** isovalue alias */
typedef        double real  ;

//-----------------------------------------------------------------------------
// Vertex structure
/** \struct Vertex "MarchingCubes.h" MarchingCubes
 * Position && normal of a vertex
 * \brief vertex structure
 * \param x X coordinate
 * \param y Y coordinate
 * \param z Z coordinate
 * \param nx X component of the normal
 * \param ny Y component of the normal
 * \param nz Z component of the normal
 */
typedef struct
{
  real  x,  y,  z ;  /**< Vertex coordinates */
  real nx, ny, nz ;  /**< Vertex normal */
} Vertex ;

//-----------------------------------------------------------------------------
// Triangle structure
/** \struct Triangle "MarchingCubes.h" MarchingCubes
 * Indices of the oriented triange vertices
 * \brief triangle structure
 * \param v1 First vertex index
 * \param v2 Second vertex index
 * \param v3 Third vertex index
 */
typedef struct
{
  int v1,v2,v3 ;  /**< Triangle vertices */
} Triangle ;
//_____________________________________________________________________________



//_____________________________________________________________________________
/** Marching Cubes algorithm wrapper */
/** \class MarchingCubes
  * \brief Marching Cubes algorithm.
  */
class MarchingCubes
//-----------------------------------------------------------------------------
{
// Constructors
public :
  /**
   * Main && default constructor
   * \brief constructor
   * \param size_x width  of the grid
   * \param size_y depth  of the grid
   * \param size_z height of the grid
   */
  MarchingCubes ( const int size_x = -1, const int size_y = -1, const int size_z = -1 ) ;
  /** Destructor */
  ~MarchingCubes() ;

//-----------------------------------------------------------------------------
// Accessors
public :
  inline const float xmin () const {return _xmin;}
  inline const float ymin () const {return _ymin;}
  inline const float zmin () const {return _zmin;}
  inline const float xmax () const {return _xmax;}
  inline const float ymax () const {return _ymax;}
  inline const float zmax () const {return _zmax;}


  /** accesses the number of vertices of the generated mesh */
  inline const int nverts() const { return _nverts ; }
  /** accesses the number of triangles of the generated mesh */
  inline const int ntrigs() const { return _ntrigs ; }
  /** accesses a specific vertex of the generated mesh */
  inline Vertex   * vert( const int i ) const { if( i < 0  || i >= _nverts ) return ( Vertex *)NULL ; return _vertices  + i ; }
  /** accesses a specific triangle of the generated mesh */
  inline Triangle * trig( const int i ) const { if( i < 0  || i >= _ntrigs ) return (Triangle*)NULL ; return _triangles + i ; }

  /** accesses the vertex buffer of the generated mesh */
  inline Vertex   *vertices () { return _vertices  ; }
  /** accesses the triangle buffer of the generated mesh */
  inline Triangle *triangles() { return _triangles ; }

  /**  accesses the width  of the grid */
  inline const int size_x() const { return _size_x ; }
  /**  accesses the depth  of the grid */
  inline const int size_y() const { return _size_y ; }
  /**  accesses the height of the grid */
  inline const int size_z() const { return _size_z ; }

  inline const float to_real_x(float i) const { return _xmin+i*(_xmax-_xmin)/_size_x ; }
  inline const float to_real_y(float j) const { return _ymin+j*(_ymax-_ymin)/_size_y ; }
  inline const float to_real_z(float k) const { return _zmin+k*(_zmax-_zmin)/_size_z ; }

  inline const int to_cube_i(float x) const { return (int) ((x-_xmin)/(_xmax-_xmin)*_size_x);}
  inline const int to_cube_j(float y) const { return (int) ((y-_ymin)/(_ymax-_ymin)*_size_y);}
  inline const int to_cube_k(float z) const { return (int) ((z-_zmin)/(_zmax-_zmin)*_size_z);}
  /**
   * changes the size of the grid
   * \param size_x width  of the grid
   * \param size_y depth  of the grid
   * \param size_z height of the grid
   */
  inline void set_limits (const float xmin, const float ymin, const float zmin, const float xmax, const float ymax, const float zmax) {_xmin=xmin; _ymin=ymin; _zmin=zmin; _xmax=xmax; _ymax=ymax; _zmax=zmax;};

  inline void set_resolution( const int size_x, const int size_y, const int size_z ) { _size_x = size_x ;  _size_y = size_y ;  _size_z = size_z ; }
  /**
   * selects wether the algorithm will use the enhanced topologically controlled lookup table || the original MarchingCubes
   * \param originalMC true for the original Marching Cubes
   */
  inline void set_method    ( const bool originalMC = false ) { _originalMC = originalMC ; }
  /**
   * selects to use data from another class
   * \param data is the pointer to the external data, allocated as a size_x*size_y*size_z vector running in x first
   */
  inline void set_ext_data  ( real *data )
  { if( !_ext_data ) delete [] _data ;  _ext_data = data != NULL ;  if( _ext_data ) _data = data ; }
  /**
   * selects to allocate data
   */
  inline void set_int_data  () { _ext_data = false ;  _data = NULL ; }

  // Data access
  /**
   * accesses a specific cube of the grid
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline const real get_data  ( const int i, const int j, const int k ) const { return _data[ i + j*_size_x + k*_size_x*_size_y] ; }
  /**
   * sets a specific cube of the grid
   * \param val new value for the cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
	
	real get_interpolated_data (const real x, const real y, const real z) {
		int i1, i2;
		i1 = to_cube_i (x);
		i2 = i1 + 1;
		if (i1 < 0) {
			i1 = 0;
			i2 = 0;
		}	
		else if (i2 > _size_x -1) {
			i2 = _size_x -1 ;
		}
		
		int j1, j2;
		j1 = to_cube_j (y);
		j2 = j1 + 1;
		if (j1 < 0) {
			j1 = 0;
			j2 = 0;
		}	
		else if (j2 > _size_y -1) {
			j2 = _size_y -1 ;
		}
		
		int k1, k2;
		k1 = to_cube_k (z);
		k2 = k1 + 1;
		if (k1 < 0) {
			k1 = 0;
			k2 = 0;
		}	
		else if (k2 > _size_z -1) {
			k2 = _size_z -1 ;
		}
		
		real V000 = get_data (i1, j1, k1);
		real V001 = get_data (i1, j1, k2);
		real V010 = get_data (i1, j2, k1);
		real V011 = get_data (i1, j2, k2);
		real V100 = get_data (i2, j1, k1);
		real V101 = get_data (i2, j1, k2);
		real V110 = get_data (i2, j2, k1);
		real V111 = get_data (i2, j2, k2);
		real x_space = (_xmax - _xmin) / _size_x;
		real y_space = (_ymax - _ymin) / _size_y;
		real z_space = (_zmax - _zmin) / _size_z;
		
		real xx = (x - to_real_x(i1)) / x_space;
		real yy = (y - to_real_y(j1)) / y_space;
		real zz = (z - to_real_z(k1)) / z_space;
	//	cerr <<j1<<" "<<j2<<" "<< y <<" "<< to_real_y(j1)<<" "<<y_space<<" "<<xx<<" "<<yy<<" "<<zz<<endl;
		return 
		V000 *(1 - xx) *(1 - yy) *(1 - zz) +
		V100 * xx *(1 - yy) *(1 - zz) +
		V010 *(1 - xx) *yy *(1 - zz) +
		V001 *(1 - xx) *(1 - yy)* zz +
		V101 *xx *(1 - yy)* zz +
		V011 *(1 - xx) *yy *zz +
		V110 *xx* yy *(1 - zz) +
		V111 *xx *yy* zz; 
		
		
		
		
	}
	
	
  inline void  set_data  ( const real val, const int i, const int j, const int k ) { _data[ i + j*_size_x + k*_size_x*_size_y] = val ; }

  // Data initialization
  /** inits temporary structures (must set sizes before call) : the grid && the vertex index per cube */
  void init_temps () ;
  /** inits all structures (must set sizes before call) : the temporary structures && the mesh buffers */
  void init_all   () ;
  /** clears temporary structures : the grid && the main */
  void clean_temps() ;
  /** clears all structures : the temporary structures && the mesh buffers */
  void clean_all  () ;

  void clean_mesh () ;
  void init_mesh ();

//-----------------------------------------------------------------------------
// Exportation
public :
  /**
   * PLY exportation of the generated mesh
   * \param fn  name of the PLY file to create
   * \param bin if true, the PLY will be written in binary mode
   */
  void writePLY( const char *fn, bool bin = false ) ;

  /**
   * VRML / Open Inventor exportation of the generated mesh
   * \param fn  name of the IV file to create
   */
  void writeIV ( const char *fn ) ;

  /**
   * ISO exportation of the input grid
   * \param fn  name of the ISO file to create
   */
  void writeISO( const char *fn ) ;


//-----------------------------------------------------------------------------
// Algorithm
public :
  /**
   * Main algorithm : must be called after init_all
   * \param iso isovalue
   */
  void run( real iso = (real)0.0, int *prog = 0, int *max = 0 ) ;

protected :
  /** tesselates one cube */
  void process_cube ()             ;
  /** tests if the components of the tesselation of the cube should be connected by the interior of an ambiguous face */
  bool test_face    ( schar face ) ;
  /** tests if the components of the tesselation of the cube should be connected through the interior of the cube */
  bool test_interior( schar s )    ;


//-----------------------------------------------------------------------------
// Operations
protected :
  /**
   * computes almost all the vertices of the mesh by interpolation along the cubes edges
   * \param iso isovalue
   */
  void compute_intersection_points( real iso ) ;

  /**
   * routine to add a triangle to the mesh
   * \param trig the code for the triangle as a sequence of edges index
   * \param n    the number of triangles to produce
   * \param v12  the index of the interior vertex to use, if necessary
   */
  void add_triangle ( const char* trig, char n, int v12 = -1 ) ;

  /** tests && eventually doubles the vertex buffer capacity for a new vertex insertion */
  void test_vertex_addition() ;
  /** adds a vertex on the current horizontal edge */
  int add_x_vertex() ;
  /** adds a vertex on the current longitudinal edge */
  int add_y_vertex() ;
  /** adds a vertex on the current vertical edge */
  int add_z_vertex() ;
  /** adds a vertex inside the current cube */
  int add_c_vertex() ;

  /**
   * interpolates the horizontal gradient of the implicit function at the lower vertex of the specified cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */


  real get_x_grad( const int i, const int j, const int k ) const ;
  /**
   * interpolates the longitudinal gradient of the implicit function at the lower vertex of the specified cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  real get_y_grad( const int i, const int j, const int k ) const ;
  /**
   * interpolates the vertical gradient of the implicit function at the lower vertex of the specified cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  real get_z_grad( const int i, const int j, const int k ) const ;

  /**
   * accesses the pre-computed vertex index on the lower horizontal edge of a specific cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */



  inline int   get_x_vert( const int i, const int j, const int k ) const { return _x_verts[ i + j*_size_x + k*_size_x*_size_y] ; }
  /**
   * accesses the pre-computed vertex index on the lower longitudinal edge of a specific cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline int   get_y_vert( const int i, const int j, const int k ) const { return _y_verts[ i + j*_size_x + k*_size_x*_size_y] ; }
  /**
   * accesses the pre-computed vertex index on the lower vertical edge of a specific cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline int   get_z_vert( const int i, const int j, const int k ) const { return _z_verts[ i + j*_size_x + k*_size_x*_size_y] ; }

  /**
   * sets the pre-computed vertex index on the lower horizontal edge of a specific cube
   * \param val the index of the new vertex
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_x_vert( const int val, const int i, const int j, const int k ) { _x_verts[ i + j*_size_x + k*_size_x*_size_y] = val ; }
  /**
   * sets the pre-computed vertex index on the lower longitudinal edge of a specific cube
   * \param val the index of the new vertex
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_y_vert( const int val, const int i, const int j, const int k ) { _y_verts[ i + j*_size_x + k*_size_x*_size_y] = val ; }
  /**
   * sets the pre-computed vertex index on the lower vertical edge of a specific cube
   * \param val the index of the new vertex
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_z_vert( const int val, const int i, const int j, const int k ) { _z_verts[ i + j*_size_x + k*_size_x*_size_y] = val ; }

  /** prints cube for debug */
  void print_cube() ;

//-----------------------------------------------------------------------------
// Elements
protected :
  bool      _originalMC ;   /**< selects wether the algorithm will use the enhanced topologically controlled lookup table || the original MarchingCubes */
  bool      _ext_data   ;   /**< selects wether to allocate data || use data from another class */

  int       _size_x     ;  /**< width  of the grid */
  int       _size_y     ;  /**< depth  of the grid */
  int       _size_z     ;  /**< height of the grid */
  real     *_data       ;  /**< implicit function values sampled on the grid */

  int      *_x_verts    ;  /**< pre-computed vertex indices on the lower horizontal   edge of each cube */
  int      *_y_verts    ;  /**< pre-computed vertex indices on the lower longitudinal edge of each cube */
  int      *_z_verts    ;  /**< pre-computed vertex indices on the lower vertical     edge of each cube */

  int       _nverts     ;  /**< number of allocated vertices  in the vertex   buffer */
  int       _ntrigs     ;  /**< number of allocated triangles in the triangle buffer */
  int       _Nverts     ;  /**< size of the vertex   buffer */
  int       _Ntrigs     ;  /**< size of the triangle buffer */
  Vertex   *_vertices   ;  /**< vertex   buffer */
  Triangle *_triangles  ;  /**< triangle buffer */

  int       _i          ;  /**< abscisse of the active cube */
  int       _j          ;  /**< height of the active cube */
  int       _k          ;  /**< ordinate of the active cube */

  float     _xmin       ;  /**< min abscisse of the rendered cube */
  float     _ymin       ;  /**< min height of the rendered cube */
  float     _zmin       ;  /**< min ordinate of the rendered cube */
  float     _xmax       ;  /**< min abscisse of the rendered cube */
  float     _ymax       ;  /**< min abscisse of the rendered cube */
  float     _zmax       ;  /**< min abscisse of the rendered cube */

  real      _cube[8]    ;  /**< values of the implicit function on the active cube */
  uchar     _lut_entry  ;  /**< cube sign representation in [0..255] */
  uchar     _case       ;  /**< case of the active cube in [0..15] */
  uchar     _config     ;  /**< configuration of the active cube */
  uchar     _subconfig  ;  /**< subconfiguration of the active cube */
};
//_____________________________________________________________________________


#endif // _MARCHINGCUBES_H_
