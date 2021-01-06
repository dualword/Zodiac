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

#include "thread.h"
#ifdef HAPTICS
#include <HD/hd.h>
#endif// HAPTICS


float anti_gravity = 2.0;


Thread::Thread (QObject *parent) : QThread (parent), is_running (true), _name ("Thread"), _doing ("running"), _min_step (0), _max_step (0), _current_step (0){}

void Thread::stop () {
	is_running = false;
}


MapThread::MapThread (QObject *parent, Map *mp, DDWin *ddw): Thread (parent) {
    ddwin = ddw;
    map = mp;
	_name = "Computing map values";
}


void MapThread::run () {
    compute_map_data ();
}


void MapThread::compute_map_data () {
	_doing = "computing map";
	MarchingCubes *cube;
	if (!map ->cube) {
		map ->cube = new MarchingCubes (); 
		cube = map ->cube;
		float res = map ->resolution;
		
		
		
		Chemscore *chemscore = new Chemscore;
		ZNMolecule *mol = new ZNMolecule;
		chemscore ->load_environment (ddwin ->molecules, mol);
		
		
		
		
		
		int xm, ym, zm;
		
		
		float x, y, z;
		
		float xmin, ymin, zmin, xmax, ymax, zmax;
		
		
		//  find_limits (map ->molecule);
		
		
		vect c = map ->site_center;
		float r = map ->site_radius;
		
		xmin = c.x() -r ;
		ymin = c.y() -r ;
		zmin = c.z() -r ;
		
		xmax = c.x() +r ;
		ymax = c.y() +r ;
		zmax = c.z() +r ;
		
		
		
		
		
		if ((xmin < xmax) && (ymin < ymax) && (zmin < zmax)) {		
			/*     xmin -= distance;
			 ymin -= distance;
			 zmin -= distance;
			 xmax += distance;
			 ymax += distance;
			 zmax += distance;*/
			
			xm = (int)((xmax-xmin)* res);
			ym = (int)((ymax-ymin)* res);
			zm = (int)((zmax-zmin)* res);
			//   cout <<xm<<" "<<ym<<" "<<zm<<endl;
			
			
			cube->set_resolution( xm, ym, zm) ;
			cube->set_method (false); //use original MC algo?
			cube->set_limits (xmin, ymin, zmin, xmax, ymax, zmax);
			cube->init_all() ;
			_max_step = zm -1;
			for (unsigned int k=0; k<zm; k++) {
				_current_step = k;
				z =  cube->to_real_z (k);
				for (unsigned int j=0; j<ym; j++) {
					y =  cube->to_real_y (j);
					for (unsigned int i=0; i<xm; i++) {
						
						x = cube->to_real_x (i);
						vect p (x, y, z);
						float value;
						switch (map ->type) {
							case 0:
								value = chemscore -> Clvalue (p) + chemscore ->Livalue (p);
								break;
							case 1:
								value = chemscore -> Clvalue (p) + chemscore ->Acceptorvalue (p);
								break;
							case 2:
								value = chemscore -> Clvalue (p) + chemscore ->Donorvalue (p);
								break;
							case 3:
								value =  chemscore ->Electrostatic_potential (p);
								break;
							default:
								value = chemscore -> Clvalue (p) + chemscore ->Livalue (p);
								break;
						}
						//	cerr << value << endl;
						cube->set_data(value, i, j, k );
					}    
				}
			}
		}
	}
	else {	cube = map ->cube;}
	cube->run(map ->threshold, &_current_step, &_max_step) ;
	
	
	_doing = "finalising";
	_max_step = 0;
	int count = 0;
	map -> clean ();
	for (unsigned int i=0; i<cube->ntrigs (); i++) {
		int n1, n2, n3;
		n1 = cube->trig(i)->v1;
		n2 = cube->trig(i)->v2;
		n3 = cube->trig(i)->v3;
		
		
		bool accept = true;
		float x1 = cube->to_real_x (cube->vert(n1)->x);
		float y1 = cube->to_real_y (cube->vert(n1)->y); 
		float z1 = cube->to_real_z (cube->vert(n1)->z);
		
		float x2 = cube->to_real_x (cube->vert(n2)->x);
		float y2 = cube->to_real_y (cube->vert(n2)->y); 
		float z2 = cube->to_real_z (cube->vert(n2)->z);
		
		float x3 = cube->to_real_x (cube->vert(n3)->x);
		float y3 = cube->to_real_y (cube->vert(n3)->y); 
		float z3 = cube->to_real_z (cube->vert(n3)->z);		
		
		
		
		
		if (accept) {
			SurfFace *face = new SurfFace;
			
			SurfVertex *vert1 = new SurfVertex;
			vert1->n = count*3;
			vert1->normal.x() = cube->vert (n1)->nx;        
			vert1->normal.y() = cube->vert (n1)->ny;
			vert1->normal.z() = cube->vert (n1)->nz;
			vert1->coordinates.x() = x1;       
			vert1->coordinates.y() = y1;       
			vert1->coordinates.z() = z1;
			vert1 ->col = map ->solid_color;
			map->vertices.push_back (vert1);		
			
			SurfVertex *vert2 = new SurfVertex;
			vert2->n = count*3+1;
			vert2->normal.x() = cube->vert (n2)->nx;        
			vert2->normal.y() = cube->vert (n2)->ny;
			vert2->normal.z() = cube->vert (n2)->nz;
			vert2->coordinates.x() = x2;       
			vert2->coordinates.y() = y2;       
			vert2->coordinates.z() = z2;
			vert2 ->col = map ->solid_color;
			map->vertices.push_back (vert2);	
			
			SurfVertex *vert3 = new SurfVertex;
			vert3->n = count*3+2;
			vert3->normal.x() = cube->vert (n3)->nx;        
			vert3->normal.y() = cube->vert (n3)->ny;
			vert3->normal.z() = cube->vert (n3)->nz;
			vert3->coordinates.x() = x3;       
			vert3->coordinates.y() = y3;       
			vert3->coordinates.z() = z3;
			vert3 ->col = map ->solid_color;
			map->vertices.push_back (vert3);	
			
			
			face->v1 = vert1;
			face->v2 = vert2;
			face->v3 = vert3;
			
			map->faces.push_back (face);
			count ++;
		}
	}
	//   surface->color_by_atom (alph);
	
	
	
	cube ->clean_mesh ();
	cube ->init_mesh ();
	
	
}



SurfaceThread::SurfaceThread (QObject *parent, Surface *surf, DDWin *ddw): Thread (parent)  {
	_name = "Computing molecular surface";
    ddwin = ddw;
    surface = surf;
    alph = 1.f;
    res = 1.f;
}




void SurfaceThread::run()    {
    compute_surface_data ();	

}




void SurfaceThread::compute_surface_data () {
	_doing = "computing surface";
	bool near_bool = true;
	if ((! surface ->near_to )|| (surface ->near_to == surface ->molecule)) near_bool = false;
	cutoffGrid<Atom*> *grid = surface -> grid;
	MarchingCubes *cube = new MarchingCubes (); 
	
	
	
	int xm, ym, zm;
	
	
	float x, y, z;
	float rad = 1.4f;
	float threshold = rad;
	float xmin, ymin, zmin, xmax, ymax, zmax;
	
	float distance = 4.f;
	if (surface ->near_to_dist > distance) distance = surface ->near_to_dist;
	float ratio = 1; //grid resolution ratio between cube2 && cube
    find_limits (surface ->molecule);
	
	
	if (near_bool) {
		find_limits (surface ->near_to);
		vect min_c = get_min_corner(surface -> molecule);
		vect min_c1 = get_min_corner(surface -> near_to);		
		if (min_c.x() > min_c1.x()) xmin = min_c.x();
		else xmin = min_c1.x();
		
		if (min_c.y() > min_c1.y()) ymin = min_c.y();
		else ymin = min_c1.y();
		
		if (min_c.z() > min_c1.z()) zmin = min_c.z();
		else zmin = min_c1.z();
		
		
		vect max_c = get_max_corner(surface -> molecule);
		vect max_c1 = get_max_corner(surface -> near_to);	
		
		if (max_c.x() < max_c1.x()) xmax = max_c.x();
		else xmax = max_c1.x();
		
		if (max_c.y() < max_c1.y()) ymax = max_c.y();
		else ymax = max_c1.y();
		
		if (max_c.z() < max_c1.z()) zmax = max_c.z();
		else zmax = max_c1.z();

		
		
	}
	else {
		vect min_c = get_min_corner(surface ->molecule);
		vect max_c = get_max_corner(surface ->molecule);
		xmin = min_c.x();
		ymin = min_c.y();
		zmin = min_c.z();
		
		xmax = max_c.x();
		ymax = max_c.y();
		zmax = max_c.z();
		
	}
	
	xmin -= distance;
	ymin -= distance;
	zmin -= distance;
	xmax += distance;
	ymax += distance;
	zmax += distance;
	
	if ((xmin < xmax) && (ymin < ymax) && (zmin < zmax)) {
		

		
        xm = (int)((xmax-xmin)* res/ratio);
        ym = (int)((ymax-ymin)* res/ratio);
        zm = (int)((zmax-zmin)* res/ratio);
		//   cout <<xm<<" "<<ym<<" "<<zm<<endl;
		
		
        cube->set_resolution( xm, ym, zm) ;
        cube->set_method (false); //use original MC algo?
        cube->set_limits (xmin, ymin, zmin, xmax, ymax, zmax);
        cube->init_all() ;
        for (unsigned int k=0; k<zm; k++) {
            z =  cube->to_real_z (k);
            for (unsigned int j=0; j<ym; j++) {
                y =  cube->to_real_y (j);
                for (unsigned int i=0; i<xm; i++) {
                    float value = 1000;
                    x = cube->to_real_x (i);
                    vect p (x, y, z);
                    objectList<Atom*>* nbAtoms = grid->getNeighborObjects(p);
                    if (nbAtoms) {
                        vector <Atom *> neighbours = nbAtoms->objects;
                        for (unsigned int a =0; a<neighbours.size (); a++) {
                            double vdw = etab.GetVdwRad (neighbours[a] -> GetAtomicNum ());
                            float dis = subtract (get_coordinates (neighbours[a]), p).module () - vdw;
                            if ( dis < value) {
                                value = dis;
                            }
                        }
                    }
                    cube->set_data(value, i, j, k );
                }    
            }
        }

        cube->run(threshold, &_current_step, &_max_step) ;
		
		

		_max_step = 0;
		
        MarchingCubes *cube2 = new MarchingCubes (); 
        vector<SurfVertex*> myverts;
        for (unsigned int i=0; i<cube->nverts (); i++) {
            SurfVertex *vert = new SurfVertex;
            vert->coordinates.x() = cube->to_real_x (cube->vert (i)->x);
            vert->coordinates.y() = cube->to_real_y (cube->vert (i)->y);
            vert->coordinates.z() = cube->to_real_z (cube->vert (i)->z);;
            myverts.push_back (vert);
        }
		//   cube->clean_all() ;
        cutoffGrid<SurfVertex*>* vgrid = new cutoffGrid<SurfVertex*>(myverts, 2);
		
		
        int xm2= (int)((xmax-xmin)* res);
        int ym2 = (int)((ymax-ymin)* res);
        int zm2 = (int)((zmax-zmin)* res);
		
        cube2->set_resolution( xm2, ym2, zm2 ) ;
        cube2->set_method (false); //use original MC algo?
        cube2->set_limits (xmin, ymin, zmin, xmax, ymax, zmax);
        cube2->init_all() ;
		_max_step = zm2 -1;
        for (unsigned int k=0; k<zm2; k++) {
			_current_step = k;
            z =  cube2->to_real_z (k);
            for (unsigned int j=0; j<ym2; j++) {
                y =  cube2->to_real_y (j);
                for (unsigned int i=0; i<xm2; i++) {
                    x =  cube2->to_real_x (i);
                    float value = 1000;
					
                    vect p (x, y, z);
                    objectList<SurfVertex*>* nbAtoms = vgrid->getNeighborObjects(p);
                    if (cube->get_data (cube->to_cube_i(x),cube->to_cube_j(y), cube->to_cube_k(z))>threshold+0.1) value =0.;
                    else {
                        if (nbAtoms) {
                            vector <SurfVertex *> neighbours = nbAtoms->objects;
                            for (unsigned int a =0; a<neighbours.size (); a++) {
                                float dis = subtract (get_coordinates (neighbours[a]), p).module ();
                                if ( dis < value) {
                                    value = dis;
                                }
                            }
                        }
                    }
					//          if (value<1.4) cout << value<<endl;
                    cube2->set_data(value, i, j, k );
                }    
            }
        }

        cube2->run(threshold,  &_current_step, &_max_step) ;
		_doing = "finalising";
		int count = 0;
        for (unsigned int i=0; i<cube2->ntrigs (); i++) {
			int n1, n2, n3;
			n1 = cube2->trig(i)->v1;
			n2 = cube2->trig(i)->v2;
			n3 = cube2->trig(i)->v3;
			
			float dist_f = surface -> near_to_dist;
			bool accept = false;
			float x1 = cube2->to_real_x (cube2->vert(n1)->x);
			float y1 = cube2->to_real_y (cube2->vert(n1)->y); 
			float z1 = cube2->to_real_z (cube2->vert(n1)->z);
			
			float x2 = cube2->to_real_x (cube2->vert(n2)->x);
			float y2 = cube2->to_real_y (cube2->vert(n2)->y); 
			float z2 = cube2->to_real_z (cube2->vert(n2)->z);
			
			float x3 = cube2->to_real_x (cube2->vert(n3)->x);
			float y3 = cube2->to_real_y (cube2->vert(n3)->y); 
			float z3 = cube2->to_real_z (cube2->vert(n3)->z);		
			
			if (!near_bool) {accept = true;}
			else {
				
				
				float xm = (x1 +x2+x3) / 3;
				float ym = (y1 +y2+y3) / 3;
				float zm = (z1 +z2+z3) / 3;
				vect p (xm, ym, zm);
				
				objectList<Atom*>* nbAtoms = surface -> near_to_grid->getNeighborObjects(p); 
				if (nbAtoms) {
					vector <Atom *> neighbours = nbAtoms->objects;
					for (unsigned int a =0; a<neighbours.size (); a++) {
						vect v1 = get_coordinates (neighbours[a]);
						double dis_c = dist (v1, p);
						//			cerr << dis_c <<" "<<v1<<" "<<p<< endl;
						if ( dis_c < dist_f) {
							accept = true;
							break;
						}
					}
				}
				else cerr << "no neighbours" << endl;
			}
			if (accept) {
				SurfFace *face = new SurfFace;
				
				SurfVertex *vert1 = new SurfVertex;
				vert1->n = count*3;
				vert1->normal.x() = cube2->vert (n1)->nx;        
				vert1->normal.y() = cube2->vert (n1)->ny;
				vert1->normal.z() = cube2->vert (n1)->nz;
				vert1->coordinates.x() = x1;       
				vert1->coordinates.y() = y1;       
				vert1->coordinates.z() = z1;
				surface->vertices.push_back (vert1);		
				
				SurfVertex *vert2 = new SurfVertex;
				vert2->n = count*3+1;
				vert2->normal.x() = cube2->vert (n2)->nx;        
				vert2->normal.y() = cube2->vert (n2)->ny;
				vert2->normal.z() = cube2->vert (n2)->nz;
				vert2->coordinates.x() = x2;       
				vert2->coordinates.y() = y2;       
				vert2->coordinates.z() = z2;
				surface->vertices.push_back (vert2);	
				
				SurfVertex *vert3 = new SurfVertex;
				vert3->n = count*3+2;
				vert3->normal.x() = cube2->vert (n3)->nx;        
				vert3->normal.y() = cube2->vert (n3)->ny;
				vert3->normal.z() = cube2->vert (n3)->nz;
				vert3->coordinates.x() = x3;       
				vert3->coordinates.y() = y3;       
				vert3->coordinates.z() = z3;
				surface->vertices.push_back (vert3);	
				
				//     cout << cube2->trig(i)->v1<<endl;
				
				face->v1 = vert1;
				face->v2 = vert2;
				face->v3 = vert3;
				
				surface->faces.push_back (face);
				count ++;
			}
        }
		color c (1., 1., 1., alph);
        surface->color_by_color (c);
        for (unsigned int ii=0; ii<myverts.size (); ii++) {
            delete myverts[ii];
        }
        cube->clean_all() ;
        cube2->clean_all();
		//	cerr << "faces " << surface->faces.size () << endl;
	}
		else cerr << "molecules do not overlap" << endl;
	
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SystematicConformationThread::SystematicConformationThread (QObject *parent, ZNMolecule *mol, Database *dat, DDWin *ddw) : Thread (parent), ddwin (ddw), molecule (mol), database (dat) {
	_name = "systematic conformational search";
}

void SystematicConformationThread::run () {
	int nstep = 10;
	vector <double> dihedrals = get_dihedrals (molecule);
	if (dihedrals.size ()) generate_conformation (dihedrals, 0);
}



void SystematicConformationThread::generate_conformation (vector <double> &dofs, int level) {
	int nstep = 10;
	for (unsigned int i = 0; i < nstep; i ++) {
		dofs[level] = M_PI/nstep*2*i;
		if (level < dofs.size ()-1) {
			generate_conformation (dofs, level+1);
		}
		else {
			set_dihedrals (molecule, dofs);
			Database_molecule *new_mol = new Database_molecule (*molecule);
			database ->safe_add_mol (new_mol);
		}
	}
	
}



StochasticConformationThread::StochasticConformationThread (QObject *parent, ZNMolecule *mol, Database *dat, DDWin *ddw) : Thread (parent), ddwin (ddw), database (dat), molecule (mol),  _results_number (10), _time_limit (60) {
	_name = "stochastic conformational search";
	_doing = "exploring conformational space";
}

void StochasticConformationThread::run () {
	build_kinematic_chain(molecule);
	Minimize *min = new Minimize (ddwin ->data);
	min ->set_molecule (molecule);
	min ->internal_ff ->initialize_internal (molecule, ddwin ->molecules);

	min ->interaction_ffs[0] ->initialize_interaction (molecule, ddwin ->molecules);

	ScoreMolecule *function = new ScoreMolecule (min, false);

	PSO *optimiser = new PSO (function);
	//		ILS *optimiser = new ILS (function);
	optimiser ->set_time_limit (_time_limit);

	optimiser ->run ();

	vector <conformation *> confs = optimiser ->get_results ();
	int limit = _results_number;
	if (confs.size () < _results_number) limit = confs.size ();
	for (unsigned int i = 0; i < limit; i++) {
		build_molecule_from_dofs (molecule, confs[i] ->state);
		Database_molecule *new_mol = new Database_molecule (*molecule);
		database ->safe_add_mol (new_mol);
	}
	
}

DockingThread::DockingThread (QObject *parent, ZNMolecule *mol, Database *dat, DDWin *ddw) : Thread (parent), ddwin (ddw), database (dat), molecule (mol), _results_number (10), _time_limit (60),
_bindingsite_radius (10), _bindingsite_center (vect (0., 0., 0.)){
	_name = "molecular docking";
}

void DockingThread::run () {
	MMFF *mmff = new MMFF;
	Chemscore *chemscore = new Chemscore;
	build_kinematic_chain(molecule);
	Minimize *min = new Minimize (ddwin ->data, mmff, chemscore);
	min ->set_molecule (molecule);
	min ->internal_ff ->initialize_internal (molecule, ddwin ->molecules);
	min ->interaction_ffs[0] ->initialize_interaction (molecule, ddwin ->molecules, _bindingsite_center, _bindingsite_radius);
	ScoreMolecule *function = new ScoreMolecule (min, true, _bindingsite_center, _bindingsite_radius);
	//		ILS *optimiser = new ILS (function);
	PSO *optimiser = new PSO (function);
	optimiser ->set_time_limit (_time_limit);
	optimiser ->set_iteration_limit(0);
	msleep (10);
	optimiser ->run ();
	vector <conformation *> confs = optimiser ->get_results ();
	int limit = _results_number;
	if (confs.size () < _results_number) limit = confs.size ();
	for (unsigned int i = 0; i < limit; i++) {
		build_molecule_from_dofs (molecule, confs[i] ->state, true);
		Database_molecule *new_mol = new Database_molecule (*molecule);
		ddwin ->set_lists (new_mol);
		database ->safe_add_mol (new_mol);
	}
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



MinimiseThread::MinimiseThread (QObject *parent, DDWin *ddw) : Thread (parent){
	_name = "energy minimisation";
    ddwin = ddw;
	
    molecule = NULL;
	
	
}

void MinimiseThread::run () {
	//  cerr << "minimise" << endl;
    assert (molecule);
	
    OBForceField *ff;
    ff = OBForceField::FindForceField( "MMFF94s" );
	assert (ff);
    ff -> Setup (*molecule);
    ff -> SteepestDescentInitialize ();
    while (ff -> SteepestDescentTakeNSteps(10)) {
		//   ff -> ConjugateGradientsInitialize ();
		//  while (ff -> ConjugateGradientsTakeNSteps(1)) {
		
        ff -> GetCoordinates (*molecule);
		ddwin -> gl -> draw_molecule (molecule);
    }
	/*
	 MMFF ff;
	 ff.clear_internal_interactions ();
	 ff.clear_nonbonded_interactions (); 
	 ff.initialize (molecule, ddwin -> molecules);  
	 
	 bool converged = false;
	 int iterations = 0;
	 float last_E = 0.f;
	 while (!converged && iterations < MAX_ITERATIONS) {
	 iterations ++;
	 ff.update ();
	 ff.compute_forces ();
	 FOR_ATOMS_OF_MOL(a, molecule) {
	 flush_forces (&*a);
	 }
	 ddwin -> data-> minimize -> apply_forces (molecule, 500.f);
	 
	 float this_E = ff.compute_total_energy ();
	 
	 if (last_E-this_E<MIN_ENERGY) converged = true;
	 if (iterations < 5) {
	 converged = false; 
	 }
	 if (last_E-this_E<0) {
	 converged = false;
	 }
	 last_E = this_E;
	 
	 set_needs_redraw (molecule, true);
	 //   if (converged || iterations > MAX_ITERATIONS) {
	 //            ddwin -> data -> minimize ->deinitialise_minimisation ();
	 //   }
	 }
	 */
}


////////////////////////////////////////////////////////////////////////////////



DatabaseMinimiseThread::DatabaseMinimiseThread (QObject *parent, Database *dat, DDWin *ddw) : Thread (parent) {
	_name = "database energy minimisation";
	ddwin = ddw;
	database = dat;
}


void DatabaseMinimiseThread::run () {
	ddwin -> data -> undo_stack -> beginMacro ("Database Minimisation");
	_max_step = database -> count_entries () -1;
	for (unsigned int i=0; i<database -> count_entries (); i++) {
		_current_step = i;
		ZNMolecule *mol = database -> get_molecule (i);
		ddwin -> data -> undo_stack -> beginMacro ("Energy Minimisation");
		MoveAtomsCommand *command = new MoveAtomsCommand (ddwin -> gl, 1);
		FOR_ATOMS_OF_MOL(a, mol) {
			command -> add (&*a, get_coordinates(&*a));
		}
		ddwin -> execute (command);
		
		MinimiseThread *thread = new MinimiseThread (0, ddwin);
		thread -> set_molecule (mol);
		thread -> start ();
		thread -> wait ();
		delete thread;
		ddwin -> end_minimisation ();
	}
	ddwin -> data -> undo_stack -> endMacro ();
}

/////////////////////////////////////////////////////////////////////////////////

HapticForceThread::HapticForceThread (QObject *parent, DDWin *ddw) : Thread (parent) {
	ddwin = ddw;

}

void HapticForceThread::run () {

	while (is_running) {

		haptic_callback ((void *) ddwin ->data);
	msleep (1);

	}
}


/////////////////////////////////////////////////////////////////////////////////



HapticThread::HapticThread (QObject *parent, DDWin *ddw) : Thread (parent){
	last_is_user_selected = false;
	_name = "continuous minimisation";
	number_of_threads = 1;
    ddwin = ddw;
    is_running = true;
	save_result = false;
	saving = false;
	is_worse = false;
	user_save_result = false;
	save_E_threshold = 1000.f;
	mult = 1.f;
	
    connect (this, SIGNAL (ask_color_by_score (ZNMolecule *)), ddwin, SLOT (recolor_by_score (ZNMolecule *)));
	
    haptic_dof_mode = 1;
    molecule = NULL;
	results = NULL;
    automove = true;
    color_by_score = false;
	next_atom_to_check = 0;
	number_of_atoms = 0;
	internal_interaction_counter = 0;
	//need to load internal interactions first
}


void HapticThread::initialise (Minimize *min) {
	ncycles = 0;
	worker_threads.clear ();
	internal_interactions.clear ();
	
	ddwin -> vertex_list.clear ();
	Hbonds.clear ();
	atoms.clear ();
	number_of_threads = min -> haptic_number_of_threads;
    haptic_dof_mode = min -> haptic_dof_mode;
    set_molecule (min -> haptic_molecule);
	results = new Database;
	vector <double> v;
	results -> safe_add_field ("Energy", v);
	FOR_ATOMS_OF_MOL (a, molecule) {
		atoms.push_back (&*a);
	}
	number_of_atoms = atoms.size ();
}


void HapticThread::stop () {
    for (unsigned int i=0; i< worker_threads.size (); i++) {
		worker_threads [i] -> stop ();
	}
	terminate ();
}

void HapticThread::run () {
	
	minimise = ddwin -> data -> minimize;
	if (haptic_dof_mode != 0) {
		minimise -> internal_ff ->load_internal_interactions (&internal_interactions);
	}
	if (number_of_threads < 1) number_of_threads = 1;
	if (number_of_threads == 1) 	{
		HapticWorkerThread *single_thread = new HapticWorkerThread (this, 1, HapticWorkerThread::SINGLE);
		single_thread -> start ();
		worker_threads.push_back (single_thread);
	}
	else {
		HapticWorkerThread *internal_thread = new HapticWorkerThread (this, 1, HapticWorkerThread::INTERNAL);
		internal_thread -> start ();
		worker_threads.push_back (internal_thread);
		for (int i = 1; i < number_of_threads; i++) {
			HapticWorkerThread *worker_thread = new HapticWorkerThread (this, i, HapticWorkerThread::INTERACTION);
			worker_thread -> start ();
			worker_threads.push_back (worker_thread);
		}
	}
	worker_threads[0] ->wait ();

}




HapticWorkerThread::HapticWorkerThread (HapticThread *ht, int i, HapticWorkerType t) : Thread (0)  {
	master = ht;
	number = i;
	thread_type = t;
}

void HapticWorkerThread::run () {
	while (is_running) {
		master ->automove = master ->ddwin -> haptic_menu ->automove_b;
		master ->color_by_score = master ->ddwin ->haptic_menu ->color_by_score_b;
		master ->mult = master ->ddwin ->haptic_menu ->mult;
		if (thread_type == INTERNAL) {
			for (unsigned int i = 0; i < master -> internal_interactions.size (); i++) {
				master -> compute_internal_interaction (i);
			}
			//	master ->compute_one_internal_interaction (); 
		}
		else if (thread_type == INTERACTION) {
			for (unsigned int a = 0; a < master ->atoms.size (); a++) {
				master ->compute_nonbonded_interactions_for_atom (a);
			}
			master -> end_of_cycle ();				
		}
		else if (thread_type == SINGLE) {

//#pragma omp parallel sections


{
			if (master -> haptic_dof_mode != 0) {
//#pragma omp for
				for (unsigned int i = 0; i < master -> internal_interactions.size (); i++) {
					master -> compute_internal_interaction (i);
				}
			}
//#pragma omp section
{
//#pragma omp for
			for (unsigned int a = 0; a < master ->atoms.size (); a++) {
				master ->compute_nonbonded_interactions_for_atom (a);
			}
}
//#pragma omp section
{
			master ->ddwin ->haptic_menu ->restrain_lock ->lockForRead ();
			for (unsigned int r = 0; r < master ->ddwin ->haptic_menu ->restrains.size (); r++) {
				master ->compute_restrain (r);
			}			
			master ->ddwin ->haptic_menu ->restrain_lock ->unlock ();
}
}
			master -> end_of_cycle ();
		}
		
		
		/*
		 if (master -> is_end_of_cycle ()) {
		 //		cout << "Thread "<<number<<" end of cycle" << endl;
		 master ->end_of_cycle ();	
		 
		 }
		 else if (master ->have_to_work_on_internal_interactions ()) {
		 //		cout << "Thread "<<number<<" one internal interaction" << endl;
		 master ->compute_one_internal_interaction (); 
		 
		 }
		 else if (master ->have_to_work_on_nonbonded_interactions ()) {
		 //		cout << "Thread "<<number<<" one nonbonded interaction" << endl;
		 master ->compute_one_nonbonded_interaction (); 
		 
		 }
		 else {
		 //		cout << "Thread "<<number<<" new atom" << endl;
		 master ->load_interactions_for_new_atom ();
		 
		 
		 }
		 */
	}
}


bool HapticThread::is_end_of_cycle () {
	if (have_to_work_on_nonbonded_interactions () || !are_atoms_finished ()) return false;
	else return true;
}

bool HapticThread::have_to_work_on_nonbonded_interactions () {
	if (nonbonded_interactions.size () == 0) return false;
	else return true;
}

bool HapticThread::have_to_work_on_internal_interactions () {
	if (internal_interactions.size () <= internal_interaction_counter) return false;
	else return true;
}

bool HapticThread::are_atoms_finished () {
	//	cerr << next_atom_to_check<<"  "<<number_of_atoms<<endl;
	if (next_atom_to_check >= number_of_atoms) return true;
	else return false;
}


void HapticThread::end_of_cycle () {
	ncycles++;
	next_atom_to_check = 0;
	internal_interaction_counter = 0;
	saving = ddwin ->haptic_menu ->saving;
	if (haptic_dof_mode == 1) {
		if (!automove) {
			vector<Atom*> alist;
			find_center (molecule);
			vect cent = get_center (molecule);
			FOR_ATOMS_OF_MOL (a, molecule) {
				vect vec = get_coordinates (&*a);
				vec = subtract (vec, cent);
				vec = sum (vec, last_center);
				set_coordinates (&*a, vec);
			}
			
			
		}
		else {
			last_center = get_center (molecule);
		}	
	}
	else if (haptic_dof_mode == 0) {
		
		if (automove) {
			vect total_force (0., 0., 0.);
			FOR_ATOMS_OF_MOL (a, molecule) {
				vect force = get_force (&*a);
				total_force = sum (total_force, force);
			}
			total_force.trunc_at(200);
			total_force.multiply(0.001f);
			ddwin -> gl ->translate_molecule(molecule, total_force);
		}
		
	} 
	if (color_by_score) {
		vector <color_mask> masks;
		color_mask mask;
		mask.intensity = 1.0f;
		mask.only_to = 0;
		mask.excluding = 0;
		mask.type =  2; //score //see menu.cc
		masks.push_back (mask);
		if (!(molecule -> needs_recolor)) {
			molecule -> needs_recolor = true;
			ask_color_by_score (molecule);
		}
	}
	ddwin -> Hbonds = Hbonds;
	double energy_threshold = 0.;
	FOR_ATOMS_OF_MOL (a, molecule) {
		double score = get_score (&*a);
		
		total_energy += score;
	}
	ddwin ->data ->total_energy_haptic = total_energy;
	if (save_result && saving) {
		bool save = true;
		if (results ->count_entries()) {
			if (!last_is_user_selected) {
				ZNMolecule *last_mol = results ->get_molecule(results ->count_entries()-1);
				float RMSD = very_fast_RMSD(molecule, last_mol);
				if (RMSD < ddwin ->haptic_menu ->cluster_RMSD) {
					if (is_worse) {
						save = false;
					}
					else results -> delete_entry (results ->count_entries()-1);
				}
			}
		}
		if (save) {
			last_is_user_selected = false;
			Database_molecule *new_mol = new Database_molecule (*molecule);
			ddwin ->set_lists (new_mol);
			results ->add_mol (new_mol);
			results ->set_double (1, results ->count_entries () -1, total_energy);
		}
	}
	if (user_save_result) {
		Database_molecule *new_mol = new Database_molecule (*molecule);
		ddwin ->set_lists (new_mol);
		new_mol ->SetTitle("USER mol");
		results ->add_mol (new_mol);
		results ->set_double (1, results ->count_entries () -1, total_energy);
		last_is_user_selected = true;
	}
	save_result = false;
	is_worse = false;
	user_save_result = false;
	//if (total_energy < energy_threshold) {
	//	if (ddwin ->vertex_list.size ()) {
	//		vect c = get_center(molecule);
	//		if (c != ddwin ->vertex_list[ddwin ->vertex_list.size ()-1])
	//			ddwin -> vertex_list.push_back (c);
	//	}
	//	else 				ddwin -> vertex_list.push_back (get_center (molecule));
	//}
	Hbonds.clear ();
	total_energy = 0.;
	
	if (haptic_dof_mode == 0) {
		quaternion quat = axis_angle_to_quaternion(total_torque, total_torque.module()*0.001);
		rotate_molecule(molecule, quat, get_center (molecule));
		total_torque = vect (0., 0., 0.);
		
	}
	
#ifdef HAPTICS
	vect haptic_force (0., 0., 0.), actual_force (0., 0., 0.);
	FOR_ATOMS_OF_MOL (a, molecule) {
		
		vect force = get_force (&*a);
		haptic_force = sum (haptic_force, force);
	}
	//ddwin->gl->world_to_haptic_coordinates (haptic_force, actual_force, -5, 5, -5, 5, -5, 5);
	haptic_force = ddwin ->gl ->deapply_world_rotation (haptic_force);
	//cerr << haptic_force << haptic_force.module ()<< endl;
	haptic_force.trunc_at (100.);
	/*
	 viscous_force_x.push_back (haptic_force.x());
	 viscous_force_y.push_back (haptic_force.y());
	 viscous_force_z.push_back (haptic_force.z());
	 int viscosity = 6;
	 while (viscous_force_x.size () > viscosity) {
	 viscous_force_x.erase (viscous_force_x.begin ());
	 viscous_force_y.erase (viscous_force_y.begin ());
	 viscous_force_z.erase (viscous_force_z.begin ());
	 }
	 float x = 0.f;
	 float y = 0.f;
	 float z = 0.f;
	 for (unsigned int i = 0; i < viscous_force_x.size (); i++) {
	 x += viscous_force_x[i];
	 y += viscous_force_y[i];
	 z += viscous_force_z[i];
	 }
	 x /= viscous_force_x.size();
	 y /= viscous_force_y.size();
	 z /= viscous_force_z.size();
	 */
	
	//temporary
	
#ifdef IANS_HACK
	ddwin ->data ->current_force_x=x / 20.0f *mult;
	ddwin ->data ->current_force_y=y / 20.0f *mult;
	ddwin ->data ->current_force_z=z / 20.0f *mult;
#else // IANS_HACK
	double dBlendFactor = 0.1;
	ddwin->data->current_force_x = ((1.0-dBlendFactor)*ddwin->data->current_force_x) + (dBlendFactor*(haptic_force.x()/20.0f*mult));
	ddwin->data->current_force_y = ((1.0-dBlendFactor)*ddwin->data->current_force_y) + (dBlendFactor*(haptic_force.y()/20.0f*mult));
	ddwin->data->current_force_z = ((1.0-dBlendFactor)*ddwin->data->current_force_z) + (dBlendFactor*(haptic_force.z()/20.0f*mult));
#endif // IANS_HACK
	
	//	cerr << "current force "<<current_force[0]<<" "<<current_force[1]<<" "<<current_force[2]<<endl;
	
	// hduVector3Dd force (haptic_force.x() / 5.0f, haptic_force.y() / 5.0f, haptic_force.z() / 5.0f );
	//     hdSetDoublev(HD_CURRENT_FORCE, force);
	
	/*	gHaptics.setForce(
	 haptic_force.x() / 5.0f, 
	 haptic_force.y() / 5.0f + anti_gravity,
	 haptic_force.z() / 5.0f
	 
	 
	 =======
	 gHaptics.setForce(
	 haptic_force.x() / 10.0f, 
	 haptic_force.y() / 10.0f + anti_gravity,
	 haptic_force.z() / 10.0f 
	 >>>>>>> 1.23
	 );
	 */ 
	minimise -> update_molecule_position_with_haptic_pointer (molecule);
	
#endif //HAPTICS
	
	QTime curr_time = QTime::currentTime ();
	stringstream ss;
	int msecs = last_time.msecsTo (curr_time);
	if (msecs > 10000) {
		cerr << ncycles / 10 << " fps;  ";
		cerr << ddwin -> gl ->redraw_counter  / 10 << " graphics updates;  ";

		cerr << endl;
		last_time = curr_time;
		ncycles = 0;
		ddwin -> gl ->redraw_counter = 0;

		
	}
	set_needs_redraw(molecule, true);
}

void HapticThread::end_of_calculation_for_atom (Atom *at) {
	const float maxforce = 200.0f;
	flush_forces (at);
	flush_scores (at);
	if (haptic_dof_mode == 0) {
		vect tor;
		tor = torque (get_force(at), get_coordinates(at), get_center (molecule));
		total_torque = sum (tor, total_torque);
	}
	else minimise ->apply_force_to_atom (at, maxforce); 
}

void HapticThread::compute_restrain (int i) {
	ddwin ->haptic_menu ->restrains[i] ->set_forces ();
}

void HapticThread::compute_internal_interaction (int i) {
	internal_interactions[i] ->set_forces ();
}

void HapticThread::compute_nonbonded_interactions_for_atom (int a) {
	queue <ForceFieldInteraction *> interactions;
	Atom *atom = atoms [a];
	for (unsigned int i = 0; i < minimise ->interaction_ffs.size (); i++) {
		minimise ->interaction_ffs[i] ->load_nonbonded_interactions_for_atom (atom, &interactions); 
		while (interactions.size()) {
			ForceFieldInteraction *interaction = interactions.front ();
			interactions.pop ();
			//		if (interaction ->isElectrostatic ()) interaction ->set_forces (true);
			//	else
			

			interaction ->set_forces (true);

			if (interaction -> isHbond ()) {
				vect v1 = get_coordinates (interaction ->at1);
				vect v2 = get_coordinates (interaction ->at2);
				float dis = dist (v1, v2); 
				float ref_dis = 2.2f;
				//			float ref_dis = 1.85f;
				float ref_delta = 0.4f;
				if ( dis < (ref_delta + ref_dis)) { 
					double perc = 0.;
					if (dis < ref_dis) perc = 1.;
					else {
						float dd = dis - ref_dis;
						perc = 1- (dd/0.4);
					}
					HBond hb;
					hb.v1 = v1;
					hb.v2 = v2;
					
					hb.perc = perc;
					Hbonds.push_back (hb);
				}
			}
			
			delete interaction;
		}
	}
	end_of_calculation_for_atom (atom);
}


void HapticThread::compute_one_internal_interaction () {
	ForceFieldInteraction *interaction = 0;
	//	internal_interactions_mutex.lock ();
	if (internal_interaction_counter < internal_interactions.size ()) {
		interaction = internal_interactions[internal_interaction_counter];
		internal_interaction_counter ++;
	}
	//	internal_interactions_mutex.unlock ();
	if (interaction) interaction ->set_forces ();
}

void HapticThread::compute_one_nonbonded_interaction () {
	nonbonded_interactions_mutex.lock ();
	ForceFieldInteraction *interaction = nonbonded_interactions.front ();
	nonbonded_interactions.pop ();
	nonbonded_interactions_mutex.unlock ();
	if (interaction ->isElectrostatic ()) interaction ->set_forces (true);
	else interaction ->set_forces (true);
	if (interaction -> isHbond ()) {
		vect v1 = get_coordinates (interaction ->at1);
		vect v2 = get_coordinates (interaction ->at2);
		float dis = dist (v1, v2); 
		float ref_dis = 2.2f;
		//			float ref_dis = 1.85f;
		float ref_delta = 0.4f;
		if ( dis < (ref_delta + ref_dis)) { 
			double perc = 0.;
			if (dis < ref_dis) perc = 1.;
			else {
				float dd = dis - ref_dis;
				perc = 1- (dd/0.4);
			}
			HBond hb;
			hb.v1 = v1;
			hb.v2 = v2;
			
			hb.perc = perc;
			Hbonds.push_back (hb);
		}
	}
	
	//	delete interaction;
	
	
}


void HapticThread::load_interactions_for_new_atom () {
	atoms_mutex.lock ();
	int n = next_atom_to_check;
	if (n < number_of_atoms) {
		Atom *atom = atoms [n];
		next_atom_to_check++;
		atoms_mutex.unlock ();
		for (unsigned int i = 0; i < minimise ->interaction_ffs.size (); i++) {
			minimise ->interaction_ffs[i] ->load_nonbonded_interactions_for_atom (atom, &nonbonded_interactions); 
		}
		if (n) end_of_calculation_for_atom (atoms [n-1]);
	}
	else atoms_mutex.unlock ();
}




GridThread::GridThread (QObject *parent, Grid *gri, DDWin *ddw) : Thread (parent) {
	_name = "grid generation";
    ddwin = ddw;
    grid = gri;
	
}


void GridThread::run () {
    grid -> load ();
}


/////////////////////////////////////////////////////////////////////////////////

HeadTrackingThread::HeadTrackingThread (QObject *parent, DDWin *ddw) : Thread (parent) {
    ddwin = ddw;
	//  lastx = lasty =  500;
	connect (this, SIGNAL (head_moved (int, int)), ddwin -> gl, SLOT (head_tracking_update (int, int)));
	
}


void HeadTrackingThread::run () {
	
	
    wiimote = init_wiimote (wiimote);
	
	
    while (true) {
        vector <int> x, y;
        get_IR_data (wiimote, x, y);
		
		//        if (x.size () == 1) {
		//           int currx = x[0];
		//	    int curry = y[0];
		//	    if (!ddwin -> gl -> needs_GL_update) {
		//	      ddwin -> gl -> needs_GL_update = true;
		//	      head_moved (currx-500, curry-500);
		//	    }
		//	}
        /*else*/ if (x.size () == 2) {
            int currx = (x[0] + x[1]) * 0.5;
			int curry = (y[0] + y[1]) * 0.5;
			if (!ddwin -> gl -> needs_GL_update) {
				ddwin -> gl -> needs_GL_update = true;
				head_moved (currx-500, curry-500);
			}
		}
		
		msleep (10);
    }
	
    if (is_wiimote_closed (wiimote)) {
		fprintf(stderr, "Error on wiimote disconnect\n");
    }
}


//////////////////////////////////////////////////////////////////////
WiimoteTrackingThread::WiimoteTrackingThread (QObject *parent, DDWin *ddw) : Thread (parent) {
    ddwin = ddw;
    lastx = lasty =  500;
    lastd = 400;
	connect (this, SIGNAL (move_target (float, float, float)), ddwin -> gl, SLOT (move_target (float, float, float)));	
	connect (this, SIGNAL (move_camera (float, float, float)), ddwin -> gl, SLOT (move_camera (float, float, float)));
	connect (this, SIGNAL (rotate_world (double, double, double, double, double, double)), ddwin -> gl, SLOT (map_vector_on_vector_world (double, double, double, double, double, double)));
	connect (this, SIGNAL (rotate_target (double, double, double, double, double, double)), ddwin -> gl, SLOT (map_vector_on_vector_target (double, double, double, double, double, double)));
    is_moving = false;
    is_rotating = false;
	
}


void WiimoteTrackingThread::run () {
	float scale = 0.01f;
	float scale2 = 0.05f;
	wiimote = init_wiimote (wiimote);
	
	
	//switch_to_mode (1);
	
	
	while (true) {
		poll (wiimote);
		
		vector <int> x, y;
		get_IR_data (wiimote, x, y);
		vect current_orientation = get_Acc_data (wiimote);
		
		
		
		
		bool a = is_A_pressed (wiimote);
		bool b = is_B_pressed (wiimote);
		bool one = is_1_pressed (wiimote);
		bool two = is_2_pressed (wiimote);
		
		
		if (one) {
			cerr << "one" << endl;
			switch_to_mode (1);
		}
		else if (two) {
			switch_to_mode (2);
		}
		
		if (mode == 1) {
			if (x.size () == 2) {
				if (a && !is_moving) {
					is_moving = true;
					lastx = (float (x[0] + x[1])) * 0.5;
					lasty = (float (y[0] + y[1])) * 0.5;
					lastd = (x[1] -x[0]) * (x[1] -x[0]) + (y[1] -y[0]) * (y[1] -y[0]); 
					
				}
				else if ((!a) && is_moving) {
					is_moving = false;
				}
				else if (a && is_moving) {
					float currx = (float (x[0] + x[1])) * 0.5;
					float curry = (float (y[0] + y[1])) * 0.5;
					float currd = (x[1] -x[0]) * (x[1] -x[0]) + (y[1] -y[0]) * (y[1] -y[0]); 
					if (!ddwin -> gl -> needs_GL_update) {
						ddwin -> gl -> needs_GL_update = true;
						move_camera (-(currx-lastx)*scale, -(curry-lasty)*scale, -(currd-lastd)*scale*scale2);
						lastx = currx;
						lasty = curry;
						lastd = currd;
					}
				}
				
			}
			if (b) {
				cerr << "b"<<endl;
				double x2 =  current_orientation.x();
				double z2 = -current_orientation.y();
				double y2 =  current_orientation.z();
				double x1 =  last_orientation.x();
				double z1 = -last_orientation.y();
				double y1 =  last_orientation.z();
				
				if (is_rotating) {rotate_world (x1, y1, z1, x2, y2, z2);}
				else is_rotating = true;
				last_orientation = current_orientation;
			}
			else is_rotating = false;
			
		}
		else if (mode == 2) {
			if (b) {
				double x2 =  current_orientation.x();
				double z2 = -current_orientation.y();
				double y2 =  current_orientation.z();
				double x1 =  last_orientation.x();
				double z1 = -last_orientation.y();
				double y1 =  last_orientation.z();
				
				if (is_rotating) {rotate_target (x1, y1, z1, x2, y2, z2);}
				else is_rotating = true;
				last_orientation = current_orientation;
			}
			else if (a) {
				if (x.size () == 2) {
					float currx = (float (x[0] + x[1])) * 0.5;
					float curry = (float (y[0] + y[1])) * 0.5;
					float currd = (x[1] -x[0]) * (x[1] -x[0]) + (y[1] -y[0]) * (y[1] -y[0]); 
					if (is_moving) {move_target (-(currx-lastx)*scale, -(curry-lasty)*scale, -(currd-lastd)*scale*scale2);}
					else is_moving = true;
					lastx = currx;
					lasty = curry;
					lastd = currd;
				}
			}
			else {
				is_rotating = false;
				is_moving   = false;
			}
		}
		
		
		msleep (100);
	}
	
	
    if (is_wiimote_closed(wiimote)) {
		fprintf(stderr, "Error on wiimote disconnect\n");
    }
	
}


void WiimoteTrackingThread::switch_to_mode (int i) {
	switch (i) {
		case 1:
			mode = 1;
			switch_led (wiimote, 1, true);
			switch_led  (wiimote, 2, false);
			break;
		case 2:
			mode = 2;
			switch_led (wiimote, 2, true);
			switch_led  (wiimote, 1, false);
			break;
	}
	
}

