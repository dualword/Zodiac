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

#include "haptics.h"

Thread::Thread (QObject *parent) : QThread (parent) {}

void Thread::stop () {
	is_running = false;
}

SurfaceThread::SurfaceThread (QObject *parent, Surface *surf, DDWin *ddw) : Thread (parent) {
 //   mesh = false;
    ddwin = ddw;
    surface = surf;
    a = 1.f;
    res = 1.f;
}




void SurfaceThread::run()    {
    compute_surface_data ();
}




void SurfaceThread::compute_surface_data () {
	bool near_bool = true;
	if ((! surface ->near_to )|| (surface ->near_to == surface ->molecule)) near_bool = false;
	cutoffGrid<Atom*> *grid = surface -> grid;
	MarchingCubes *cube = new MarchingCubes (); 
	
	
	
	int xm, ym, zm;
	
	
	float x, y, z;
	float rad = 1.4f;
	float treshold = rad;
	float xmin, ymin, zmin, xmax, ymax, zmax;
	
	float distance = 4.f;
	float ratio = 1; //grid resolution ratio between cube2 && cube
    surface ->molecule->find_limits ();

	
	if (near_bool) {
		surface ->near_to -> find_limits ();
		if (surface -> molecule -> min_corner.x() > surface -> near_to -> min_corner.x()) xmin = surface -> molecule -> min_corner.x();
		else xmin = surface -> near_to -> min_corner.x();
		
		if (surface -> molecule -> min_corner.y() > surface -> near_to -> min_corner.y()) ymin = surface -> molecule -> min_corner.y();
		else ymin = surface -> near_to -> min_corner.y();
		
		if (surface -> molecule -> min_corner.z() > surface -> near_to -> min_corner.z()) zmin = surface -> molecule -> min_corner.z();
		else zmin = surface -> near_to -> min_corner.z();
		
		
		if (surface -> molecule -> max_corner.x() < surface -> near_to -> max_corner.x()) xmax = surface -> molecule -> max_corner.x();
		else xmax = surface -> near_to -> max_corner.x();
		
		if (surface -> molecule -> max_corner.y() < surface -> near_to -> max_corner.y()) ymax = surface -> molecule -> max_corner.y();
		else ymax = surface -> near_to -> max_corner.y();
		
		if (surface -> molecule -> max_corner.z() < surface -> near_to -> max_corner.z()) zmax = surface -> molecule -> max_corner.z();
		else zmax = surface -> near_to -> max_corner.z();	
	//	cerr << "mol1 "<< surface -> near_to ->min_corner.x () <<"-"<<	surface -> near_to ->max_corner.x ()<<"            mol2 "<< 
		//surface -> molecule ->min_corner.x () <<" - "<<	surface -> molecule ->max_corner.x ()<<"       "<<xmin<<" - "<<xmax<<endl;
		
		//	cerr << "mol1 "<< surface -> near_to ->min_corner.y () <<"-"<<	surface -> near_to ->max_corner.y ()<<"            mol2 "<< 
	//	surface -> molecule ->min_corner.y () <<" - "<<	surface -> molecule ->max_corner.y ()<<"       "<<ymin<<" - "<<ymax<<endl;
		
			//		cerr << "mol1 "<< surface -> near_to ->min_corner.z () <<"-"<<	surface -> near_to ->max_corner.z ()<<"            mol2 "<< 
	//	surface -> molecule ->min_corner.z () <<" - "<<	surface -> molecule ->max_corner.z ()<<"       "<<zmin<<" - "<<zmax<<endl;	
		
		
	}
	else {
		xmin = surface -> molecule -> min_corner.x();
		ymin = surface -> molecule -> min_corner.y();
		zmin = surface -> molecule -> min_corner.z();
		
		xmax = surface -> molecule -> max_corner.x();
		ymax = surface -> molecule -> max_corner.y();
		zmax = surface -> molecule -> max_corner.z();
		
	}
	
	
	
	if ((xmin < xmax) && (ymin < ymax) && (zmin < zmax)) {
		
        xmin -= distance;
        ymin -= distance;
        zmin -= distance;
        xmax += distance;
        ymax += distance;
        zmax += distance;
		
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
        cube->run(treshold) ;
		
		
		
		
		
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
		
		
        for (unsigned int k=0; k<zm2; k++) {
            z =  cube2->to_real_z (k);
            for (unsigned int j=0; j<ym2; j++) {
                y =  cube2->to_real_y (j);
                for (unsigned int i=0; i<xm2; i++) {
                    x =  cube2->to_real_x (i);
                    float value = 1000;

                    vect p (x, y, z);
                    objectList<SurfVertex*>* nbAtoms = vgrid->getNeighborObjects(p);
                    if (cube->get_data (cube->to_cube_i(x),cube->to_cube_j(y), cube->to_cube_k(z))>treshold+0.1) value =0.;
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
        cube2->run(treshold) ;
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
        surface->color_by_atom (a);
        for (unsigned int ii=0; ii<myverts.size (); ii++) {
            delete myverts[ii];
        }
        cube->clean_all() ;
        cube2->clean_all();
	//	cerr << "faces " << surface->faces.size () << endl;
	}
//	else cerr << "molecules do not overlap" << endl;
	
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



MinimiseThread::MinimiseThread (QObject *parent, DDWin *ddw) : Thread (parent) {
    ddwin = ddw;
    connect (this, SIGNAL (ask_redraw (Molecule *)), ddwin, SLOT (redraw (Molecule *)));
    molecule = NULL;

}

void MinimiseThread::run () {
    assert (molecule);
/*
    assert (molecule);
    OBForceField *ff;
    ff = OBForceField::FindForceField( "Ghemical" );
	assert (ff);
    ff -> Setup (*molecule);
 //   ff -> SteepestDescentInitialize ();
 //   while (ff -> SteepestDescentTakeNSteps(1)) {
    ff -> ConjugateGradientsInitialize ();
    while (ff -> ConjugateGradientsTakeNSteps(1)) {


        ff -> UpdateCoordinates (*molecule);
        if (!molecule -> needs_redraw) {
            molecule -> needs_redraw = true;
            ask_redraw (molecule);
        }
    }
	*/
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
        if (!molecule -> needs_redraw) {
            molecule -> needs_redraw = true;
            ask_redraw (molecule);
        }
		
	 //   if (converged || iterations > MAX_ITERATIONS) {
			//            ddwin -> data -> minimize ->deinitialise_minimisation ();
	 //   }
    }
}


////////////////////////////////////////////////////////////////////////////////



DatabaseMinimiseThread::DatabaseMinimiseThread (QObject *parent, Database *dat, DDWin *ddw) : Thread (parent) {
	ddwin = ddw;
	database = dat;
}


void DatabaseMinimiseThread::run () {
	ddwin -> data -> undo_stack -> beginMacro ("Database Minimisation");
	for (unsigned int i=0; i<database -> molecules.size (); i++) {
		Molecule *mol = database -> molecules [i];
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



HapticThread::HapticThread (QObject *parent, DDWin *ddw) : Thread (parent) {
	number_of_threads = 1;
    ddwin = ddw;
    is_running = true;
    connect (this, SIGNAL (ask_redraw (Molecule *)), ddwin, SLOT (redraw (Molecule *)));
    connect (this, SIGNAL (ask_color_by_score (Molecule *)), ddwin, SLOT (recolor_by_score (Molecule *)));

    haptic_dof_mode = 2;
    molecule = NULL;
    automove = true;
    color_by_score = false;
	next_atom_to_check = 0;
	number_of_atoms = 0;
	internal_interaction_counter = 0;
	//need to load internal interactions first
}


void HapticThread::initialise (Minimize *min) {
	atoms.clear ();
	number_of_threads = min -> haptic_number_of_threads;
    haptic_dof_mode = min -> haptic_dof_mode;
    set_molecule (min -> haptic_molecule);
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
	minimise -> internal_ff ->load_internal_interactions (&internal_interactions);
	if (number_of_threads < 1) number_of_threads = 1;
	for (int i = 0; i < number_of_threads; i++) {
		HapticWorkerThread *worker_thread = new HapticWorkerThread (this, i);
		worker_thread -> start ();
		worker_threads.push_back (worker_thread);
	}
	worker_threads[0] ->wait ();
/*
    is_running = true;
    Minimize *minimise = ddwin -> data -> minimize;

	const float maxforce = 200.0f;

	Molecule * min_mol = molecule;
    assert (min_mol);
    while (is_running) {
    #ifdef HAPTICS
	    minimise -> update_molecule_position_with_haptic_pointer (min_mol);
    #endif //HAPTICS

        vect haptic_force (0.f, 0.f, 0.f);

        FOR_ATOMS_OF_MOL (a, min_mol) {
			vect v (0., 0., 0.);
			set_force (&*a, v);
			set_old_score (&*a, 0);
	    }



	    if (haptic_dof_mode ==0) {
		
		  // 6 dofs model


            vect tot_force, tot_torque;
            tot_force.null ();
            tot_torque.null ();
		    minimise -> interaction_ff -> update ();
		    minimise -> interaction_ff -> compute_forces ();




		    for (unsigned int fi=0; fi<minimise -> fragments.size (); fi++) { 

			    for (unsigned int i=0; i<minimise -> fragments[fi].atoms.size (); i++) {
            vect *force = (vect *) minimise -> fragments[fi].atoms[i].atom -> GetData ("force");
			        tot_force = sum (tot_force, *force);
				    vect to = torque (*force, minimise -> fragments[fi].atoms[i].atom->GetVector (), minimise -> fragments[fi].translation);
                    tot_torque = sum (tot_torque, to);
			    }

			    haptic_force = sum (haptic_force, tot_force);
			    if (automove) {
                    vect force = tot_force;
                    force.trunc_at (500.f);
                    force.multiply (STEP_SIZE);
				    minimise -> fragments[fi].translation = sum (minimise -> fragments[fi].translation, force); 
			    }
		        float cut = 10.f;
                tot_torque.trunc_at (cut);

			    float new_quat [4], mult_quat [4];
			    axis_angle_to_quaternion (tot_torque, tot_torque.module ()*0.001, new_quat);

			    multiply_quaternions (new_quat, minimise -> fragments[fi].rotation_quat, mult_quat);
			    normalize_quaternion (mult_quat);
			    minimise -> fragments[fi].rotation_quat[0] = mult_quat [0];
			    minimise -> fragments[fi].rotation_quat[1] = mult_quat [1];
			    minimise -> fragments[fi].rotation_quat[2] = mult_quat [2];
			    minimise -> fragments[fi].rotation_quat[3] = mult_quat [3];
			    minimise -> update_fragment_position (minimise -> fragments[fi]);
		    }

	    }



	    else if (haptic_dof_mode ==1) { // 6+r dofs model
	    }
	    else { //3N model
		    vect lastCenter;
			min_mol -> find_center ();
            lastCenter = min_mol -> center;


		    minimise -> internal_ff -> compute_forces ();
		    minimise -> interaction_ff -> update ();
		    minimise -> interaction_ff -> compute_forces ();
            minimise -> apply_forces (min_mol, maxforce);

            FOR_ATOMS_OF_MOL (a, min_mol) {
                vect force = get_force (&*a);

        		haptic_force = sum (haptic_force, force);
		    }




		    if (!automove) {
                vector<Atom*> alist;
				FOR_ATOMS_OF_MOL (a, min_mol) {
					alist.push_back (&*a);
				}
				vect cent = find_mass_center (alist);
                FOR_ATOMS_OF_MOL (a, min_mol) {
					vect vec = (vect &) a -> GetVector ();
	        		vec -= cent;
				    vec += lastCenter;
					a -> SetVector (vec);

			    }
		    }
        //	total_interaction_E = interaction_ff->total_energy;
        //	total_internal_E = internal_ff->total_energy;
	    //    total_E = total_interaction_E + total_internal_E;

    //	    data->ddwin->haptic_menu->update ();

            #ifdef HAPTICS


	        // HAPTIC FUNCTION needs some force cut off and scaling 
	        vect actual_force;
	        ddwin->gl->world_to_haptic_coordinates (haptic_force, actual_force, -5, 5, -5, 5, -5, 5);
	        gHaptics.setForce(
		    haptic_force.x() / maxforce, /////////////why not actual force???
		    haptic_force.y() / maxforce,
		    haptic_force.z() / maxforce
		    );
            #endif //HAPTICS

        }

        FOR_ATOMS_OF_MOL (a, molecule) {
			set_score (&*a, get_old_score (&*a));
	    }



        if (color_by_score) {
        //  vector <color_mask> masks;
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
        if (!(molecule -> needs_redraw)) {
            molecule -> needs_redraw = true;
            ask_redraw (molecule);
        }
    }
*/
}




HapticWorkerThread::HapticWorkerThread (HapticThread *ht, int i) : Thread (ht) {
	master = ht;
	number = i;
}

void HapticWorkerThread::run () {
	while (is_running) {
		if (master -> is_end_of_cycle ()) {
		//	cout << "Thread "<<number<<" end of cycle" << endl;
			master ->end_of_cycle ();	

		}
		else if (master ->have_to_work_on_internal_interactions ()) {
				//	cout << "Thread "<<number<<" one internal interaction" << endl;
			master ->compute_one_internal_interaction (); 

		}
		else if (master ->have_to_work_on_nonbonded_interactions ()) {
	//		cout << "Thread "<<number<<" one nonbonded interaction" << endl;
			master ->compute_one_nonbonded_interaction (); 

		}
		else {
		//	cout << "Thread "<<number<<" new atom" << endl;
			master ->load_interactions_for_new_atom ();

		}
	}
}


bool HapticThread::is_end_of_cycle () {
	if (have_to_work_on_nonbonded_interactions () or !are_atoms_finished ()) return false;
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
	end_of_calculation_for_atom (atoms [number_of_atoms-1]);
	next_atom_to_check = 0;
	internal_interaction_counter = 0;
	if (!automove) {
		vector<Atom*> alist;
		molecule -> find_center ();
		vect cent = molecule -> center;
		FOR_ATOMS_OF_MOL (a, molecule) {
			vect vec = get_coordinates (&*a);
			vec -= cent;
			vec += last_center;
			set_coordinates (&*a, vec);
		}
		last_center = cent;
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
	
	
	if (!(molecule -> needs_redraw)) {
		molecule -> needs_redraw = true;
		ask_redraw (molecule);
	}
 
}

void HapticThread::end_of_calculation_for_atom (Atom *at) {
	const float maxforce = 200.0f;
	flush_forces (at);
	flush_scores (at);
	minimise ->apply_force_to_atom (at, maxforce); 
}

void HapticThread::compute_one_internal_interaction () {
	ForceFieldInteraction *interaction = 0;
	internal_interactions_mutex.lock ();
	if (internal_interaction_counter < internal_interactions.size ()) {
		interaction = internal_interactions[internal_interaction_counter];
		internal_interaction_counter ++;
	}
	internal_interactions_mutex.unlock ();
	if (interaction) interaction ->set_forces ();
}

void HapticThread::compute_one_nonbonded_interaction () {
	nonbonded_interactions_mutex.lock ();
	ForceFieldInteraction *interaction = nonbonded_interactions.front ();
	nonbonded_interactions.pop ();
	interaction ->set_forces ();
	delete interaction;
	nonbonded_interactions_mutex.unlock ();
}


void HapticThread::load_interactions_for_new_atom () {
	atoms_mutex.lock ();
	int n = next_atom_to_check;
	if (n < number_of_atoms) {
		Atom *atom = atoms [n];
		next_atom_to_check++;
		atoms_mutex.unlock ();
		minimise ->interaction_ff ->load_nonbonded_interactions_for_atom (atom, &nonbonded_interactions); 
		if (n) end_of_calculation_for_atom (atoms [n-1]);
	}
	else atoms_mutex.unlock ();
}




GridThread::GridThread (QObject *parent, Grid *gri, DDWin *ddw) : Thread (parent) {
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

    set_led_state (wiimote, 3);

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

    if (cwiid_close(wiimote)) {
	fprintf(stderr, "Error on wiimote disconnect\n");
    }
}


//////////////////////////////////////////////////////////////////////
WiimoteTrackingThread::WiimoteTrackingThread (QObject *parent, DDWin *ddw) : Thread (parent) {
    ddwin = ddw;
    lastx = lasty =  500;
    lastd = 400;
	connect (this, SIGNAL (move_camera (float, float, float)), ddwin -> gl, SLOT (move_camera (float, float, float)));
    is_moving = false;

}


void WiimoteTrackingThread::run () {
	float scale = 0.06f;
	float scale2 = 0.01f;

    wiimote = init_wiimote (wiimote);

    set_led_state (wiimote, 2);

    while (true) {
        vector <int> x, y;
        get_IR_data (wiimote, x, y);
	bool a = is_A_pressed (wiimote);

//        if (x.size () == 1) {
 //           int currx = x[0];
//	    int curry = y[0];
//	    if (!ddwin -> gl -> needs_GL_update) {
//	      ddwin -> gl -> needs_GL_update = true;
//	      head_moved (currx-500, curry-500);
//	    }
//	}
        /*else*/ if (x.size () == 2) {
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

	msleep (10);
    }

    if (cwiid_close(wiimote)) {
	fprintf(stderr, "Error on wiimote disconnect\n");
    }
}

