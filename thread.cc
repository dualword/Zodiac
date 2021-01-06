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

Thread::Thread (QObject *parent) : QThread (parent) {}

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
        cutoffGrid<Atom*> *grid = surface -> grid;
        MarchingCubes *cube = new MarchingCubes (); 



        int xm, ym, zm;


        float x, y, z;
        float rad = 1.4f;
        float treshold = rad;
        float xmin, ymin, zmin, xmax, ymax, zmax;

        float distance = 4.f;
        float ratio = 1; //grid resolution ratio between cube2 && cube
    

        xmin = surface -> molecule -> min_corner.x() - distance;
        ymin = surface -> molecule ->min_corner.y() - distance;
        zmin = surface -> molecule->min_corner.z() - distance;
        xmax = surface -> molecule->max_corner.x() + distance;
        ymax = surface -> molecule->max_corner.y() + distance;
        zmax = surface -> molecule->max_corner.z() + distance;

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
       //             cout <<x<<"  "<<y<<"  "<<z<<endl;
      //              cout<<i<<"  "<<j<<"  "<<k<<endl;
                    vect p (x, y, z);
                    objectList<Atom*>* nbAtoms = grid->getNeighborObjects(p);
                    if (nbAtoms) {
                        vector <Atom *> neighbours = nbAtoms->objects;
                        for (unsigned int a =0; a<neighbours.size (); a++) {
                            double vdw = etab.GetVdwRad (neighbours[a] -> GetAtomicNum ());
                            float dis = subtract ((vect&) neighbours[a] -> GetVector (), p).module () - vdw;
                            if ( dis < value) {
                                value = dis;
                            }
                        }
                    }
       //             cout << value<<endl;
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
        cout <<xm2<<" "<<ym2<<" "<<zm2<<endl;

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
       //             cout <<x<<"  "<<y<<"  "<<z<<endl;
      //              cout<<i<<"  "<<j<<"  "<<k<<endl;
                    vect p (x, y, z);
                    objectList<SurfVertex*>* nbAtoms = vgrid->getNeighborObjects(p);
                    if (cube->get_data (cube->to_cube_i(x),cube->to_cube_j(y), cube->to_cube_k(z))>treshold+0.1) value =0.;
                    else {
                        if (nbAtoms) {
                            vector <SurfVertex *> neighbours = nbAtoms->objects;
                            for (unsigned int a =0; a<neighbours.size (); a++) {
                                float dis = subtract (neighbours[a] -> GetVector (), p).module ();
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
    //    cout << cube->nverts ()<<endl;
      
      //  glNewList(mol->surface_list,GL_COMPILE);
     //   draw_list (cube2); 
     //   glEndList ();


        for (unsigned int i=0; i<cube2->nverts (); i++) {
            SurfVertex *vert = new SurfVertex;
            vert->n = i;
            vert->normal.x() = cube2->vert (i)->nx;        
            vert->normal.y() = cube2->vert (i)->ny;
            vert->normal.z() = cube2->vert (i)->nz;
            vert->coordinates.x() = cube2->to_real_x (cube2->vert (i)->x);       
            vert->coordinates.y() = cube2->to_real_y (cube2->vert (i)->y);       
            vert->coordinates.z() = cube2->to_real_z (cube2->vert (i)->z);
            surface->vertices.push_back (vert);
        }

     //   cout <<surface->vertices.size ();

        for (unsigned int i=0; i<cube2->ntrigs (); i++) {
            SurfFace *face = new SurfFace;
       //     cout << cube2->trig(i)->v1<<endl;

            face->v1 = surface->vertices [cube2->trig(i)->v1];
            face->v2 = surface->vertices [cube2->trig(i)->v2];
            face->v3 = surface->vertices [cube2->trig(i)->v3];

            surface->faces.push_back (face);

        }
        surface->color_by_atom (a);
        for (unsigned int ii=0; ii<myverts.size (); ii++) {
            delete myverts[ii];
        }
        cube->clean_all() ;
        cube2->clean_all();


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
        FOR_ATOMS_OF_MOL (a, molecule) {
            vect nul = (0., 0.,0.);
			set_force (&*a, nul);
	    }
	    iterations ++;
	    ff.update ();
	    ff.compute_forces ();
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
			command -> add (&*a, (vect &) a -> GetVector ());
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
    ddwin = ddw;
    is_running = true;
    connect (this, SIGNAL (ask_redraw (Molecule *)), ddwin, SLOT (redraw (Molecule *)));
    connect (this, SIGNAL (ask_color_by_score (Molecule *)), ddwin, SLOT (recolor_by_score (Molecule *)));

    haptic_dof_mode = 2;
    molecule = NULL;
    automove = true;
    color_by_score = false;
}


void HapticThread::initialise (Minimize *min) {
    haptic_dof_mode = min -> haptic_dof_mode;
  //  automove = min -> automove;
 //   color_by_score = min -> color_by_score;
    set_molecule (min -> haptic_molecule);
}
void HapticThread::stop () {
    is_running = false;
}

void HapticThread::run () {


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
        //    cout << "color_bu"<<endl;
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

}




GridThread::GridThread (QObject *parent, Grid *gri, DDWin *ddw) : Thread (parent) {
    ddwin = ddw;
    grid = gri;

}


void GridThread::run () {
    grid -> load ();
}


