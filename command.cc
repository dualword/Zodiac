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


#include "command.h"



Command::Command () : QUndoCommand (0)
{
    first_time = true;
 //   get_time ();
}

void Command::get_time () {
 //   commandtime = time (NULL);
   
}


bool Command::same_time (Command *other) {

  //  double diff = difftime (commandtime, other -> commandtime);
 //   assert (diff >= 0);
 //   if (diff < 1) return true;
  //  else return false;
    return false;
}

/*
bool Command::mergeWith (const QUndoCommand *other) {
    cout <<" merge"<<endl;
    return true;
}
*/

void Command::name (string str) {
    setText (QString (str.c_str ()));
}

ColorAtomCommand::ColorAtomCommand (MyGl *gl_point) : Command ()
{
    redraw = false;
    gl = gl_point;
    setText ( QString (string("color:").c_str ()));
   
}





ColorAtomCommand::ColorAtomCommand (Atom *at, color pres_col, MyGl *gl_point) : Command ()
{
    redraw = false;
    gl = gl_point;
    add (at, pres_col);
    setText ( QString (string("color: 1  atom").c_str ()));
 
}

void ColorAtomCommand::redo () {
    for (unsigned int i=0; i < target_atoms.size (); i++) {
        set_color (target_atoms[i], present_colors [i]);
    }
    if (redraw) {
        vector <Molecule *> drawn_mols;
        for (unsigned int i =0; i < target_atoms.size (); i++) {
            Molecule * mol = (Molecule *) target_atoms[i] -> GetParent ();
            bool redr = true;
            for (unsigned int j =0; j < drawn_mols.size (); j++) {    
                if (mol == drawn_mols[j]) {
                    redr = false;
                    break;
                }
            }
            if (redr)  {
                gl -> draw_molecule (mol);
                drawn_mols.push_back (mol);
            }
        }
    }
    redraw = true;
}

void ColorAtomCommand::add (Atom *at, color col) {
    target_atoms.push_back (at);
    color c = get_color (at);
    present_colors.push_back ( col);
    previous_colors.push_back ( c );  
    
}


void ColorAtomCommand::undo () {
    for (unsigned int i=0; i < target_atoms.size (); i++) {
        set_color (target_atoms[i], previous_colors [i]);
    }
    if (redraw) {
        vector <Molecule *> drawn_mols;
        for (unsigned int i =0; i < target_atoms.size (); i++) {
            Molecule * mol = (Molecule *) target_atoms[i] -> GetParent ();
            bool redr = true;
            for (unsigned int j =0; j < drawn_mols.size (); j++) {    
                if (mol == drawn_mols[j]) {
                    redr = false;
                    break;
                }
            }

            if (redr)  {
                gl -> draw_molecule (mol);
                drawn_mols.push_back (mol);
            }
        }
    }
}

int ColorAtomCommand::id () const { 
    return 1; 
}

void ColorAtomCommand::set_name () {
    stringstream ss;
    ss << target_atoms.size ();
    setText ( QString          (       (    string("color: ")+ ss.str ()+ string (" atoms") ).c_str ()));
}

bool ColorAtomCommand::mergeWith (const QUndoCommand *other) {
    if (other->id() != id()) 
            return false;
    const ColorAtomCommand *oth = static_cast<const ColorAtomCommand*>(other);
 //   if (! same_time (oth)) return false;
    target_atoms.push_back (oth ->target_atoms[0]);
    previous_colors.push_back (oth ->previous_colors[0]);
    present_colors.push_back (oth ->present_colors[0]);
    stringstream ss;
    ss << target_atoms.size ();
    setText ( QString          (       (    string("color: ")+ ss.str ()+ string (" atoms") ).c_str ()));

    return true;
}

//////////////////////////////////////

ChangeDisplayStyleCommand::ChangeDisplayStyleCommand (MyGl *gl_point) : Command ()
{
    redraw = false;
    gl = gl_point;
    setText ( QString (string("change display style: 1  bond").c_str ()));
   
}


ChangeDisplayStyleCommand::ChangeDisplayStyleCommand (Atom *at, int disp_style, MyGl *gl_point) : Command ()
{
    redraw = false;
    gl = gl_point;
    add (at, disp_style);
    setText ( QString (string("change display style: 1  atom").c_str ()));

}

ChangeDisplayStyleCommand::ChangeDisplayStyleCommand (Bond *bo, int disp_style, MyGl *gl_point) : Command ()
{
    gl = gl_point;
    add (bo, disp_style);
    setText ( QString (string("change display style: 1  bond").c_str ()));
  
}

void ChangeDisplayStyleCommand::add (Molecule *mol, int at_di, int bo_di) {
	target_molecules.push_back (mol);
	present_styles_molecule_atoms.push_back (at_di);
	present_styles_molecule_bonds.push_back (bo_di);
	previous_styles_molecule_atoms.push_back (mol -> atoms_default_ds);
	previous_styles_molecule_bonds.push_back (mol -> bonds_default_ds);
}

void ChangeDisplayStyleCommand::add (Bond *bo, int disp_style) {
    target_bonds.push_back (bo);
    present_styles_bonds.push_back (disp_style);
    previous_styles_bonds.push_back (get_ds (bo));  
}


void ChangeDisplayStyleCommand::add (Atom *at, int disp_style) {
    target_atoms.push_back (at);
    present_styles_atoms.push_back (disp_style);
    previous_styles_atoms.push_back (get_ds (at));  
}




void ChangeDisplayStyleCommand::set_name () {

    string out = "change display style: ";
    if (target_atoms.size ()) {
        stringstream ss;
        ss << target_atoms.size ();
        string plur = "s";
        if (target_atoms.size () == 1)    plur = "";
        out += ss.str();
        out += " atom"+plur;
        out += " ";
    }
    if (target_bonds.size ()) {
        stringstream ss;
        ss << target_bonds.size ();
        string plur = "s";
        if (target_bonds.size () == 1)    plur = "";
        out += ss.str();
        out += " bond"+plur;
        out += " ";
    }



    setText ( QString (out.c_str ()));
}


void ChangeDisplayStyleCommand::redo () {
    for (unsigned int i=0; i < target_atoms.size (); i++) {
        set_ds (target_atoms[i],present_styles_atoms [i]);
    }
    for (unsigned int i=0; i < target_bonds.size (); i++) {
        set_ds (target_bonds[i], present_styles_bonds [i]);
    }
    for (unsigned int i=0; i < target_molecules.size (); i++) {
        target_molecules [i] -> atoms_default_ds = present_styles_molecule_atoms [i];
        target_molecules [i] -> bonds_default_ds = present_styles_molecule_bonds [i];
    }
    if (redraw) {
        vector <Molecule *> drawn_mols;
        for (unsigned int i =0; i < target_atoms.size (); i++) {
            Molecule * mol = (Molecule *) target_atoms[i] -> GetParent ();
            bool redr = true;
            for (unsigned int j =0; j < drawn_mols.size (); j++) {    
                if (mol == drawn_mols[j]) {
                    redr = false;
                    break;
                }
            }
            if (redr)  {
                drawn_mols.push_back (mol);
                gl -> draw_molecule (mol);
            }
        }
        for (unsigned int i =0; i < target_bonds.size (); i++) {
            Molecule * mol = (Molecule *) target_bonds[i] -> GetBeginAtom () -> GetParent ();
            bool redr = true;
            for (unsigned int j =0; j < drawn_mols.size (); j++) {    
                if (mol == drawn_mols[j]) {
                    redr = false;
                    break;
                }
            }

            if (redr)  {
                drawn_mols.push_back (mol);
                gl -> draw_molecule (mol);

            }
        }
    }
    redraw = true;
}

void ChangeDisplayStyleCommand::undo () {
    for (unsigned int i=0; i < target_atoms.size (); i++) {
        set_ds (target_atoms [i], previous_styles_atoms [i]);
    }
    for (unsigned int i=0; i < target_bonds.size (); i++) {
        set_ds (target_bonds [i], previous_styles_bonds [i]);
    }
    
    for (unsigned int i=0; i < target_molecules.size (); i++) {
        target_molecules [i] -> atoms_default_ds = previous_styles_molecule_atoms [i];
        target_molecules [i] -> bonds_default_ds = previous_styles_molecule_bonds [i];
    }
    
    
    vector <Molecule *> drawn_mols;
    for (unsigned int i =0; i < target_atoms.size (); i++) {
        Molecule * mol = (Molecule *) target_atoms[i] -> GetParent ();
        bool redr = true;
        for (unsigned int j =0; j < drawn_mols.size (); j++) {    
            if (mol == drawn_mols[j]) {
                redr = false;
                break;
            }
        }
        if (redr)  {
            drawn_mols.push_back (mol);
            gl -> draw_molecule (mol);
        }

        for (unsigned int i =0; i < target_bonds.size (); i++) {
            Molecule * mol = (Molecule *) target_bonds[i] -> GetBeginAtom () -> GetParent ();
            bool redr = true;
            for (unsigned int j =0; j < drawn_mols.size (); j++) {    
                if (mol == drawn_mols[j]) {
                    redr = false;
                    break;
                }
            }

            if (redr)  {
                drawn_mols.push_back (mol);
                gl -> draw_molecule (mol);

            }
        }
    }

}
/*
int ChangeDisplayStyleCommand::id () const { 
    return 2; 
}


bool ChangeDisplayStyleCommand::mergeWith (const QUndoCommand *other) {
    if (other->id() != id()) 
            return false;
    const ChangeDisplayStyleCommand *oth = static_cast<const ChangeDisplayStyleCommand*>(other);
  //  if (! same_time (oth)) return false;
    if (oth -> target_atoms.size ()) {
        target_atoms.push_back (oth ->target_atoms[0]);
        assert (oth -> previous_styles_atoms.size ());
        assert (oth -> present_styles_atoms.size ());
        previous_styles_atoms.push_back (oth -> previous_styles_atoms[0]); 
        present_styles_atoms.push_back (oth -> present_styles_atoms[0]);
    }
    if (oth -> target_bonds.size ()) {
        target_bonds.push_back (oth ->target_bonds[0]);
        assert (oth -> previous_styles_bonds.size ());
        assert (oth -> present_styles_bonds.size ());
        previous_styles_bonds.push_back (oth -> previous_styles_bonds[0]); 
        present_styles_bonds.push_back (oth -> present_styles_bonds[0]);
    }

    set_name ();
    return true;
}
*/
/////////////////////////////////////////////////////////////////////////////

ChangeVectorCommand::ChangeVectorCommand (vect &vec, vect val, MyGl *gl_point, string name) : Command ()
{
    gl = gl_point;
    setText ( QString (name.c_str ()));
    v = &vec;
    previous_value = vec;
    present_value = val;
  
}




void ChangeVectorCommand::undo () {
    *v = previous_value;
    gl -> updateGL ();
}


void ChangeVectorCommand::redo () {
    *v = present_value;
    gl -> updateGL ();
}

///////////////////////////////////////////////////////////////////////////////////

/*

ChangeMoleculeCommand::ChangeMoleculeCommand (MyGl *gl_point, string name, bool first_d) : Command ()
    
{
    first_do = first_d;
    gl = gl_point;
    setText ( QString (name.c_str ()));
    mol = NULL;
    present_value = NULL;
    previous_value = NULL;

}


void ChangeMoleculeCommand::add (Molecule *m, Molecule *prev_val) {
    cerr << "adding"<<endl;
    mol = m;
    if (previous_value) delete previous_value;
    if (present_value) delete present_value;
    previous_value = new Molecule ();
    previous_value -> copy_from (prev_val);
    present_value = new Molecule ();
    present_value -> copy_from (mol);
    cerr << "added"<<endl;
}

void ChangeMoleculeCommand::undo () {
 //   for (unsigned int i=0; i < mols.size (); i++) {
 //       Molecule *mol = mols [i];
 //       *mol = previous_values [i];
   //     gl -> draw_molecule (mol);


  //  }

}


void ChangeMoleculeCommand::redo () {
   // for (unsigned int i=0; i < mols.size (); i++) {
     //   Molecule *mol = mols [i];
      //  *mol = present_values [i];
//        gl -> draw_molecule (mol);

    //}
}
////////////////////////////////////////////////
*/

MoveAtomsCommand::MoveAtomsCommand (MyGl *gl_point, int mod) : Command ()
{
    mode = mod;
    redraw = true;
    gl = gl_point;
    setText ( QString (string("move atoms").c_str ()));
   
}


void MoveAtomsCommand::redo () {
    if (mode != 1) {
        for (unsigned int i=0; i < target_atoms.size (); i++) {
            set_coordinates (target_atoms [i], present_coordinates [i]);
        }
        if (redraw) {
            vector <Molecule *> drawn_mols;
            for (unsigned int i =0; i < target_atoms.size (); i++) {
                Molecule * mol = (Molecule *) target_atoms[i] -> GetParent ();
                bool redr = true;
                for (unsigned int j =0; j < drawn_mols.size (); j++) {    
                    if (mol == drawn_mols[j]) {
                        redr = false;
                        break;
                    }
                }
                if (redr)  {
                    gl -> draw_molecule (mol);
                    drawn_mols.push_back (mol);
                }
            }
        }
        redraw = true;
    }
}
void MoveAtomsCommand::add (Atom *at, vect coo) {
    target_atoms.push_back (at);
    previous_coordinates.push_back ( coo);
    present_coordinates.push_back ( get_coordinates(at));  
    
}


void MoveAtomsCommand::undo () {
    if (mode != 0) {
        for (unsigned int i=0; i < target_atoms.size (); i++) {
            set_coordinates (target_atoms [i] ,previous_coordinates [i]);
        }
        if (redraw) {
            vector <Molecule *> drawn_mols;
            for (unsigned int i =0; i < target_atoms.size (); i++) {
                Molecule * mol = (Molecule *) target_atoms[i] -> GetParent ();
                bool redr = true;
                for (unsigned int j =0; j < drawn_mols.size (); j++) {    
                    if (mol == drawn_mols[j]) {
                        redr = false;
                        break;
                    }
                }

                if (redr)  {
                    gl -> draw_molecule (mol);
                    drawn_mols.push_back (mol);
                }
            }
    }
    }
}


void MoveAtomsCommand::set_name () {
    stringstream ss;
    ss << target_atoms.size ();
    setText ( QString          (       (    string("move: ")+ ss.str ()+ string (" atoms") ).c_str ()));
}

//////////////////////////////////////////////////////////////////////

CreateGraphicalObjectCommand::CreateGraphicalObjectCommand (GraphicalObject *obj, DDWin *ddwin_point): Command ()
 {
    rank = -1;
    ddwin = ddwin_point;
    setText ( QString ("Create Graphical Object"));
    object = obj;
}



void CreateGraphicalObjectCommand::redo () {

    if (rank == -1) ddwin -> graphical_objects.push_back (object);
    else ddwin -> graphical_objects.insert (ddwin -> graphical_objects.begin () + rank, object);
    ddwin -> graphical_objects_menu -> update_slot ();
    ddwin -> gl -> updateGL ();
}

void CreateGraphicalObjectCommand::undo () {
    for (unsigned int i =0; i< ddwin -> graphical_objects.size (); i++) {
        if (ddwin -> graphical_objects[i] == object) {
            ddwin -> graphical_objects.erase (ddwin -> graphical_objects.begin ()+i);
            ddwin -> graphical_objects_menu -> update_slot ();
            rank = i;
            ddwin -> gl -> updateGL ();
            break;
        }
    }
}
//////////////////////////////////////////////////////////////////

DeleteGraphicalObjectCommand::DeleteGraphicalObjectCommand (GraphicalObject *obj, DDWin *ddwin_point) : Command ()
{
    rank = -1;
    ddwin = ddwin_point;
    setText ( QString ("Delete Graphical Object"));

    object = obj;
}



void DeleteGraphicalObjectCommand::undo () {

    if (rank == -1) ddwin -> graphical_objects.push_back (object);
    else ddwin -> graphical_objects.insert (ddwin -> graphical_objects.begin () + rank, object);
    ddwin -> graphical_objects_menu -> update_slot ();
    ddwin -> gl -> updateGL ();
}

void DeleteGraphicalObjectCommand::redo () {
    for (unsigned int i =0; i< ddwin -> graphical_objects.size (); i++) {
        if (ddwin -> graphical_objects[i] == object) {
            ddwin -> graphical_objects.erase (ddwin -> graphical_objects.begin ()+i);
            ddwin -> graphical_objects_menu -> update_slot ();
            rank = i;
            ddwin -> gl -> updateGL ();
            break;
        }
    }
}









////////////////////////////////////////////////////////////////////

CreateMoleculeCommand::CreateMoleculeCommand ( Molecule * molecule, DDWin *ddwin_point) : Command ()
{
    ddwin = ddwin_point;
    builder = ddwin -> builder;
    setText ( QString ("Create ")+QString (molecule -> GetTitle ()));
    mol = molecule;
    mol -> ZNinit ();
 //   ddwin->gl->draw_molecule (ddwin->target_molecule); 
}

void CreateMoleculeCommand::redo () {
    ddwin -> add_molecule (mol);
  //  cerr << "conformers "<<mol -> NumConformers () << endl;
    ddwin -> set_current_target (-1);
    ddwin -> gl -> draw_molecule (mol);
    ddwin -> gl -> updateGL ();


}


void CreateMoleculeCommand::undo () {
    ddwin -> remove_molecule (mol);
    ddwin -> set_current_target (0);
    ddwin -> gl -> draw_molecule (mol);
    ddwin -> gl -> updateGL ();
}
///////////////////////////////////////////////////

DeleteMoleculeCommand::DeleteMoleculeCommand (Molecule *molecule, DDWin *ddwin_point): CreateMoleculeCommand (molecule, ddwin_point) {
    setText ( QString ("Delete ")+QString (molecule -> GetTitle ()));  
}


void DeleteMoleculeCommand::redo () {
    CreateMoleculeCommand::undo ();
}


void DeleteMoleculeCommand::undo () {
    CreateMoleculeCommand::redo ();
}


////////////////////////////////////////////////
MutateAtomCommand::MutateAtomCommand ( Atom *at, int atn, DDWin *ddwin_point): Command ()
 {
    ddwin = ddwin_point;
    builder = ddwin -> builder;
    setText ( QString ("Mutate Atom"));
    currentn = atn;
    precn = at-> GetAtomicNum ();
    atom = at;


    Molecule *mol = (Molecule *) atom -> GetParent ();
 //   first_time = true;

    builder -> save_Hs (atom, prev_H_atoms, prev_H_bonds);
    builder -> delete_Hs (at);
    atom -> SetAtomicNum (currentn);
    mol -> ZNAddHydrogens (atom);    
    builder -> save_Hs (atom, curr_H_atoms, curr_H_bonds);
    mol -> ZNinit ();

}

void MutateAtomCommand::redo () {
    Molecule *mol = (Molecule *) atom -> GetParent ();
    if (!first_time) {
        builder -> delete_Hs (atom);
        atom -> SetAtomicNum (currentn);
        builder -> add_Hs (atom, curr_H_atoms, curr_H_bonds);

    }
    first_time = false;
    ddwin -> gl -> draw_molecule (mol);


}

void MutateAtomCommand::undo () {
    Molecule *mol = (Molecule *) atom -> GetParent ();
    builder -> delete_Hs (atom);
    atom -> SetAtomicNum (precn);
    mol -> set_color_mw (atom) ;
    builder -> add_Hs (atom, prev_H_atoms, prev_H_bonds);
    ddwin -> gl -> draw_molecule (mol);
}
///////////////////////////////////////////////////////////

DeleteAtomCommand::DeleteAtomCommand (Atom *at, DDWin *ddwin_point) : Command ()
{
 /*   ddwin = ddwin_point;
    builder = ddwin -> builder;
    setText ( QString ("Remove atom"));
    atom = at;
    for (unsigned int i = 0; i < at -> bound.size (); i++) {
        if (atom -> bound [i] -> isHydrogen ()) {
            bound.push_back (atom -> bound[i]);
            bonds.push_back (atom -> bonds[i]);
            vector <Atom *> prevH;
            vector <Atom *> postH;
            vector <Bond *> prevb;
            vector <Bond *> postb;
            prev_targets_H_atoms.push_back (prevH);
            curr_targets_H_atoms.push_back (postH);
            prev_targets_H_bonds.push_back (prevb);
            curr_targets_H_bonds.push_back (postb);

        }
    }

*/
}


void DeleteAtomCommand::redo () {
/*    Molecule *mol = atom -> GetParent ();
    for (unsigned int i = 0; i < bonds.size (); i++) {
        mol -> remove (bonds[i]);
    }
    builder -> delete_Hs (atom, prev_H_atoms, prev_H_bonds);
    mol -> remove (atom);
    for (unsigned int i = 0; i < bonds.size (); i++) {
        builder -> delete_Hs (bound[i], prev_targets_H_atoms[i], prev_targets_H_bonds[i]);
        builder -> add_H (bound[i], curr_targets_H_atoms[i], curr_targets_H_bonds[i]);
    }

    builder -> redefine_mol (mol);
    ddwin -> gl -> draw_molecule (atom->residue->molecule);
*/
}



void DeleteAtomCommand::undo () {
/*    Molecule *mol = atom -> GetParent ();
    mol -> add (atom);
    for (unsigned int i = 0; i < bonds.size (); i++) {
        mol -> add (bonds[i]);
    }
    builder -> delete_Hs (atom, curr_H_atoms, curr_H_bonds);
    builder -> add_H (atom, prev_H_atoms, prev_H_bonds);
    for (unsigned int i = 0; i < bonds.size (); i++) {
        builder -> delete_Hs (bound[i], curr_targets_H_atoms[i], curr_targets_H_bonds[i]);
        builder -> add_H (bound[i], prev_targets_H_atoms[i], prev_targets_H_bonds[i]);
    }
    
    builder -> redefine_mol (mol);
    ddwin -> gl -> draw_molecule (atom->residue->molecule);
*/
}





//////////////////////////////////////////////////////////////



AddAtomCommand::AddAtomCommand ( Atom *at, Bond *bo, Atom *part,  DDWin *ddwin_point) : Command ()
{
    Molecule *mol = (Molecule *) part -> GetParent ();
 //   first_time = true;
    ddwin = ddwin_point;
    builder = ddwin -> builder;
    setText ( QString ("Add atom"));
    atom = at;
    partner = part;
    bond = bo;
    builder -> save_Hs (partner, prev_H_atoms, prev_H_bonds);
  //  builder -> save_Hs (at, prev_target_H_atoms, prev_target_H_bonds);
    builder -> delete_Hs (partner);
    mol -> add_atom_bonded_to (atom, bond, partner);
    mol -> ZNAddHydrogens (partner);
    mol -> ZNAddHydrogens (at);    
    builder -> save_Hs (partner, curr_H_atoms, curr_H_bonds);
    builder -> save_Hs (at, curr_target_H_atoms, curr_target_H_bonds);
    mol -> ZNinit ();
}


void AddAtomCommand::redo () {
    if (!first_time) {

        Molecule *mol = (Molecule *) partner -> GetParent ();
        builder -> delete_Hs (partner);
        mol -> add_atom_bonded_to (atom, bond, partner);

      
    //    builder -> redefine_mol (mol);

  //      builder -> delete_Hs (atom);
        builder -> add_Hs (partner, curr_H_atoms, curr_H_bonds);
        builder -> add_Hs (atom, curr_target_H_atoms, curr_target_H_bonds);

    //    builder -> redefine_mol (mol);

        ddwin -> gl -> draw_molecule (mol);
    }
    first_time = false;
}

void AddAtomCommand::undo () {
    Molecule *mol = (Molecule *) atom -> GetParent ();
    builder -> delete_Hs (atom);
    mol -> RemoveAtom (atom);
  //  mol -> RemoveBond (bond);
    builder -> delete_Hs (partner);
    builder -> add_Hs (partner, prev_H_atoms, prev_H_bonds);
    ddwin -> gl -> draw_molecule (mol);



}


///////////////////////////////////////////////////////////

AddBondCommand::AddBondCommand (Bond *bo,  DDWin *ddwin_point): Command ()
 {
    ddwin = ddwin_point;
    builder = ddwin -> builder;
    setText ( QString ("Add bond"));
    bond = bo;
    Atom *at1 = bond -> GetBeginAtom ();
    Atom *at2 = bond -> GetEndAtom ();
    Molecule *mol = (Molecule *) at1 -> GetParent ();

    builder -> save_Hs (at1, prev_H_atoms0, prev_H_bonds0);
    builder -> save_Hs (at2, prev_H_atoms1, prev_H_bonds1);


    builder -> delete_Hs (at1);
    builder -> delete_Hs (at2);
    mol -> ZNAddBond (bond);
    mol -> ZNAddHydrogens (at1);
    mol -> ZNAddHydrogens (at2);    
    builder -> save_Hs (at1, curr_H_atoms0, curr_H_bonds0);
    builder -> save_Hs (at2, curr_H_atoms1, curr_H_bonds1);
    mol -> ZNinit ();

}


void AddBondCommand::redo () {
    Atom *at1 = bond -> GetBeginAtom ();
    Atom *at2 = bond -> GetEndAtom ();
    Molecule *mol = (Molecule *) at1 -> GetParent ();
    if (!first_time) {
        builder -> delete_Hs (at1);
        builder -> delete_Hs (at2);
        mol -> ZNAddBond (bond);
        builder -> add_Hs (at1, curr_H_atoms0, curr_H_bonds0);
        builder -> add_Hs (at2, curr_H_atoms1, curr_H_bonds1);    
    }
    ddwin -> gl -> draw_molecule (mol);
    first_time = false;

}

void AddBondCommand::undo () {
    Atom *at1 = bond -> GetBeginAtom ();
    Atom *at2 = bond -> GetEndAtom ();
    Molecule *mol = (Molecule *) at1 -> GetParent ();
    builder -> delete_Hs (at1);
    builder -> delete_Hs (at2);
    mol -> RemoveBond (bond);
    builder -> add_Hs (at1, prev_H_atoms0, prev_H_bonds0);
    builder -> add_Hs (at2, prev_H_atoms1, prev_H_bonds1);    
    ddwin -> gl -> draw_molecule (mol);
}


DeleteBondCommand::DeleteBondCommand (Bond *bo,  DDWin *ddwin_point) : AddBondCommand (bo, ddwin_point) {
    setText ( QString ("Delete bond"));

}

void DeleteBondCommand::undo () {
    AddBondCommand::redo ();
}


void DeleteBondCommand::redo () {
    AddBondCommand::undo ();
}
///////////////////////////////////////////////////////////

ModifyBondCommand::ModifyBondCommand (Bond *bo, int new_o,  DDWin *ddwin_point, bool red) : Command ()
{
    first_time = true;
    redefine = red;
    new_order = new_o;
    old_order = bo -> GetBondOrder ();
    ddwin = ddwin_point;
    builder = ddwin -> builder;
    setText ( QString ("Change Bond Order"));
    bond = bo;

    Atom *at1 = bond -> GetBeginAtom ();
    Atom *at2 = bond -> GetEndAtom ();
    Molecule *mol = (Molecule *) at1 -> GetParent ();
    builder -> save_Hs (at1, prev_H_atoms0, prev_H_bonds0);
    builder -> save_Hs (at2, prev_H_atoms1, prev_H_bonds1);
  //  builder -> save_Hs (at, prev_target_H_atoms, prev_target_H_bonds);
    builder -> delete_Hs (at1);
    builder -> delete_Hs (at2);

    bond -> SetBondOrder (new_order);
    mol -> ZNAddHydrogens (at1);
    mol -> ZNAddHydrogens (at2);    
    builder -> save_Hs (at1, curr_H_atoms0, curr_H_bonds0);
    builder -> save_Hs (at2, curr_H_atoms1, curr_H_bonds1);
    mol -> ZNinit ();

}


void ModifyBondCommand::redo () {
        Atom *at1 = bond -> GetBeginAtom ();
        Atom *at2 = bond -> GetEndAtom ();
        Molecule *mol = (Molecule *) at1 -> GetParent ();
    if (!first_time) {

        builder -> delete_Hs (at1);
        builder -> delete_Hs (at2);
        bond -> SetBondOrder (new_order);
        builder -> add_Hs (at1, curr_H_atoms0, curr_H_bonds0);
        builder -> add_Hs (at2, curr_H_atoms1, curr_H_bonds1);    

    }
        ddwin -> gl -> draw_molecule (mol);
    first_time = false;
}


void ModifyBondCommand::undo () {
    Atom *at1 = bond -> GetBeginAtom ();
    Atom *at2 = bond -> GetEndAtom ();
    Molecule *mol = (Molecule *) at1 -> GetParent ();
    builder -> delete_Hs (at1);
    builder -> delete_Hs (at2);
    bond -> SetBondOrder (old_order);
    builder -> add_Hs (at1, prev_H_atoms0, prev_H_bonds0);
    builder -> add_Hs (at2, prev_H_atoms1, prev_H_bonds1);      
    ddwin -> gl -> draw_molecule (mol);

}


///////////////////////////////////////////////////////////

RedefineMoleculeCommand::RedefineMoleculeCommand (Molecule *mo, Builder *bu, int mod) : Command ()
{
/*    mode = mod;
    builder = bu;
    setText ( QString ("Redefine Mol"));
    mol = mo;
*/
}


void RedefineMoleculeCommand::redo () {
/*    if (mode !=1 ) {
        builder -> redefine_mol (mol);
        builder -> ddwin -> gl -> draw_molecule (mol);
    }
*/
}

void RedefineMoleculeCommand::undo () {
/*    if (mode != 0) {
        builder -> redefine_mol (mol);
        builder -> ddwin -> gl -> draw_molecule (mol);
    }
*/
}


///////////////////////////////////////////////////////////////////////


ReprotonateCommand::ReprotonateCommand (Molecule *molecule, DDWin *ddw) {
    setText ( QString ("Reprotonate"));
	ddwin = ddw;
	mol = molecule;
    first_time = true;	
	FOR_ATOMS_OF_MOL (a, mol) {
		if (a ->IsHydrogen ()) {
			prev_H_atoms.push_back (&*a);
			OBBondIterator i;
			Atom *b = a -> BeginNbrAtom (i);
			prev_partners.push_back (&*b);
			prev_H_bonds.push_back (mol -> GetBond (&*a, b));
	//		prec_bonds.push_back ()
		}
	}
	
        ddwin ->builder -> delete_Hs (mol);
        mol -> ZNAddHydrogens ();
        ddwin ->gl -> draw_molecule (mol);
	
	FOR_ATOMS_OF_MOL (a, mol) {
		if (a ->IsHydrogen ()) {
			curr_H_atoms.push_back (&*a);
			OBBondIterator i;
			Atom *b = a -> BeginNbrAtom (i);
			curr_partners.push_back (&*b);
			curr_H_bonds.push_back (mol -> GetBond (&*a, b));
	//		prec_bonds.push_back ()
		}
	}
	
}


void ReprotonateCommand::redo () {
	if (!first_time) {
		ddwin ->builder -> delete_Hs (mol);
		for (unsigned int i = 0; i< curr_H_atoms.size (); i++) {
            mol -> add_atom_bonded_to (curr_H_atoms[i], curr_H_bonds[i], curr_partners[i]);
        }

	}
    ddwin -> gl -> draw_molecule (mol);
    first_time = false;
}


void ReprotonateCommand::undo () {
	ddwin ->builder -> delete_Hs (mol);
	for (unsigned int i = 0; i< prev_H_atoms.size (); i++) {
		mol -> add_atom_bonded_to (prev_H_atoms[i], prev_H_bonds[i], prev_partners[i]);
	}
	ddwin -> gl -> draw_molecule (mol);
}
