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


#ifndef COMMAND_H
#define COMMAND_H

#include <QUndoCommand>
#include "ZNmolecule.h"
#include "constants.h"
#include "ddwin.h"
#include <time.h>

class MyGl; class DDWin; class Builder;

class Command : public QUndoCommand
{

public:
Command ();
//    virtual bool mergeWith (const QUndoCommand *other);
void get_time ();
//time_t commandtime;
bool same_time (Command *other);
void name (string str);
bool first_time;

};



class ColorAtomCommand : public Command {
    public:
    ColorAtomCommand (MyGl *gl);
    ColorAtomCommand (Atom *at, color col, MyGl *gl);
    void add (Atom *at, color col);
    vector <color> previous_colors;
    vector <color> present_colors;
    vector <Atom *>target_atoms;
    void undo ();
    void redo ();
    int id () const;
    bool mergeWith (const QUndoCommand *command);
    void set_name ();
    private:
    MyGl *gl;
};


class ChangeDisplayStyleCommand : public Command {
    public:
    bool redraw;
    ChangeDisplayStyleCommand (MyGl *gl);
    ChangeDisplayStyleCommand (Atom *at, int dis, MyGl *gl);
    ChangeDisplayStyleCommand (ZNBond *bo, int dis, MyGl *gl);
    void add (Atom *at, int dis);
    void add (ZNBond *bo, int dis);
    void add (ZNMolecule *mol, int at_ds, int bo_ds, int backbone_ds);
    void set_name ();
    vector <int> previous_styles_atoms;
    vector <int> present_styles_atoms;
    vector <int> previous_styles_bonds;
    vector <int> present_styles_bonds;
    vector <int> present_styles_molecule_bonds;
    vector <int> present_styles_molecule_atoms;
	vector <int> present_styles_molecule_backbone;
    vector <int> previous_styles_molecule_atoms;
    vector <int> previous_styles_molecule_bonds;
	vector <int> previous_styles_molecule_backbone;

    vector <Atom *>target_atoms;
    vector <ZNBond *>target_bonds;
    vector <ZNMolecule *>target_molecules;

    void undo ();
    void redo ();
 //   int id () const;
 //   bool mergeWith (const QUndoCommand *command);
    private:
    MyGl *gl;
};


class ChangeFloatCommand : public Command {
    public:
    ChangeFloatCommand (float &f, float new_val, MyGl *gl, string name = "change value");
    float *f;
    float present_value;
    float previous_value;

    void undo ();
    void redo ();
//    int id () const;
//    bool mergeWith (const QUndoCommand *command);
    private:
    MyGl *gl;	
};

class ChangeVectorCommand : public Command {
    public:
    ChangeVectorCommand (vect &v, vect new_val, MyGl *gl, string name = "change vector");
    vect *v;
    vect present_value;
    vect previous_value;

    void undo ();
    void redo ();
//    int id () const;
//    bool mergeWith (const QUndoCommand *command);
    private:
    MyGl *gl;
};


/*
class ChangeZNMoleculeCommand : public Command {
    public:
    ChangeZNMoleculeCommand (MyGl *gl, string name = "change molecules", bool first_do = FALSE);

    bool first_do;
    ZNMolecule * mol;
    ZNMolecule*  present_value;
    ZNMolecule*  previous_value;
    void add (ZNMolecule *mol, ZNMolecule *previous_val);

    void undo ();
    void redo ();
//    int id () const;
//    bool mergeWith (const QUndoCommand *command);
    private:
    MyGl *gl;
};
*/

class MoveAtomsCommand : public Command {
    public:
    bool redraw;
    MoveAtomsCommand (MyGl *gl, int mod = 2); //0 only on redo, 1 only on undo, 2 both
    void add (Atom *at, vect vec);
    int mode;
    vector <vect> previous_coordinates;
    vector <vect> present_coordinates;
    vector <Atom *>target_atoms;
    void undo ();
    void redo ();
 //   int id () const;
 //   bool mergeWith (const QUndoCommand *command);
    void set_name ();
    private:
    MyGl *gl;
};


class CreateGraphicalObjectCommand : public Command {
    public:
    CreateGraphicalObjectCommand (GraphicalObject *obj, DDWin *ddwin_pt);
    int rank;
    GraphicalObject *object;
    void undo ();
    void redo ();

   // void set_name ();
    private:
    DDWin *ddwin;
};

class DeleteGraphicalObjectCommand : public Command {
    public:
    DeleteGraphicalObjectCommand (GraphicalObject *obj, DDWin *ddwin_pt);
    int rank;
    GraphicalObject *object;
    void undo ();
    void redo ();

  //  void set_name ();
    private:
    DDWin *ddwin;
};


class CreateZNMoleculeCommand : public Command {
    public:
    CreateZNMoleculeCommand (ZNMolecule *mol, DDWin *ddwin_pt);
    ZNMolecule *mol;
    virtual void undo ();
    virtual void redo ();
    private:
    DDWin *ddwin;
    Builder *builder;
};


class DeleteZNMoleculeCommand : public CreateZNMoleculeCommand {
    public:
    DeleteZNMoleculeCommand (ZNMolecule *mol, DDWin *ddwin_pt);
    void undo ();
    void redo ();

};


class MutateAtomCommand : public Command {
    public:
    MutateAtomCommand (Atom *at, int atn , DDWin *ddwin_pt);
    Atom *atom;
    int precn;
    int currentn;
    void undo ();
    void redo ();
    vector <Atom *> prev_H_atoms;
    vector <ZNBond *> prev_H_bonds;
    vector <Atom *> curr_H_atoms;
    vector <ZNBond *> curr_H_bonds;

    DDWin *ddwin;
    Builder *builder;
};


class AddAtomCommand : public Command {
    public:

    vector <Atom *> prev_H_atoms;
    vector <ZNBond *> prev_H_bonds;
    vector <Atom *> curr_H_atoms;
    vector <ZNBond *> curr_H_bonds;

    vector <Atom *> curr_target_H_atoms;
    vector <ZNBond *> curr_target_H_bonds;

    vector <Atom *> prev_target_H_atoms;
    vector <ZNBond *> prev_target_H_bonds;


    AddAtomCommand (Atom *to_add, ZNBond *bond, Atom *bond_partner, DDWin *ddwin_pt);
    Atom *atom;
    ZNBond *bond;
    Atom *partner;
    void undo ();
    void redo ();
    DDWin *ddwin;
    Builder *builder;
};


class DeleteAtomCommand : public Command {
    public:
    vector <Atom *> prev_H_atoms;
    vector <ZNBond *> prev_H_bonds;
    vector <Atom *> curr_H_atoms;
    vector <ZNBond *> curr_H_bonds;
    vector <Atom *> neighs;
    vector <vector <Atom *> > curr_targets_H_atoms;
    vector <vector <ZNBond *> > curr_targets_H_bonds;

    vector <vector <Atom *> > prev_targets_H_atoms;
    vector <vector <ZNBond *> > prev_targets_H_bonds;


    DeleteAtomCommand (Atom *to_remove, DDWin *ddwin_pt);
    Atom *atom;
    vector <ZNBond *> bonds;

    void undo ();
    void redo ();
    DDWin *ddwin;
    Builder *builder;
};




class AddBondCommand : public Command {
    public:
    AddBondCommand (ZNBond *bo, DDWin *ddwin);
    vector <Atom *> prev_H_atoms0;
    vector <ZNBond *> prev_H_bonds0;
    vector <Atom *> curr_H_atoms0;
    vector <ZNBond *> curr_H_bonds0;

    vector <Atom *> prev_H_atoms1;
    vector <ZNBond *> prev_H_bonds1;
    vector <Atom *> curr_H_atoms1;
    vector <ZNBond *> curr_H_bonds1;


    ZNBond *bond;
    virtual void undo ();
    virtual void redo ();
    DDWin *ddwin;
    Builder *builder;
};


class DeleteBondCommand : public AddBondCommand {
    public:
    DeleteBondCommand (ZNBond *bo, DDWin *ddwin);
    void undo ();
    void redo ();
};

class ModifyBondCommand : public Command {
    public:
    ModifyBondCommand (ZNBond *bo, int ord, DDWin *ddwin, bool redefine = true);
    bool redefine;
    int new_order, old_order;
    vector <Atom *> prev_H_atoms0;
    vector <ZNBond *> prev_H_bonds0;
    vector <Atom *> curr_H_atoms0;
    vector <ZNBond *> curr_H_bonds0;

    vector <Atom *> prev_H_atoms1;
    vector <ZNBond *> prev_H_bonds1;
    vector <Atom *> curr_H_atoms1;
    vector <ZNBond *> curr_H_bonds1;


    ZNBond *bond;
    void undo ();
    void redo ();
    DDWin *ddwin;
    Builder *builder;
};


class RedefineZNMoleculeCommand : public Command {
    public:
    RedefineZNMoleculeCommand (ZNMolecule *mol, Builder *bu, int mod = 2); //0 only on redo, 1 only on undo, 2 both
    int mode;
    ZNMolecule *mol;
    void undo ();
    void redo ();
    Builder *builder;
};


class ReprotonateCommand : public Command {
	public:
	ReprotonateCommand (ZNMolecule *mol, DDWin *ddw, double ph = 7.4);
	ZNMolecule *mol;
	void undo ();
	void redo ();
	DDWin *ddwin;
	bool first_time;
	
	vector <Atom *> prev_H_atoms;
    vector <ZNBond *> prev_H_bonds;
    vector <Atom *> curr_H_atoms;
    vector <ZNBond *> curr_H_bonds;
	vector <Atom *> prev_partners;
	vector <Atom *> curr_partners;
	 
};

#endif
