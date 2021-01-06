/***************************************************************************
 *   Copyright (C) 2007 by Harm van Eersel                                 *
 *   devsciurus@xs4all.nl                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/** @file
 * This file is part of molsKetch and contains the command classes used in
 * the Qt Redo/Undo system. Each class represents an action on the scene and
 * contains the code to redo/undo that action.
 *
 * @author Harm van Eersel <devsciurus@xs4all.nl>
 * @since Deuterium
 */


#ifndef commands_H
#define commands_H

#include <QUndoCommand>
#include <QPointF>

class MsKAtom;
class Molecule;
class MsKBond;
class MolScene;
class QGraphicsItem;
class QGraphicsScene;
class QTransform;

namespace Commands
{

// MsKAtom command

/**
 * Command to add an atom
 *
 * @author Harm van Eersel
 */
class AddAtom : public QUndoCommand
  {
  public:
    /** 
     * Creates a new AddAtom command.
     *
     * @param newAtom pointer to the atom that should be added
     * @param molecule pointer to the molecule for the new atom
     * @param text a description of the command
     */
    AddAtom(MsKAtom* newAtom, Molecule* molecule, const QString & text = "");
    /**
     * Destructor 
     *
     * Deletes m_atom if command is destructed in an undone state.
     */
    ~AddAtom();
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** MsKAtom of this command. */
    MsKAtom* m_atom;
    /** Molecule of this command. */
    Molecule* m_molecule;
  };


/**
 * Command to change the element of an atom
 *
 * @author Harm van Eersel
 */
class ChangeElement : public QUndoCommand
  {
  public:
    /** 
     * Creates a new ChangeElement command.
     *
     * @param changeMsKAtom the atom which element symbol should be changed
     * @param newElementSymbol the new element symbol
     * @param text a description of the command
     */
    ChangeElement(MsKAtom* changeMsKAtom, const QString & newElementSymbol, const QString & text = "");
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** Old element symbol of the atom. */
    QString m_oldName;
    /** New element symbol of the atom. */
    QString m_newName;
    /** MsKAtom of this command. */
   MsKAtom* m_atom;
  };
  
/**
 * Command to decrease the charge of an atom
 *
 * @author Harm van Eersel
 */
class DecCharge : public QUndoCommand
  {
  public:
    /** 
     * Constructor
     *
     * @param atom the atom to decrease the charge of
     * @param text a description of the command
     */
    DecCharge(MsKAtom* atom, const QString & text = "");
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** MsKAtom of this command. */
    MsKAtom* m_atom;
  };
  
 /**
 * Command to increase the charge of an atom
 *
 * @author Harm van Eersel
 */
class IncCharge : public QUndoCommand
  {
  public:
    /** 
     * Constructor
     *
     * @param atom the atom to increase the charge of
     * @param text a description of the command
     */
    IncCharge(MsKAtom* atom, const QString & text = "");
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** MsKAtom of this command. */
    MsKAtom* m_atom;
  };
  
/**
 * Command to add an implicit hydrogen
 *
 * @author Harm van Eersel
 */
class AddImplicitHydrogen : public QUndoCommand
  {
  public:
    /** 
     * Constructor
     *
     * @param atom the atom to decrease the charge of
     * @param text a description of the command
     */
    AddImplicitHydrogen(MsKAtom* atom, const QString & text = "");
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** MsKAtom of this command. */
    MsKAtom* m_atom;
  };
  
 /**
 * Command to remove an implicit hydrogen
 *
 * @author Harm van Eersel
 */
class RemoveImplicitHydrogen : public QUndoCommand
  {
  public:
    /** 
     * Constructor
     *
     * @param atom the atom to increase the charge of
     * @param text a description of the command
     */
    RemoveImplicitHydrogen(MsKAtom* atom, const QString & text = "");
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** MsKAtom of this command. */
    MsKAtom* m_atom;
  };

/**
 * Command to delete an atom
 *
 * @author Harm van Eersel
 */
class DelAtom : public QUndoCommand
  {
  public:
    /**
     * Creates a new DelAtom command.
     *
     * @param delAtom the atom to be removed
     * @param text a description of the command
     */
    DelAtom(MsKAtom* delAtom, const QString & text = "");
    /**
     * Destructor
     *
     * Deletes m_atom if the command is destructed in a done state.
     */
    virtual ~DelAtom();
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** MsKAtom of this command. */
    MsKAtom* m_atom;
    /** Molecule of this command. */
    Molecule* m_molecule;
    /** The list of bonds that were connected to m_atom. */
    QList<MsKBond*> m_bondList;
  };

// MsKBond commands

/**
 * Command to add a bond
 *
 * @author Harm van Eersel
 */
class AddBond : public QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param newBond the new bond to add
     * @param text a description of the command
     */
    
    AddBond(MsKBond* newBond, const QString & text = "");
    /**
     * Destructor
     *
     * Deletes m_bond if the command is in an undone state.
     */
    ~AddBond();
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The bond of this command. */
    MsKBond* m_bond;
    /** Molecule of this command. */
    Molecule* m_mol;
  };


/**
 * Command to remove a bond
 *
 * @author Harm van Eersel
 */
class DelBond : public  QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param delBond the bond that should be removed
     * @param text a description of the command
     */
    DelBond(MsKBond* delBond, const QString & text = "");
    /**
     * Destructor
     *
     * Deletes m_bond if the command is destructed in a done state.
     */
    virtual ~DelBond();
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The bond of this command. */
    MsKBond* m_bond;
    /** Molecule of this command. */
    Molecule* m_mol;
  };


/**
 * Command to increase the bond type
 *
 * @author Harm van Eersel
 */
class IncType : public  QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param incBond bond of which the type should be changed
     * @param text a description of the command
     */
    IncType(MsKBond* incBond, const QString & text = "");
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The bond of this command. */
    MsKBond* m_bond;
  };


/**
 * Command to increase the bond order
 *
 * @author Harm van Eersel
 */
class IncOrder : public  QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param incBond bond to increase the order of 
     * @param text a description of the command
     */
    IncOrder(MsKBond* incBond, const QString & text = "");
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The bond of this command. */
    MsKBond* m_bond;
  };


// Molecule commands


/**
 * Command to merge two molecules
 *
 * @author Harm van Eersel
 */
class MergeMol : public QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param oldMolA pointer to the first of the two molecules that should be merged
     * @param oldMolB pointer to the second of the two molecules that should be merged
     * @param newMol pointer to the new merged molecule, passed as reference
     * @param text a description of the command
     */
    MergeMol(Molecule* oldMolA, Molecule* oldMolB, Molecule*& newMol, const QString & text = "");
    /**
     * Destructor
     *
     * Deletes the two unmerged molecules if the command is destucted in a done state
     * and the merged molecule if destructed in an undone state.
     */
    virtual ~MergeMol();
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The first of the two molecules that should be merged. */
    Molecule* m_molA;
    /** The second of the two molecules that should be merged. */
    Molecule* m_molB;
    /** The merged molecule. */
    Molecule* m_molC;
    /** The scene of this command. */
    MolScene* m_scene;
  };


/**
 * Command to split a molecule
 *
 * @author Harm van Eersel
 */
class SplitMol : public QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param molecule the molecule to be split
     * @param text a description of the command
     */
    SplitMol(Molecule* molecule, const QString & text = "");
    /**
     * Destructor
     *
     * Deletes the original molecule if destructed in a done state
     * and the submolecules if destructed in an undone state.
     */
    ~SplitMol();
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The molecule before the split. */
    Molecule* m_oldMol;
    /** The list of molecules after the split. */
    QList<Molecule*> m_newMolList;
    /** The scene of this command. */
    MolScene* m_scene;
  };


// Generic item commands

/**
 * Command to add an item to the scene
 *
 * @author Harm van Eersel
 */
class AddItem : public QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param newItem the item that should be added to the scene
     * @param addScene the scene for the new item
     * @param text a description of the command
     */
    AddItem(QGraphicsItem* newItem, MolScene* addScene, const QString & text = "");
    /**
     * Destructor
     *
     * Deletes m_item if the command is destructed in an undone state.
     */
    ~AddItem();
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The item of this command. */
    QGraphicsItem* m_item;
    /** The scene of this command. */
    MolScene* m_scene;
  };


/**
 * Command to remove an item from the scene
 *
 * @author Harm van Eersel
 */
class DelItem : public QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param delItem item to be deleted
     * @param text a description of the command
     */
    DelItem(QGraphicsItem* delItem, const QString & text = "");
    /**
     * Destructor
     *
     * Deletes m_item if the command is destructed in a done state.
     */
    ~DelItem();
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The item of this command. */
    QGraphicsItem* m_item;
    /** The scene of this command. */
    MolScene* m_scene;
  };


/**
 * Command to move an item on the scene
 *
 * @author Harm van Eersel
 */
class MoveItem : public QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param moveItem the item to be moved
     * @param moveVector the vector representation of the move
     * @param text a description of the command
     */
    MoveItem(QGraphicsItem* moveItem, const QPointF & moveVector, const QString & text = "");
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The position of the item before the move. */
    QPointF m_oldPos;
    /** The position of the item after the move. */
    QPointF m_newPos;
    /** The item of this command. */
    QGraphicsItem* m_item;
  };

/**
 * Command to rotate an item on the scene
 *
 * @author Harm van Eersel
 */
class RotateItem : public QUndoCommand
  {
  public:
    /**
     * Constructor
     *
     * @param rotateItem the item to be rotated
     * @param transform the matrix representation of the rotation
     * @param text a description of the command
     */
    RotateItem(QGraphicsItem* rotateItem, const QTransform & transform, const QString & text = "");
    /** Undo this command. */
    virtual void undo();
    /** Redo this command. */
    virtual void redo();
  private:
    /** Undo state of the command. */
    bool m_undone;
    /** The item of this command. */
    QGraphicsItem* m_item;
    /** The position of the item before the move. */
    QTransform m_transform;

  };

}

#endif
