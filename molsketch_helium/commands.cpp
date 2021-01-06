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

#include <QTransform>

#include "commands.h"

#include "atom.h"
#include "bond.h"
#include "molecule.h"

using namespace Commands;

// AddAtom

AddAtom::AddAtom(MsKAtom * newAtom, Molecule * newMol, const QString & text) : QUndoCommand(text), m_atom(newAtom), m_molecule(newMol)
{}

AddAtom::~AddAtom()
{
  if (m_undone) delete m_atom;
}

void AddAtom::undo()
{
  m_molecule->delAtom(m_atom);
  m_undone = true;
}

void AddAtom::redo()
{
  m_molecule->addMsKAtom(m_atom);
  m_atom->setFlag(QGraphicsItem::ItemIsSelectable, m_molecule->scene()->editMode()==MolScene::MoveMode);
  m_undone = false;
}

// Change element
ChangeElement::ChangeElement(MsKAtom* changeMsKAtom, const QString &newEl, const QString & text) : QUndoCommand(text), m_oldName(changeMsKAtom->element()), m_newName(newEl), m_atom(changeMsKAtom)
{}

void ChangeElement::undo()
{
  m_atom->setElement(m_oldName);
  m_undone = true;
}

void ChangeElement::redo()
{
  m_atom->setElement(m_newName);
  m_undone = false;
}

// IncCharge
IncCharge::IncCharge(MsKAtom* atom, const QString & text) : QUndoCommand(text), m_atom(atom)
{};
void IncCharge::undo()
{
  m_atom->setOxidationState(m_atom->oxidationState() - 1);
  if (m_atom->scene()) m_atom->scene()->update();
  m_undone = true;
}
void IncCharge::redo()
{
  m_atom->setOxidationState(m_atom->oxidationState() + 1);
  if (m_atom->scene()) m_atom->scene()->update();
  m_undone = false;
}

// DecCharge
DecCharge::DecCharge(MsKAtom* atom, const QString & text) : QUndoCommand(text), m_atom(atom)
{};
void DecCharge::undo()
{
  m_atom->setOxidationState(m_atom->oxidationState() + 1);
  if (m_atom->scene()) m_atom->scene()->update();
  m_undone = true;
}
void DecCharge::redo()
{
  m_atom->setOxidationState(m_atom->oxidationState() - 1);
  if (m_atom->scene()) m_atom->scene()->update();
  m_undone = false;
}

// AddImplicitHydrogen
AddImplicitHydrogen::AddImplicitHydrogen(MsKAtom* atom, const QString & text) : QUndoCommand(text), m_atom(atom)
{};
void AddImplicitHydrogen::undo()
{
  m_atom->setNumberOfImplicitHydrogens(m_atom->numberOfImplicitHydrogens() - 1);
  if (m_atom->scene()) m_atom->scene()->update();
  m_undone = true;
}
void AddImplicitHydrogen::redo()
{
  m_atom->setNumberOfImplicitHydrogens(m_atom->numberOfImplicitHydrogens() + 1);
  if (m_atom->scene()) m_atom->scene()->update();
  m_undone = false;
}

// RemoveImplicitHydrogen
RemoveImplicitHydrogen::RemoveImplicitHydrogen(MsKAtom* atom, const QString & text) : QUndoCommand(text), m_atom(atom)
{};
void RemoveImplicitHydrogen::undo()
{
  m_atom->setNumberOfImplicitHydrogens(m_atom->numberOfImplicitHydrogens() + 1);
  if (m_atom->scene()) m_atom->scene()->update();
  m_undone = true;
}
void RemoveImplicitHydrogen::redo()
{
  m_atom->setNumberOfImplicitHydrogens(m_atom->numberOfImplicitHydrogens() - 1);
  if (m_atom->scene()) m_atom->scene()->update();
  m_undone = false;
}

// DelAtom
DelAtom::DelAtom(MsKAtom* delAtom, const QString & text) : QUndoCommand(text), m_atom(delAtom), m_molecule(delAtom->molecule())
{};
DelAtom::~DelAtom()
{
  if (!m_undone)
  {
    foreach(MsKBond* bond, m_bondList) delete bond;
    delete m_atom;
  }
}
void DelAtom::undo()
{
  m_molecule->addMsKAtom(m_atom);
  m_atom->setFlag(QGraphicsItem::ItemIsSelectable, m_molecule->scene()->editMode()==MolScene::MoveMode);
  for (int i = 0; i < m_bondList.size(); i++) m_molecule->addBond(m_bondList.at(i));
  m_undone = true;
}
void DelAtom::redo()
{
  m_bondList = m_molecule->delAtom(m_atom);
  m_undone = false;
}

// MsKBond commands

AddBond::AddBond(MsKBond* newBond, const QString & text) : QUndoCommand(text), m_bond(newBond), m_mol(newBond->firstMsKAtom()->molecule())
{}

AddBond::~AddBond()
{
  if (m_undone) delete m_bond;
}

void AddBond::undo()
{
  m_mol->delBond(m_bond);
  m_undone = true;
}
void AddBond::redo()
{
  m_mol->addBond(m_bond);
  m_undone = false;
}


DelBond::DelBond(MsKBond* delBond, const QString & text) : QUndoCommand(text), m_bond(delBond), m_mol(delBond->molecule())
{}
DelBond::~DelBond()
{
  if(!m_undone) delete m_bond;
}
void DelBond::undo()
{
  m_mol->addBond(m_bond);
  m_undone = true;
}
void DelBond::redo()
{
  m_mol->delBond(m_bond);
  m_undone = false;
}

IncType::IncType(MsKBond* incBond, const QString & text) : QUndoCommand(text), m_bond(incBond)
{}
void IncType::undo()
{
  m_bond->decType();
  m_undone = true;
}
void IncType::redo()
{
  m_bond->incType();
  m_undone = false;
}

IncOrder::IncOrder(MsKBond* incBond, const QString & text) : QUndoCommand(text), m_bond(incBond)
{}
void IncOrder::undo()
{
  m_bond->decOrder();
  m_undone = true;
}
void IncOrder::redo()
{
  m_bond->incOrder();
  m_undone = false;
}

MergeMol::MergeMol(Molecule* moleculeA, Molecule* moleculeB, Molecule*& mergedMolecule, const QString & text) :  QUndoCommand(text), m_molA(moleculeA), m_molB(moleculeB), m_molC(mergedMolecule), m_scene(moleculeA->scene())
{
  //pre: molA.scene = molB.scene
  Q_ASSERT(moleculeA->scene() == moleculeB->scene());

  mergedMolecule = m_molC = m_scene->merge(m_molA,m_molB);
}
MergeMol::~MergeMol()
{
  if (!m_undone) delete m_molA;
  if (!m_undone) delete m_molB;
  if (m_undone) delete m_molC;
}
void MergeMol::undo()
{
  m_scene->removeItem(m_molC);
  m_scene->addItem(m_molA);
  m_molA->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  foreach(MsKAtom* atom, m_molC->atoms()) 
    atom->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  m_scene->addItem(m_molB);
  m_molB->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  foreach(MsKAtom* atom, m_molC->atoms()) 
    atom->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  m_undone = true;
}
void MergeMol::redo()
{
  m_scene->removeItem(m_molA);
  m_scene->removeItem(m_molB);
  m_scene->addItem(m_molC);
  m_molC->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);  
  foreach(MsKAtom* atom, m_molC->atoms()) 
    atom->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  m_undone = false;
}

SplitMol::SplitMol(Molecule* mol, const QString & text) :  QUndoCommand(text), m_oldMol(mol), m_scene(mol->scene())
{
  m_newMolList = m_oldMol->split();
}
SplitMol::~SplitMol()
{
  if (!m_undone) delete m_oldMol;
  if (m_undone) foreach(Molecule* mol,m_newMolList) delete mol;
}
void SplitMol::undo()
{
  foreach(Molecule* mol,m_newMolList) m_scene->removeItem(mol);
  m_scene->addItem(m_oldMol);
  m_oldMol->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  foreach(MsKAtom* atom, m_oldMol->atoms()) 
    atom->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  m_undone = true;
}
void SplitMol::redo()
{
  m_scene->removeItem(m_oldMol);
  foreach(Molecule* mol,m_newMolList) 
  {
    m_scene->addItem(mol);
	mol->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
	foreach(MsKAtom* atom, mol->atoms()) 
      atom->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  }
  m_undone = false;
}

// Generic item commands

AddItem::AddItem(QGraphicsItem* newItem, MolScene* addScene, const QString & text) : QUndoCommand(text), m_item(newItem), m_scene(addScene)
{}

AddItem::~AddItem()
{
  if (m_undone) delete m_item;
}

void AddItem::undo()
{
  m_scene->removeItem(m_item);
  m_undone = true;
}
void AddItem::redo()
{
  m_scene->addItem(m_item);
  m_item->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  if (m_item->type() == Molecule::Type) foreach(MsKAtom* atom, dynamic_cast<Molecule*>(m_item)->atoms()) 
    atom->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  m_scene->update();
  m_undone = false;
}

DelItem::DelItem(QGraphicsItem* delItem, const QString & text) : QUndoCommand(text), m_item(delItem)
{
  m_scene = dynamic_cast<MolScene*>(delItem->scene());
}

DelItem::~DelItem()
{
  if (!m_undone) delete m_item;
}

void DelItem::undo()
{
  m_scene->addItem(m_item);
  m_item->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  if (m_item->type() == Molecule::Type) foreach(MsKAtom* atom, dynamic_cast<Molecule*>(m_item)->atoms()) 
    atom->setFlag(QGraphicsItem::ItemIsSelectable, m_scene->editMode()==MolScene::MoveMode);
  m_scene->update();
  m_undone = true;
}
void DelItem::redo()
{
  m_scene->removeItem(m_item);
  m_undone = false;
}

MoveItem::MoveItem(QGraphicsItem* moveItem, const QPointF & moveVector, const QString & text) : QUndoCommand(text), m_oldPos(moveItem->pos()), m_newPos(moveItem->pos() + moveVector), m_item(moveItem)
{}
void MoveItem::undo()
{
  m_item -> setPos(m_oldPos);
  if (m_item->type()==MsKAtom::Type) dynamic_cast<MsKAtom*>(m_item)->molecule()->rebuild();
  m_undone = true;
}
void MoveItem::redo()
{
  m_item -> setPos(m_newPos);
  if (m_item->type()==MsKAtom::Type) dynamic_cast<MsKAtom*>(m_item)->molecule()->rebuild();
  m_undone = false;
}


RotateItem::RotateItem(QGraphicsItem* rotateItem, const QTransform & transform, const QString & text) : QUndoCommand(text), m_item(rotateItem), m_transform( transform )
{}
void RotateItem::undo()
{
  m_item -> setTransform( m_transform.inverted(), true );
  m_undone = true;
}
void RotateItem::redo()
{
  m_item -> setTransform( m_transform, true );
  m_undone = false;
}
