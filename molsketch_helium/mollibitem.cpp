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

#include <QtGui>

#include "molecule.h"
#include "mollibitem.h"
#include "fileio.h"

MolLibItem::MolLibItem( Molecule* molecule, const QString & name )
{
  //pre: molecule is a valid molecule
  Q_CHECK_PTR(molecule);
//   Q_ASSERT(!molecule->atoms().isEmpty());

  // Copying the molecule
  m_molecule = new Molecule(molecule);
  m_molecule->setPos(0,0);

  // Creating pixmap
  MolScene renderScene;
  renderScene.addItem(m_molecule);

  QBitmap bitmap(renderScene.sceneRect().width(),renderScene.sceneRect().height());
  bitmap.clear();

  // Creating and setting the painter
  QPainter painter(&bitmap);
  painter.setRenderHint(QPainter::Antialiasing);
  renderScene.render(&painter);

  setIcon(QIcon(bitmap));

  // Checking dir
  QDir dir;
  if (!dir.exists(QDir::homePath() + "/.molsketch/library/custom/"))
    {
      dir.mkpath(QDir::homePath() + "/.molsketch/library/custom/");
    }

  // Setting name
  m_fileName = name;
  if (!m_fileName.exists())
    {
      m_fileName.setFile(QDir::homePath() + "/.molsketch/library/custom/" + name + ".mol");
      molsKetch::saveFile(m_fileName.filePath(),&renderScene);
    }

  setText(m_fileName.baseName());

  // Remove the m_molecule before destroying the scene
  renderScene.removeItem(m_molecule);

}

MolLibItem::~ MolLibItem( )
{
  // Delete all bonds and atoms and finally the molecule
  foreach(MsKBond* bond, m_molecule->bonds()) delete bond;
  foreach(MsKAtom* atom, m_molecule->atoms()) delete atom;
  delete m_molecule;
}

Molecule* MolLibItem::getMolecule( )
{
  // Return a copy of the m_molecule
//     Molecule* mol = new Molecule(m_molecule);
  return new Molecule(m_molecule);
}

void MolLibItem::setMolecule( Molecule* mol )
{
  // Delete the old molecule
  foreach(MsKBond* bond, m_molecule->bonds()) delete bond;
  foreach(MsKAtom* atom, m_molecule->atoms()) delete atom;
  delete m_molecule;
  
  // Set a copy of mol as the new molecule
  m_molecule = new Molecule(mol);
}

