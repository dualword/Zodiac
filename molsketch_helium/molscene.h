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
* This file is part of molsKetch and contains the MolScene class.
*
* @author Harm van Eersel <devsciurus@xs4all.nl
*/
 
 
#ifndef molscene_H
#define molscene_H

#include <QGraphicsScene>
#include <QUndoCommand>

#include "mollibitem.h"

class QString;
class QImage;
class QListWidgetItem;
class QTableWidgetItem;
class QUndoStack;
class Molecule;
class MsKAtom;
class MsKBond;


/**
 * This class is a subclassed version of QGraphicsScene with addition of some extra methods and
 * event handling for molecules. In can be used in the same fashion as QGraphicsScene.
 *
 * @author Harm van Eersel
 * @since Hydrogen
 */
class MolScene : public QGraphicsScene
{
  Q_OBJECT

public:
  // Constructor & destructor
  /** Creates a new MolScene with @p parent */
  MolScene(QObject* parent = 0);
  /** Cleans up the UndoStack and deletes the molscene. */
  ~MolScene();

  // Queries
  /** Returns the EditMode of the scene. */
  int editMode() const;
  /** Returns @c true  if neutral carbons are visible, return @c false otherwise. */
  bool carbonVisible() const;
  /** Returns @c true if hydrogens are visible, return @c false otherwise. */
  bool hydrogenVisible() const;
  /** Returns @c true if atom valency is visible, return @c false otherwise. */
  bool chargeVisible() const { return m_chargeVisible; };
  /** Returns @c true if hydrogens are automaticly added, return @c false otherwise. */
  bool autoAddHydrogen() const { return m_autoAddHydrogen; };
  /** Returns the current bond length. */
  qreal bondLength() const;
  /** Returns the current bond width */
  qreal bondWidth() const;
  /** Returns the current bond angle. */
  int bondAngle() const { return m_bondAngle; };
  /** Returns the current bond order. */
  int bondOrder() const { return m_bondOrder; };
  /** Returns the current bond type. */
  int bondType() const { return m_bondType; };
  /** Returns the current atom size. */
  qreal atomSize() const;
  /** Returns the current atomsymbol font */
  QFont atomSymbolFont() const;
  /** Returns the current dynamic grid radius. */
  int dynamicGridSize() const;

  // Commands  
  /** Renders the @p rect on the scene in a image. */
  QImage renderImage(const QRectF &rect);
  /** Sets whether neutral carbons are drawn. */
  void setCarbonVisible(bool value);
  /** Sets whether hydrogens are drawn. */
  void setHydrogenVisible(bool value);
  /** Sets whether atom charges are drawn. */
  void setChargeVisible(bool value) { m_chargeVisible = value; };
  /** Sets whether hydrogens are automaticly added. */
  void setAutoAddHydrogen(bool value) { m_autoAddHydrogen = value; };
  /** Sets the atomsymbol font */
  void setAtomSymbolFont(const QFont & font);
  /** Sets the atom size. */
  void setMsKAtomSize(qreal size);
  /** Sets the number of hint points in the dynamic grid. */
  void setHintPointSize(int size);

  /** Merges @p molA and @p molB. Used two merge two molecules when connected with a bond. */
  Molecule* merge(const Molecule* molA, const Molecule* molB);

  /** Access to the stack */
  QUndoStack * stack();

  // Advanced queries
  /** Returns the first molecule at position @p pos or NULL if none. */
  Molecule* moleculeAt(const QPointF &pos);
  /** Returns the first atom at position @p pos or NULL if none. */
  MsKAtom* atomAt(const QPointF &pos);
  /** Returns the first bond at position @p pos or NULL if none. */
  MsKBond* bondAt(const QPointF &pos);

  /** Enum for the different edit modes. */
  enum editModes {
  MoveMode, /**< MsKAtoms and molecules are movable. */
  AddMode, /**< Mode to add atoms and bonds. */
  RemoveMode, /**< Mode to remove atoms, bonds and molecules. */
  RotateMode /**< Mode to rotate molecules. */
  };

signals:
  /** Signal emitted if copy becomes available. */
  void copyAvailable(bool);
  /** Signal emitted if paste becomes available. */
  void pasteAvailable(bool);
  /** Signal emitted if the edit mode changes. */
  void editModeChange(int);
  /** Signal emitted if the contents of the scene change. */
  void documentChange( );
  /** Signal emitted if the selection on the scene changes. */
  void selectionChange( );
//  /** Signal emitted if a new molecule is added to the scene. */
//   void newMolecule(QPointF,QString);

public slots:
  /** Slot to cut the current selection to the clipboard. */
  void cut();
  /** Slot to copy the current selection to the clipboard. */
  void copy();
  /** Slot to paste the current clipboard contents. */
  void paste();
  /** Slot to clear the scene. */
  void clear();
  /** Slot to set the edit mode to @p mode. */
  void setEditMode(int mode);
  /** Slot to set the edit mode to increase charge mode */
  void setIncChargeMode();
  /** Slot to set the edit mode to decrease charge mode */
  void setDecChargeMode();
  /** Slot to set the edit mode to increase hydrogen mode */
  void setIncHydrogenMode();
  /** Slot to set the edit mode to decrease hydrogen mode */
  void setDecHydrogenMode();
  /** Slot to select all contents of the scene. */
  void selectAll();
  /** Slot to set the current element to the element of @p item. */
  void setElement(QListWidgetItem* item);
  /** Slot to set the current element to the element of @p item. */
  void setElement(QTableWidgetItem* item);
  /** Slot to set the current element to the element of @p symbol. */
  void setElement(const QString & symbol);
  /** Sets the bond angle. */
  void setBondAngle(double angle);
  /** Sets the bond length. */
  void setBondLength(double length);
  /** Sets the bond width */
  void setBondWidth(double width);
  /** Sets the bond order. */
  void setBondOrder(int order);
  /** Sets the bond type. */
  void setBondType(int type);
  /** Slot to add a copy of molecule @p mol. */
  void addMolecule(QListWidgetItem* mol);
  /** Slot to add a copy of molecule @p mol. */
  void addMolecule(Molecule* mol);
  /** Slot to align the molecules of the scene to the grid. */
  void alignToGrid();

protected:
  /** Generic event handler. Reimplementation for sceneChanged signals. */
  bool event (QEvent* event);

  // Event handlers
  /** Event handler for mouse press events. */
  virtual void mousePressEvent(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse move events. */
  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
  /** Event handler fot mouse release events. */
  virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);
  /** Event handler for double click events. */
  virtual void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event);
  /** Event handler for keyboard interaction */
  virtual void keyPressEvent(QKeyEvent* keyEvent);

  //void contextMenuEvent(Reason reason, );

private:

  // Global properties
  /** Contains the current element */
  QString m_currentElementSymbol;

  // Hinting methods and properties
  /** Shows a hintline from @p startPos to @p curPos. */
  void setHintLine(const QPointF &startPos,const QPointF &curPos);
  /** Hides the hintline. */
  void removeHintLine();
  /** Draws the dynamic grid around point @p curPos. */
  void drawHintPoints(const QPointF &curPos );
  /** Hides the dynamic grid. */
  void removeHintPoints();
  /** Shows a highlight rectangle around @p item. */
  void setHoverRect( QGraphicsItem* item );
  /** Stores the rotating angle during rotation*/
  qreal m_tempRotateAngle;

  // Auxillary methods
  /** Calculates the nearest magnetic point around @p curPos. */
  QPointF nearestPoint(const QPointF &curPos);
  /** Returns the nearest grid point, starting from @p position. */
  QPointF toGrid(const QPointF &position);
  /** Returns the nearest rounded angle, starting from @p angle. */
  int toRoundedAngle(int angle);

  // Scene properties
  /** Stores the edit mode of the scene as an integer. */
  int m_editMode;
  /** Stores the current bond length. */
  qreal m_bondLength;
  /** Stores the current bond width. */
  qreal m_bondWidth;
  /** Stores the current bond order. */
  int m_bondOrder;
  /** Stores the current bond type. */
  int m_bondType;
  /** Strores the bond angle. */
  int m_bondAngle;
  /** Stores the current atom size. */
  qreal m_atomSize;
  /** Stores the current atomsymbol font */
  QFont m_atomSymbolFont;
  /** Stores whether hydrogens should be visible.*/
  bool m_hydrogenVisible;
  /** Stores whether neutral carbons are to be shown.*/
  bool m_carbonVisible;
  /** Stores whether the charge of the atoms is to be shown.*/
  bool m_chargeVisible;
  /** Stores whether hydrogens are to be added automaticly.*/
  bool m_autoAddHydrogen;

  // Internal clipboard
  /** Internal clipboard for scene items.*/
  QList<Molecule*> m_clipItems;

  // Undo stack
  /** The undo stack of the commands used to edit the scene. */
  QUndoStack * m_stack;

  // Event handlers
  /** Event handler for mouse presses in add mode. */
  void addModePress(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse moves in add mode.*/
  void addModeMove(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse releases in add mode.*/
  void addModeRelease(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse doubleclicks in add mode.*/
  void addModeDoubleRelease(QGraphicsSceneMouseEvent* event);

  /** Event handler for mouse presses in move mode. */
  void moveModePress(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse moves in move mode.*/
  void moveModeMove(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse releases in move mode.*/
  void moveModeRelease(QGraphicsSceneMouseEvent* event);

  /** Event handler for mouse presses in rotate mode. */
  void rotateModePress(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse moves in rotate mode.*/
  void rotateModeMove(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse releases in rotate mode.*/
  void rotateModeRelease(QGraphicsSceneMouseEvent* event);

  /** Event handler for mouse presses in delete mode.*/
  void delModePress(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse moves in delete mode.*/
  void delModeMove(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse releases in delete mode.*/
  void delModeRelease(QGraphicsSceneMouseEvent* event);
  /** Event handler for mouse double clicks in delete mode.*/
  void delModeDoubleClick(QGraphicsSceneMouseEvent* event);

  // Auxillary feedback items
  /** The highlight rectangle. */
  QGraphicsPathItem* m_hoverRect;
  /** The hint line. */
  QGraphicsLineItem* m_hintLine;
  /** The list of dynamic grid points. */
  QList<QGraphicsItem*> m_hintPoints;
  /** The group of dynamic grid points. */
  QGraphicsItemGroup* m_hintPointsGroup;
  /** Method to initialize the hinting.*/
  void initHintItems();

};

#endif
