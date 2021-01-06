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

#define ATOM_SIZE 30

#include <QPainter>

#include "atom.h"
#include "element.h"
#include "molscene.h"
#include "molecule.h"

// Constructor

MsKAtom::MsKAtom(const QPointF &position, const QString &element, bool implicitHydrogens, QGraphicsItem* parent, QGraphicsScene* scene) : QGraphicsItem(parent,scene)
{
    //pre: position is a valid position in scene coordinates
    // Setting initial parameters
    setPos(position);
    setZValue(3);
    // setFlag(QGraphicsItem::ItemIsMovable);
//     setFlag(QGraphicsItem::ItemIgnoresTransformations);
    if (dynamic_cast<MolScene*>(scene))
        setFlag(QGraphicsItem::ItemIsSelectable, dynamic_cast<MolScene*>(scene)->editMode()==MolScene::MoveMode);

    // Enabling hovereffects
    setAcceptsHoverEvents(true);
    setAcceptedMouseButtons(Qt::LeftButton);

    // Setting private fields
    m_elementSymbol = element;
    m_hidden = true;
    m_drawn = false;

//     m_hidden = invisible;// || (m_elementSymbolement == "C" && ((MolScene*)parentItem()->scene())->getShowCarbon()) || (m_elementSymbolement == "H" && ((MolScene*)parentItem()->scene())->getShowHydrogen());
//     m_valency = 0;
    m_charge = 0; // The initial additional charge is zero
    m_oxidationState = molsKetch::oxidationStateOfElement(molsKetch::symbol2number(m_elementSymbol));
    m_valency = molsKetch::valencyOfElement(molsKetch::symbol2number(m_elementSymbol));
    m_weight = molsKetch::weightOfElement(molsKetch::symbol2number(m_elementSymbol));

    setImplicitHydrogens(implicitHydrogens);
    m_numberOfImplicitHydrogens =  m_implicitHydrogens ? molsKetch::valencyOfElement(molsKetch::symbol2number(m_elementSymbol)) : 0;
//     m_numberOfImplicitHydrogens = 0;

    m_numberOfBonds = m_numberOfImplicitHydrogens;
    m_oxidationState -= m_numberOfBonds;

    // 	QPainterPath path;
    // 	QFont font;
    // 	path.addEllipse(QRectF(position.x(),position.y(),20,20));
    // 	path.addText(position + QPointF(-5,5),font,m_elementSymbolement);
    //
    // 	setPath(path);
}


// Inherited methods

QRectF MsKAtom::boundingRect() const
{
    // 	qreal penWidth = 1;
    // 	return QRectF(-10 - penWidth / 2, -10 - penWidth / 2,
    //                        20 + penWidth / 2, 20 + penWidth / 2);
    return QRectF(-ATOM_SIZE/2, -ATOM_SIZE/2, ATOM_SIZE, ATOM_SIZE);
}

QPainterPath MsKAtom::shape() const
{
	QPainterPath path;
	path.addRect(QRectF(-10,-10,20,20));
	return path;
}

void MsKAtom::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
    // Save the original painter state
//     painter->save();

    // 	painter->drawRect(boundingRect());
    QRectF rect = QRectF(-12,-12,24,24);

    // Check the scene
    MolScene* molScene = dynamic_cast<MolScene*>(scene());
    Q_CHECK_PTR(molScene);

    // If element is m_hidden, don't draw the atoms
    if (m_hidden && !isSelected()) {
        if (m_elementSymbol == "H" && !molScene->autoAddHydrogen() && (m_oxidationState == 1 || !molScene->chargeVisible())) {
            m_drawn = false;
            return;
        }
        if (m_elementSymbol == "C" && !molScene->carbonVisible() && ((m_numberOfBonds - m_numberOfImplicitHydrogens) > 2 && charge() == 0 || !molScene->chargeVisible())) {
            m_drawn = false;
             return;
        }
    }

    m_drawn = true;

    // Drawing background
    QBrush brush(Qt::white);

    // Use a different color if selected
    if (this->isSelected()) painter->setPen(Qt::blue);

    painter->save();
    painter->setPen(Qt::NoPen);
    painter->setBrush(brush);
    painter->drawEllipse(rect);
    painter->restore();

    painter->setFont(molScene->atomSymbolFont());

    // Check to add implicit hydrogens
    int hydrogenCount = numberOfImplicitHydrogens();
    if (hydrogenCount > 0) {
        // Check if hydrogen
        if (m_elementSymbol =="H") {
            painter->drawText(rect, Qt::AlignCenter, "H");
	    hydrogenCount += 1;
	} else
            // Draw element + hydrogen
            painter->drawText(rect, Qt::AlignCenter, m_elementSymbol + "H");
	
        if (hydrogenCount > 1) {
            painter->setFont( QFont("",7) );
            painter->drawText( rect.translated(3,0), Qt::AlignBottom|Qt::AlignRight, QString::number(hydrogenCount));
        }
    } else {
        // Draw element
        painter->drawText(rect, Qt::AlignCenter, m_elementSymbol);
    }

    // Draw charge
    if (molScene->chargeVisible()) {
        QString chargeId = chargeID();
        painter->setFont( QFont("",7) );
        painter->drawText(rect.translated(5,0), Qt::AlignTop|Qt::AlignRight, chargeId);
    }


    // Draw unbound electrons
    if (0) /*molScene->chargeVisible()*/ {
        int unboundElectrons = 8;
        QList<QRectF> layoutList;

        // Loading different layouts
        layoutList << QRectF(-3,-10,2,2);
        layoutList << QRectF(3,-10,2,2);
        layoutList << QRectF(-3,10,2,2);
        layoutList << QRectF(3,10,2,2);
        layoutList << QRectF(10,-3,2,2);
        layoutList << QRectF(10,3,2,2);
        layoutList << QRectF(-10,-3,2,2);
        layoutList << QRectF(-10,3,2,2);

        painter->save();
        painter->setBrush(Qt::black);

        for (int i = 0; i < unboundElectrons; i++)
            painter->drawEllipse(layoutList[i]);

        painter->restore();
    }

    // Restore the original painter state
//     painter->restore();
}

QVariant MsKAtom::itemChange(GraphicsItemChange change, const QVariant &value)
{
    if (change == ItemPositionChange && parentItem())
    {
//         setTransform(parentItem()->transform().transposed());
        parentItem()->update();
        dynamic_cast<Molecule*>(parentItem())->rebuild();

// 	    setGroup(dynamic_cast<Molecule*>(parentItem()));
    };
    if (change == ItemSelectedChange && molecule())
    {
//       molecule()->setSm_elementSymbolected(isSm_elementSymbolected());
      molecule()->setFlag(ItemIsSelectable,isSelected());
    }
    
    return QGraphicsItem::itemChange(change, value);
}

// Commands

void MsKAtom::setElement(const QString &element)
{
    m_elementSymbol = element;
    update();
}

void MsKAtom::setValency( int valency )
{
    // Adjusting the charge
//     m_charge = m_charge - (valency - m_valency);

    // Setting the new valency
    m_valency = valency;
}

void MsKAtom::setOxidationState(int state)
{
  m_oxidationState = state;
}

void MsKAtom::setNumberOfBonds(int number)
{
  Q_ASSERT (number >= 0);

  int deltaNoB = number - m_numberOfBonds;
  m_oxidationState -= deltaNoB;

  if (m_implicitHydrogens) {
    int oldNoIH = m_numberOfImplicitHydrogens;
    int newNoIH = oldNoIH - deltaNoB;
    int deltaNoIH = (newNoIH - oldNoIH > -oldNoIH) ? newNoIH - oldNoIH : -oldNoIH;
    m_numberOfImplicitHydrogens += deltaNoIH;
    m_oxidationState -= deltaNoIH;
  }

  m_numberOfBonds = number;
}

void MsKAtom::setNumberOfImplicitHydrogens(int number)
{
  Q_ASSERT (number >= 0);

  m_implicitHydrogens = true;

  int deltaNoIH = number - m_numberOfImplicitHydrogens;
  int newNoB = m_numberOfBonds - deltaNoIH;
  m_numberOfBonds = (newNoB < 0) ? 0 : newNoB;

  m_oxidationState += deltaNoIH;

  m_numberOfImplicitHydrogens = number;
}


// Query methods

int MsKAtom::numberOfBonds() const
{
  return m_numberOfBonds;
}

int MsKAtom::numberOfImplicitHydrogens() const
{
  return m_numberOfImplicitHydrogens;
}

QString MsKAtom::element() const
{
    return m_elementSymbol;
}

int MsKAtom::oxidationState() const
{
  return m_oxidationState;
}

int MsKAtom::valency( )  const // return the number of expected bonds
{
    return m_valency;
}

int MsKAtom::charge( )  const // Return the additional charge of the atom
{
    return m_oxidationState;
}

// int MsKAtom::hydrogenNeeded( )  const // Returs the number of hydrogens needed to make the atom neutral
// {
//     int needed = molsKetch::valencyOfElement(molsKetch::symbol2number(m_elementSymbol)) - m_valency;
//     if (molsKetch::symbol2number(m_elementSymbol) > 10) needed = 0;
// //     if (needed > -molsKetch::valencyOfElement(molsKetch::symbol2number(m_elementSymbol))) needed = -molsKetch::valencyOfElement(molsKetch::symbol2number(m_elementSymbol));
//     if (needed < 0) needed = 0;
//     if (scene()) return static_cast<MolScene*>(scene())->autoAddHydrogen()?needed:0;
//       return 0;
// }

QString MsKAtom::chargeID( ) const
{
    // Get the charge
    int c = charge();

    // Drawing text
    QString chargeId;
    chargeId.setNum(c);
    if (c < -1) chargeId =  chargeId.remove(0,1) + "-";
    if (c == -1) chargeId = "-";
    if (c == 0) chargeId = "";
    if (c == 1) chargeId = "+";
    if (c > 1) chargeId = chargeId + "+";

    // Return the charge ID string
    return chargeId;
}

qreal MsKAtom::weight( ) const
{
    return m_weight + m_numberOfImplicitHydrogens*molsKetch::weightOfElement(1);
}

Molecule * MsKAtom::molecule() const
{
    return dynamic_cast<Molecule*>(this->parentItem());
}

// Event handlers

void MsKAtom::mousePressEvent( QGraphicsSceneMouseEvent* event )
{
    // Execute default behavior
    QGraphicsItem::mousePressEvent( event );
}

void MsKAtom::hoverEnterEvent( QGraphicsSceneHoverEvent * event )
{
//   show();
    m_hidden = false;
    // Execute default behavior
    QGraphicsItem::hoverEnterEvent( event );
}

void MsKAtom::hoverLeaveEvent( QGraphicsSceneHoverEvent * event )
{
    m_hidden = true;
//   if (m_hidden) hide();
    // Execute default behavior
    QGraphicsItem::hoverLeaveEvent( event );
}

bool MsKAtom::implicitHydrogens() const
{
  return m_implicitHydrogens;
}

bool MsKAtom::isDrawn() const
{
  return m_drawn;
}

bool MsKAtom::isHidden() const
{
  return m_hidden;
}

void MsKAtom::setImplicitHydrogens(bool enabled)
{
  m_implicitHydrogens = enabled && (m_elementSymbol == "C" || m_elementSymbol == "N" || m_elementSymbol == "O");
}

