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

#include <QtGlobal>
#include <QPainter>

#include "bond.h"
#include "element.h"

// Constructor

MsKBond::MsKBond(MsKAtom* atomA, MsKAtom* atomB, int order, int type, QGraphicsItem* parent,QGraphicsScene* scene) : QGraphicsItem(parent,scene)
{
  Q_CHECK_PTR(atomA);
  Q_CHECK_PTR(atomB);

  m_bondType = type;
  m_bondOrder = order;
  m_firstMsKAtom = atomA;
  m_lastMsKAtom = atomB;
  /*
    // Increasing the valency of the atoms
    m_firstMsKAtom->setValency(m_firstMsKAtom->valency() + m_bondOrder);
    m_lastMsKAtom->setValency(m_lastMsKAtom->valency() + m_bondOrder);*/

  setPos(m_firstMsKAtom->scenePos());
//   setFlag(QGraphicsItem::ItemIsSelectable);
//   setAcceptedMouseButtons(Qt::LeftButton);
}

void MsKBond::redoValency()
{
  // Check if the atoms still exist
  if (!m_firstMsKAtom || !m_lastMsKAtom) return;
  Q_CHECK_PTR(m_firstMsKAtom);
  Q_CHECK_PTR(m_lastMsKAtom);

  // Set the new number of bonds of the atoms
  m_firstMsKAtom->setNumberOfBonds(m_firstMsKAtom->numberOfBonds() + bondOrder());
  m_lastMsKAtom->setNumberOfBonds(m_lastMsKAtom->numberOfBonds() + bondOrder());

/*  // Setting the new valency of the atoms
  int e1 = molsKetch::valencyOfElement(molsKetch::symbol2number(m_firstMsKAtom->element()));
  int e2 = molsKetch::valencyOfElement(molsKetch::symbol2number(m_lastMsKAtom->element()));
  if (e1 < 0 && e2 < 0 )
    {
      m_firstMsKAtom->setValency(m_firstMsKAtom->valency() - m_bondOrder);
      m_lastMsKAtom->setValency(m_lastMsKAtom->valency() - m_bondOrder);
      return;
    }

  if (e1 > 0 && e2 > 0 )
    {
      m_firstMsKAtom->setValency(m_firstMsKAtom->valency() + m_bondOrder);
      m_lastMsKAtom->setValency(m_lastMsKAtom->valency() + m_bondOrder);
      return;
    }

  if ( e1 > e2 )
    {
      m_firstMsKAtom->setValency(m_firstMsKAtom->valency() + m_bondOrder);
      m_lastMsKAtom->setValency(m_lastMsKAtom->valency() - m_bondOrder);
      return;
    }
  if ( e1 < e2 )
    {
      m_firstMsKAtom->setValency(m_firstMsKAtom->valency() - m_bondOrder);
      m_lastMsKAtom->setValency(m_lastMsKAtom->valency() + m_bondOrder);
      return;
    }

  m_firstMsKAtom->setValency(m_firstMsKAtom->valency() - m_bondOrder);
  m_lastMsKAtom->setValency(m_lastMsKAtom->valency() - m_bondOrder);*/
}

void MsKBond::undoValency()
{
  // Check if the atoms still exist
  if (!m_firstMsKAtom || !m_lastMsKAtom) return;
  Q_CHECK_PTR(m_firstMsKAtom);
  Q_CHECK_PTR(m_lastMsKAtom);

  // Set the new number of bonds of the atoms
  m_firstMsKAtom->setNumberOfBonds(m_firstMsKAtom->numberOfBonds() - bondOrder());
  m_lastMsKAtom->setNumberOfBonds(m_lastMsKAtom->numberOfBonds() - bondOrder());

/*  // Setting the new valency of the atoms
  int e1 = -molsKetch::valencyOfElement(molsKetch::symbol2number(m_firstMsKAtom->element()));
  int e2 = -molsKetch::valencyOfElement(molsKetch::symbol2number(m_lastMsKAtom->element()));
  if (e1 < 0 && e2 < 0 )
    {
      m_firstMsKAtom->setValency(m_firstMsKAtom->valency() - m_bondOrder);
      m_lastMsKAtom->setValency(m_lastMsKAtom->valency() - m_bondOrder);
      return;
    }

  if (e1 > 0 && e2 > 0 )
    {
      m_firstMsKAtom->setValency(m_firstMsKAtom->valency() + m_bondOrder);
      m_lastMsKAtom->setValency(m_lastMsKAtom->valency() + m_bondOrder);
      return;
    }

  if ( e1 > e2 )
    {
      m_firstMsKAtom->setValency(m_firstMsKAtom->valency() + m_bondOrder);
      m_lastMsKAtom->setValency(m_lastMsKAtom->valency() - m_bondOrder);
      return;
    }
  if ( e1 < e2 )
    {
      m_firstMsKAtom->setValency(m_firstMsKAtom->valency() - m_bondOrder);
      m_lastMsKAtom->setValency(m_lastMsKAtom->valency() + m_bondOrder);
      return;
    }

  m_firstMsKAtom->setValency(m_firstMsKAtom->valency() - m_bondOrder);
  m_lastMsKAtom->setValency(m_lastMsKAtom->valency() - m_bondOrder);*/
}

// Inherited methods

QRectF MsKBond::boundingRect() const
  {
    Q_CHECK_PTR(m_firstMsKAtom);
    Q_CHECK_PTR(m_lastMsKAtom);

    qreal w = m_lastMsKAtom->x()-m_firstMsKAtom->x();
    qreal h = m_lastMsKAtom->y()-m_firstMsKAtom->y();

    // 	qreal x = qMax(m_firstMsKAtom->pos().x(),m_lastMsKAtom->pos().x());
    // 	qreal y = qMax(m_firstMsKAtom->pos().y(),m_firstMsKAtom->pos().y());
    return QRectF(mapFromParent(m_firstMsKAtom->pos()) - QPointF(5,5),QSizeF(w+10,h+10));
  }

void MsKBond::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
  // Check the scene
  MolScene* molScene = dynamic_cast<MolScene*>(scene());
  Q_CHECK_PTR(molScene);

  // 	painter->drawRect(boundingRect());
  QPointF points[4] =
    {
      mapFromParent(m_firstMsKAtom->pos()),
      shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),3).p2(),
      shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),-3).p2(),
      mapFromParent(m_firstMsKAtom->pos())
    };
  QPointF points2[4] =
    {
      mapFromParent(m_lastMsKAtom->pos()),
      shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),3).p1(),
      shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),-3).p1(),
      mapFromParent(m_lastMsKAtom->pos())
    };

  // Set painter defaults
  painter->save();
  QPen pen;
  pen.setWidthF(molScene->bondWidth());
  painter->setPen(pen);

  // Create dash pattern for dot
  QVector<qreal> dash;
  dash << 2 << 5;

  // Create a gradient for down
  QRadialGradient radialGrad(mapFromScene(m_firstMsKAtom->scenePos()), 5);
  radialGrad.setColorAt(0, Qt::white);
  radialGrad.setColorAt(0.3,Qt::white);
  radialGrad.setColorAt(0.5,Qt::black);
  radialGrad.setColorAt(0.7, Qt::white);
  radialGrad.setColorAt(1, Qt::white);
  radialGrad.setSpread(QGradient::RepeatSpread);

  switch ( m_bondType )
    {
    case MsKBond::Down:
      painter->setPen( Qt::NoPen );
      painter->setBrush( radialGrad );
      painter->drawConvexPolygon( points2, 4);
      break;
    case MsKBond::DownR:
      painter->setPen( Qt::NoPen );
      painter->setBrush( radialGrad );
      painter->drawConvexPolygon( points, 4);
      break;
    case MsKBond::Dot:
      pen.setDashPattern(dash);
      painter->setPen(pen);
      painter->drawLine(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())));
      break;
    case MsKBond::Up:
      painter->setBrush( QBrush(Qt::black) );
      painter->drawConvexPolygon( points, 4);
      break;
    case MsKBond::UpR:
      painter->setBrush( QBrush(Qt::black) );
      painter->drawConvexPolygon( points2, 4);
      break;

    default:
      switch ( m_bondOrder )
        {
        case MsKBond::Single:
          painter->drawLine(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())));
          break;
        case MsKBond::Double:
          painter->drawLine(shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),2));
          painter->drawLine(shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),-2));
          break;
        case MsKBond::Triple:
          painter->drawLine(shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),3));
          painter->drawLine(shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),-3));
          painter->drawLine(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())));
          break;
        default:
          painter->drawLine(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())));
        }
    }

  // Restore old painter
  painter->restore();
}

QVariant MsKBond::itemChange(GraphicsItemChange change, const QVariant &value)
{
  if (change == ItemPositionChange && parentItem()) parentItem()->update();
  return QGraphicsItem::itemChange(change, value);
}

QPainterPath MsKBond::shape() const
  {
    QPolygonF polygon;
    polygon << shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),10).p1()
    << shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),10).p2()
    << shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),-10).p2() << shiftVector(QLineF(mapFromParent(m_firstMsKAtom->pos()),mapFromParent(m_lastMsKAtom->pos())),-10).p1();

    QPainterPath path(mapFromParent(m_firstMsKAtom->pos()));
//     path.quadTo(QPointF(),m_lastMsKAtom->pos());
    path.addPolygon( polygon );
    path.closeSubpath();

    return path;
  }

// Manipulation methods

void MsKBond::setOrder(int order)
{
  //pre: order>0
  //post: m_bondOrder=order
  Q_ASSERT( order > 0 );
  undoValency();
  m_bondOrder = order;
  redoValency();
  update();
}

void MsKBond::incOrder()
{
  //pre: true
  //post: m_bondOrder = oldBondOrder % 3 + 1

  // Calculating the new order
  setOrder( m_bondOrder % 3 + 1 );

}

void MsKBond::decOrder()
{
  //pre: true
  //post: m_bondOrder = (oldBondOrder + 1) % 3 + 1

  // Calculating the new order
  setOrder( (m_bondOrder + 1) % 3 + 1 );

}

void MsKBond::setType(int t)
{
  //pre: 0 <= t < 6
  //post: bondType = t
  Q_ASSERT(0 <= t && t < 6);

  m_bondType = t;
  update();
}

void MsKBond::incType()
{
  //pre: true
  //post: bondType = bondType % 6 + 1
  setType((m_bondType + 1) % 6);
}

void MsKBond::decType()
{
  //pre: true
  //post: bondType = (bondType + 5) % 6
  setType( (m_bondType + 5) % 6);
}


// Query methods

int MsKBond::bondOrder() const
  {
    return m_bondOrder;
  }

int MsKBond::bondType() const
  {
    return m_bondType;
  }

MsKAtom* MsKBond::firstMsKAtom() const
  {
    return m_firstMsKAtom;
  }

MsKAtom* MsKBond::lastMsKAtom() const
  {
    return m_lastMsKAtom;
  }

bool MsKBond::hasMsKAtom(MsKAtom* atom) const
  {
    return m_firstMsKAtom == atom || m_lastMsKAtom == atom;
  }
  
Molecule* MsKBond::molecule() const
  {
    return dynamic_cast<Molecule*>(this->parentItem());
  }

// Auxilary methods

QLineF MsKBond::shiftVector(const QLineF &vector, qreal shift) // Shifts a vector on the perpendicular axis
{
  //pre: true
  //ret:shifted vector

  // Calculating the new coordinates
  qreal rx1 = vector.x1() + shift*(vector.unitVector().y2()-vector.unitVector().y1());
  qreal ry1 = vector.y1() + shift*-(vector.unitVector().x2()-vector.unitVector().x1());
  qreal rx2 = vector.x2() + shift*(vector.unitVector().y2()-vector.unitVector().y1());
  qreal ry2 = vector.y2() + shift*-(vector.unitVector().x2()-vector.unitVector().x1());

  // Returning the new vector
  return QLineF(rx1,ry1,rx2,ry2);
}
