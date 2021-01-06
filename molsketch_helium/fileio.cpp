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
#include <QtSvg/QtSvg>

//#include <openbabel/obiter.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include "fileio.h"
#include "molecule.h"
#include "element.h"
#include "molscene.h"

QString svgLine(QLineF l, bool dashed = false, QString stroke = "#000000", qreal width = 1)
{
  // Create the output object
  QString line;
  QTextStream out(&line);

  // Code a SVG line
  out << "<line x1=\"" << l.x1() << "\" y1=\"" << l.y1() << "\" x2=\"" << l.x2() << "\" y2=\"" << l.y2() << "\" stroke=\"" << stroke << "\" stroke-width=\"" << width << "\" ";

  // Check if the line is dotted
  if (dashed) out << "stroke-dasharray = \"9, 5\" ";

  // End and return the SVG tag
  out << "/>" << endl;
  return line;
}

QString svgPolygon(QLineF l, bool dashed = false)
{
  // Create the output object
  QString polygon;
  QTextStream out(&polygon);

  // Get the needed positions
  QLineF lup = MsKBond::shiftVector(l,3);
  QLineF ldown = MsKBond::shiftVector(l,-3);

  // Code the SVG polygon
  out << "<polygon points=\"" <<  l.x1() << "," << l.y1() << " " << lup.x2() << "," << lup.y2() << " " << ldown.x2() << "," << ldown.y2() << "\" ";

  // Code the correct fill method
  if (dashed)
    out << "fill=\"url(#StripedGradient)\" stroke=\"none\" ";
  else
    out << "fill=\"black\" stroke=\"black\" ";

  // End and return the SVG tag
  out << "/>" << endl;
  return polygon;
}

QString svgDefs()
{
  // Create the output object
  QString defs;
  QTextStream out(&defs);

  // Begin tag
  out << "<defs>" << endl;

  // Code the stripe pattern
  out << "<pattern id=\"StripedPattern\" x=\"0\" y=\"0\" width=\"100\" height=\"10\" patternUnits=\"userSpaceOnUse\" >" << endl;
  out << svgLine(QLineF(0,0,100,0));
  out << "</pattern>" << endl;

  // Code the stripe gradient
  ///TODO finetuning
//     out << "<radialGradient id=\"StripedGradient\" r=\"5%\" spreadMethod=\"repeat\" >" << endl;
  out << "<linearGradient id=\"StripedGradient\" x2=\"5%\" spreadMethod=\"repeat\" >"  << endl;
  out << "<stop offset=\"0\" stop-color=\"white\" />" << endl;
  out << "<stop offset=\"0.3\" stop-color=\"white\" />" << endl;
  out << "<stop offset=\"0.5\" stop-color=\"black\" />" << endl;
  out << "<stop offset=\"0.7\" stop-color=\"white\" />" << endl;
  out << "<stop offset=\"1\" stop-color=\"white\" />" << endl;
  out << "</linearGradient>" << endl;
//     out << "</radialGradient>" << endl;

  // End tag
  out << "</defs>" << endl;

  return defs;
}
/*
bool molsKetch::saveToSVG( const QString & fileName, MolScene * scene )
{
  QSvgGenerator svgGen;
  svgGen.setFileName( fileName );
  QPainter painter(&svgGen);
    
  // Clear selection
  QList<QGraphicsItem*> selList(scene->selectedItems());
  scene->clearSelection();
  
  scene->render( &painter );
  
  // Restore selection
  foreach(QGraphicsItem* item, selList) item->setSelected(true);
  
}*/

bool molsKetch::saveToSVG( const QString & fileName, MolScene * scene )
{
  // Trying to open a file with filename
  QFile file(fileName);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return false;

  // Making the writeobject
  QTextStream out(&file);

  // Writing the header
  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  out<< "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << endl;
  out << "<svg version=\"1.1\"  width=\"" << scene->itemsBoundingRect().width() << "\" height=\"" << scene->itemsBoundingRect().height() << "\" viewBox=\"" <<  scene->itemsBoundingRect().left() << " " <<  scene->itemsBoundingRect().top() << " " << scene->itemsBoundingRect().width() << " " << scene->itemsBoundingRect().height() << "\">" << endl;

  out << svgDefs();

  // Writing the molecule
  for (int i = 0; i < scene->items().count(); i++)
    {
      if (scene->items().at(i)->type() == Molecule::Type)
        {
          Molecule* mol = (Molecule*)(scene->items().at(i));

          foreach (MsKBond* bond, mol->bonds())
            {
              // Get the bond and its attributes
//               MsKBond* bond = mol->bond(j);
              QPointF p1 = bond->firstMsKAtom()->scenePos();
              QPointF p2 = bond->lastMsKAtom()->scenePos();
              QLineF line = QLineF(p1,p2);

              // Check the type and order
              switch (bond->bondType())
                {
                case MsKBond::Up:
                  out << svgPolygon(line);
                  break;
                case MsKBond::UpR:
                  out << svgPolygon(line);
                  break;
                case MsKBond::Down:
                  out << svgPolygon(line,true);
                  break;
                case MsKBond::DownR:
                  out << svgPolygon(line,true);
                  break;
                case MsKBond::Dot:
                  out << svgLine(line,true);
                  break;
                default:
                  switch (bond->bondOrder())
                    {
                    case MsKBond::Single:
                      out << svgLine(line);
                      break;
                    case MsKBond::Double:
                      out << svgLine(MsKBond::shiftVector(line,2));
                      out << svgLine(MsKBond::shiftVector(line,-2));
                      break;
                    case MsKBond::Triple:
                      out << svgLine(MsKBond::shiftVector(line,3));
                      out << svgLine(MsKBond::shiftVector(line,-3));
                      out << svgLine(line);
                      break;
                    }
                }
            }

          foreach (MsKAtom* atom, mol->atoms())
            {
//               MsKAtom* atom = mol->atom(j);

              // Check if atom is visible and draw
              if (atom->isDrawn())
                {
                  // Mask the bonds under the atom
//                     out << "<mask>" << endl;
                  out << "<circle cx=\"" << atom->scenePos().x() << "\" cy=\"" << atom->scenePos().y() << "\" r=\"" << 13 << "\" fill=\"#FFFFFF\"/>" << endl;
//                     out << "</mask>" << endl;

                  // Draw the element
                  out << "<text x=\"" << atom->scenePos().x() << "\" y=\"" << atom->scenePos().y() + 5 << "\" font-family=\"Helvetica\" font-size=\"12\" fill=\"#000000\">";
                  out << "<tspan text-anchor=\"middle\">" << atom->element();
                  if (atom->numberOfImplicitHydrogens()>0) out << "H";
                  if (atom->numberOfImplicitHydrogens()>1) out << "<tspan dy=\"4\" font-size=\"8\" >" << atom->numberOfImplicitHydrogens() << "</tspan>"; // baseline-shift=\"sub\"
                  out << "<tspan font-size=\"8\" dy=\"-10\" >" << atom->chargeID() << "</tspan>";
                  out << "</tspan>"; //baseline-shift=\"super\"
                  out << "</text>" << endl;
                }
            }

        }
    }

  // Writing the final tag and closing the file
  out << "</svg>" << endl;
  file.close();

  return true;
}



bool molsKetch::saveFile(const QString &fileName, QGraphicsScene* scene)
{
//     QMessageBox::warning(this,tr(PROGRAM_NAME),tr("Saving is only partially implemented. You may lose data if you overwrite an existing file."),QMessageBox::Ok,QMessageBox::Ok);
  using namespace OpenBabel;
  OBConversion * conversion = new OBConversion;

  if (conversion->SetOutFormat(QFileInfo(fileName).suffix().toAscii()))
    {
      // Create the output molecule
      OBMol* obmol = new OBMol;
      QHash<MsKAtom*,OpenBabel::OBAtom*> hash;

      // Add all molecules on the scene
      foreach(QGraphicsItem* item, scene->items())
      {
        if (item->type() == Molecule::Type)
          {
            Molecule* mol = dynamic_cast<Molecule*>(item);

            hash.clear();

            obmol->BeginModify();
//                 obmol->ReserveMsKAtoms(mol->countMsKAtoms());
            foreach (MsKAtom* atom, mol->atoms())
              {
                OpenBabel::OBAtom* obatom = obmol->NewAtom();
//                 MsKAtom* atom = mol->atom(j);
                obatom->SetVector(atom->scenePos().x()/40,atom->scenePos().y()/40,0);
                std::string element = atom->element().toStdString();
//                 obatom->SetType(element);
                obatom->SetAtomicNum(molsKetch::symbol2number(atom->element()));
//                 obmol->AddAtom(*obatom);
                hash.insert(atom,obatom);
//                 cerr << hash.count() << "\n";
              }
            foreach (MsKBond* bond, mol->bonds())
              {
//                 MsKBond* bond = mol->bonds[j];
                MsKAtom* a1 = bond->firstMsKAtom();
                MsKAtom* a2 = bond->lastMsKAtom();

                OBAtom* oba1 = hash.value(a1);
                OBAtom* oba2 = hash.value(a2);

                OpenBabel::OBBond* obbond = new OpenBabel::OBBond();
//                 OBBond* obbond = obmol->NewBond();

                // Set identifier
//                 obbond->SetIdx(j);

                // Set bondorder
                obbond->SetBO(bond->bondOrder());

                // Setting bondtype
                switch (bond->bondType())
                  {
                  case MsKBond::Up:
                    obbond->SetWedge();
                  case MsKBond::UpR:
                    obbond->SetBegin(oba2);
                    obbond->SetEnd(oba1);
                    obbond->SetWedge();
                    break;
                  case MsKBond::Down:
                    obbond->SetHash();
                  case MsKBond::DownR:
                    obbond->SetBegin(oba2);
                    obbond->SetEnd(oba1);
                    obbond->SetHash();
                    break;
                  default:
                    obbond->SetBegin(oba1);
                    obbond->SetEnd(oba2);
                  }

                // Adding the bond
//                 obmol->AddBond(oba1->GetIdx(),oba2->GetIdx(),bond->getOrder());

                obmol->AddBond(*obbond);

              }
            obmol->EndModify();
          }
      }

      // Checking if the file exists and making a backup
      if (QFile::exists(fileName)) 
      {
        QFile::remove(fileName + "~");
        QFile::copy(fileName,fileName + "~");
      }

      // Writing the final result to the file
      conversion->WriteFile(obmol,fileName.toStdString());
//       qDebug << "File saved: " << fileName.toStdString() << "\n";
    }
  else
    {
      return false;
    }

  return true;
}


Molecule* molsKetch::loadFile(const QString &fileName)
{
  // Creating and setting conversion classes
  using namespace OpenBabel;
  OBConversion * conversion = new OBConversion;
  conversion->SetInFormat(conversion->FormatFromExt(fileName.toAscii()));
  OBMol obmol;

  // Try to load a file
  if (conversion->ReadFile(&obmol, fileName.toStdString()))
    {
      // Create a new molecule
      Molecule* mol = new Molecule();
      mol->setPos(QPointF(0,0));
//         cerr << "mol " << "x: " << mol->pos().x() << " y: " << mol->pos().y() << "\n";

      // Initialize normalization factor
//       qreal factor = 1;
      // molfile.GetInternalCoord(0,0,0);

      // Add atoms one-by-ons
      for (unsigned int i = 1; i <= obmol.NumAtoms();i++)
// 	FOR_ATOMS_OF_MOL(obatom,obmol)
        {
          OpenBabel::OBAtom *obatom = obmol.GetAtom(i);
          //  			scene->addRect(QRectF(atom->GetX(),atom->GetY(),5,5));
          MsKAtom* atom = mol->addMsKAtom(molsKetch::number2symbol(obatom->GetAtomicNum()),QPointF(obatom->x()*40,obatom->y()*40), false);

        }

      // Add bonds one-by-one
      /// Mind the numbering!
      for (unsigned int i = 0; i < obmol.NumBonds();i++)
// 	FOR_BONDS_OF_MOL(obbond,obmol)
        {
          // Loading the OpenBabel objects
          OpenBabel::OBBond *obbond = obmol.GetBond(i);
          OBAtom *a1 = obbond->GetBeginAtom();
          OBAtom *a2 = obbond->GetEndAtom();

          // Creating their internal counterparts
          MsKAtom* atomA = mol->atomAt(QPointF(a1->x()*40,a1->y()*40));
          MsKAtom* atomB = mol->atomAt(QPointF(a2->x()*40,a2->y()*40));
          MsKBond* bond  = mol->addBond(atomA,atomB,obbond->GetBondOrder());

          // Set special bond types
          if (obbond->IsWedge()) bond->setType( MsKBond::Up );
          if (obbond->IsHash()) bond->setType( MsKBond::Down );
//             if (obbond->IsUp()) bond->setType( MsKBond::Up );
//             if (obbond->IsDown()) bond->setType( MsKBond::Down );
          // if (obbond->IsHash()) bond->setType( MsKBond::Dot );

          // Normalizing
//             factor = scene->getBondLength()/obbond->GetLength();
        }

      // // Normalizing molecule
      // mol->scale(factor,factor);
      // mol->setMsKAtomSize(LABEL_SIZE/factor);

      return mol;
    }
  else
    {
      return 0;
    }

}

// Molecule* smiles(QString formula)
// {
//     Molecule* mol = new Molecule();
// //   QGraphicsScene scene;
// //   scene.addItem(mol);
//
//     for (int i = 0; i < formula.length();i++)
//     {
//         if (formula.at( i ).isLetter())
//             mol->addMsKAtom( QString(formula.at( i )), QPoint(0,0) );
//     }
//
//     return mol;
// }

bool molsKetch::exportFile(const QString &fileName, MolScene * scene)
{
  // Clear selection
  QList<QGraphicsItem*> selList(scene->selectedItems());
  scene->clearSelection();
  
  QImage image = scene->renderImage(scene->itemsBoundingRect());
  
  // Restore selection
  foreach(QGraphicsItem* item, selList) item->setSelected(true);
  
  return image.save(fileName);
}


bool molsKetch::printFile(QPrinter &printer, MolScene * scene)
{
  // Creating the painter
  QPainter painter(&printer);
  
  // Clear selection
  QList<QGraphicsItem*> selList(scene->selectedItems());
  scene->clearSelection();
  
  // Rendering on the printer
  QRectF rect(scene->itemsBoundingRect());
  scene->render(&painter,printer.pageRect(),rect);
  
  // Restore selection
  foreach(QGraphicsItem* item, selList) item->setSelected(true);
  
  return true;
}


/////////////////////////////////////////////////////////////////ZN CODE//////////////////////////////////////

OpenBabel::OBMol *molsKetch::toOBMolecule (QGraphicsItem *item) {
	using namespace OpenBabel;
	OpenBabel::OBMol *obmol = new OpenBabel::OBMol ();
	if (item->type() == Molecule::Type)
	{

		Molecule* mol = dynamic_cast<Molecule*>(item);
		QHash<MsKAtom*,OpenBabel::OBAtom*> hash;
		hash.clear();
		
		obmol->BeginModify();
		//                 obmol->ReserveMsKAtoms(mol->countMsKAtoms());
		foreach (MsKAtom* atom, mol->atoms())
		{
			OpenBabel::OBAtom* obatom = obmol->NewAtom();
			//                 MsKAtom* atom = mol->atom(j);
			obatom->SetVector(atom->scenePos().x()/40*1.5,-atom->scenePos().y()/40*1.5,0);
			std::string element = atom->element().toStdString();
			//                 obatom->SetType(element);
			obatom->SetAtomicNum(molsKetch::symbol2number(atom->element()));
			//                 obmol->AddAtom(*obatom);
			hash.insert(atom,obatom);
			//                 cerr << hash.count() << "\n";
		}
		foreach (MsKBond* bond, mol->bonds())
		{
			//                 MsKBond* bond = mol->bonds[j];
			MsKAtom* a1 = bond->firstMsKAtom();
			MsKAtom* a2 = bond->lastMsKAtom();
			
			OBAtom* oba1 = hash.value(a1);
			OBAtom* oba2 = hash.value(a2);
			
			OpenBabel::OBBond* obbond = new OpenBabel::OBBond();
			//                 OBBond* obbond = obmol->NewBond();
			
			// Set identifier
			//                 obbond->SetIdx(j);
			
			// Set bondorder
			obbond->SetBO(bond->bondOrder());
			
			// Setting bondtype
			switch (bond->bondType())
			{
				case MsKBond::Up:
                    obbond->SetWedge();
				case MsKBond::UpR:
                    obbond->SetBegin(oba2);
                    obbond->SetEnd(oba1);
                    obbond->SetWedge();
                    break;
				case MsKBond::Down:
                    obbond->SetHash();
				case MsKBond::DownR:
                    obbond->SetBegin(oba2);
                    obbond->SetEnd(oba1);
                    obbond->SetHash();
                    break;
				default:
                    obbond->SetBegin(oba1);
                    obbond->SetEnd(oba2);
			}
			
			// Adding the bond
			//                 obmol->AddBond(oba1->GetIdx(),oba2->GetIdx(),bond->getOrder());
			
			obmol->AddBond(*obbond);
			
		}
		obmol->EndModify();
	}
	return obmol;
}

