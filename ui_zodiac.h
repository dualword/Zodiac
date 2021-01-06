/********************************************************************************
** Form generated from reading ui file 'zodiac.ui'
**
** Created: Tue May 6 09:18:05 2008
**      by: Qt User Interface Compiler version 4.3.4
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_ZODIAC_H
#define UI_ZODIAC_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>

class Ui_zodiacClass
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *zodiacClass)
    {
    if (zodiacClass->objectName().isEmpty())
        zodiacClass->setObjectName(QString::fromUtf8("zodiacClass"));
    zodiacClass->resize(600, 400);
    menuBar = new QMenuBar(zodiacClass);
    menuBar->setObjectName(QString::fromUtf8("menuBar"));
    zodiacClass->setMenuBar(menuBar);
    mainToolBar = new QToolBar(zodiacClass);
    mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
    zodiacClass->addToolBar(mainToolBar);
    centralWidget = new QWidget(zodiacClass);
    centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
    zodiacClass->setCentralWidget(centralWidget);
    statusBar = new QStatusBar(zodiacClass);
    statusBar->setObjectName(QString::fromUtf8("statusBar"));
    zodiacClass->setStatusBar(statusBar);

    retranslateUi(zodiacClass);

    QMetaObject::connectSlotsByName(zodiacClass);
    } // setupUi

    void retranslateUi(QMainWindow *zodiacClass)
    {
    zodiacClass->setWindowTitle(QApplication::translate("zodiacClass", "zodiac", 0, QApplication::UnicodeUTF8));
    Q_UNUSED(zodiacClass);
    } // retranslateUi

};

namespace Ui {
    class zodiacClass: public Ui_zodiacClass {};
} // namespace Ui

#endif // UI_ZODIAC_H
