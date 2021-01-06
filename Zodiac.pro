######################################################################
# Automatically generated by qmake (2.01a) Tue 30. Sep 18:26:14 2008
######################################################################

TEMPLATE = app
TARGET = Zodiac
DEPENDPATH += . \
              debug \
              forcefields \
              headers \
              molsketch_helium \
              molsketch_helium\i18n \
              molsketch_helium\part
INCLUDEPATH += . headers molsketch_helium 
CONFIG += qt assistant
QT += qt3support opengl svg
unix:LIBS += -L/usr/local/lib -L/usr/lib -lopenbabel
win32:LIBS += -L"C:\Program Files (x86)\openbabel\lib" 
win32:LIBS += -L"C:\Program Files\openbabel\lib" -lopenbabel-2
win32 {
INCLUDEPATH += "C:\Program Files (x86)\openbabel\include\openbabel-2.0"
INCLUDEPATH += "C:\Program Files\openbabel\include\openbabel-2.0"
}
unix: INCLUDEPATH += /usr/local/include/openbabel-2.0 /usr/include/openbabel-2.0
unix: RC_FILE = icons/zeden.icns
win32: RC_FILE = win_rc.rc

unix:exe.path = /usr/local/bin
unix:exe.files = Zodiac

unix:INSTALLS +=exe

# Input
HEADERS += obabel_includes.h \
           forcefields/forcefieldghemical.h \
           forcefields/forcefieldmmff94.h \
           forcefields/forcefielduff.h \
           headers/actions.h \
           headers/arcball.h \
           headers/builder.h \
           headers/chemscore.h \
           headers/cmat.h \
           headers/command.h \
           headers/constants.h \
           headers/cutoffGrid.h \
           headers/database.h \
           headers/datagrid.h \
           headers/ddwin.h \
           headers/FF.h \
           headers/function.h \
           headers/graphical_object.h \
           headers/ils.h \
           headers/iodevice.h \
           headers/LookUpTable.h \
           headers/MarchingCubes.h \
           headers/maths.h \
           headers/menu.h \
           headers/minimize.h \
           headers/MMFF.h \
           headers/myline.h \
           headers/nms.h \
           headers/optimiser.h \
           headers/plants.h \
           headers/plants_mainwin.h \
           headers/plantsconfig.h \
           headers/PLP.h \
           headers/ply.h \
           headers/pso.h \
           headers/subscrpt.h \
           headers/thread.h \
           headers/vec.h \
           headers/wiimote.h \
           headers/ZNdata.h \
           headers/ZNmolecule.h \
           molsketch_helium/atom.h \
           molsketch_helium/bond.h \
           molsketch_helium/commands.h \
           molsketch_helium/element.h \
           molsketch_helium/fileio.h \
           molsketch_helium/mainwindow.h \
           molsketch_helium/molecule.h \
           molsketch_helium/mollibitem.h \
           molsketch_helium/molscene.h \
           molsketch_helium/molview.h \
           molsketch_helium/periodictablewidget.h \
           molsketch_helium/settings.h \
           molsketch_helium/part/molsketch_factory.h \
           molsketch_helium/part/molsketchpart.h \
           molsketch_helium/part/molsketchpart_shell.h
FORMS += zodiac.ui molsketch_helium/settings.ui
SOURCES += actions.cc \
           arcball.cc \
           builder.cc \
           chemscore.cc \
           command.cc \
           constants.cc \
           database.cc \
           datagrid.cc \
           ddwin.cc \
           FF.cc \
           graphical_object.cc \
           iodevice.cc \
           MarchingCubes.cc \
           maths.cc \
           menu.cc \
           minimize.cc \
           MMFF.cc \
           myline.cc \
           plants.cc \
           plants_mainwin.cc \
           plantsconfig.cc \
           PLP.cc \
           ply.c \
           povray.cc \
           thread.cc \
           wiimote.cc \
           ZNdata.cc \
           ZNmolecule.cc \
           Zodiac.cc \
           forcefields/forcefieldghemical.cpp \
           forcefields/forcefieldmmff94.cpp \
           forcefields/forcefielduff.cpp \
           molsketch_helium/atom.cpp \
           molsketch_helium/bond.cpp \
           molsketch_helium/commands.cpp \
           molsketch_helium/element.cpp \
           molsketch_helium/fileio.cpp \
           molsketch_helium/mainwindow.cpp \
           molsketch_helium/molecule.cpp \
           molsketch_helium/mollibitem.cpp \
           molsketch_helium/molscene.cpp \
           molsketch_helium/molview.cpp \
           molsketch_helium/periodictablewidget.cpp \
           molsketch_helium/settings.cpp \
           molsketch_helium/part/molsketch_factory.cpp \
           molsketch_helium/part/molsketchpart.cpp \
           molsketch_helium/part/molsketchpart_shell.cpp
RESOURCES += resources.qrc \
             zodiac.qrc \
             molsketch_helium/molsketch.qrc \
             molsketch_helium/part/molsketchpart_shell.qrc
TRANSLATIONS += molsketch_helium/i18n/molsketch_cs.ts \
                molsketch_helium/i18n/molsketch_nl.ts \
                molsketch_helium/i18n/molsketch_pt_BR.ts
