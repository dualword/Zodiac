/****************************************************************************
** Meta object code from reading C++ file 'menu.h'
**
** Created: Tue May 6 09:18:58 2008
**      by: The Qt Meta Object Compiler version 59 (Qt 4.3.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "menu.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'menu.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.3.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_SurfaceMenu[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x08,
      27,   12,   12,   12, 0x08,
      42,   12,   12,   12, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_SurfaceMenu[] = {
    "SurfaceMenu\0\0add_surface()\0draw_surface()\0"
    "update_near_to()\0"
};

const QMetaObject SurfaceMenu::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_SurfaceMenu,
      qt_meta_data_SurfaceMenu, 0 }
};

const QMetaObject *SurfaceMenu::metaObject() const
{
    return &staticMetaObject;
}

void *SurfaceMenu::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_SurfaceMenu))
	return static_cast<void*>(const_cast< SurfaceMenu*>(this));
    return QWidget::qt_metacast(_clname);
}

int SurfaceMenu::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: add_surface(); break;
        case 1: draw_surface(); break;
        case 2: update_near_to(); break;
        }
        _id -= 3;
    }
    return _id;
}
static const uint qt_meta_data_SphereMenu[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_SphereMenu[] = {
    "SphereMenu\0\0draw_sphere()\0"
};

const QMetaObject SphereMenu::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_SphereMenu,
      qt_meta_data_SphereMenu, 0 }
};

const QMetaObject *SphereMenu::metaObject() const
{
    return &staticMetaObject;
}

void *SphereMenu::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_SphereMenu))
	return static_cast<void*>(const_cast< SphereMenu*>(this));
    return QWidget::qt_metacast(_clname);
}

int SphereMenu::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: draw_sphere(); break;
        }
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_BuilderMenu[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x08,
      26,   12,   12,   12, 0x08,
      34,   12,   12,   12, 0x08,
      42,   12,   12,   12, 0x08,
      50,   12,   12,   12, 0x08,
      64,   12,   12,   12, 0x08,
      78,   12,   12,   12, 0x08,
      92,   12,   12,   12, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_BuilderMenu[] = {
    "BuilderMenu\0\0add_smiles()\0add_C()\0"
    "add_N()\0add_O()\0single_bond()\0"
    "double_bond()\0triple_bond()\0no_bond()\0"
};

const QMetaObject BuilderMenu::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_BuilderMenu,
      qt_meta_data_BuilderMenu, 0 }
};

const QMetaObject *BuilderMenu::metaObject() const
{
    return &staticMetaObject;
}

void *BuilderMenu::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_BuilderMenu))
	return static_cast<void*>(const_cast< BuilderMenu*>(this));
    return QWidget::qt_metacast(_clname);
}

int BuilderMenu::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: add_smiles(); break;
        case 1: add_C(); break;
        case 2: add_N(); break;
        case 3: add_O(); break;
        case 4: single_bond(); break;
        case 5: double_bond(); break;
        case 6: triple_bond(); break;
        case 7: no_bond(); break;
        }
        _id -= 8;
    }
    return _id;
}
static const uint qt_meta_data_HapticMenu[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,
      17,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_HapticMenu[] = {
    "HapticMenu\0\0Ok()\0end()\0"
};

const QMetaObject HapticMenu::staticMetaObject = {
    { &QTabWidget::staticMetaObject, qt_meta_stringdata_HapticMenu,
      qt_meta_data_HapticMenu, 0 }
};

const QMetaObject *HapticMenu::metaObject() const
{
    return &staticMetaObject;
}

void *HapticMenu::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_HapticMenu))
	return static_cast<void*>(const_cast< HapticMenu*>(this));
    return QTabWidget::qt_metacast(_clname);
}

int HapticMenu::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QTabWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: Ok(); break;
        case 1: end(); break;
        }
        _id -= 2;
    }
    return _id;
}
static const uint qt_meta_data_BrowserMenu[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x08,
      26,   12,   12,   12, 0x08,
      38,   12,   12,   12, 0x08,
      50,   12,   12,   12, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_BrowserMenu[] = {
    "BrowserMenu\0\0first_slot()\0prev_slot()\0"
    "next_slot()\0last_slot()\0"
};

const QMetaObject BrowserMenu::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_BrowserMenu,
      qt_meta_data_BrowserMenu, 0 }
};

const QMetaObject *BrowserMenu::metaObject() const
{
    return &staticMetaObject;
}

void *BrowserMenu::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_BrowserMenu))
	return static_cast<void*>(const_cast< BrowserMenu*>(this));
    return QWidget::qt_metacast(_clname);
}

int BrowserMenu::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: first_slot(); break;
        case 1: prev_slot(); break;
        case 2: next_slot(); break;
        case 3: last_slot(); break;
        }
        _id -= 4;
    }
    return _id;
}
static const uint qt_meta_data_Clicked_atomMenu[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      26,   18,   17,   17, 0x08,
      51,   18,   17,   17, 0x08,
      77,   18,   17,   17, 0x08,
     105,   17,   17,   17, 0x08,
     114,   17,   17,   17, 0x08,
     151,   17,   17,   17, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_Clicked_atomMenu[] = {
    "Clicked_atomMenu\0\0lab,val\0"
    "set_value(QLabel*,float)\0"
    "set_value(QLabel*,string)\0"
    "set_value(QLineEdit*,float)\0add_Hs()\0"
    "set_clicked_atom_as_center_of_view()\0"
    "set_clicked_atom_as_center_of_rotation()\0"
};

const QMetaObject Clicked_atomMenu::staticMetaObject = {
    { &QTabWidget::staticMetaObject, qt_meta_stringdata_Clicked_atomMenu,
      qt_meta_data_Clicked_atomMenu, 0 }
};

const QMetaObject *Clicked_atomMenu::metaObject() const
{
    return &staticMetaObject;
}

void *Clicked_atomMenu::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Clicked_atomMenu))
	return static_cast<void*>(const_cast< Clicked_atomMenu*>(this));
    return QTabWidget::qt_metacast(_clname);
}

int Clicked_atomMenu::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QTabWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: set_value((*reinterpret_cast< QLabel*(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2]))); break;
        case 1: set_value((*reinterpret_cast< QLabel*(*)>(_a[1])),(*reinterpret_cast< string(*)>(_a[2]))); break;
        case 2: set_value((*reinterpret_cast< QLineEdit*(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2]))); break;
        case 3: add_Hs(); break;
        case 4: set_clicked_atom_as_center_of_view(); break;
        case 5: set_clicked_atom_as_center_of_rotation(); break;
        }
        _id -= 6;
    }
    return _id;
}
static const uint qt_meta_data_ColorMenu[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      11,   10,   10,   10, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_ColorMenu[] = {
    "ColorMenu\0\0ok_slot()\0"
};

const QMetaObject ColorMenu::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_ColorMenu,
      qt_meta_data_ColorMenu, 0 }
};

const QMetaObject *ColorMenu::metaObject() const
{
    return &staticMetaObject;
}

void *ColorMenu::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ColorMenu))
	return static_cast<void*>(const_cast< ColorMenu*>(this));
    return QWidget::qt_metacast(_clname);
}

int ColorMenu::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ok_slot(); break;
        }
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_DDSettingsMenu[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      16,   15,   15,   15, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_DDSettingsMenu[] = {
    "DDSettingsMenu\0\0ok_slot()\0"
};

const QMetaObject DDSettingsMenu::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_DDSettingsMenu,
      qt_meta_data_DDSettingsMenu, 0 }
};

const QMetaObject *DDSettingsMenu::metaObject() const
{
    return &staticMetaObject;
}

void *DDSettingsMenu::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DDSettingsMenu))
	return static_cast<void*>(const_cast< DDSettingsMenu*>(this));
    return QWidget::qt_metacast(_clname);
}

int DDSettingsMenu::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ok_slot(); break;
        }
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_GraphicalObjectsMenu[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      22,   21,   21,   21, 0x0a,
      36,   21,   21,   21, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_GraphicalObjectsMenu[] = {
    "GraphicalObjectsMenu\0\0update_slot()\0"
    "delete_selected_slot()\0"
};

const QMetaObject GraphicalObjectsMenu::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_GraphicalObjectsMenu,
      qt_meta_data_GraphicalObjectsMenu, 0 }
};

const QMetaObject *GraphicalObjectsMenu::metaObject() const
{
    return &staticMetaObject;
}

void *GraphicalObjectsMenu::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GraphicalObjectsMenu))
	return static_cast<void*>(const_cast< GraphicalObjectsMenu*>(this));
    return QWidget::qt_metacast(_clname);
}

int GraphicalObjectsMenu::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: update_slot(); break;
        case 1: delete_selected_slot(); break;
        }
        _id -= 2;
    }
    return _id;
}
static const uint qt_meta_data_DatabaseGrid[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
      10,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      16,   14,   13,   13, 0x08,
      32,   14,   13,   13, 0x08,
      54,   50,   13,   13, 0x08,
      81,   79,   13,   13, 0x08,
     111,  107,   13,   13, 0x08,
     142,   50,   13,   13, 0x08,
     171,   13,   13,   13, 0x08,
     184,   13,   13,   13, 0x08,
     196,   13,   13,   13, 0x08,
     208,   13,   13,   13, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_DatabaseGrid[] = {
    "DatabaseGrid\0\0r\0select_row(int)\0"
    "deselect_row(int)\0r,c\0set_row_color(int,color)\0"
    "i\0set_current_molecule(int)\0str\0"
    "manage_number_changed(QString)\0"
    "manage_double_click(int,int)\0first_slot()\0"
    "prev_slot()\0next_slot()\0last_slot()\0"
};

const QMetaObject DatabaseGrid::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_DatabaseGrid,
      qt_meta_data_DatabaseGrid, 0 }
};

const QMetaObject *DatabaseGrid::metaObject() const
{
    return &staticMetaObject;
}

void *DatabaseGrid::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DatabaseGrid))
	return static_cast<void*>(const_cast< DatabaseGrid*>(this));
    return QWidget::qt_metacast(_clname);
}

int DatabaseGrid::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: select_row((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: deselect_row((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: set_row_color((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< color(*)>(_a[2]))); break;
        case 3: set_current_molecule((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: manage_number_changed((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: manage_double_click((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 6: first_slot(); break;
        case 7: prev_slot(); break;
        case 8: next_slot(); break;
        case 9: last_slot(); break;
        }
        _id -= 10;
    }
    return _id;
}
static const uint qt_meta_data_MySlider[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      12,   10,    9,    9, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_MySlider[] = {
    "MySlider\0\0s\0setValue(QString)\0"
};

const QMetaObject MySlider::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_MySlider,
      qt_meta_data_MySlider, 0 }
};

const QMetaObject *MySlider::metaObject() const
{
    return &staticMetaObject;
}

void *MySlider::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MySlider))
	return static_cast<void*>(const_cast< MySlider*>(this));
    return QWidget::qt_metacast(_clname);
}

int MySlider::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: setValue((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        }
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_MyFloatEditLine[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      21,   17,   16,   16, 0x0a,
      36,   16,   16,   16, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_MyFloatEditLine[] = {
    "MyFloatEditLine\0\0val\0set_value(int)\0"
    "set(QString)\0"
};

const QMetaObject MyFloatEditLine::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_MyFloatEditLine,
      qt_meta_data_MyFloatEditLine, 0 }
};

const QMetaObject *MyFloatEditLine::metaObject() const
{
    return &staticMetaObject;
}

void *MyFloatEditLine::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyFloatEditLine))
	return static_cast<void*>(const_cast< MyFloatEditLine*>(this));
    return QWidget::qt_metacast(_clname);
}

int MyFloatEditLine::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: set_value((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: set((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        }
        _id -= 2;
    }
    return _id;
}
static const uint qt_meta_data_MyIntegerEditLine[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      23,   19,   18,   18, 0x0a,
      38,   18,   18,   18, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_MyIntegerEditLine[] = {
    "MyIntegerEditLine\0\0val\0set_value(int)\0"
    "set(QString)\0"
};

const QMetaObject MyIntegerEditLine::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_MyIntegerEditLine,
      qt_meta_data_MyIntegerEditLine, 0 }
};

const QMetaObject *MyIntegerEditLine::metaObject() const
{
    return &staticMetaObject;
}

void *MyIntegerEditLine::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyIntegerEditLine))
	return static_cast<void*>(const_cast< MyIntegerEditLine*>(this));
    return QWidget::qt_metacast(_clname);
}

int MyIntegerEditLine::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: set_value((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: set((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        }
        _id -= 2;
    }
    return _id;
}
static const uint qt_meta_data_MyCompleteColorSettings[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      25,   24,   24,   24, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MyCompleteColorSettings[] = {
    "MyCompleteColorSettings\0\0setAlpha(QString)\0"
};

const QMetaObject MyCompleteColorSettings::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_MyCompleteColorSettings,
      qt_meta_data_MyCompleteColorSettings, 0 }
};

const QMetaObject *MyCompleteColorSettings::metaObject() const
{
    return &staticMetaObject;
}

void *MyCompleteColorSettings::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyCompleteColorSettings))
	return static_cast<void*>(const_cast< MyCompleteColorSettings*>(this));
    return QWidget::qt_metacast(_clname);
}

int MyCompleteColorSettings::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: setAlpha((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        }
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_MyLabelf[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_MyLabelf[] = {
    "MyLabelf\0"
};

const QMetaObject MyLabelf::staticMetaObject = {
    { &Q3HBox::staticMetaObject, qt_meta_stringdata_MyLabelf,
      qt_meta_data_MyLabelf, 0 }
};

const QMetaObject *MyLabelf::metaObject() const
{
    return &staticMetaObject;
}

void *MyLabelf::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyLabelf))
	return static_cast<void*>(const_cast< MyLabelf*>(this));
    return Q3HBox::qt_metacast(_clname);
}

int MyLabelf::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Q3HBox::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_MyColorButton[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MyColorButton[] = {
    "MyColorButton\0\0my_clicked()\0"
};

const QMetaObject MyColorButton::staticMetaObject = {
    { &QPushButton::staticMetaObject, qt_meta_stringdata_MyColorButton,
      qt_meta_data_MyColorButton, 0 }
};

const QMetaObject *MyColorButton::metaObject() const
{
    return &staticMetaObject;
}

void *MyColorButton::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyColorButton))
	return static_cast<void*>(const_cast< MyColorButton*>(this));
    return QPushButton::qt_metacast(_clname);
}

int MyColorButton::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QPushButton::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: my_clicked(); break;
        }
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_MyCheckBox[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,
      18,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MyCheckBox[] = {
    "MyCheckBox\0\0get()\0set()\0"
};

const QMetaObject MyCheckBox::staticMetaObject = {
    { &QCheckBox::staticMetaObject, qt_meta_stringdata_MyCheckBox,
      qt_meta_data_MyCheckBox, 0 }
};

const QMetaObject *MyCheckBox::metaObject() const
{
    return &staticMetaObject;
}

void *MyCheckBox::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyCheckBox))
	return static_cast<void*>(const_cast< MyCheckBox*>(this));
    return QCheckBox::qt_metacast(_clname);
}

int MyCheckBox::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QCheckBox::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: get(); break;
        case 1: set(); break;
        }
        _id -= 2;
    }
    return _id;
}