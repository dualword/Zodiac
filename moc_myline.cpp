/****************************************************************************
** Meta object code from reading C++ file 'myline.h'
**
** Created: Tue May 6 09:19:00 2008
**      by: The Qt Meta Object Compiler version 59 (Qt 4.3.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "myline.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'myline.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.3.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_MyLine[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_MyLine[] = {
    "MyLine\0"
};

const QMetaObject MyLine::staticMetaObject = {
    { &Q3HBox::staticMetaObject, qt_meta_stringdata_MyLine,
      qt_meta_data_MyLine, 0 }
};

const QMetaObject *MyLine::metaObject() const
{
    return &staticMetaObject;
}

void *MyLine::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyLine))
	return static_cast<void*>(const_cast< MyLine*>(this));
    return Q3HBox::qt_metacast(_clname);
}

int MyLine::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Q3HBox::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_MyLineF[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      17,    8,    9,    8, 0x09,
      28,    8,    8,    8, 0x09,
      39,    8,    8,    8, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_MyLineF[] = {
    "MyLineF\0\0QString\0ask_file()\0set_file()\0"
    "set_line()\0"
};

const QMetaObject MyLineF::staticMetaObject = {
    { &Q3HBox::staticMetaObject, qt_meta_stringdata_MyLineF,
      qt_meta_data_MyLineF, 0 }
};

const QMetaObject *MyLineF::metaObject() const
{
    return &staticMetaObject;
}

void *MyLineF::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyLineF))
	return static_cast<void*>(const_cast< MyLineF*>(this));
    return Q3HBox::qt_metacast(_clname);
}

int MyLineF::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Q3HBox::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: { QString _r = ask_file();
            if (_a[0]) *reinterpret_cast< QString*>(_a[0]) = _r; }  break;
        case 1: set_file(); break;
        case 2: set_line(); break;
        }
        _id -= 3;
    }
    return _id;
}
