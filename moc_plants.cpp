/****************************************************************************
** Meta object code from reading C++ file 'plants.h'
**
** Created: Tue May 6 09:19:01 2008
**      by: The Qt Meta Object Compiler version 59 (Qt 4.3.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "plants.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'plants.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.3.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_Plants[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
       8,    7,    7,    7, 0x08,
      16,    7,    7,    7, 0x08,
      23,    7,    7,    7, 0x0a,
      31,    7,    7,    7, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Plants[] = {
    "Plants\0\0reset()\0save()\0about()\0load()\0"
};

const QMetaObject Plants::staticMetaObject = {
    { &Q3VBox::staticMetaObject, qt_meta_stringdata_Plants,
      qt_meta_data_Plants, 0 }
};

const QMetaObject *Plants::metaObject() const
{
    return &staticMetaObject;
}

void *Plants::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Plants))
	return static_cast<void*>(const_cast< Plants*>(this));
    return Q3VBox::qt_metacast(_clname);
}

int Plants::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Q3VBox::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: reset(); break;
        case 1: save(); break;
        case 2: about(); break;
        case 3: load(); break;
        }
        _id -= 4;
    }
    return _id;
}