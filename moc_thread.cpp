/****************************************************************************
** Meta object code from reading C++ file 'thread.h'
**
** Created: Tue May 6 09:19:05 2008
**      by: The Qt Meta Object Compiler version 59 (Qt 4.3.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "thread.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'thread.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.3.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_Thread[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_Thread[] = {
    "Thread\0"
};

const QMetaObject Thread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_Thread,
      qt_meta_data_Thread, 0 }
};

const QMetaObject *Thread::metaObject() const
{
    return &staticMetaObject;
}

void *Thread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Thread))
	return static_cast<void*>(const_cast< Thread*>(this));
    return QThread::qt_metacast(_clname);
}

int Thread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_SurfaceThread[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_SurfaceThread[] = {
    "SurfaceThread\0"
};

const QMetaObject SurfaceThread::staticMetaObject = {
    { &Thread::staticMetaObject, qt_meta_stringdata_SurfaceThread,
      qt_meta_data_SurfaceThread, 0 }
};

const QMetaObject *SurfaceThread::metaObject() const
{
    return &staticMetaObject;
}

void *SurfaceThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_SurfaceThread))
	return static_cast<void*>(const_cast< SurfaceThread*>(this));
    return Thread::qt_metacast(_clname);
}

int SurfaceThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Thread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_MinimiseThread[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      20,   16,   15,   15, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_MinimiseThread[] = {
    "MinimiseThread\0\0mol\0ask_redraw(Molecule*)\0"
};

const QMetaObject MinimiseThread::staticMetaObject = {
    { &Thread::staticMetaObject, qt_meta_stringdata_MinimiseThread,
      qt_meta_data_MinimiseThread, 0 }
};

const QMetaObject *MinimiseThread::metaObject() const
{
    return &staticMetaObject;
}

void *MinimiseThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MinimiseThread))
	return static_cast<void*>(const_cast< MinimiseThread*>(this));
    return Thread::qt_metacast(_clname);
}

int MinimiseThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Thread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ask_redraw((*reinterpret_cast< Molecule*(*)>(_a[1]))); break;
        }
        _id -= 1;
    }
    return _id;
}

// SIGNAL 0
void MinimiseThread::ask_redraw(Molecule * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_DatabaseMinimiseThread[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_DatabaseMinimiseThread[] = {
    "DatabaseMinimiseThread\0"
};

const QMetaObject DatabaseMinimiseThread::staticMetaObject = {
    { &Thread::staticMetaObject, qt_meta_stringdata_DatabaseMinimiseThread,
      qt_meta_data_DatabaseMinimiseThread, 0 }
};

const QMetaObject *DatabaseMinimiseThread::metaObject() const
{
    return &staticMetaObject;
}

void *DatabaseMinimiseThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DatabaseMinimiseThread))
	return static_cast<void*>(const_cast< DatabaseMinimiseThread*>(this));
    return Thread::qt_metacast(_clname);
}

int DatabaseMinimiseThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Thread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_HapticThread[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      18,   14,   13,   13, 0x05,
      40,   14,   13,   13, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_HapticThread[] = {
    "HapticThread\0\0mol\0ask_redraw(Molecule*)\0"
    "ask_color_by_score(Molecule*)\0"
};

const QMetaObject HapticThread::staticMetaObject = {
    { &Thread::staticMetaObject, qt_meta_stringdata_HapticThread,
      qt_meta_data_HapticThread, 0 }
};

const QMetaObject *HapticThread::metaObject() const
{
    return &staticMetaObject;
}

void *HapticThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_HapticThread))
	return static_cast<void*>(const_cast< HapticThread*>(this));
    return Thread::qt_metacast(_clname);
}

int HapticThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Thread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ask_redraw((*reinterpret_cast< Molecule*(*)>(_a[1]))); break;
        case 1: ask_color_by_score((*reinterpret_cast< Molecule*(*)>(_a[1]))); break;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void HapticThread::ask_redraw(Molecule * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void HapticThread::ask_color_by_score(Molecule * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
static const uint qt_meta_data_GridThread[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_GridThread[] = {
    "GridThread\0"
};

const QMetaObject GridThread::staticMetaObject = {
    { &Thread::staticMetaObject, qt_meta_stringdata_GridThread,
      qt_meta_data_GridThread, 0 }
};

const QMetaObject *GridThread::metaObject() const
{
    return &staticMetaObject;
}

void *GridThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GridThread))
	return static_cast<void*>(const_cast< GridThread*>(this));
    return Thread::qt_metacast(_clname);
}

int GridThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Thread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_HeadTrackingThread[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      24,   20,   19,   19, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_HeadTrackingThread[] = {
    "HeadTrackingThread\0\0x,y\0head_moved(int,int)\0"
};

const QMetaObject HeadTrackingThread::staticMetaObject = {
    { &Thread::staticMetaObject, qt_meta_stringdata_HeadTrackingThread,
      qt_meta_data_HeadTrackingThread, 0 }
};

const QMetaObject *HeadTrackingThread::metaObject() const
{
    return &staticMetaObject;
}

void *HeadTrackingThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_HeadTrackingThread))
	return static_cast<void*>(const_cast< HeadTrackingThread*>(this));
    return Thread::qt_metacast(_clname);
}

int HeadTrackingThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Thread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: head_moved((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        }
        _id -= 1;
    }
    return _id;
}

// SIGNAL 0
void HeadTrackingThread::head_moved(int _t1, int _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_WiimoteTrackingThread[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      29,   23,   22,   22, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_WiimoteTrackingThread[] = {
    "WiimoteTrackingThread\0\0x,y,z\0"
    "move_camera(float,float,float)\0"
};

const QMetaObject WiimoteTrackingThread::staticMetaObject = {
    { &Thread::staticMetaObject, qt_meta_stringdata_WiimoteTrackingThread,
      qt_meta_data_WiimoteTrackingThread, 0 }
};

const QMetaObject *WiimoteTrackingThread::metaObject() const
{
    return &staticMetaObject;
}

void *WiimoteTrackingThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_WiimoteTrackingThread))
	return static_cast<void*>(const_cast< WiimoteTrackingThread*>(this));
    return Thread::qt_metacast(_clname);
}

int WiimoteTrackingThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = Thread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: move_camera((*reinterpret_cast< float(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2])),(*reinterpret_cast< float(*)>(_a[3]))); break;
        }
        _id -= 1;
    }
    return _id;
}

// SIGNAL 0
void WiimoteTrackingThread::move_camera(float _t1, float _t2, float _t3)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)), const_cast<void*>(reinterpret_cast<const void*>(&_t3)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}