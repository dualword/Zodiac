/****************************************************************************
** Meta object code from reading C++ file 'ddwin.h'
**
** Created: Tue May 6 09:18:56 2008
**      by: The Qt Meta Object Compiler version 59 (Qt 4.3.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "ddwin.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ddwin.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.3.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_DDWin[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
      48,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
       7,    6,    6,    6, 0x05,

 // slots: signature, parameters, type, tag, flags
      25,    6,    6,    6, 0x08,
      52,    6,    6,    6, 0x08,
      69,    6,    6,    6, 0x08,
      84,    6,    6,    6, 0x08,
     102,    6,    6,    6, 0x08,
     115,    6,    6,    6, 0x08,
     138,    6,    6,    6, 0x08,
     153,    6,    6,    6, 0x08,
     168,    6,    6,    6, 0x08,
     183,    6,    6,    6, 0x08,
     199,    6,    6,    6, 0x08,
     221,    6,    6,    6, 0x08,
     252,    6,    6,    6, 0x08,
     274,    6,    6,    6, 0x08,
     296,    6,    6,    6, 0x08,
     320,    6,    6,    6, 0x08,
     339,    6,    6,    6, 0x08,
     352,    6,    6,    6, 0x08,
     367,    6,    6,    6, 0x08,
     381,    6,    6,    6, 0x08,
     406,    6,    6,    6, 0x08,
     420,    6,    6,    6, 0x08,
     430,    6,    6,    6, 0x08,
     445,    6,    6,    6, 0x08,
     458,    6,    6,    6, 0x08,
     481,    6,    6,    6, 0x08,
     508,    6,    6,    6, 0x08,
     530,    6,    6,    6, 0x08,
     542,    6,    6,    6, 0x08,
     565,    6,    6,    6, 0x08,
     581,    6,    6,    6, 0x08,
     595,    6,    6,    6, 0x08,
     615,  609,    6,    6, 0x08,
     644,    6,    6,    6, 0x08,
     658,    6,    6,    6, 0x08,
     672,    6,    6,    6, 0x08,
     684,    6,    6,    6, 0x08,
     696,    6,    6,    6, 0x08,
     708,    6,    6,    6, 0x08,
     720,    6,    6,    6, 0x08,
     732,    6,    6,    6, 0x08,
     744,    6,    6,    6, 0x08,
     756,    6,    6,    6, 0x0a,
     788,  784,    6,    6, 0x0a,
     806,  784,    6,    6, 0x0a,
     834,    6,    6,    6, 0x0a,
     853,    6,    6,    6, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_DDWin[] = {
    "DDWin\0\0targets_updated()\0"
    "resizeEvent(QResizeEvent*)\0open_file_slot()\0"
    "save_as_slot()\0screenshot_slot()\0"
    "movie_slot()\0raytraced_movie_slot()\0"
    "builder_slot()\0history_slot()\0"
    "wiimote_slot()\0wiimote2_slot()\0"
    "hide_hydrogens_slot()\0"
    "hide_nonpolar_hydrogens_slot()\0"
    "show_all_atoms_slot()\0hide_all_atoms_slot()\0"
    "display_settings_slot()\0DD_settings_slot()\0"
    "color_slot()\0surface_slot()\0sphere_slot()\0"
    "graphical_objects_slot()\0add_Hs_slot()\0"
    "disp_ok()\0atdebug_slot()\0about_slot()\0"
    "partial_charges_slot()\0"
    "scores_from_charges_slot()\0"
    "compute_energy_slot()\0logP_slot()\0"
    "minimise_energy_slot()\0not_impl_slot()\0"
    "haptic_slot()\0plants_slot()\0index\0"
    "set_current_target_slot(int)\0del_pressed()\0"
    "esc_pressed()\0a_pressed()\0b_pressed()\0"
    "c_pressed()\0m_pressed()\0n_pressed()\0"
    "o_pressed()\0s_pressed()\0"
    "raytraced_screenshot_slot()\0mol\0"
    "redraw(Molecule*)\0recolor_by_score(Molecule*)\0"
    "end_minimisation()\0emit_targets_updated()\0"
};

const QMetaObject DDWin::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_DDWin,
      qt_meta_data_DDWin, 0 }
};

const QMetaObject *DDWin::metaObject() const
{
    return &staticMetaObject;
}

void *DDWin::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DDWin))
	return static_cast<void*>(const_cast< DDWin*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int DDWin::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: targets_updated(); break;
        case 1: resizeEvent((*reinterpret_cast< QResizeEvent*(*)>(_a[1]))); break;
        case 2: open_file_slot(); break;
        case 3: save_as_slot(); break;
        case 4: screenshot_slot(); break;
        case 5: movie_slot(); break;
        case 6: raytraced_movie_slot(); break;
        case 7: builder_slot(); break;
        case 8: history_slot(); break;
        case 9: wiimote_slot(); break;
        case 10: wiimote2_slot(); break;
        case 11: hide_hydrogens_slot(); break;
        case 12: hide_nonpolar_hydrogens_slot(); break;
        case 13: show_all_atoms_slot(); break;
        case 14: hide_all_atoms_slot(); break;
        case 15: display_settings_slot(); break;
        case 16: DD_settings_slot(); break;
        case 17: color_slot(); break;
        case 18: surface_slot(); break;
        case 19: sphere_slot(); break;
        case 20: graphical_objects_slot(); break;
        case 21: add_Hs_slot(); break;
        case 22: disp_ok(); break;
        case 23: atdebug_slot(); break;
        case 24: about_slot(); break;
        case 25: partial_charges_slot(); break;
        case 26: scores_from_charges_slot(); break;
        case 27: compute_energy_slot(); break;
        case 28: logP_slot(); break;
        case 29: minimise_energy_slot(); break;
        case 30: not_impl_slot(); break;
        case 31: haptic_slot(); break;
        case 32: plants_slot(); break;
        case 33: set_current_target_slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 34: del_pressed(); break;
        case 35: esc_pressed(); break;
        case 36: a_pressed(); break;
        case 37: b_pressed(); break;
        case 38: c_pressed(); break;
        case 39: m_pressed(); break;
        case 40: n_pressed(); break;
        case 41: o_pressed(); break;
        case 42: s_pressed(); break;
        case 43: raytraced_screenshot_slot(); break;
        case 44: redraw((*reinterpret_cast< Molecule*(*)>(_a[1]))); break;
        case 45: recolor_by_score((*reinterpret_cast< Molecule*(*)>(_a[1]))); break;
        case 46: end_minimisation(); break;
        case 47: emit_targets_updated(); break;
        }
        _id -= 48;
    }
    return _id;
}

// SIGNAL 0
void DDWin::targets_updated()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
static const uint qt_meta_data_MyGl[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
       6,    5,    5,    5, 0x08,
      22,    5,    5,    5, 0x08,
      34,    5,    5,    5, 0x08,
      44,    5,    5,    5, 0x08,
      62,    5,    5,    5, 0x08,
      85,   81,    5,    5, 0x08,
     121,  115,    5,    5, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MyGl[] = {
    "MyGl\0\0set_wireframe()\0set_stick()\0"
    "set_cpk()\0set_ballandline()\0"
    "set_ballandstick()\0x,y\0"
    "head_tracking_update(int,int)\0x,y,z\0"
    "move_camera(float,float,float)\0"
};

const QMetaObject MyGl::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_MyGl,
      qt_meta_data_MyGl, 0 }
};

const QMetaObject *MyGl::metaObject() const
{
    return &staticMetaObject;
}

void *MyGl::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyGl))
	return static_cast<void*>(const_cast< MyGl*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int MyGl::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: set_wireframe(); break;
        case 1: set_stick(); break;
        case 2: set_cpk(); break;
        case 3: set_ballandline(); break;
        case 4: set_ballandstick(); break;
        case 5: head_tracking_update((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 6: move_camera((*reinterpret_cast< float(*)>(_a[1])),(*reinterpret_cast< float(*)>(_a[2])),(*reinterpret_cast< float(*)>(_a[3]))); break;
        }
        _id -= 7;
    }
    return _id;
}