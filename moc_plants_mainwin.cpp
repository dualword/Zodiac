/****************************************************************************
** Meta object code from reading C++ file 'plants_mainwin.h'
**
** Created: Tue May 6 09:19:03 2008
**      by: The Qt Meta Object Compiler version 59 (Qt 4.3.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "plants_mainwin.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'plants_mainwin.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.3.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_PlantsMainWin[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
      26,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x08,
      32,   14,   14,   14, 0x08,
      42,   14,   14,   14, 0x08,
      51,   14,   14,   14, 0x08,
      60,   14,   14,   14, 0x08,
      69,   14,   14,   14, 0x08,
      78,   14,   14,   14, 0x08,
      94,   14,   14,   14, 0x08,
     110,   14,   14,   14, 0x08,
     125,   14,   14,   14, 0x08,
     140,   14,   14,   14, 0x08,
     155,   14,   14,   14, 0x08,
     170,   14,   14,   14, 0x08,
     184,   14,   14,   14, 0x08,
     198,   14,   14,   14, 0x08,
     212,   14,   14,   14, 0x08,
     228,   14,   14,   14, 0x08,
     261,  255,   14,   14, 0x08,
     282,   14,   14,   14, 0x08,
     297,   14,   14,   14, 0x08,
     311,   14,   14,   14, 0x08,
     326,   14,   14,   14, 0x08,
     346,   14,   14,   14, 0x08,
     358,   14,   14,   14, 0x08,
     371,   14,   14,   14, 0x08,
     382,   14,   14,   14, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_PlantsMainWin[] = {
    "PlantsMainWin\0\0add_water_slot()\0"
    "addphbc()\0addshc()\0addsdc()\0addidc()\0"
    "addldc()\0add_flsc_slot()\0add_phbc_slot()\0"
    "add_shc_slot()\0add_sdc_slot()\0"
    "add_idc_slot()\0add_ldc_slot()\0"
    "add_lf_slot()\0add_ll_slot()\0set_sf_slot()\0"
    "add_fixp_slot()\0load_bs_from_ligand_slot()\0"
    "index\0ligand_selected(int)\0load_protein()\0"
    "load_ligand()\0check_outdir()\0"
    "check_bindingsite()\0del_ifile()\0"
    "del_constr()\0del_flex()\0del_water()\0"
};

const QMetaObject PlantsMainWin::staticMetaObject = {
    { &QTabWidget::staticMetaObject, qt_meta_stringdata_PlantsMainWin,
      qt_meta_data_PlantsMainWin, 0 }
};

const QMetaObject *PlantsMainWin::metaObject() const
{
    return &staticMetaObject;
}

void *PlantsMainWin::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PlantsMainWin))
	return static_cast<void*>(const_cast< PlantsMainWin*>(this));
    return QTabWidget::qt_metacast(_clname);
}

int PlantsMainWin::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QTabWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: add_water_slot(); break;
        case 1: addphbc(); break;
        case 2: addshc(); break;
        case 3: addsdc(); break;
        case 4: addidc(); break;
        case 5: addldc(); break;
        case 6: add_flsc_slot(); break;
        case 7: add_phbc_slot(); break;
        case 8: add_shc_slot(); break;
        case 9: add_sdc_slot(); break;
        case 10: add_idc_slot(); break;
        case 11: add_ldc_slot(); break;
        case 12: add_lf_slot(); break;
        case 13: add_ll_slot(); break;
        case 14: set_sf_slot(); break;
        case 15: add_fixp_slot(); break;
        case 16: load_bs_from_ligand_slot(); break;
        case 17: ligand_selected((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 18: load_protein(); break;
        case 19: load_ligand(); break;
        case 20: check_outdir(); break;
        case 21: check_bindingsite(); break;
        case 22: del_ifile(); break;
        case 23: del_constr(); break;
        case 24: del_flex(); break;
        case 25: del_water(); break;
        }
        _id -= 26;
    }
    return _id;
}
