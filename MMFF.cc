/*   Zodiac - molecular modelling environment. www.zeden.org
    Copyright (C) 2008  Nicola Zonta (nicola.zonta(at)zeden.org)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#include <iostream>
#include "MMFF.h"


#include <limits>   // To enable NaN usage
#ifdef _OPENMP

#include <omp.h>
#endif

#define ACCEPTOR 1
#define DONOR 2
#define NEITHER 0

using namespace std;



#ifdef WIN32
#define isnan(x) (x != x)
#endif // WIN32


/*
MMFFabInteraction::MMFFabInteraction () {
    at3 = NULL;
}
*/







float MMFFbsInteraction::value () {
    assert (at1);
    assert (at2);
    assert (at1 != at2);
    double r = dist (get_coordinates (at1), get_coordinates (at2));
    double dr = r - r0;
    double cs = -2.;
	double out = 143.9325*kb*0.5*(dr*dr)*(1.+cs*dr+(7./12)*(cs*cs*dr*dr));
  //   cout <<"BOND STRETCHINGS "<<at1->GetIdx ()<<" "<<at1->GetAtomicNum ()<<" "<< at1->GetVector ()<<" "<<at2->GetVector ()<<" "<<get_MMFFtype(at1)<<" "<< get_MMFFtype(at1)<<" "<< type<<" "<< r<<" "<<r0<<" "<<dr<<" "<<out<<" "<<kb<<endl;

//	cerr << out<<"out"<<endl;
    assert (!isnan (out));
    return (float) out;

}




float MMFFabInteraction::value () {
    assert (at1);
    assert (at2);
    assert (at3);
    assert (at1 != at3);
    double theta = angle (get_coordinates (at1), get_coordinates (at2), get_coordinates (at3));
    double cb = -0.4*PI/180;
    double dtheta = theta - theta0;

    if (linear) {  
 //       cout << "linear"<<endl;
        double out = 143.9325*ka* (1+ cos (theta/180*PI)); 
    assert (!isnan (out));

        return (float) out;
    }
    else {
		double out = 0.043844*ka*0.5*dtheta*dtheta*(1+cb*dtheta);
		assert (!isnan (out));
	//	cout <<"ANGLE BENDINGS "<<at1->GetIdx ()<<" "<<at2->GetIdx ()<<" "<<at3->GetIdx ()<<" "<< at1->GetVector ()<<" "<<at2->GetVector ()<<" "<<at3->GetVector ()<<" "<<theta0<<" "<<theta<<" "<< out <<" "<<endl;
        return (float) out;
    }
    return 0.f;
}



float MMFFsbInteraction::value () {
    assert (at1);
    assert (at2);
    assert (at3);
//	cerr <<at1 -> GetVector () << at2 -> GetVector ()<< at3 -> GetVector ()<<angle ((vect &)at1-> GetVector (), (vect &)at2-> GetVector (), (vect &)at3-> GetVector ())<<" "<<theta0<<endl;
    double theta = angle (get_coordinates (at1), get_coordinates (at2), get_coordinates (at3));
    double rij = dist (get_coordinates (at1), get_coordinates (at2));
    double rkj = dist (get_coordinates (at3), get_coordinates (at2));
    double drij = rij - r0ij;
    double drkj = rkj - r0kj;
    double dtheta = theta - theta0;   
    double out = 2.5121 *(kijk*drij + kkji*drkj) * dtheta;
//	cerr <<"STRETCH BENDINGS "<< out << " "<< at1 -> GetVector () << at2 -> GetVector ()<< at3 -> GetVector ()<<drij<<" "<<drkj<<dtheta<<endl;

    assert (!isnan (out));
    return (float) out;
}



float MMFFopInteraction::value () {
    assert (at1);
    assert (at2);
    assert (at3);
    assert (at4);
    double chi = 180.* wilson (get_coordinates (at1), get_coordinates (at2), get_coordinates (at3), get_coordinates (at4)) / PI;
  //     cout <<"OUT OF PLANE "<< opint->at1->MMFFtype<<" "<< opint->at2->MMFFtype<<" "<<opint->at3->MMFFtype<<" "<< opint->at4->MMFFtype<<" "<<chi*180/PI<<" "<<k<<" "<< 0.043844*0.5*k*chi*chi <<endl;
    double out = 0.043844*0.5*koop*chi*chi;
    assert (!isnan (out));
    return (float) out;

}


float MMFFtoInteraction::value () {
    assert (at1);
    assert (at2);
    assert (at3);
    assert (at4);
    double phi = dihedral (get_coordinates (at1), get_coordinates (at2), get_coordinates (at3), get_coordinates (at4));
    phi = phi*PI/180; //cos takes rad
  //     cout <<"TORSION INTERACTIONS "<< at1->MMFFtype<<" "<< at2->MMFFtype<<" "<<at3->MMFFtype<<" "<< at4->MMFFtype<<" " <<type<<" "<< phi*180/PI<<" "<<0.5*(v1 * (1 + cos (phi)) + v2 * (1 - cos (2*phi)) + v3 * (1 + cos (3*phi)))<<" "<<v1<<" "<<v2<<" "<<v3<<endl;
    double out = 0.5*(v1 * (1.0 + cos (phi)) + v2 * (1.0 - cos (2.0*phi)) + v3 * (1.0 + cos (3.0*phi)));
    assert (!isnan (out));
    return (float) out;
}



float MMFFvwInteraction::value () {
    assert (at1);
    assert (at2);
    assert (at1 != at2);
#if 1
    double r = dist(get_coordinates (at1), get_coordinates (at2));
    double out = e*pow((1.07*r0/(r+0.07*r0)),7)*((1.12*pow(r0,7)/(pow(r,7)+0.12*pow(r0,7)))-2.);
#else
    double r2 = square_distance(get_coordinates (at1), get_coordinates (at2));
    double r = sqrt (r2);
	
	double rpow7 = r2 * r2 * r2 * r;
     //  cout <<"VDW "<< at1->ID<<" "<< at2->ID<<" "<< r<<" "<<r0<<" "<<e<<" "<<e*pow((1.07*r0/(r+0.07*r0)),7)*((1.12*pow(r0,7)/(pow(r,7)+0.12*pow(r0,7)))-2)<<endl;ble 
	
	double r0pow7 = (r0 * r0); // pow(r0, 7);
	r0pow7 = r0pow7 * r0pow7 * r0pow7 * r0;
	
	double rX = (1.07*r0/(r+0.07*r0));
	double rXpow7 = rX * rX;
	rXpow7 = rXpow7 * rXpow7 * rXpow7 * rX;
	
    double out = e*rXpow7*((1.12*r0pow7/(rpow7+0.12*r0pow7))-2.0);
#endif
    assert (!isnan (out));
    return (float) out;
}


float MMFFelInteraction::value () {
    assert (at1);
    assert (at2);

    double qi = at1-> GetPartialCharge ();
    double qj = at2-> GetPartialCharge ();
    double delta = 0.05;
    double r = dist (get_coordinates (at1), get_coordinates (at2));
  //  double D = 1.0; D assumed to be 1
   // int n = 1; n assumed to be 1
    double out = 332.0716 * qi * qj / (r + delta) * scale;
    assert (!isnan (out));
    return (float) out;
}

bool MMFFelInteraction::isHbond () {
	int at1n = at1 ->GetAtomicNum ();
	int at2n = at2 ->GetAtomicNum ();
	if (at1n == 1) {
		Atom *r = root_at (at1);
		if (r) {
			if (is_polar (r)) {
				if (is_polar (at2)) return true;
			};
		}
	}
	else if (at2n == 1) {
		Atom *r = root_at (at2);
		if (r) {
			if (is_polar (r)) {
				if (is_polar (at1)) return true;
			};
		}
	}
	
	return false;
}





void MMFFabInteraction::set_forces (bool score, double multbm) {
    assert (at1);
    assert (at2);
    assert (at3);
    float force_1x = derive (at1);
    float force_1y = derive_y (at1);
    float force_1z = derive_z (at1);
    vect force1 (-force_1x, -force_1y, -force_1z);

    float force_2x = derive (at2);
    float force_2y = derive_y (at2);
    float force_2z = derive_z (at2);
    vect force2 (-force_2x, -force_2y, -force_2z);

    float force_3x = derive (at3);
    float force_3y = derive_y (at3);
    float force_3z = derive_z (at3);
    vect force3 (-force_3x, -force_3y, -force_3z);
	assert (!isnan (force_1x));	assert (!isnan (force_1x));	assert (!isnan (force_1x));
	assert (!isnan (force_1y));	assert (!isnan (force_1y));	assert (!isnan (force_1y));
	assert (!isnan (force_1z));	assert (!isnan (force_1z));	assert (!isnan (force_1z));



//	cerr << force1<<force2<<force3<<endl;
	vect pforce1 = get_back_force (at1);
    vect pforce2 = get_back_force (at2);
    vect pforce3 = get_back_force (at3);
	//cerr << "force was " << pforce1;
	pforce1 = sum (pforce1, force1);
	pforce2 = sum (pforce2, force2);
	pforce3 = sum (pforce3, force3);
//	cerr << " and now after "<<force1<< " is "<<pforce1<<endl;

	set_back_force (at1, pforce1);
	set_back_force (at2, pforce2);
	set_back_force (at3, pforce3);
//	cerr << "AB interaction "<< pforce1<<pforce2<<pforce3<<endl;
}



void MMFFsbInteraction::set_forces (bool score, double mult) {
    assert (at1);
    assert (at2);
    assert (at3);
    float force_1x = derive (at1);
    float force_1y = derive_y (at1);
    float force_1z = derive_z (at1);
    vect force1 (-force_1x, -force_1y, -force_1z);

    float force_2x = derive (at2);
    float force_2y = derive_y (at2);
    float force_2z = derive_z (at2);
    vect force2 (-force_2x, -force_2y, -force_2z);

    float force_3x = derive (at3);
    float force_3y = derive_y (at3);
    float force_3z = derive_z (at3);
    vect force3 (-force_3x, -force_3y, -force_3z);
	assert (!isnan (force_1x));	assert (!isnan (force_1x));	assert (!isnan (force_1x));
	assert (!isnan (force_1y));	assert (!isnan (force_1y));	assert (!isnan (force_1y));
	assert (!isnan (force_1z));	assert (!isnan (force_1z));	assert (!isnan (force_1z));



//	cerr << force1<<force2<<force3<<endl;
	vect pforce1 = get_back_force (at1);
    vect pforce2 = get_back_force (at2);
    vect pforce3 = get_back_force (at3);
	//cerr << "force was " << pforce1;
	pforce1 = sum (pforce1, force1);
	pforce2 = sum (pforce2, force2);
	pforce3 = sum (pforce3, force3);
//	cerr << " and now after "<<force1<< " is "<<pforce1<<endl;

	set_back_force (at1, pforce1);
	set_back_force (at2, pforce2);
	set_back_force (at3, pforce3);
//	cerr << "SB interaction "<< pforce1<<pforce2<<pforce3<<endl;
}




void MMFFopInteraction::set_forces (bool score, double mult) {

    assert (at1);
    assert (at2);
    assert (at3);
    assert (at4);
    float force_1x = derive (at1);
    float force_1y = derive_y (at1);
    float force_1z = derive_z (at1);
    vect force1 (-force_1x, -force_1y, -force_1z);

    float force_2x = derive (at2);
    float force_2y = derive_y (at2);
    float force_2z = derive_z (at2);
    vect force2 (-force_2x, -force_2y, -force_2z);

    float force_3x = derive (at3);
    float force_3y = derive_y (at3);
    float force_3z = derive_z (at3);
    vect force3 (-force_3x, -force_3y, -force_3z);

    float force_4x = derive (at4);
    float force_4y = derive_y (at4);
    float force_4z = derive_z (at4);
    vect force4 (-force_4x, -force_4y, -force_4z);

    vect pforce1 = get_back_force (at1);
    vect pforce2 = get_back_force (at2);
    vect pforce3 = get_back_force (at3);
    vect pforce4 = get_back_force (at4);

    pforce1 += force1;
    pforce2 += force2;
    pforce3 += force3;
    pforce4 += force4;

	set_back_force (at1, pforce1);
	set_back_force (at2, pforce2);
	set_back_force (at3, pforce3);
	set_back_force (at4, pforce4); 
}


void MMFFtoInteraction::set_forces (bool score, double mult) {

    assert (at1);
    assert (at2);
    assert (at3);
    assert (at4);
    float force_1x = derive (at1);
    float force_1y = derive_y (at1);
    float force_1z = derive_z (at1);
    vect force1 (-force_1x, -force_1y, -force_1z);

    float force_2x = derive (at2);
    float force_2y = derive_y (at2);
    float force_2z = derive_z (at2);
    vect force2 (-force_2x, -force_2y, -force_2z);

    float force_3x = derive (at3);
    float force_3y = derive_y (at3);
    float force_3z = derive_z (at3);
    vect force3 (-force_3x, -force_3y, -force_3z);

    float force_4x = derive (at4);
    float force_4y = derive_y (at4);
    float force_4z = derive_z (at4);
    vect force4 (-force_4x, -force_4y, -force_4z);

    vect pforce1 = get_back_force (at1);
    vect pforce2 = get_back_force (at2);
    vect pforce3 = get_back_force (at3);
    vect pforce4 = get_back_force (at4);

    pforce1 += force1;
    pforce2 += force2;
    pforce3 += force3;
    pforce4 += force4;
	
	set_back_force (at1, pforce1);
	set_back_force (at2, pforce2);
	set_back_force (at3, pforce3);
	set_back_force (at4, pforce4); 
}






MMFF::MMFF () {
    is_initialised = false;
    is_initialised = load_parameters ();
    
}

int MMFF::load_parameters () {
    int out = 1;
    out *= load_aromatic_parameters ();
    out *= load_ci_parameters ();
    out *= load_H_parameters ();
    out *= load_atype_parameters ();
    out *= load_atom_parameters ();
    out *= load_ab_parameters ();
    out *= load_bs_parameters ();
    out *= load_sb_parameters ();
    out *= load_sb2_parameters ();
    out *= load_to_parameters ();
    out *= load_vw_parameters ();
    out *= load_op_parameters ();
    return out;
}


double MMFF::compute_total_energy () {

    double Einternal, Eb, Ea, Eba, Eoop, Et, Evdw, Eq;

    
    Eb = compute_bond_stretchings ();
    Ea = compute_angle_bendings ();
    Eba = compute_stretch_bend_interactions ();
    Eoop = compute_out_of_plane_bendings ();
    Et = compute_torsion_interactions ();
    Evdw = compute_van_der_waals_interactions ();
    Eq = compute_electrostatic_interactions ();
    Einternal = Eb + Ea + Eba + Eoop + Et + Evdw + Eq;
	
	double Einteraction = compute_interaction_energy ();
    return Einternal + Einteraction;

}

double MMFF::compute_interaction_energy () {
	double Eq = compute_nonbonded_electrostatic_interactions ();
	double Evdw = compute_nonbonded_van_der_waals_interactions ();
	return Eq + Evdw;
}
    

void MMFF::compute_forces () {

    compute_bond_stretching_forces ();
	compute_angle_bending_forces ();
    compute_stretch_bend_interaction_forces ();
 //   compute_out_of_plane_bending_forces ();
    compute_torsion_forces ();
    compute_van_der_waals_forces ();
    compute_electrostatic_forces ();
    compute_nonbonded_van_der_waals_forces ();
    compute_nonbonded_electrostatic_forces ();

}



void MMFF::compute_bond_stretching_forces () {
	//cerr << "bs " << bsInteractions.size () << endl;
/*#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
*/
    for (int i=0; i<bsInteractions.size (); i++) {
        bsInteractions[i] -> set_forces ();
    }
	
}

void MMFF::compute_angle_bending_forces () {
//	cerr << "ab " << abInteractions.size () << endl;
/*#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
*/
    for (int i=0; i<abInteractions.size (); i++) {
        abInteractions[i] -> set_forces ();
    }
}

void MMFF::compute_stretch_bend_interaction_forces () {
/*#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
*/
    for (int i=0; i<sbInteractions.size (); i++) {
        sbInteractions[i] -> set_forces ();
    }
}

void MMFF::compute_out_of_plane_bending_forces () {
/*#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
*/
    for (int i=0; i<opInteractions.size (); i++) {
        opInteractions[i] -> set_forces ();
    }
}

void MMFF::compute_torsion_forces () {
/*#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
*/
    for (int i=0; i<toInteractions.size (); i++) {
        toInteractions[i] -> set_forces ();
    }
}

void MMFF::compute_electrostatic_forces () {
/*#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
*/
    for (int i=0; i<elInteractions.size (); i++) {
        if (elInteractions[i]->at1->GetPartialCharge () && elInteractions[i]->at2->GetPartialCharge ()) {
            elInteractions[i] -> set_forces ();
        }
    }
}

void MMFF::compute_nonbonded_electrostatic_forces () {
/*#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
*/
    for (int i=0; i<elNBInteractions.size (); i++) {
        if (elNBInteractions[i]->at1->GetPartialCharge () && elNBInteractions[i]->at2->GetPartialCharge ()) {
            elNBInteractions[i] -> set_forces (true); //bool enables scores to be written to atoms
        }
    }
}

void MMFF::compute_van_der_waals_forces () {
/*#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
*/
    for (int i=0; i<vwInteractions.size (); i++) {
        vwInteractions[i] -> set_forces (); 
    }
}

void MMFF::compute_nonbonded_van_der_waals_forces () {
/*#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
*/
    for (int i=0; i<vwNBInteractions.size (); i++) {
        vwNBInteractions[i] -> set_forces (true); //bool enables scores to be written to atoms
    }
}




double MMFF::compute_nonbonded_electrostatic_interactions () {
    double Eq = 0.f;
    for (unsigned int i=0; i<elNBInteractions.size (); i++) {
        Eq += elNBInteractions[i] -> value ();

    }
    return Eq;
}


double MMFF::compute_nonbonded_van_der_waals_interactions () {
    double Evdw = 0.f;
    for (unsigned int i=0; i<vwNBInteractions.size (); i++) {
        Evdw += vwNBInteractions[i] -> value ();

    }
    return Evdw;
}


double MMFF::compute_bond_stretchings () {
    double Eb = 0.f;
    for (unsigned int i=0; i<bsInteractions.size (); i++) {
        Eb += bsInteractions[i] -> value ();

    }
    return Eb;
}






double MMFF::compute_electrostatic_interactions () {
    double Eq = 0.f;
    for (unsigned int i=0; i<elInteractions.size (); i++) {
        Eq += elInteractions[i] -> value ();
    }
    return Eq;
}


double MMFF::compute_angle_bendings () {

    double Ea = 0.f;
    for (unsigned int i=0; i<abInteractions.size (); i++) {
        Ea += abInteractions[i] -> value ();
    }
    return Ea;
}

double MMFF::compute_stretch_bend_interactions () {
    double Eba = 0.f;
    for (unsigned int i=0; i<sbInteractions.size (); i++) {
        Eba += sbInteractions[i] -> value ();
    }
    return Eba;
}

double MMFF::compute_torsion_interactions () {
    double Et = 0.f;
    for (unsigned int i=0; i<toInteractions.size (); i++) {
        Et += toInteractions[i] -> value ();
    } 
    return Et;
}

double MMFF::compute_van_der_waals_interactions () {
    double Evdw = 0.f;
    for (unsigned int i=0; i<vwInteractions.size (); i++) {
        Evdw += vwInteractions[i] -> value ();
    }
    return Evdw;
}

double MMFF::compute_out_of_plane_bendings () {
    double Eoop = 0.f;
    for (unsigned int i=0; i<opInteractions.size (); i++) {
        Eoop += opInteractions[i] -> value ();
    }
    return Eoop;
}





/*

float MMFF::compute_total_electrostatic_energy () {
    float Eq = 0.;
    for (unsigned int i=0; i<elNBInteractions.size (); i++) {
        float E = compute_electrostatic_interaction (elNBInteractions[i]);
        Eq += E;
        elNBInteractions[i]->at1->score += E;
    }
    return Eq;
}

float MMFF::compute_total_van_der_waals_energy () {
    float Evdw = 0.;
    for (unsigned int i=0; i<vwNBInteractions.size (); i++) {
        float E = compute_van_der_waals_interaction (vwNBInteractions[i]);
        Evdw += E;
        vwNBInteractions[i]->at1->score += E;
    }
    return Evdw;
}



*/


















///////////////////////////////////////////FORCES////////////////////////////////////////////////////////////////////////


/*


void MMFF::compute_bond_stretching_force (MMFFbsInteraction *bsint) {
    float E, E2;
    for (unsigned int i=0; i<3; i++) {
        float old_coordinate = bsint->at1-> GetVector ()[i];
        bsint->at1-> GetVector ()[i] = old_coordinate - DX;
        E = compute_bond_stretching (bsint);
        bsint->at1-> GetVector ()[i] = old_coordinate + DX;
        E2 = compute_bond_stretching (bsint);
        bsint->at1-> GetVector ()[i] = old_coordinate;
        bsint->at1->force[i]-= (E2-E)/(2*DX);
        bsint->at2->force[i]+= (E2-E)/(2*DX);


    }

        total_energy += (E2+E)/2;
}

void MMFF::compute_angle_bending_force (MMFFabInteraction *abint) {
    float E, E2;
    for (unsigned int i=0; i<3; i++) {
        abint->at1-> GetVector ()[i]-= DX;
        E = compute_angle_bending (abint);
        abint->at1-> GetVector ()[i]+= 2*DX;
        E2 = compute_angle_bending (abint);
        abint->at1-> GetVector ()[i]-= DX;
        abint->at1->force[i]-= (E2-E)/(2*DX);

        abint->at3-> GetVector ()[i]-= DX;
        E = compute_angle_bending (abint);
        abint->at3-> GetVector ()[i]+= 2*DX;
        E2 = compute_angle_bending (abint);
        abint->at3-> GetVector ()[i]-= DX;
        abint->at3->force[i]-= (E2-E)/(2*DX);

        abint->at2-> GetVector ()[i]-= DX;
        E = compute_angle_bending (abint);
        abint->at2-> GetVector ()[i]+= 2*DX;
        E2 = compute_angle_bending (abint);
        abint->at2-> GetVector ()[i]-= DX;
        abint->at2->force[i]-= (E2-E)/(2*DX);

    }

        total_energy += (E2+E)/2;
}

void MMFF::compute_stretch_bend_interaction_force (MMFFsbInteraction *sbint) {
    float E, E2;
    for (unsigned int i=0; i<3; i++) {
        sbint->at1-> GetVector ()[i]-= DX;
        E = compute_stretch_bend_interaction (sbint);
        sbint->at1-> GetVector ()[i]+= 2*DX;
        E2 = compute_stretch_bend_interaction (sbint);
        sbint->at1-> GetVector ()[i]-= DX;
        sbint->at1->force[i]-= (E2-E)/(2*DX);

        sbint->at3-> GetVector ()[i]-= DX;
        E = compute_stretch_bend_interaction (sbint);
        sbint->at3-> GetVector ()[i]+= 2*DX;
        E2 = compute_stretch_bend_interaction (sbint);
        sbint->at3-> GetVector ()[i]-= DX;
        sbint->at3->force[i]-= (E2-E)/(2*DX);

        sbint->at2-> GetVector ()[i]-= DX;
        E = compute_stretch_bend_interaction (sbint);
        sbint->at2-> GetVector ()[i]+= 2*DX;
        E2 = compute_stretch_bend_interaction (sbint);
        sbint->at2-> GetVector ()[i]-= DX;
        sbint->at2->force[i]-= (E2-E)/(2*DX);


    }

        total_energy += (E2+E)/2;
}

void MMFF::compute_out_of_plane_bending_force (MMFFopInteraction *opint) {
    float E, E2;
    for (unsigned int i=0; i<3; i++) {
        opint->at4-> GetVector ()[i]-= DX;
        E = compute_out_of_plane_bending (opint);
        opint->at4-> GetVector ()[i]+= 2*DX;
        E2 = compute_out_of_plane_bending (opint);
        opint->at4-> GetVector ()[i]-= DX;
        opint->at4->force[i]-= (E2-E)/(2*DX);

        opint->at1-> GetVector ()[i]-= DX;
        E = compute_out_of_plane_bending (opint);
        opint->at1-> GetVector ()[i]+= 2*DX;
         E2 = compute_out_of_plane_bending (opint);
        opint->at1-> GetVector ()[i]-= DX;
        opint->at1->force[i]-= (E2-E)/(2*DX);

        opint->at2-> GetVector ()[i]-= DX;
        E = compute_out_of_plane_bending (opint);
        opint->at2-> GetVector ()[i]+= 2*DX;
        E2 = compute_out_of_plane_bending (opint);
        opint->at2-> GetVector ()[i]-= DX;
        opint->at2->force[i]-= (E2-E)/(2*DX);

        opint->at3-> GetVector ()[i]-= DX;
        E = compute_out_of_plane_bending (opint);
        opint->at3-> GetVector ()[i]+= 2*DX;
        E2 = compute_out_of_plane_bending (opint);
        opint->at3-> GetVector ()[i]-= DX;
        opint->at3->force[i]-= (E2-E)/(2*DX);


    }

        total_energy += (E2+E)/2;
}


void MMFF::compute_torsion_force (MMFFtoInteraction *toint) {
    float E, E2;
    for (unsigned int i=0; i<3; i++) {
        toint->at1-> GetVector ()[i]-= DX;
        E = compute_torsion_interaction (toint);
        toint->at1-> GetVector ()[i]+= 2*DX;
        E2 = compute_torsion_interaction (toint);
        toint->at1-> GetVector ()[i]-= DX;
        toint->at1->force[i]-= (E2-E)/(2*DX);

        toint->at4-> GetVector ()[i]-= DX;
        E = compute_torsion_interaction (toint);
        toint->at4-> GetVector ()[i]+= 2*DX;
        E2 = compute_torsion_interaction (toint);
        toint->at4-> GetVector ()[i]-= DX;
        toint->at4->force[i]-= (E2-E)/(2*DX);

        toint->at2-> GetVector ()[i]-= DX;
        E = compute_torsion_interaction (toint);
        toint->at2-> GetVector ()[i]+= 2*DX;
        E2 = compute_torsion_interaction (toint);
        toint->at2-> GetVector ()[i]-= DX;
        toint->at2->force[i]-= (E2-E)/(2*DX);

        toint->at3-> GetVector ()[i]-= DX;
        E = compute_torsion_interaction (toint);
        toint->at3-> GetVector ()[i]+= 2*DX;
        E2 = compute_torsion_interaction (toint);
        toint->at3-> GetVector ()[i]-= DX;
        toint->at3->force[i]-= (E2-E)/(2*DX);





        toint->at4->force[i]-= (E2-E)/(2*DX);

    }

        total_energy += (E2+E)/2;
}


void MMFF::compute_electrostatic_force (MMFFelInteraction *elint) {
    float E, E2;
    for (unsigned int i=0; i<3; i++) {
        elint->at1-> GetVector ()[i]-= DX;
        E = compute_electrostatic_interaction (elint);
        elint->at1-> GetVector ()[i]+= 2*DX;
        E2 = compute_electrostatic_interaction (elint);
        elint->at1-> GetVector ()[i]-= DX;
        elint->at1->force[i]-= (E2-E)/(2*DX);
        elint->at2->force[i]+= (E2-E)/(2*DX);

    }

        total_energy += (E2+E)/2;
}

void MMFF::compute_van_der_waals_force (MMFFvwInteraction *vwint) {
    float E, E2;
    for (unsigned int i=0; i<3; i++) {
        vwint->at1-> GetVector ()[i]-= DX;
        E = compute_van_der_waals_interaction (vwint);
        vwint->at1-> GetVector ()[i]+= 2*DX;
        E2 = compute_van_der_waals_interaction (vwint);
        vwint->at1-> GetVector ()[i]-= DX;
        vwint->at1->force[i]-= (E2-E)/(2*DX);
        vwint->at2->force[i]+= (E2-E)/(2*DX);

    }

            total_energy += (E2+E)/2;
}



vector<float> MMFF::compute_van_der_waals_force_vector (MMFFvwInteraction *vwint) {

    vector<float> out (3, 0);
    double r = distance (vwint->at1-> GetVector (), vwint->at2-> GetVector ());
    float force = compute_van_der_waals_force_module (vwint);
    for (unsigned int i=0; i<3; i++) {
        out[i] = vwint->at2-> GetVector ()[i]-vwint->at1-> GetVector ()[i];
        out[i]*=force/r;
    }
    return out;
}





vector<float> MMFF::compute_electrostatic_force_vector (MMFFelInteraction *elint) {
    float r = distance (elint->at1-> GetVector (), elint->at2-> GetVector ());
    float force = compute_electrostatic_force_module (elint);
    vector<float> out (3, 0);
    for (unsigned int i=0; i<3; i++) {
        out[i] = elint->at2-> GetVector ()[i]-elint->at1-> GetVector ()[i];
        out[i]*=force/r;
    }
    return out;
}





float MMFF::compute_van_der_waals_force_module (MMFFvwInteraction *vwint) {

    
    double e = get_vdw_e (vwint->at1, vwint->at2);
    double r = distance (vwint->at1-> GetVector (), vwint->at2-> GetVector ());
    double r0 = get_vdw_r0 (vwint->at1, vwint->at2);

    return e*pow(1.07*r0,7)*    (     (pow((1/(r+0.07*r0)),7))       *      (1.12*pow(r0,7)*7*(-1)*pow (r0,6)/((pow(r,7)+
    0.12*pow(r0,7))*(pow(r,7)+0.12*pow(r0,7))))  +    ((1.12*pow(r0,7)/(pow(r,7)+0.12 *pow(r0,7)))-2) *   (-7*pow ((r+0.07*r0),-8)));   
}



float MMFF::compute_electrostatic_force_module (MMFFelInteraction *elint) {

    double qi = elint->at1-> charge;
    double qj = elint->at2-> charge;
    double delta = 0.05;
    double r = distance (elint->at1-> GetVector (), elint->at2-> GetVector ());

    return elint->scale*-332.0716*qi*qj/((r+delta)*(r+delta));
}






float MMFF::compute_bond_stretching_force_module (MMFFbsInteraction *bsint) {
    float kb = bsint->kb;
    float r = distance (bsint->at1-> GetVector (), bsint->at2-> GetVector ());
    float r0 = bsint->r0;
    float dr = r-r0;
    float cs = -2.;
    float out= 143.9325*0.5*kb*(2*dr+cs*3*dr*dr+(7./12)*cs*cs*4*dr*dr*dr);
    return out;

}


float MMFF::compute_angle_bending_force_module (MMFFabInteraction *abint) {

    float theta = angle (abint->at1-> GetVector (), abint->at2-> GetVector (), abint->at3-> GetVector ());
    float ka = abint->ka;
    float cb = -0.4*PI/180;
    float theta0 = abint->theta0;
    float dtheta = theta - theta0;
    bool linear = abint->linear;
    if (linear) {  
        return -143.9325*ka* (sin (theta/180*PI));
    }
    else {
        return 0.043844*ka*0.5*(2*dtheta+cb*3*dtheta*dtheta);
    }
    return 0;
}



float MMFF::compute_out_of_plane_bending_force_module (MMFFopInteraction *opint) {
    double k = opint->koop;
    double chi = 180.*wilson (opint->at1-> GetVector (), opint->at2-> GetVector (),opint->at3-> GetVector (),opint->at4-> GetVector ())/PI;
    return 0.043844*k*chi;
}



float MMFF::compute_torsion_force_module (MMFFtoInteraction *toint) {
    int I = toint->at1->MMFFtype;
    int J = toint->at2->MMFFtype;
    int K = toint->at3->MMFFtype;
    int L = toint->at4->MMFFtype;


//    int pl = toint->pl;
//get_to_pl (toint->type, I, J, K, L); //parameters line
    float v1 = toint->v1;
    float v2 = toint->v2;
    float v3 = toint->v3;
    float phi = dihedral (toint->at1-> GetVector (), toint->at2-> GetVector (),toint-> at3-> GetVector (), toint->at4-> GetVector ());
    phi = phi*PI/180; //cos takes rads
    return 0.5*(v1 * (-sin(phi)) + v2 * (2*sin (2.*phi)) + v3 * (-3*sin (3.*phi)));
}

*/

///////////////////////////////////////////ENERGIES////////////////////////////////////////////////////////////////////////


/*
double MMFF::compute_bond_stretching (MMFFbsInteraction *bsint) {


    double kb = bsint->kb;
    double r = distance (bsint->at1-> GetVector (), bsint->at2-> GetVector ());
    double r0 = bsint->r0;
    double dr = r-r0;
    double cs = -2.;
   //     cout <<"BOND STRETCHINGS "<< bsint->at1->MMFFtype<<" "<< bsint->at2->MMFFtype<<" "<< bsint->type<<" "<< r<<" "<<r0<<" "<<dr<<" "<<143.9325*kb*0.5*dr*dr*(1+cs*dr+7/12*cs*cs*dr*dr)<<" "<<kb<<endl;
    long double out = 143.9325*kb*0.5*(dr*dr)*(1.+cs*dr+(7./12)*(cs*cs*dr*dr));
  //  if (out<0) cout<< "dr "<<dr<<"   cs "<<cs<<"   primo pezzo "<< 143.9325*kb*0.5*(dr*dr)<<"   secondo pezzo "<<(1+cs*dr+(7/12)*(cs*cs*dr*dr))<<endl;
    return out;

}





double MMFF::compute_angle_bending (MMFFabInteraction *abint) {

   double theta = angle (abint->at1-> GetVector (), abint->at2-> GetVector (), abint->at3-> GetVector ());
    double ka = abint->ka;
    double cb = -0.4*PI/180;
    double theta0 = abint->theta0;
    double dtheta = theta - theta0;
//  cout << theta0<< " "<<dtheta<<endl;
    
    bool linear = abint->linear;
  //  cout <<"ANGLE BENDINGS "<< abint->at1->MMFFtype<<" "<< abint->at2->MMFFtype<<" "<< abint->type<<" "<<theta0<<" "<<theta<<" "<< 0.043844*ka*0.5*dtheta*dtheta*(1+cb*dtheta) <<" "<<endl;
    if (linear) {  
        cout << "linear"<<endl;
        return 143.9325*ka* (1+ cos (theta/180*PI)); 
    }
    else {
        return 0.043844*ka*0.5*dtheta*dtheta*(1+cb*dtheta);
    }
    return 0;
}









double MMFF::compute_stretch_bend_interaction (MMFFsbInteraction *sbint) {

    double theta = angle (sbint->at1-> GetVector (), sbint->at2-> GetVector (), sbint->at3-> GetVector ());
    double theta0 = sbint->theta0;
    double kijk = sbint->kijk;
    double kkji = sbint->kkji;
    double rij = distance (sbint->at1-> GetVector (), sbint->at2-> GetVector ());
    double rkj = distance (sbint->at3-> GetVector (), sbint->at2-> GetVector ());
    double r0ij = sbint->r0ij;
    double r0kj = sbint->r0kj;
    double drij = rij - r0ij;
    double drkj = rkj - r0kj;
    double dtheta = theta - theta0;   
//       cout <<"STRETCH BENDINGS "<< sbint->at1->MMFFtype<<" "<< sbint->at2->MMFFtype<<" "<<sbint->at3->MMFFtype<<" "<< theta<<" "<<dtheta<<" "<<2.5121 *(kijk*drij) * dtheta<<" "<<2.5121 *(kkji*drkj) * dtheta<<" "<<kijk<<" "<<kkji<<endl;
    return 2.5121 *(kijk*drij + kkji*drkj) * dtheta;
}



double MMFF::compute_out_of_plane_bending (MMFFopInteraction *opint) {
    double k = opint->koop;
    double chi = 180.*wilson (opint->at1-> GetVector (), opint->at2-> GetVector (),opint->at3-> GetVector (),opint->at4-> GetVector ())/PI;
  //     cout <<"OUT OF PLANE "<< opint->at1->MMFFtype<<" "<< opint->at2->MMFFtype<<" "<<opint->at3->MMFFtype<<" "<< opint->at4->MMFFtype<<" "<<chi*180/PI<<" "<<k<<" "<< 0.043844*0.5*k*chi*chi <<endl;
    return 0.043844*0.5*k*chi*chi;

}






double MMFF::compute_torsion_interaction (MMFFtoInteraction *toint) {
    int I = toint->at1->MMFFtype;
    int J = toint->at2->MMFFtype;
    int K = toint->at3->MMFFtype;
    int L = toint->at4->MMFFtype;


    double v1 = toint->v1;
    double v2 = toint->v2;
    double v3 = toint->v3;
    double phi = dihedral (toint->at1-> GetVector (), toint->at2-> GetVector (),toint-> at3-> GetVector (), toint->at4-> GetVector ());
    phi = phi*PI/180; //cos takes rad
  //     cout <<"TORSION INTERACTIONS "<< toint->at1->MMFFtype<<" "<< toint->at2->MMFFtype<<" "<<toint->at3->MMFFtype<<" "<< toint->at4->MMFFtype<<" " <<toint->type<<" "<< phi*180/PI<<" "<<0.5*(v1 * (1 + cos (phi)) + v2 * (1 - cos (2*phi)) + v3 * (1 + cos (3*phi)))<<" "<<v1<<" "<<v2<<" "<<v3<<endl;
    return 0.5*(v1 * (1.0 + cos (phi)) + v2 * (1.0 - cos (2.0*phi)) + v3 * (1.0 + cos (3.0*phi)));
}



double MMFF::compute_van_der_waals_interaction (MMFFvwInteraction *vwint) {
    
    double e = get_vdw_e (vwint->at1, vwint->at2);
    double r = distance (vwint->at1-> GetVector (), vwint->at2-> GetVector ());
    double r0 = get_vdw_r0 (vwint->at1, vwint->at2);
   //    cout <<"VDW "<< vwint->at1->ID<<" "<< vwint->at2->ID<<" "<< r<<" "<<r0<<" "<<e<<" "<<e*pow((1.07*r0/(r+0.07*r0)),7)*((1.12*pow(r0,7)/(pow(r,7)+0.12*pow(r0,7)))-2)<<endl;
    return e*pow((1.07*r0/(r+0.07*r0)),7)*((1.12*pow(r0,7)/(pow(r,7)+0.12*pow(r0,7)))-2.);
}


double MMFF::compute_electrostatic_interaction (MMFFelInteraction *elint) {

    double qi = elint->at1-> charge;
    double qj = elint->at2-> charge;
    double delta = 0.05;
    double r = distance (elint->at1-> GetVector (), elint->at2-> GetVector ());
  //  double D = 1.0; D assumed to be 1
   // int n = 1; n assumed to be 1
    return 332.0716*qi*qj/(r+delta)*elint->scale;
}

*/


//______________________________________________________________________________________________________________

int MMFF::get_atype_line (int I) 
{
    for (unsigned int i = 0; i< atomParameters.size (); i++) 
    {
        if (atomParameters[i]->type == I) 
        {
            return i;
        }
    }

//    cerr <<i<<"non parametrized atom"<<endl;
    cerr << I << "non parametrized atom" << endl;
    return 0;
}

int MMFF::get_bond_type (int I, int J) {
    int iline = 0, jline =0;
    for (unsigned int i = 0; i< atomParameters.size (); i++) {
        if (atomParameters[i]->type == I) {
            iline = i;
            break;
        }
    }
    for (unsigned int i = 0; i< atomParameters.size (); i++) {
        if (atomParameters[i]->type == J) {
            jline = i;
            break;
        }    
    }
 //   cout <<iline<<" "<<jline<<endl;
 //   cout << atomParameters[iline]->multipleBond <<atomParameters[jline]->multipleBond<<endl;
//    if ((atomParameters[iline]->aromatic || atomParameters[iline]->multipleBond) && (atomParameters[jline]->aromatic || atomParameters[jline]->multipleBond) ) return 1;

    if (atomParameters[iline]->multipleBond && atomParameters[jline]->multipleBond ) return 1;
    else return 0;
}

double MMFF::get_charge_increment (int I, int J) {
    if (I==J) return 0;
    int correct = 1;
    if (J<I) {
        correct = -1;
        int swap = J;
        J = I;
        I = swap;
    }
    for (unsigned int n=0; n<ciParameters.size (); n++) {
        if (ciParameters[n]->I == I && ciParameters[n]->J == J ){
            return ciParameters[n]->increment*correct;
        }
    }
    cout <<"could not retrieve ci informations for atom types "<<I<<" && "<<J<<endl;
    return 0;
}


string MMFF::get_aromatic_string (Atom *at) {
    return "not implemented";
/*
    Ring *ring=NULL;
    int last_size = 100;
    for (unsigned int r=0; r<at->in_ring.size (); r++) {
        if (at->in_ring[r]->IsAromatic () && at->in_ring[r]->Size ()<last_size) {
            ring = at->in_ring[r];
            last_size = ring->Size ();
        }
    }
    bool N5ANION, IMCAT;
    IMCAT = ring->IM_CAT;
    N5ANION = ring->N5ANION;
    int L5;
    if (last_size ==5) {
        if (ring->alpha_atom == NULL) L5 =4;
        else if (ring->alpha_atom == at) L5 =1;
        else if (at-> GetBond (ring->alpha_atom)) L5=2;
        else L5=3;
    }
    else {
        L5 =0;
    }

//    cout <<at->MMFFstring<<" "<<at->GetAtomicNum ()<<" "<<last_size<<" "<<L5<<" "<<IMCAT<<" "<<N5ANION<<endl;
    for (unsigned int n=0; n<aromaticParameters.size (); n++) {
        if (aromaticParameters[n]->old_type==at->MMFFstring && aromaticParameters[n]->at_number==at->GetAtomicNum () && aromaticParameters[n]->ring_size==last_size && aromaticParameters[n]->L5==L5 && aromaticParameters[n]->IMCAT==IMCAT && aromaticParameters[n]->N5ANION==N5ANION){
            return aromaticParameters[n]->arom_type;
        }
    }
    string atype;
    if (at->GetAtomicNum () ==6) atype = "C*";
    else if (at->GetAtomicNum () ==7) atype = "N*";
    else if (at->GetAtomicNum () ==8) atype = "O*";
    else if (at->GetAtomicNum () ==16) atype = "S*";
    else atype ="*";
    for (unsigned int n=0; n<aromaticParameters.size (); n++) {
        if (aromaticParameters[n]->old_type==atype && aromaticParameters[n]->at_number==at->GetAtomicNum () && aromaticParameters[n]->ring_size==last_size && aromaticParameters[n]->L5==L5){
            return aromaticParameters[n]->arom_type;
        }
    }

    return "unknown aromatic atom";
*/
}

int MMFF::get_bs_pl (int type, int I, int J) {
    for (unsigned int n=0; n<bsParameters.size (); n++) {
        if (bsParameters[n]->I == I && bsParameters[n]->J == J && bsParameters[n]->type == type){
            return n;
        }
    }
    cout <<"could not retrieve bs informations for atom types "<<I<<" , "<<J<<" && bond type "<<type<<endl;
    return 0;
}

int MMFF::get_ab_pl (int type, int I, int J, int K) {
    for (unsigned int n=0; n<abParameters.size (); n++) {
        if (abParameters[n]->I == I && abParameters[n]->J == J && abParameters[n]->K == K && abParameters[n]->type == type){
            return n;
        }
    }
    cout <<"could not retrieve ab informations for atom types "<<I<<" , "<<J<<" , "<<K<<" && angle type "<<type<<endl;
    return 0;
}

int MMFF::get_sb_pl (int I, int J, int K) {
    for (unsigned int n=0; n<sbParameters.size (); n++) {
        if (sbParameters[n]->I == I && sbParameters[n]->J == J && sbParameters[n]->K == K){
            return n;
        }
    }
//    cout <<"could not retrieve sb informations for atom types "<<I<<" , "<<J<<" && "<<K<<endl;
    return -1;
}

int MMFF::get_sb2_pl (int RI, int RJ, int RK) {
  //  cout<<RI<<" "<<RJ<<" "<<RK<<endl;
    for (unsigned int n=0; n<sb2Parameters.size (); n++) {
        if (sb2Parameters[n]->I == RI && sb2Parameters[n]->J == RJ && sb2Parameters[n]->K == RK){
            return n;
        }
    }
    for (unsigned int n=0; n<sb2Parameters.size (); n++) {
        if (sb2Parameters[n]->I == RK && sb2Parameters[n]->J == RJ && sb2Parameters[n]->I == RI){
            return n;
        }
    }

    cout <<"could not retrieve sb2 informations for rows "<<RI<<" , "<<RJ<<" && "<<RK<<endl;
    return -1;
}


int MMFF::get_to_pl (int type, int I, int J, int K, int L) {
    for (unsigned int n=0; n<toParameters.size (); n++) {
        if (toParameters[n]->I == I && toParameters[n]->J == J && toParameters[n]->K == K && toParameters[n]->L == L && toParameters[n]->type==type){
            return n;
        }
    }
    for (unsigned int n=0; n<toParameters.size (); n++) {
        if (toParameters[n]->I == 0 && toParameters[n]->J == J && toParameters[n]->K == K && toParameters[n]->L == L && toParameters[n]->type==type){
            return n;
        }
    }
    for (unsigned int n=0; n<toParameters.size (); n++) {
        if (toParameters[n]->I == I && toParameters[n]->J == J && toParameters[n]->K == K && toParameters[n]->L == 0 && toParameters[n]->type==type){
            return n;
        }
    }
    for (unsigned int n=0; n<toParameters.size (); n++) {
        if (toParameters[n]->I == 0 && toParameters[n]->J == J && toParameters[n]->K == K && toParameters[n]->L == 0 && toParameters[n]->type==type){
            return n;
        }
    }
 //   cout <<"could not retrieve to informations for atom types "<<I<<" , "<<J<<" , "<<K<<" && "<<L<<endl;
    return -1;
}

int MMFF::get_op_pl (int I, int J, int K, int L) {
    for (unsigned int n=0; n<opParameters.size (); n++) {
        if (opParameters[n]->I == I && opParameters[n]->J == J && opParameters[n]->K == K && opParameters[n]->L == L){
            return n;
        }
    }
    for (unsigned int n=0; n<opParameters.size (); n++) {
        if (opParameters[n]->I == 0 && opParameters[n]->J == J && opParameters[n]->K == 0 && opParameters[n]->L == 0){
            return n;
        }
    }
    cout <<"could not retrieve op informations for atom types "<<I<<" , "<<J<<" , "<<K<<" && "<<L<<endl;
    return 0;
}


int MMFF::get_vw_pl (int I) {
    for (unsigned int n=0; n<vwParameters.size (); n++) {
        if (vwParameters[n]->I == I){
            return n;
        }
    }
    cout <<"could not retrieve vw informations for atom type "<<I<<endl;
    return 0;
}




bool MMFF::get_linear (int J) {
        for (unsigned int i=0; i< atomParameters.size (); i++) {
            if (atomParameters[i]->type == J) return atomParameters[i]->linear;
        }
        return false;
    }

bool MMFF::get_sbmb (int J) {
        for (unsigned int i=0; i< atomParameters.size (); i++) {
            if (atomParameters[i]->type == J) return atomParameters[i]->multipleBond;
        }
        return false;
    }


bool MMFF::get_arom (int J) {
        for (unsigned int i=0; i< atomParameters.size (); i++) {
            if (atomParameters[i]->type == J) return atomParameters[i]->aromatic;
        }
        return false;
    }




/*
double MMFF::get_bs_kb (int I, int J) {
    for (unsigned int n=0; n<bsParameters.size (); n++) {
//        cout<< bsParameters[n]->J <<" "<<bsParameters[n]->J<<endl;
        if (bsParameters[n]->I == I && bsParameters[n]->J == J) {
            return bsParameters[n]->kb;
        }
    }
    cout <<"could not retrieve bs Kb for atom types "<<I<<" && "<<J<<endl;
    return 0;
}

double MMFF::get_bs_r0 (int I, int J) {
    for (unsigned int n=0; n<bsParameters.size (); n++) {
        if (bsParameters[n]->I != I) continue;
        if (bsParameters[n]->J != J) continue;
        return bsParameters[n]->r0;
    }
}
*/


/*
double MMFF::get_ab_theta0 (int I, int J, int K) {
    for (unsigned int n=0; n<abParameters.size (); n++) {
        if (abParameters[n]->I == I && abParameters[n]->J == J && abParameters[n]->K != K) {
            return abParameters[n]->theta0;
        }
    }
    cout <<"could not retrieve ab theta0 for atom types "<<I<<" , "<<J<<" && "<<K<<endl;
    return 0;
}

double MMFF::get_ab_ka (int I, int J, int K) {
    for (unsigned int n=0; n<abParameters.size (); n++) {
        if (abParameters[n]->I == I && abParameters[n]->J == J && abParameters[n]->K != K) {
            return abParameters[n]->ka;
        }
    }
    cout <<"could not retrieve ab ka for atom types "<<I<<" , "<<J<<" && "<<K<<endl;
    return 0;
}

*/


void MMFF::empiric_to_parameters (Atom *at1, Atom *at2, Atom *at3, Atom *at4, float& v1, float &v2, float &v3) {
/*
    MMFFatomParameter *an1, *an2, *an3, *an4;
    
    an1 = atomParameters[get_atype_line (at1->MMFFtype)];
    an2 = atomParameters[get_atype_line (at2->MMFFtype)];
    an3 = atomParameters[get_atype_line (at3->MMFFtype)];
    an4 = atomParameters[get_atype_line (at4->MMFFtype)];
    ZNBond *cb = at2 -> GetBond (at3);
    bool arom;
    if (!cb) arom = false;
    else (arom = cb->IsAromatic ());
    bool bo;
    if (!cb) bo = 0;
    else (bo = cb->GetBO ());


    if (an2->linear || an3->linear) {
        v1 = 0;
        v2 = 0;
        v3 = 0;
    }
    else if ((an2->aromatic && an3->aromatic && arom) || (bo ==2)) {
        float b, p, Uj, Uk;
        if (!an2->lonePairs && !an3->lonePairs) p = 0.5;
        else p = 0.3;
        if (bo==2) p = 0.4;
        if (bo==2 && an2->multipleBond==2 && an3->multipleBond==2) p = 1;

        if ((an3->valence ==3 && an2->valence ==4) || (an2->valence ==3 && an3->valence ==4)) b = 3;
        else b=6;
        if (bo==2) b = 6;
        if (at2->GetAtomicNum () == 6 || at2->GetAtomicNum () == 7 || at2->GetAtomicNum () == 8) Uj = 2.0;
        else Uj = 1.25;
        if (at3->GetAtomicNum () == 6 || at3->GetAtomicNum () == 7 || at3->GetAtomicNum () == 8) Uj = 2.0;
        else Uk = 1.25;
        v1 = 0;
        v2 =  b*p*sqrt (Uj*Uk);
        v3 = 0;
    }
    else if (an2->crd ==4 || an3->crd == 4) {
        int n = 9;
        float Vj, Vk;
        bool no_v3 = false;
        if (an2->crd == 4 && an3->crd == 3 && ((an3->valence ==4 || an3->valence==34) || an3->multipleBond)) no_v3=true;
        else if (an3->crd == 4 && an2->crd == 3 && ((an2->valence ==4 || an2->valence==34) || an2->multipleBond)) no_v3=true;
        else if (an2->crd ==4 && an3->crd==2 && (an3->valence==3 || an3->multipleBond)) no_v3=true;
        else if (an3->crd ==4 && an2->crd==2 && (an2->valence==3 || an2->multipleBond)) no_v3=true;

        if (at2->GetAtomicNum () ==6) Vj = 2.12;
        else if (at2->GetAtomicNum () ==7) Vj = 1.5;
        else if (at2->GetAtomicNum () ==8) Vj = 0.2;
        else if (at2->GetAtomicNum () ==14) Vj = 1.22;
        else if (at2->GetAtomicNum () ==15) Vj = 2.4;
        else  Vj = 0.49;

        if (at3->GetAtomicNum () ==6) Vk = 2.12;
        else if (at3->GetAtomicNum () ==7) Vk = 1.5;
        else if (at3->GetAtomicNum () ==8) Vk = 0.2;
        else if (at3->GetAtomicNum () ==14) Vk = 1.22;
        else if (at3->GetAtomicNum () ==15) Vk = 2.4;
        else  Vk = 0.49;
        v1=0;
        v2 = 0;
        if (!no_v3) v3 = sqrt (Vj*Vk)/n;
        else v3 = 0;
    }
    else if (bo ==1 && ((an2->multipleBond && an3->multipleBond) || (an2->multipleBond && an3->lonePairs) || (an3->multipleBond && an2->lonePairs))) {
        bool no_v = false;
        float b, p, Uj, Uk;
        if (an2->lonePairs && an3->lonePairs) no_v = true;
        else if ((an2->lonePairs && an3->multipleBond) || (an3->lonePairs && an2->multipleBond) ) {
            b=6;
            if (an2->multipleBond && an3->multipleBond) p = 0.5;
            else if (at2->GetAtomicNum ()>8 || at3->GetAtomicNum ()>8) p = 0.15;
            else p=0.3;
        }

        else if ((an2->multipleBond || an3->multipleBond) && (at1->GetAtomicNum ()!=6 || at2->GetAtomicNum ()!=6)) {
            b = 6;
            p = 0.4;
        }
        else{
            b =6;
            p = 0.15;
        }
        if (at2->GetAtomicNum () == 6 || at2->GetAtomicNum () == 7 || at2->GetAtomicNum () == 8) Uj = 2.0;
        else Uj = 1.25;
        if (at3->GetAtomicNum () == 6 || at3->GetAtomicNum () == 7 || at3->GetAtomicNum () == 8) Uj = 2.0;
        else Uk = 1.25;
        v1 = 0;
        if (!no_v) v2 = b*p*sqrt (Uj*Uk);
        else v2 =0;
        v3 = 0;
                
    }
    else {
        if ((at2->GetAtomicNum () == 8 || at2->GetAtomicNum () == 16) && (at3->GetAtomicNum () == 8 || at3->GetAtomicNum () == 16)) {
            float Wj, Wk;
            Wj = 8.;
            Wk = 8.;
            if (at2->GetAtomicNum () ==8) Wj = 2.;
            if (at3->GetAtomicNum () ==8) Wk = 2.;
            v1 =0;
            v2=-sqrt(Wj*Wk);
            v3 = 0;
        }
        else {
            float Vj, Vk;
            int n = (an2->crd-1)*(an3->crd-1);
            if (at2->GetAtomicNum () ==6) Vj = 2.12;
            else if (at2->GetAtomicNum () ==7) Vj = 1.5;
            else if (at2->GetAtomicNum () ==8) Vj = 0.2;
            else if (at2->GetAtomicNum () ==14) Vj = 1.22;
            else if (at2->GetAtomicNum () ==15) Vj = 2.4;
            else  Vj = 0.49;
    
            if (at3->GetAtomicNum () ==6) Vk = 2.12;
            else if (at3->GetAtomicNum () ==7) Vk = 1.5;
            else if (at3->GetAtomicNum () ==8) Vk = 0.2;
            else if (at3->GetAtomicNum () ==14) Vk = 1.22;
            else if (at3->GetAtomicNum () ==15) Vk = 2.4;
            else  Vk = 0.49;
            v1=0;
            v2 = 0;
            v3 = sqrt (Vj*Vk)/n;
            
        }
    }
*/
}







double MMFF::get_vdw_e (Atom *at1, Atom *at2) {

    int I = get_MMFFtype (at1);
    int J = get_MMFFtype (at2);
    int pli = get_vw_pl (I);
    int plj = get_vw_pl (J);

    double Gi =  vwParameters [pli] -> G;
    double Gj =  vwParameters [plj] -> G;
    double ai =  vwParameters [pli] -> alpha;
    double aj =  vwParameters [plj] -> alpha;
    double Ni =  vwParameters [pli] -> N;
    double Nj =  vwParameters [plj] -> N;
    double r0 = get_vdw_r0 (at1, at2);
    return (181.16*Gi*Gj*ai*aj) / ((pow((ai/Ni),0.5)+pow((aj/Nj),0.5))*pow(r0,6));

}


double MMFF::get_vdw_r0 (Atom *at1, Atom *at2) {

    double DARAD = 0.8;
    double B = 0.2;
    double beta = 12.;
    int I = get_MMFFtype (at1);
    int J = get_MMFFtype (at2);
    int pli = get_vw_pl (I);
    if ((at1->IsPolarHydrogen ()) || (at2->IsPolarHydrogen ())) B = 0; 
    double Ai = vwParameters [pli] -> A;
    double ai = vwParameters [pli] -> alpha;
    double r0ii = Ai*pow(ai,0.25);
    double scale = 1;

    if (I == J) return r0ii;
    else {
        int plj = get_vw_pl (J);
        if ((vwParameters[pli]->DA ==ACCEPTOR && vwParameters[plj]->DA ==DONOR) ||  (vwParameters[plj]->DA ==ACCEPTOR && vwParameters[pli]->DA ==DONOR)) scale = DARAD;
        double Aj = vwParameters [plj] -> A;
        double aj = vwParameters [plj] -> alpha;
        double r0jj = Aj*pow (aj,0.25);
        double gamma = (r0ii-r0jj) / (r0ii + r0jj);
        return 0.5*(r0ii + r0jj)*(1 + B* (1- exp (- beta * gamma*gamma)))*scale;
    } 
}

int MMFF::get_pt_row (int an) {
    if (an<2) return 0;
    else if (an < 12) return 1;
    else if (an <20) return 2;
    else if (an <38) return 3;
    else return 4;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MMFF::load_internal_interactions (vector <ForceFieldInteraction *> *vect = 0) {
    clear_internal_interactions ();
    ZNMolecule *mol = target_mol;

//bond stretchings
    FOR_BONDS_OF_MOL (bond, mol) {
        MMFFbsInteraction *bsint = new MMFFbsInteraction;
        int I = get_MMFFtype (bond->GetBeginAtom ());
        int J = get_MMFFtype (bond->GetEndAtom ());
        if (J<I) {
            bsint->at1 = bond->GetEndAtom ();
            bsint->at2 = bond->GetBeginAtom ();        
        }
        else {
            bsint->at1 = bond->GetBeginAtom ();
            bsint->at2 = bond->GetEndAtom ();   
        }
        I = get_MMFFtype (bsint->at1);
        J = get_MMFFtype (bsint->at2);
		bsint->type = getBondType (bsint ->at1, bsint ->at2);
        int pl = get_bs_pl (bsint->type, I, J);
        bsint->kb = bsParameters[pl] -> kb;
        bsint->r0 = bsParameters[pl] -> r0;
		if (vect) vect -> push_back (bsint);

        else bsInteractions.push_back (bsint);
					//				cerr << "adding BS interaction "<<endl;
        
    }

    

//angle interactions
    FOR_ATOMS_OF_MOL (a, mol) {
        if (CountBonds (&*a) >1) {
            FOR_NBORS_OF_ATOM (n1, &*a) {
                FOR_NBORS_OF_ATOM (n2, &*a) {
					if (n1 -> GetIdx () < n2 -> GetIdx ()) {
						MMFFabInteraction *abint = new MMFFabInteraction;
						abint->at2 = &*a;
						ZNBond *bondij, *bondjk;
						if (get_MMFFtype (&*n1)  < get_MMFFtype (&*n2)) {
							abint->at1 = &*n1;
							abint->at3 = &*n2;
							bondij = a -> GetBond (&*n1);   
							bondjk = a -> GetBond (&*n2);
						}
                    else {
							abint->at1 = &*n2;
							abint->at3 = &*n1;
							bondij = a -> GetBond (&*n2);   
							bondjk = a -> GetBond (&*n1);
                    }
                    abint->type = getAngleType (abint ->at1, abint ->at2, abint ->at3);


					int I = get_MMFFtype (abint->at1);
					int J = get_MMFFtype (abint->at2);
					int K = get_MMFFtype (abint->at3);
					
					
                    int abpl = get_ab_pl (abint->type, I, J, K);
                    abint->ka = abParameters[abpl]->ka;
                    abint->theta0 = abParameters[abpl]->theta0;                   
                    abint->linear = get_linear (get_MMFFtype (abint->at2));
					
					if (vect) vect -> push_back (abint);
                    else abInteractions.push_back (abint);
		
//stretch-bend					
					MMFFsbInteraction *sbint = new MMFFsbInteraction;
					sbint->at1 =abint->at1;
					sbint->at2 =abint->at2;
					sbint->at3 =abint->at3;

					sbint->theta0 = abint->theta0;
			//		cerr << "adding AB interaction "<<endl;


					
                    int sbpl = get_sb_pl (I, J, K); 
                    if (sbpl != -1) { 
                        sbint->kijk = sbParameters[sbpl]->kijk;
                        sbint->kkji = sbParameters[sbpl]->kkji;

                    }
                    else {
                        int sb2pl = get_sb2_pl (get_pt_row (abint->at1->GetAtomicNum ()), get_pt_row (abint->at2->GetAtomicNum ()),get_pt_row (abint->at3->GetAtomicNum ()));
						if (sb2pl == -1) sb2pl = 1;
                        sbint->kijk = sb2Parameters[sb2pl]->kijk;
                        sbint->kkji = sb2Parameters[sb2pl]->kkji;
                    }
                    int typeij, typejk;
                    if (bondij->GetBO () ==1) typeij = get_bond_type  (I, J);
                    else typeij = 0;
                    if (bondjk->GetBO () ==1) typejk = get_bond_type (K, J);
                    else typejk = 0;

                    int ijbspl=0, kjbspl=0;
                    if (J < I) ijbspl = get_bs_pl (typeij, J, I);
                    else ijbspl = get_bs_pl (typeij, I, J);

                    if (K < J) kjbspl = get_bs_pl (typeij, K, J);
                    else kjbspl = get_bs_pl (typejk, J, K);
                      
                    sbint->r0ij = bsParameters[ijbspl]->r0;
                    sbint->r0kj = bsParameters[kjbspl]->r0;
					
					if (vect) vect -> push_back (sbint);
                    else sbInteractions.push_back (sbint);                    
                }
				}
            }
        }
    } 
	
//OOP interactions
    Atom *at1, *at2, *at3, *at4;
	FOR_ATOMS_OF_MOL (a, mol) {
        if (CountBonds (&*a) == 3) {
			Atom *n1= NULL, *n2= NULL, *n3= NULL;
			int nn =0;
			FOR_NBORS_OF_ATOM (n, &*a) {
				nn++;
				if (nn == 1) n1 = &*n;
				else if (nn ==2) n2 = &*n;
				else n3 = &*n;
			} 
			assert (n1); assert (n2); assert (n3);
			int t1 = get_MMFFtype (n1);
			int t2 = get_MMFFtype (n2);
			int t3 = get_MMFFtype (n3);
			int ta = get_MMFFtype (&*a);
            if (t1 >= t2 && t1 >= t3) {
                if (t2 > t3) {
                    at1 = n3;          
                    at2 = &*a;
                    at3 = n2;
                    at4 = n1; 
                }
                else {
                    at1 = n2;          
                    at2 = &*a;
                    at3 = n3;
                    at4 = n1;                 
                }
            }
            else if (t2 >= t1 && t2 >= t3) {
                if (t1 > t3) {
                    at1 = n3;          
                    at2 = &*a;
                    at3 = n1;
                    at4 = n2; 
                }
                else {
                    at1 = n1;          
                    at2 = &*a;
                    at3 = n3;
                    at4 = n2;                 
                }

            }
            else {
                if (t1 > t2) {
                    at1 = n2;          
                    at2 = &*a;
                    at3 = n1;
                    at4 = n3; 
                }
                else {
                    at1 = n1;          
                    at2 = &*a;
                    at3 = n2;
                    at4 = n3;                 
                }

            }
            int pl = get_op_pl (get_MMFFtype (at1), get_MMFFtype (at2),get_MMFFtype (at3),get_MMFFtype (at4));
            double koop = opParameters[pl]->koop;
            MMFFopInteraction *opint = new MMFFopInteraction;
            opint->at1 = at1;          
            opint->at2 = at2;
            opint->at3 = at3;
            opint->at4 = at4;
            opint->koop = koop;
			if (vect) vect -> push_back (opint);
			else opInteractions.push_back (opint);
            
            opint = new MMFFopInteraction;
            opint->at1 = at1;          
            opint->at2 = at2;
            opint->at3 = at4;
            opint->at4 = at3;
            opint->koop = koop;
			if (vect) vect -> push_back (opint);
			else opInteractions.push_back (opint);
            
            opint = new MMFFopInteraction;
            opint->at1 = at3;          
            opint->at2 = at2;
            opint->at3 = at4;
            opint->at4 = at1;
            opint->koop = koop;
			if (vect) vect -> push_back (opint);
			else opInteractions.push_back (opint);
        }
    } 
//torsions 


    FOR_TORSIONS_OF_MOL(t, mol) {
		MMFFtoInteraction *toint = new MMFFtoInteraction;
		Atom *a, *b, *c, *d;
        a = mol -> GetAtom((*t)[0] + 1);
        b = mol -> GetAtom((*t)[1] + 1);
        c = mol -> GetAtom((*t)[2] + 1);
        d = mol -> GetAtom((*t)[3] + 1);
		int atn1 = get_MMFFtype (a);
		int atn2 = get_MMFFtype (b);
		int atn3 = get_MMFFtype (c);
		int atn4 = get_MMFFtype (d);

		if (atn3==atn2) {
			if (atn4 < atn1) {
				toint->at1 = d;
				toint->at4 = a;
				toint->at2 = c;
				toint->at3 = b;
			}
			else {
				toint->at1 = a;
				toint->at4 = d;
				toint->at2 = b;
				toint->at3 = c;
			}
		}
		else if (atn2<atn3) {
			toint->at1 = a;
			toint->at4 = d;
			toint->at2 = b;
			toint->at3 = c;
		}
		else {
			toint->at1 = d;
			toint->at4 = a;
			toint->at2 = c;
			toint->at3 = b;                            
		}
		int type = 0;
		//mancano tutti i torsion type
		toint->type = type;
		int pl = get_to_pl (toint->type, get_MMFFtype (toint->at1),  get_MMFFtype (toint->at2),  get_MMFFtype (toint->at3),  get_MMFFtype (toint->at4)); //parameters line
		float V1, V2, V3;
		if (pl!=-1) {
			V1 = toParameters[pl]->V1;
			V2 = toParameters[pl]->V2;
			V3 = toParameters[pl]->V3;
		}
		else {
			empiric_to_parameters (toint->at1,toint->at2,toint->at3,toint->at4,V1, V2, V3);
		}
		toint->v1 = V1;
		toint->v2 = V2;
		toint->v3 = V3;
		if (vect) vect -> push_back (toint);
		else toInteractions.push_back (toint);
	
    } 

//nonbonded interactions
    FOR_ATOMS_OF_MOL (a1, mol) {
        FOR_ATOMS_OF_MOL (a2, mol) {
			if (a1->GetIdx () < a2 -> GetIdx ()) {

				unsigned int distance = 5;
				Atom *at1 = &*a1;
				Atom *at2 = &*a2;

				//chech for 1-2 interactions
				FOR_NBORS_OF_ATOM (n1, at1)  {
					if (&*n1 == at2) {
						distance = 2;
					}
					//check for 1-3 interactions
					FOR_NBORS_OF_ATOM (n2, &*n1) {
						if (&*n2 == at2 && distance>3) {
							distance = 3;
						}
						FOR_NBORS_OF_ATOM (n3, &*n2)  {
							if (&*n3 ==at2 && distance>4) {
								distance = 4;
							}               
						}
					}       
				}

				if (distance >3) {
					MMFFvwInteraction *vwint = new MMFFvwInteraction;
					vwint -> at1 = at1;
					vwint -> at2 = at2;
					vwint -> e = get_vdw_e (vwint->at1, vwint->at2);
					vwint -> r0 = get_vdw_r0 (vwint->at1, vwint->at2);
					if (vect) vect -> push_back (vwint);
					else vwInteractions.push_back (vwint);

					MMFFelInteraction *elint = new MMFFelInteraction;
					elint->at1 = at1;
					elint->at2 = at2;
					if (distance ==4) elint->scale = 0.75;
					else elint->scale =1;
					if (vect) vect -> push_back (elint);
					else elInteractions.push_back (elint);              
			
			     }  
			}
        }
    }

}

int MMFF::load_aromatic_parameters () {

    QFile file (":parameters/MMFFAROM.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read aromatic parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));
	    MMFFaromaticParameter *param = new MMFFaromaticParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv    >> param->old_type             ;
	        liness  >> param->arom_type            ;
	        liness  >> param->at_number            ;
	        liness  >> param->ring_size            ;
	        liness  >> param->L5                   ;
	        liness  >> param->IMCAT                ;
	        liness  >> param->N5ANION              ;
	        aromaticParameters.push_back(param)    ;
        }
    }
    if (aromaticParameters.size ()) return 1;
    return 0;
}



int MMFF::load_ci_parameters () {

    QFile file (":parameters/MMFFCHG.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read ci parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));
	    MMFFciParameter *param = new MMFFciParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv >> param->type            ;
	        liness  >> param->I            ;
	        liness  >> param->J            ;
	        liness  >> param->increment            ;
	        ciParameters.push_back(param);
        }
    }
    if (ciParameters.size ()) return 1;
    return 0;
}

int MMFF::load_H_parameters () {

    QFile file (":parameters/MMFFHDEF.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read H parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));
	    MMFFHParameter *param = new MMFFHParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv >> param->parent           ;
	        liness  >> param->strin ;
  	        HParameters.push_back(param);

        }
    }
    if (HParameters.size ()) return 1;
    return 0;
}


int MMFF::load_atype_parameters () {

    QFile file (":parameters/MMFFSYMB.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read atype parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));
	    MMFFatypeParameter *param = new MMFFatypeParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv >> param->strin           ;
	        liness  >> param->number ;
            param->description = liness.str();
  	        atypeParameters.push_back(param);

        }
    }
    if (atypeParameters.size ()) return 1;
    return 0; 
}




int MMFF::load_atom_parameters () {

    QFile file (":parameters/MMFFPROP.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read atom parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));
	    MMFFatomParameter *param = new MMFFatomParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv >> param->type            ;
	        liness  >> param->atomicNumber ;
	        liness  >> param->crd          ;
	        liness  >> param->valence      ;
	        liness  >> param->lonePairs    ;
	        liness  >> param->mltb         ;
	        liness  >> param->aromatic     ;
	        liness  >> param->linear       ;
	        liness  >> param->multipleBond ;
	        atomParameters.push_back(param);

        }
    }
    if (atomParameters.size ()) return 1;
    return 0;
}


int MMFF::load_ab_parameters () {



    QFile file (":parameters/MMFFANG.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read ab parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));
	    MMFFabParameter *param = new MMFFabParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv    >> param->type       ;
	        liness  >> param->I          ;
	        liness  >> param->J          ;
	        liness  >> param->K          ;
	        liness  >> param->ka         ;
	        liness  >> param->theta0     ;
	        abParameters.push_back(param);
        }
    }
    if (abParameters.size ()) return 1;
    return 0;
}


int MMFF::load_bs_parameters () {


    QFile file (":parameters/MMFFBOND.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read bs parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));



	    MMFFbsParameter *param = new MMFFbsParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv >> param->type            ;
	        liness  >> param->I            ;
	        liness  >> param->J            ;
	        liness  >> param->kb           ;
	        liness  >> param->r0           ;
	        bsParameters.push_back(param)  ;
        }
    }
    if (bsParameters.size ()) return 1;
    return 0;
}



int MMFF::load_sb_parameters () {



    QFile file (":parameters/MMFFSTBN.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read sb parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));
	    MMFFsbParameter *param = new MMFFsbParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        liness  >> param->I            ;
	        liness  >> param->J            ;
            liness  >> param->K            ;
	        liness  >> param->kijk         ;
	        liness  >> param->kkji         ;
	        sbParameters.push_back(param)  ;
        }
    }
    if (sbParameters.size ()) return 1;
    return 0;
}

int MMFF::load_sb2_parameters () {


    QFile file (":parameters/MMFFDFSB.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read sb2 parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));

	    MMFFsbParameter *param = new MMFFsbParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv    >> param->I            ;
	        liness  >> param->J            ;
            liness  >> param->K            ;
	        liness  >> param->kijk         ;
	        liness  >> param->kkji         ;
	        sb2Parameters.push_back(param)  ;
        }
    }
    if (sb2Parameters.size ()) return 1;
    return 0;
}

int MMFF::load_to_parameters () {

    QFile file (":parameters/MMFFTOR.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read to parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));


	    MMFFtoParameter *param = new MMFFtoParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv >> param->type            ;
	        liness  >> param->I            ;
	        liness  >> param->J            ;
            liness  >> param->K            ;
	        liness  >> param->L            ;
	        liness  >> param->V1           ;
	        liness  >> param->V2           ;
	        liness  >> param->V3           ;
	        toParameters.push_back(param)  ;
        }
    }
    if (toParameters.size ()) return 1;
    return 0;

}


int MMFF::load_vw_parameters () {

    QFile file (":parameters/MMFFVDW.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read vw parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));
	    MMFFvwParameter *param = new MMFFvwParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv >> param->I               ;
	        liness  >> param->alpha        ;
	        liness  >> param->N            ;
            liness  >> param->A            ;
	        liness  >> param->G            ;
            char da;
            liness  >> da;
            if (da == 'A') param->DA = ACCEPTOR;
            else if (da == 'D') param->DA = DONOR;
            else param->DA = NEITHER;
	        vwParameters.push_back(param)  ;
        }
    }
    if (vwParameters.size ()) return 1;
    return 0;
}

int MMFF::load_op_parameters () {

    QFile file (":parameters/MMFFOOP.PAR");
    if (!file.open(QFile::ReadOnly | QIODevice::Text)) {
        cerr << "unable to read op parameters"<<endl;
        return 0;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine();
    	istringstream liness (string (line.toStdString ()));
	    MMFFopParameter *param = new MMFFopParameter;
        string first;
        liness >> first;
        if (first != "*") {
            istringstream conv (first);
	        conv    >> param -> I            ;
	        liness  >> param -> J            ;
	        liness  >> param -> K            ;
            liness  >> param -> L            ;
	        liness  >> param -> koop         ;
	        opParameters.push_back(param)    ;
        }
    }
    if (opParameters.size ()) return 1;
    return 0;
}

void MMFF::clear_internal_interactions () {
    for (unsigned int i=0; i<vwInteractions.size(); i++) delete vwInteractions[i];
    for (unsigned int i=0; i<elInteractions.size(); i++) delete elInteractions[i];
    for (unsigned int i=0; i<bsInteractions.size(); i++) delete bsInteractions[i];
    for (unsigned int i=0; i<abInteractions.size(); i++) delete abInteractions[i];
    for (unsigned int i=0; i<toInteractions.size(); i++) delete toInteractions[i];
    for (unsigned int i=0; i<sbInteractions.size(); i++) delete sbInteractions[i];
    for (unsigned int i=0; i<opInteractions.size(); i++) delete opInteractions[i];

    bsInteractions.clear ();
    abInteractions.clear ();
    toInteractions.clear ();
    vwInteractions.clear ();
    opInteractions.clear ();
    elInteractions.clear ();
    sbInteractions.clear ();

}



void MMFF::clear_nonbonded_interactions () {
    for (unsigned int i=0; i<vwNBInteractions.size (); i++) delete vwNBInteractions[i];
    for (unsigned int i=0; i<elNBInteractions.size (); i++) delete elNBInteractions[i];
    vwNBInteractions.clear ();
    elNBInteractions.clear ();
}



void MMFF::load_nonbonded_interactions_for_atom (Atom *a, queue <ForceFieldInteraction *> *queue = 0) {
	objectList<Atom*>* nbAtoms = far_grid ->getNeighborObjects(get_coordinates (&*a));
	if (nbAtoms) {
		vector <Atom *> neighbours = nbAtoms -> objects;      
		for (unsigned int j=0; j<neighbours.size (); j++) {
			assert (a != neighbours[j]);
			MMFFelInteraction *elint = new MMFFelInteraction;
			elint -> at1   = a;
			elint -> at2   = neighbours[j];
			elint -> scale = 1.f;
			if (queue) queue -> push (elint);
			else elNBInteractions.push_back (elint);    
		}        
	}
	objectList<Atom*>* nbAtoms2 = near_grid->getNeighborObjects(get_coordinates (&*a));
	if (nbAtoms2) {
		vector <Atom *> neighbours2 = nbAtoms2 -> objects;      
		for (unsigned int j=0; j<neighbours2.size (); j++) {
			assert (a != neighbours2[j]);
			MMFFvwInteraction *vwint = new MMFFvwInteraction;
			vwint->at1 = a;
			vwint->at2 = neighbours2[j];
			vwint -> e = get_vdw_e (vwint->at1, vwint->at2);
			vwint -> r0 = get_vdw_r0 (vwint->at1, vwint->at2);
			//        vwint->scale =1;
			if (queue) queue -> push (vwint);
			else vwNBInteractions.push_back (vwint);
		}        
	}
}




void MMFF::load_nonbonded_interactions () {

    clear_nonbonded_interactions ();
    assert (target_mol);
    ZNMolecule *mol = target_mol;
    FOR_ATOMS_OF_MOL (a, target_mol) {
		load_nonbonded_interactions_for_atom (&*a);
    }
}


//_________________________________________________________________________________________________
//_________________________________________________________________________________________________

double MMFF::compute_charge (Atom *at) {

    return std::numeric_limits<double>::signaling_NaN();

/*
    double q = at->MMFFstartcharge;
    FOR_NBORS_OF_ATOM (a, at) {
        q+= get_charge_increment (a->MMFFtype, at->MMFFtype);
    }
    return q;
*/
}

void MMFF::get_strings_mol (ZNMolecule*mol) {
/*
    FOR_ATOMS_OF_MOL (mol, a)   {
        a->MMFFstring = getMMFFcarbonstring (a);
    }
    FOR_ATOMS_OF_MOL (mol, a)   {
        a->MMFFstring = getMMFFeteroatomstring (a);
    }

    for (unsigned int r=0; r<mol->rings.size (); r++ )   {
        if (mol->rings[r]->aromatic) {
            for (unsigned int i=0; i<mol->rings[r]->atoms.size (); i++ )   {
                mol->rings[r]->atoms[i]->MMFFstring = get_aromatic_string (mol->rings[r]->atoms[i]);
            }
        }
    }
    for (unsigned int i=0; i<mol->atoms.size (); i++ )   {
        if (mol->atoms[i]->GetAtomicNum ()==1) mol->atoms[i]->MMFFstring = getMMFFHstring (mol->atoms[i]);
    }
*/
}

void MMFF::compute_partial_charges (ZNMolecule *mol){
/*
    FOR_ATOMS_OF_MOL (a, mol)   {
        a -> MMFFstartcharge = getMMFFstartcharge (&*a);
    }
    FOR_ATOMS_OF_MOL (a, mol)   {
        a -> SetPartialCharge (compute_charge (&*a));
    }
*/
}



void MMFF::initialize_mol (ZNMolecule *mol) {
/*
    get_strings_mol (mol);

    FOR_ATOMS_OF_MOL (a, mol)   {
        a -> MMFFtype = getMMFFtype (&*a);
    }

    compute_partial_charges (mol);
*/
}


double MMFF::getMMFFstartcharge (Atom *at) {

    return std::numeric_limits<double>::signaling_NaN();

/*
    double q = 0;
    switch (at->MMFFtype) {
        case 32:
            if (at->MMFFstring=="O4CL") q= -0.25;
            else if (at->MMFFstring=="OSMS") q= -0.5;
            break;
        case 55:
            q = 0.5;
            break;
        case 56:
            q = 1. /3.;
            break;
        case 58:
            q = 1.;
            break;
        case 61:
            if (at->bonded_to (-2, 7)) q = 1.;
            break;
        case 62:
            q = -1.;
            break;
        case 72:
            if (at->MMFFstring=="SM") q = -1.;
            break;
        case 81:
            if (at->MMFFstring=="NIM+") q = 0.5;
            else q = 1.;
            break;
        case 90:
            q = -1.;
            break;
        case 93:
            q = 1.;
            break;
    }
    return q;
*/
}


string MMFF::getMMFFHstring (Atom *at) {
    
    return "not implemented";

 /*   if (at->GetAtomicNum ()==1 && at->CountBonds ()) {
        string st = at->bound[0]->MMFFstring;
        for (unsigned int i=0; i<HParameters.size (); i++) {
            if (st==HParameters[i]->parent) return HParameters[i]->strin; 
        }
    }
    return "unknown H";
*/
}



string MMFF::getMMFFcarbonstring (Atom *at) {

    return "not implemented";

/*    bool in_ring = at->in_ring.size()>0;
    int ring_size = 1000;
    bool aromatic=false;
    for (unsigned int r=0; r<at->in_ring.size (); r++){
        if (at->in_ring[r]->atoms.size ()< ring_size) ring_size =at->in_ring[r]->atoms.size ();
//        if (at->in_ring[r]->aromatic) aromatic = true;
    }
    if (at->GetAtomicNum () == 6) {
// 1 2 3 4 20 22 30 37 41 57 60 63 64 78 80    C
//  CR      1  ALKYL CARBON, SP3
//  C=C     2  VINYLIC CARBON, SP2
//  CSP2    2  GENERIC SP2 CARBON
//  C=O     3  GENERAL CARBONYL CARBON
//  C=N     3  SP2 CARBON IN C=N
//  CGD     3  GUANIDINE CARBON, DOUBLY BONDED TO N
//  C=OR    3  KETONE || ALDEHYDE CARBONYL CARBON
//  C=ON    3  AMIDE CARBONYL CARBON
//  CONN    3  UREA CARBONYL CARBON
//  COO     3  CARBOXYLIC ACID || ESTER CARBONYL CARBON
//  COON    3  CARBAMATE CARBONYL CARBON
//  COOO    3  C ARBONIC ACID || ESTER CARBONYL CARBON
//  C=OS    3  THIOESTER CARBONYL CARBON, DOUBLE BONDED TO O
//  C=S     3  THIOESTER CARBON, DOUBLY BONDED TO S
//  C=SN    3  THIOAMIDE, CARBON, DOUBLY BONDED TO S
//  CSO2    3  CARBON IN >C=SO2
//  CS=O    3  CARBON IN >C=S=O (SULFINYL GROUP)
//  CSS     3  THIOCARBOXYLIC ACID || ESTER CARBONYL CARBON
//  C=P     3  CARBON DOUBLE BONDED TO PHOSPHOROUS
//  CSP     4  ACETYLENIC CARBON
//  =C=     4  ALLENIC CARBON
//  CR4R   20  CARBON IN 4-MEMBERED RINGS
//  CR3R   22  CARBON IN A 3-MEMBERED RING
//  CE4R   30  OLEFINIC CARBON IN 4-MEMBERED RINGS
//  CB     37  CARBON AS IN BENZENE, PYRROLE
//  CO2M   41  CARBOXYLATE ANION CARBON
//  CS2M   41  CARBON IN THIOCARBOXYLATE ANION
//  CGD+   57  GUANIDINIUM CARBON
//  CNN+   57  C IN +N=C-N RESONANCE STRUCTURES
//  C%     60  ISONITRILE CARBON
//  C5A    63  ALPHA CARBON IN 5-MEMBERED HETEROAROMATIC RING
//  C5B    64  BETA CARBON IN 5-MEMBERED HETEROAROMATIC RING
//  C5     78  GENERAL CARBON IN 5-MEMBERED HETEROAROMATIC RING
//  CIM+   80  C IN N-C-N IN IMIDAZOLIUM ION

        int dU = 4-at->bonds.size ();
        //general carbon
        if ((at->bound.size ()==1) && at->bonded_to  (3, 7)) return "C%";
        if (at->bonded_to (3, -2)) return "CSP";
        if (at->bonded_to (2, 15)) return "C=P";

        if (at->bonded_to (-2, 7)==3 && dU ==1) {
            bool g = true;            
            for (unsigned int i=0; i<at->bound.size(); i++) {
                if (at->bound[i]->bound.size ()!=3)  g=false;
            }
            if (g) return "CGD+";
            g = true;
            for (unsigned int i=0; i<at->bound.size(); i++) {
                if (at->bound[i]->bound.size ()!=3 && at->bonds[i]->mol2Type==1 || at->bound[i]->bound.size ()!=2 && at->bonds[i]->mol2Type==2 )  g=false;
            }
            if (g) return "CGD";

        }
        bool IM = false;
        for (unsigned int r=0; r<at->in_ring.size (); r++) {
            if (at->in_ring[r]->IM_CAT) IM=true;
        }

        if (dU ==1 && at->bonded_to (2,7)==1 && at->bonded_to (1,7) ) {
            bool res = true;
            bool ar1, ar2;
            ar1 = ar2 = false;
            for (unsigned int i=0; i<at->bound.size(); i++) {
                if (at->bound[i]->GetAtomicNum ()==7 && (at->bonds[i]->kekule==2 && at->bound[i]->bound.size ()<3)) res=false;
                if (at->bound[i]->GetAtomicNum ()==7 && at->bonds[i]->kekule==2 && at->bonds[i]->mol2Type==5) ar1=true;
                if (at->bound[i]->GetAtomicNum ()==7 && at->bonds[i]->mol2Type==1 && !at->bonds[i]->in_ring.size()) ar2=true;            
            }    
            if (ar1 && ar2) res=false;
            if (res) return "CNN+";
        } 
        if (at->bonded_to (-2, 8)==2) {
            bool carboxilate = true;
            for (unsigned int i=0; i<at->bound.size (); i++) {
                if (at->bound[i]->GetAtomicNum ()==8 && at->bound[i]->bound.size()>1) carboxilate=false;
            }
            if (carboxilate) return "CO2M";
        }
        if ((at->bonded_to (-2, 16))==2) {
            bool carboxilate = true;
            for (unsigned int i=0; i<at->bound.size (); i++) {
                if ((at->bound[i]->GetAtomicNum ()==16) && at->bound[i]->bound.size()>1) carboxilate=false;
            }
            if (carboxilate) return "CS2M";
        }


        if (at->bonded_to (2, -2) == 2) return "=C=";
        if (at->bonded_to (2, 16)) {
            if (at->bonded_to (-2, 7)) return "C=SN";
            if (at->bonded_to (1, 16)) return "CSS";
            for (unsigned int i=0; i<at->bound.size(); i++) {
                if (at->bound[i]->bound.size ()==2 && at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (2, 8) && at->bound[i]->bonded_to (2, 6)) return "CS=O";                
                if (at->bound[i]->bound.size ()==3 && at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (-2, 8)==2 && at->bound[i]->bonded_to (2, 6)) return "CSO2";
            }
            if (at->bonded_to (1, 8)) return "C=S";
        }


        if (at->bonded_to (2, 8)) {
            if (at->bonded_to (1, 7) ==2 || at->bonded_to (4, 7)==2) return "CONN";
            if (at->bonded_to (-2, 7)) return "C=ON";
            if (at->bonded_to (1, 16)) return "C=OS";
            if (at->bonded_to (-2, 8)==2){
                if (at->bonded_to (-2, 7)) {
                    int termO = 0;
                    for (unsigned int i=0; i<at->bound.size (); i++) {
                        if (at->bound[i]->GetAtomicNum ()==8 && at->bound[i]->bound.size ()==1) termO+=1;
                    }
                    if (termO ==2) return "COON";

                }
                return "COO";
            }
            if (at->bonded_to (-2, 8)>2){
                return "COOO";
            }            

            if (at->bonded_to (1, 6)) return "C=OR";

            return "C=O";
        }
        if (dU ==1 && at->bonded_to (2, 7)) return "C=N";
        if (in_ring) { //general ring

            if (aromatic) { //aromatic ring
                bool boo = false;
                for (unsigned int r=0; r<at->in_ring.size(); r++) {
                    if (at->in_ring[r]->atoms.size()==5 && at->in_ring[r]->aromatic) {
 //                       Ring *ring = at->in_ring[r];
                        for (unsigned int a=0; a<at->in_ring[r]->atoms.size (); a++) {
                            if (at->in_ring[r]->atoms[a]->GetAtomicNum ()!=6 ) {
                                if ((at->in_ring[r]->atoms[a]->GetAtomicNum ()==7) && (at->in_ring[r]->atoms[a]->bound.size ()!=2)) {}
                                else boo = true;
                            }
                        }
                    }
                }
                if (boo) { //carbon in 5 membered eteroaromatic ring
                    if (at->in_ring.size() ==1 && at->bonded_to (-2, 7)==2) {
                        if (!at->in_ring[0]->count (8) && !at->in_ring[0]->count (16))  return "CIM+";      
                    }
                    
                    for (unsigned int a=0; a<at->bound.size(); a++) {
                        if (at->bound[a]->GetAtomicNum ()!=6 && at->bound[a]->in_ring.size()) {
                            if  (at->bound[a]->GetAtomicNum ()==7 && !(at->bound[a]->bound.size ()==3)) {}
                            else return "C5A";
                        }
                    }
                    return "C5B";
                    return "C5";
                }
                return "CB";
            }
            else { //non aromatic ring
                if (ring_size == 3 && dU==0) return "CR3R";
                else if (ring_size == 4) {
                    if (dU>0) return "CE4R";
                    else return "CR4R";
                }
            }
        }
        else { //non ring carbon             
        }
        //general carbons with no other type

        if (at->bonded_to (2, 6) && dU ==1) return "C=C";
        if (dU ==0) return "CR";
        if (dU == 1) return "CSP2";

    }
    else return "eteroatom";
    return "Unknown C";
*/
}

string MMFF::getMMFFeteroatomstring (Atom *at) {

    return "not implemented";

/*
//    bool in_ring = at->in_ring.size()>0;
    int ring_size = 1000;
    bool aromatic=false;
    for (unsigned int r=0; r<at->in_ring.size (); r++){
        if (at->in_ring[r]->atoms.size ()< ring_size) ring_size =at->in_ring[r]->atoms.size ();
     //   if (at->in_ring[r]->aromatic) aromatic = true;
    }
    if (at->GetAtomicNum () == 3) {// 92    Li
//  LI+    92  LITHIUM CATION
        return "LI+";
    }

    if (at->GetAtomicNum () == 6) return at->MMFFstring;
    if (at->GetAtomicNum () == 7) {// 8 9 10 34 38 39 40 42 43 45 46 47 48 53 54 55 56 58 61 62 65 66 67 68 69 76 79 81 82   N
//  NR      8  NITROGEN IN ALIPHATIC AMINES
//  N=C     9  NITROGEN IN IMINES
//  N=N     9  NITROGEN IN AZO COMPOUNDS
//  NC=O   10  NITROGEN IN AMIDES
//  NC=S   10  NITROGEN IN N-C=S, THIOAMIDE
//  NN=C   10  NITROGEN IN N-N=C
//  NN=N   10  NITROGEN IN N-N=N    
//  NR+    34  QUATERNARY NITROGEN, SP3, POSITIVELY CHARGED
//  NPYD   38  NITROGEN, AS IN PYRIDINE
//  NPYL   39  NITROGEN, AS IN PYRROLE
//  NC=C   40  NITROGEN ON N-C=C
//  NC=N   40  NITROGEN IN N-C=N
//  NC=P   40  NITROGEN IN N-C=P
//  NC%C   40  NITROGEN ATTACHED TO C-C TRIPLE BOND
//  NSP    42  NITROGEN, TRIPLE BONDED
//  NSO2   43  NITROGEN IN SULFONAMIDES
//  NSO3   43  NITROGEN IN SULFONAMIDES, THREE O'S ON S
//  NPO2   43  NITROGEN IN PHOSPHONAMIDES
//  NPO3   43  NITROGEN IN PHOSPHONAMIDES, THREE O'S ON P
//  NC%N   43  NITROGEN ATTACHED TO CYANO GROUP
//  NO2    45  NITRO GROUP NITROGEN
//  NO3    45  NITRATE GROUP NITROGEN
//  N=O    46  NITROSO NITROGEN
//  NAZT   47  TERMINAL NITROGEN IN AZIDO || DIAZO GROUP
//  NSO    48  DIVALENT NITROGEN REPLACING MONOVALENT O IN SO2 GROUP
//  =N=    53  NITROGEN IN C=N=N || -N=N=N 
//  N+=C   54  POSITIVELY CHARGED IMINIUM NITROGEN
//  N+=N   54  POSITIVELY CHARGED NITROGEN DOUBLE-BONDED TO N
//  NCN+   55  N IN +N=C-N RESONANCE STRUCTURES - FORMAL CHARGE=1/2
//  NGD+   56  GUANIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1/3
//  NPD+   58  PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
  NR%    61  ISONITRILE NITROGEN [FC = 0] || DIAZO NITROGEN [FC = 1]
//  NM     62  DEPROTONATED SULFONAMIDE N-; FORMAL CHARGE=-1
//  N5A    65  ALPHA AROM HETEROCYCLIC 5-RING  NITROGEN
//  N5B    66  BETA AROM HETEROCYCLIC 5-RING  NITROGEN
//  N2OX   67  SP2-HYDRIDIZED N-OXIDE NITROGEN
//  N3OX   68  SP3-HYDRIDIZED N-OXIDE NITROGEN
//  NPOX   69  PYRIDINE N-OXIDE NITROGEN
  N5M    76  NEGATIVELY CHARGED N IN, E.G, TRI- || TETRAZOLE ANION
  N5     79  GENERAL NITROGEN IN 5-MEMBERED HETEROCYCLIC RING
  NIM+   81  IMIDAZOLIUM-TYPE NITROGEN - FORMAL CHARGE=1/2
  N5A+   81  POSITIVE N5A NITROGEN - FORMAL CHARGE=1
  N5B+   81  POSITIVE N5B NITROGEN - FORMAL CHARGE=1
  N5+    81  POSITIVE N5 NITROGEN - FORMAL CHARGE=1
  N5AX   82  N-OXIDE NITROGEN IN 5-RING ALPHA POSITION
  N5BX   82  N-OXIDE NITROGEN IN 5-RING BETA POSITION
  N5OX   82  N-OXIDE NITROGEN IN GENERAL 5-RING POSITION 

        int dU = 3-at->bonds.size ();
        for (unsigned int i=0; i<at->bound.size(); i++) {
            if (at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (2, 7) && at->bound[i]->bonded_to (-2, 8) && at->bonds[i]->kekule==2) return "NSO";
        }
        if (dU == 1 && at->bonded_to (2, -2)==2) return "=N=";
        if (dU == 2 && at->bound[0]->bonded_to (2, -2)==2 && at->bound[0]->GetAtomicNum ()==7) return "NAZT";
        for (unsigned int r=0; r<at->in_ring.size ();r++) {
            if (at->in_ring[r]->atoms.size () ==5 && at->in_ring[r]->N5ANION) return "NM";
        }
        if (at->bonded_to ("CGD+")) return "NGD+";
        if (at->bonded_to ("CNN+") && !(at->bound.size ()==2)) return "NCN+";


        if (at->bonded_to ("C=P")) return "NC=P";


        if (dU == 1 && at->bonded_to (2, 6)) return "N=C";

        if (dU == 1 && at->bonded_to (2, 7)) return "N=N";
        if (at->bonded_to ("C=SN") || at->bonded_to ("C=S") || at->bonded_to ("CS2M")) return "NC=S";
        if (dU == 1 && at->bonded_to (2, 8)) return "N=O";
        for (unsigned int i=0; i<at->bound.size(); i++) {
            if (at->bound[i]->GetAtomicNum ()==8 && at->bonds[i]->kekule==1 && at->bound[i]->bound.size ()==1 && dU==-1) return "N3OX"; 
            if (at->bound[i]->GetAtomicNum ()==8 && at->bonds[i]->kekule==1 && at->bound[i]->bound.size ()==1 && dU==0 && (at->bonded_to (2, -2)+at->bonded_to (5, -2)) && at->bonded_to (-2, 8)<2) return "N2OX";
            if (at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (-2, 8)==2 ) {
                if (dU ==0 || at->bonded_to (2, -2)) return "NSO2";
                if (dU==1) return "NM";        
            };
        }


        if (dU == 0 && (at->bonded_to (2, 6) || (at->bonded_to (5, -2)==2) && at->is_aromatic ())) return "N+=C";


        if (dU==-1) return "NR+";
        if (dU==0 && at->bonded_to (-2, 8)==2) return "NO2";
        if (dU==0 && at->bonded_to (-2, 8)==3) return "NO3";

        
        if (at->bonded_to ("C=ON") || at->bonded_to ("CONN") && at->bound.size ()==3) return "NC=O";
        for (unsigned int i=0; i<at->bound.size(); i++) {
            if (at->bound[i]->GetAtomicNum ()==6 && (at->bound[i]->bonded_to (2, 6) || at->bound[i]->bonded_to (5, 6))) return "NC=C";
        }

        for (unsigned int i=0; i<at->bound.size(); i++) {

            if (at->bound[i]->GetAtomicNum ()==6 && !at->bound[i]->bonded_to (2, 8) && at->bound[i]->bonded_to (-2,7)>1 && ((at->bound[i]->bonded_to (2, 7) || at->bound[i]->bonded_to (5, 7)))) return "NC=N";
            if (at->bound[i]->GetAtomicNum ()==6 && at->bound[i]->bonded_to (3, 7) && at->bound[i]->bonded_to (-2, 7)>1 && !at->bonded_to (3, 6)) return "NC%N";
            if (at->bound[i]->GetAtomicNum ()==7 && at->in_ring.size ()==0 && (at->bound[i]->bonded_to (2, 6) || at->bound[i]->bonded_to (5, 6)) && !at->bound[i]->bonded_to ("CGD+") && !at->bound[i]->bonded_to ("CGD")) return "NN=C";
            if (at->bound[i]->GetAtomicNum ()==7 && (at->bound[i]->bonded_to (2, 7) && at->bound[i]->bonded_to (1,7))) return "NN=N";

            if (at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (-2, 8)==3 ) return "NSO3";
            if (at->bound[i]->GetAtomicNum ()==15 && at->bound[i]->bonded_to (-2, 8)==2 ) return "NPO2";
            if (at->bound[i]->GetAtomicNum ()==15 && at->bound[i]->bonded_to (-2, 8)==3 ) return "NPO3";
        }




        for (unsigned int r=0; r<at->in_ring.size(); r++) {
            Ring *ring = at->in_ring[r];
            if (ring->aromatic && ring->atoms.size()==5 && at->in_ring.size()==1 && at->bound.size ()==3) return "NPYL";
            if (ring->aromatic && ring->atoms.size()==6 && at->in_ring.size()==1) {
                if (at->bonded_to (-2, 8)) return "NPOX";
                if (dU==1) return "NPYD";
                if (dU==0) return "NPD+";
            } 
        }
  
        for (unsigned int r=0; r<at->in_ring.size(); r++) {
            Ring *ring = at->in_ring[r];      
            if (ring->aromatic && ring->atoms.size ()==5 && ring->count (6)<5) {
                for (unsigned int a=0; a<at->bound.size(); a++) {
                    if (at->bound[a]->GetAtomicNum ()!=6 && at->bound[a]->in_ring.size()) {
                        if  (at->bound[a]->GetAtomicNum ()==7 && !at->bound[a]->bound.size()==3) {}
                        else return "N5A";
                    }
                }
            }    
        }
        for (unsigned int r=0; r<at->in_ring.size(); r++) {
            Ring *ring = at->in_ring[r];
            if (ring->aromatic && ring->atoms.size()==5) return "N5B";
        }
        if (dU == 0 && at->bonded_to (2, 7)) return "N+=N";

        if (at->bonded_to (3, -2)) return "NSP";
        if (at->bonded_to ("CSP")) return "NC%C";
        if (dU == 0) return "NR";
    }
//____________________________________________________________________________________________________________________
    if (at->GetAtomicNum () == 8) {// 6 7 32 35 49 51 59 70   O
 //       bool in_ring = at->in_ring.size()>0;
        unsigned int ring_size = 1000;
        bool aromatic=false;
        for (unsigned int r=0; r<at->in_ring.size (); r++){
            if (at->in_ring[r]->atoms.size ()< ring_size) ring_size =at->in_ring[r]->atoms.size ();
            if (at->in_ring[r]->aromatic) aromatic = true;
        }
 // ||      6  ALCOHOL || ETHER OXYGEN
//  OC=O    6  ESTER || CARBOXYLIC ACID -O-
//  OC=C    6  ENOLIC || PHENOLIC OXYGEN
//  OC=N    6  DIVALENT OXYGEN 
//  OC=S    6  THIOESTER || THIOACID -O-
//  ONO2    6  DIVALENT NITRATE "ETHER" OXYGEN
//  ON=O    6  DIVALENT NITRITE "ETHER" OXYGEN
//  OSO3    6  DIVALENT OXYGEN ATTACHED TO SULFUR
//  OSO2    6  DIVALENT OXYGEN ATTACHED TO SULFUR
//  OSO     6  DIVALENT OXYGEN ATTACHED TO SULFUR
//  OS=O    6  DIVALENT OXYGEN ATTACHED TO SULFOXIDE SULFUR
//  -OS     6  GENERAL DIVALENT OX ATTACHED TO S
//  OPO3    6  DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
//  OPO2    6  DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
//  OPO     6  DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
//  -OP     6  DIVALENT OXYGEN ATTACHED TO PHOSPHOROUS
//  -O-     6  GENERAL DIVALENT O
//  O=C     7  GENERAL C=O
//  O=CN    7  CARBONYL OXYGEN, AMIDES
//  O=CR    7  CARBONYL OXYGEN, ALDEHYDES && KETONES
//  O=CO    7  CARBONYL OXYGEN, CARBOXYLIC ACIDS && ESTERS
//  O=N     7  NITROSO OXYGEN
//  O=S     7  O=S IN SULFOXIDES
//  O=S=    7  O=S ON SULFUR DOUBLY BONDED TO, E.G., CARBON
//  O2CM   32  OXYGEN IN CARBOXYLATE ANION
//  OXN    32  N-OXIDE OXYGEN
//  O2N    32  NITRO OXYGEN
//  O2NO   32  NITRO-GROUP OXYGEN IN NITRATE
//  O3N    32  NITRATE ANION OXYGEN
//  O-S    32  SINGLE TERMINAL OXYGEN ON TETRACOORD SULFUR
//  O2S    32  TERMINAL O-S IN SULFONES && SULFONAMIDES
//  O3S    32  TERMINAL O IN SULFONATES
//  O4S    32  TERMINAL O IN SO4(-3)
  OSMS   32  TERM O IN THIOSULFINATE ANION - FORMAL CHARGE=-0.5
//  OP     32  TERMINAL O IN PHOSPHOXIDES
//  O2P    32  TERMINAL O IN PHOSPHINATES
//  O3P    32  TERMINAL OXYGEN IN PHOSPHONATES
//  O4P    32  TERMINAL OXYGEN IN PHOSPHATES && PHOSPHODIESTERS
  O4CL   32  OXYGEN IN CLO4(-) ANION - FORMAL CHARGE=-0.25
//  OM     35  ALKOXIDE OXYGEN, NEGATIVELY CHARGED
 // OM2    35  OXIDE OXYGEN ON SP2 CARBON, NEGATIVELY CHARGED
//  O+     49  POSITIVELY CHARGED OXONIUM (TRICOORDINATE) OXYGEN
//  O=+    51  POSITIVELY CHARGED OXENIUM (DICOORDINATE) OXYGEN
//  OFUR   59  AROMATIC OXYGEN AS IN FURAN
//  OH2    70  OXYGEN ON WATER  


        int dU = 2-at->bonds.size ();
        if (at->bonded_to ("COO") && at->bonded_to (2, 6)) return "O=CO";
        if (at->bonded_to ("COO")) return "OC=O";
        for (unsigned int i=0; i <at->bound.size (); i++) {
            if (at->bound[i]->GetAtomicNum ()==6 && at->bonds[i]->kekule==1 && at->bound[i]->bonded_to (2, 8) && at->bound[i]->bonded_to (1, 8)) return "OC=O";
        }
        if (at->bonded_to ("C=OR") || at->bonded_to ("CONN")) return "O=CR";
        if (at->bonded_to ("C=ON") && dU ==1) return "O=CN";
        if (at->bonded_to ("C=S")) return "OC=S";
        if (at->bonded_to ("CO2M")) return "O2CM";
        if (aromatic && ring_size==5) return "OFUR";
        if (dU ==-1) return "O+";
        if (dU ==1 && at->bonded_to (1, 6) ) {
            if (at->bound[0]->bound.size () == 3) return "OM2";
            return "OM";
        }
        if (dU ==1 && at->bonded_to (1, 7)) {
            if (at->bound[0]->bound.size () == 2) return "OM2";
      //      return "OM";
        }
        if (dU ==0 && at->bonded_to (1, 1)==2) return "OH2"; 
        if (dU==0 && (aromatic || at->bonded_to (2, -2))) return "O=+";
        if (dU==1 && at->bonded_to (-2, 7)) {
            if (at->bound[0]->bonded_to (-2, 8) ==3) {
                int tO =0; //terminal Os on N
                for (unsigned int i=0; i<at->bound[0]->bound.size(); i++) {
                    if (at->bound[0]->bound[i]->bound.size () ==1) tO+=1;
                }
            if (tO ==3) return "O3N";
            else return "O2NO";
            }
        }


        if (dU==1 && at->bonded_to (-2, 16)) {
            if (at->bound[0]->bound.size ()==4 && at->bound[0]->bonded_to (-2, 8) ==1) return "O-S";
            if (at->bound[0]->bonded_to (-2, 8) ==2) {
                bool boo= true;
                for (unsigned int i=0; i<at->bound[0]->bound.size (); i++) {
                    if (at->bound[0]->bound[i]->GetAtomicNum ()==8 && at->bound[0]->bound[i]->bound.size ()==2) boo=false;
                }
                if (boo) return "O2S";
            }
            if (at->bound[0]->bonded_to (-2, 8) ==3) return "O3S";
            if (at->bound[0]->bonded_to (-2, 8) ==4) return "O4S";
        }
        for (unsigned int i=0; i<at->bound.size(); i++) {
            if (at->bound[i]->GetAtomicNum ()==6 && at->bound[i]->bonded_to (2, 7) || at->bound[i]->bonded_to (5, 7)) return "OC=N";
        }
        for (unsigned int i=0; i<at->bound.size(); i++) {
            if (at->bound[i]->GetAtomicNum ()==6 && (at->bound[i]->bonded_to (2, 6) || at->bound[i]->bonded_to (5, 6))) return "OC=C";
        }
        for (unsigned int i=0; i<at->bound.size(); i++) {
            if (dU==0 && at->bound[i]->GetAtomicNum ()==7 && at->bound[i]->bonded_to (-2, 8)==3) return "ONO2";
            if (dU==0 && at->bound[i]->GetAtomicNum ()==7 && at->bound[i]->bonded_to (-2, 8)==2) return "ON=O";
            if (dU==1 && at->bound[i]->GetAtomicNum ()==7 && at->bound[i]->bonded_to (-2, 8)==2) return "O2N";
            if (dU==1 && at->bound[i]->GetAtomicNum ()==7 && at->bound[i]->bonded_to (2, 8)==1 && at->bound[i]->bonds.size()==2) return "O=N";
        if (dU==1 && at->bonded_to (1, 7)) return "OXN";
        }
        for (unsigned int i=0; i<at->bound.size(); i++) {
            if (dU==0 && at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (-2, 8)==4) return "OSO3";
            if (dU==0 && at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (-2, 8)==3) return "OSO2";
            if (dU==0 && at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (-2, 8)==2 && at->bound[i]->bonded_to (2, 8)) return "OS=O";
            if (dU==0 && at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (-2, 8)==2) return "OSO";
            if (dU==1 && at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (2, 8) && at->bound[i]->bonded_to (2, -2)==2) return "O=S=";
            if (dU==1 && at->bound[i]->GetAtomicNum ()==16 && at->bound[i]->bonded_to (2, 8)) return "O=S";

        }
        if (dU==0 && at->bonded_to (-2, 16)) return "-OS";
        for (unsigned int i=0; i<at->bound.size(); i++) {
            if (dU==0 && at->bound[i]->GetAtomicNum ()==15 && at->bound[i]->bonded_to (-2, 8)==4) return "OPO3";
            if (dU==0 && at->bound[i]->GetAtomicNum ()==15 && at->bound[i]->bonded_to (-2, 8)==3) return "OPO2";
            if (dU==0 && at->bound[i]->GetAtomicNum ()==15 && at->bound[i]->bonded_to (-2, 8)==2) return "OPO";
            if (dU==1 && at->bound[i]->GetAtomicNum ()==15 && at->bound[i]->bonded_to (-2, 8)==1) return "OP";
            if (dU==1 && at->bound[i]->GetAtomicNum ()==15 && at->bound[i]->bonded_to (-2, 8)==2) return "O2P";
            if (dU==1 && at->bound[i]->GetAtomicNum ()==15 && at->bound[i]->bonded_to (-2, 8)==3) return "O3P";
            if (dU==1 && at->bound[i]->GetAtomicNum ()==15 && at->bound[i]->bonded_to (-2, 8)==4) return "O4P";

        }
        if (dU==0 && at->bonded_to (-2, 15)) return "-OP";


        if (dU==0) {
            if (at->bonded_to (1, 6) ==2 || (at->bonded_to (1, 6)==1 && at->bonded_to (1,1) ==1)) return "OR";
            return "-O-";
        }
        if (at->bonded_to (2, 6)) return "O=C";

    }
//_________________________________________________________________________________________________________________
    if (at->GetAtomicNum () == 9) { // 11 89   F
//  F      11  FLUORINE
//  F-     89  FLUORIDE ANION
            if (at->bonds.size () ==0) return "F-";
            if (at->bonds.size () ==1) return "F";
    }
    
//_________________________________________________________________________________________________    
    if (at->GetAtomicNum () == 11) {//  93   Na
//  NA+    93  SODIUM CATION
            return "NA+";
    }
//_________________________________________________________________________________________________    
    if (at->GetAtomicNum () == 12) { // 99   Mg
//  MG+2   99  DIPOSITIVE MAGNESIUM CATION
        return "MG+2";
    }
//_________________________________________________________________________________________________    
    if (at->GetAtomicNum () == 14) { // 19   Si
//  SI     19  SILICON
        return "SI";
     }   

//____________________________________________________________________________________________________________
    if (at->GetAtomicNum () == 15) {// 25 26 75   P
//  PO4    25  PHOSPHOROUS IN PHOSPHATES && PHOSPHODIESTERS
//  PO3    25  TETRACOORDINATE P WITH THREE ATTACHED OXYGENS
//  PO2    25  TETRACOORDINATE P WITH TWO ATTACHED OXYGENS
//  PO     25  TETRACOORDINATE P WITH ONE ATTACHED OXYGEN
//  PTET   25  GENERAL TETRACOORDINATE PHOSPHORUS
//  P      26  TRICOORDINATE P, AS IN PHOSPHINES
//  -P=C   75  PHOSPHOROUS DOUBLY BONDED TO CARBON
        int coord = at->bonds.size ();
        if (at->bonded_to (2, 6)) return "-P=C";
        if (coord==4 && at->bonded_to (-2, 8)==4) return "PO4";
        if (coord==4 && at->bonded_to (-2, 8)==3) return "PO3";
        if (coord==4 && at->bonded_to (-2, 8)==2) return "PO2";
        if (coord==4 && at->bonded_to (-2, 8)==1) return "PO";
        if (coord==4) return "PTET";
        if (coord==3) return "P";
    }

 //_________________________________________________________________________________________________________________   
    if (at->GetAtomicNum () ==16) { // 15 16 17 18 44 72 73 74   S
        bool aromatic=false;
  //      bool in_ring = at->in_ring.size()>0;
        int ring_size = 1000;
        for (unsigned int r=0; r<at->in_ring.size (); r++){
            if (at->in_ring[r]->atoms.size ()< ring_size) ring_size =at->in_ring[r]->atoms.size ();
            if (at->in_ring[r]->aromatic) aromatic = true;
        }
/*  
//  S      15  SULFUR IN THIOETHERS && MERCAPTANS
//  S=C    16  TERMINAL SULFUR DOUBLY BONDED TO CARBON
 // S=O    17  SULFUR IN SULFOXIDES
//  >S=N   17  SULFUR, TRICOORD, DOUBLY BONDED TO N
//  SO2    18  SULFUR IN SULFONES 
//  SO2N   18  SULFUR IN SULFONAMIDES
//  SO3    18  SULFONATE SULFUR
//  SO4    18  SULFATE SULFUR
//  =SO2   18  SULFONE SULPHER DOUBLY BONDED TO CARBON
//  SNO    18  SULFUR IN NITROGEN ANALOG OF A SULFONE
//  STHI   44  SULFUR AS IN THIOPHENE
//  S-P    72  TERMINAL SULFUR BONDED TO PHOSPHORUS
//  S2CM   72  TERMINAL SULFUR IN THIOCARBOXYLATE ANION
//  SM     72  TERMINAL SULFUR - FORMAL CHARGE=-1
  SSMO   72  TERMINAL SULFUR IN THIOSULFINATE GROUP
  SO2M   73  SULFUR IN NEGATIVELY CHARGED SULFINATE GROUP
  SSOM   73  TRICOORD SULFUR IN THIOSULFINATE GROUP
//  =S=O   74  SULFINYL SULFUR, EG. IN C=S=O

        int coord = at->bonds.size ();
    
        if (aromatic && ring_size==5) return "STHI";
        if (at->bonded_to ("CS2M")) return "S2CM";
        if (coord==2 && at->bonded_to (2, -2)==2 && at->bonded_to (2, 8)) return "=S=O";
        if (coord==3 && at->bonded_to (2, 7)) return ">S=N";
        if (coord==3 && at->bonded_to (2, 8)) return "S=O";
        if (coord==4 && at->bonded_to (-2, 8)==4) return "SO4";
        if (coord==4 && at->bonded_to (-2, 8)==3) return "SO3";
        if (coord==4 && at->bonded_to (-2, 8)==2 && at->bonded_to (-2, 7)) return "SO2N";
        if (coord==3 && at->bonded_to (-2, 8)==2 && at->bonded_to (2, 6)) return "=SO2";

        if (coord==4 && at->bonded_to (-2, 8)==2) return "SO2";
        if (coord==4 && at->bonded_to (-2, 8) && at->bonded_to (-2, 7)) return "SNO";
    
        if (coord==1 && at->bonded_to (2, 6)) return "S=C";



        if (coord==1 && at->bonded_to (-2, 15)) return "S-P"; 
        if (coord==2) return "S";
        if (coord==1) return "SM";
    }
//________________________________________________________________________________________________________________
    if (at->GetAtomicNum () ==17) { // 12 77 90 Cl
//  CL     12  CHLORINE
//  CLO4   77  CHLORINE IN PERCHLORATE ANION, CLO4(-)
//  CL-    90  CHLORIDE ANION
        if (at->bonds.size() == 0) return "CL-";
        if (at->bonded_to (-2, 8) == 4) return "CLO4";
        return "CL";
    }
//_________________________________________________________________________________________________    
    if (at->GetAtomicNum () ==19) { // 94 K
//  K+     94  POTASSIUM CATION
            return "K+";
    }
//_________________________________________________________________________________________________    
    if (at->GetAtomicNum () ==20) { // 96 Ca
//  CA+2   96  DIPOSITIVE CALCIUM 
            return "CA+2";
    }
//_________________________________________________________________________________________________    
    if (at->GetAtomicNum () ==26) { // 87 88   Fe
//  FE+2   87  IRON +2 CATION
//////////////  FE+3   88  IRON +3 CATION
        return "FE+2";  
    }
//_________________________________________________________________________________________________    
    if (at->GetAtomicNum () ==29) { // 97 98 Cu
//  CU+1   97  MONOPOSITIVE COPPER    
/////////////  CU+2   98  DIPOSITIVE COPPER     
    return "CU+1";
    }
//_________________________________________________________________________________________________        
    if (at->GetAtomicNum () ==30) { // 95 Zn
//  ZINC   95  DIPOSITIVE ZINC        
//  ZN+2   95  DIPOSITIVE ZINC   
        return "ZN+2";
    }
//_________________________________________________________________________________________________    
    if (at->GetAtomicNum () ==35) { // 13 91 Br
//  BR     13  BROMINE
//  BR-    91  BROMIDE ANION
            if (!at->bonds.size ()) return "BR-";
            return "BR";
        }
//_________________________________________________________________________________________________    
    if (at->GetAtomicNum () ==53) { // 14 I
//  I      14  IODINE
        return "I";
        }


    return "Unknown eteroatom";
*/
}

 
  int MMFF::getBondType(Atom* a, Atom* b)
  {
	ZNMolecule *mol = (ZNMolecule *) a ->GetParent ();
	ZNBond *bond = mol ->GetBond (a, b);
	int I = get_MMFFtype (a);
	int J = get_MMFFtype (b);
    if (!bond ->IsSingle())
      return 0;
    
    if (!bond ->IsAromatic())
      if (get_arom (I) && get_arom (J))
        return 1;
      
    if (get_sbmb (I) && get_sbmb (J))
      return 1;
    
    return 0;
  }

  int MMFF::getAngleType(Atom* a, Atom* b, Atom *c)
  {
    int sumbondtypes;

    sumbondtypes = getBondType(a,b) + getBondType(b, c);

    if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && IsInSameRing(a, c))
      switch (sumbondtypes) {
      case 0:
        return 3; 
      case 1:
        return 5; 
      case 2:
        return 6; 
      }
    
    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, c))
      switch (sumbondtypes) {
      case 0:
        return 4; 
      case 1:
        return 7; 
      case 2:
        return 8; 
      }
    
    return sumbondtypes;
  }



int MMFF::getMMFFtype (Atom *at) {

    return std::numeric_limits<int>::infinity();

/*
    string st = at->MMFFstring;
    for (unsigned int i=0; i<atypeParameters.size (); i++) {
    //    cout <<atypeParameters[i]->strin<<endl;
        if (st==atypeParameters[i]->strin) return atypeParameters[i]->number; 
    }
    cout << "string "<<st<<" does not name an atom type"<<endl;
    return 1;
*/
}
