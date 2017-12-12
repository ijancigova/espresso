/*
  Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _OBJECT_IN_FLUID_OIF_LOCAL_AREA_H
#define _OBJECT_IN_FLUID_OIF_LOCAL_AREA_H

/** \file oif_local_area.hpp
 *  Routines to calculate the OIF_LOCAL_AREA
 *  for three particles that form a triangle. (Dupin2007)
 *  \ref forces.cpp
 */

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"
#include "config.hpp"
#include "integrate.hpp"

// set parameters for local area
int oif_local_area_set_params(int bond_type, double A0, double kal);


/** Computes the local area forces (Dupin2007) and adds them
    to the particle forces (see \ref tclcommand_inter). 
    @param p1,p2,p3     Pointers to particles of triangle.
    @return 0
*/
inline int calc_oif_local_area(Particle *p1, Particle *p2, Particle *p3,
				 Bonded_ia_parameters *iaparams, double force[3],
				 double force2[3], double force3[3])// first-fold-then-the-same approach
{
	int i, img[3];
	double AA[3],BB[3];
	double fp1[3],fp2[3],fp3[3];
    double A,h[3],rh[3],hn;
    double fac,m1[3],m2[3],m3[3];
    double m1_length,m2_length,m3_length,t;

	// first find out which particle out of p1, p2 (possibly p3, p4) is not a ghost particle. In almost all cases it is p2, however, it might be other one. we call this particle reference particle.
	if (p2->l.ghost != 1) {
		//unfold non-ghost particle using image, because for physical particles, the structure p->l.i is correctly set
		memmove(fp2, p2->r.p, 3*sizeof(double));
		memmove(img, p2->l.i, 3*sizeof(int));
		unfold_position(fp2,img);
		// other coordinates are obtained from its relative positions to the reference particle
		get_mi_vector(AA, p1->r.p, fp2);
		get_mi_vector(BB, p3->r.p, fp2);
		for (i=0; i<3; i++) {
            fp1[i] = fp2[i] + AA[i];
            fp3[i] = fp2[i] + BB[i];
        }
	} else {
		// in case  particle p2 is a ghost particle
		if (p1->l.ghost != 1) {
			memmove(fp1, p1->r.p, 3*sizeof(double));
			memmove(img, p1->l.i, 3*sizeof(int));
			unfold_position(fp1,img);
			get_mi_vector(AA, p2->r.p, fp1);
			get_mi_vector(BB, p3->r.p, fp1);
			for (i=0; i<3; i++) {
                fp2[i] = fp1[i] + AA[i];
                fp3[i] = fp1[i] + BB[i];
            }
		} else {
			// in case the first and the second particle are ghost particles
			if (p3->l.ghost != 1) {
				memmove(fp3, p3->r.p, 3*sizeof(double));
				memmove(img, p3->l.i, 3*sizeof(int));
				unfold_position(fp3,img);
				get_mi_vector(AA, p1->r.p, fp3);
				get_mi_vector(BB, p2->r.p, fp3);
				for (i=0; i<3; i++) {
                    fp1[i] = fp3[i] + AA[i];
                    fp2[i] = fp3[i] + BB[i];
                }
			} else {
				printf("Something wrong in oif_local_forces.hpp: All particles in a bond are ghost particles, impossible to unfold the positions...\n");
				return 0;
			}
		}
	}
	
    for(i=0; i<3; i++) {
        force[i] = 0;
        force2[i] = 0;
        force3[i] = 0;
    }

    /* local area
       only 1/3 of calculated forces are added, because each triangle will enter this calculation 3 times (one time per edge)

		Proportional distribution of forces, implemented according to the article
		I.Jancigova, I.Cimrak, 
		Non-uniform force allocation for area preservation in spring network models,  
		International Journal for Numerical Methods in Biomedical Engineering, DOI: 10.1002/cnm.2757
     
    */
    if (iaparams->p.oif_local_area.kal > TINY_OIF_ELASTICITY_COEFFICIENT) {
        
        // triangle p1,p2,p3
        for(i=0; i<3; i++){            // centroid of triangle p1,p2,p3
            h[i]=1.0/3.0 *(fp1[i]+fp2[i]+fp3[i]);
        }
        A=area_triangle(fp1,fp2,fp3); 	
        t = sqrt(A/iaparams->p.oif_local_area.A0) - 1.0;
        vecsub(h,fp1,m1);
        vecsub(h,fp2,m2);
        vecsub(h,fp3,m3);					
        
        m1_length = normr(m1);
		m2_length = normr(m2);
		m3_length = normr(m3);
		
        fac = iaparams->p.oif_local_area.kal*iaparams->p.oif_local_area.A0*(2*t+t*t)/(m1_length*m1_length + m2_length*m2_length + m3_length*m3_length);
		
		for(i=0; i<3; i++) {          // local area force for p1
			force[i] += fac*m1[i]/3.0;
		}    
		for(i=0; i<3; i++) {          // local area force for p2
			force2[i] += fac*m2[i]/3.0;
		}
		for(i=0; i<3; i++) {          // local area force for p3
			force3[i] += fac*m3[i]/3.0;
		}
    }
  return 0;
}

#endif

