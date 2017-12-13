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
#ifndef _OBJECT_IN_FLUID_OIF_STRETCHING_H
#define _OBJECT_IN_FLUID_OIF_STRETCHING_H

/** \file oif_stretching.hpp
 *  Routines to calculate the OIF_STRETCHING
 *  for a particle pair (two neighboring particles with common edge). (Dupin2007)
 *  \ref forces.cpp
 */

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"
#include "config.hpp"
#include "integrate.hpp"

// set parameters for stretching
int oif_stretching_set_params(int bond_type, double r0, double ks, double kslin);

inline double ks_nonlinearity(double lambda){ // Defined by (19) from Dupin2007
    double res;
    res = (pow(lambda,0.5) + pow(lambda,-2.5))/(lambda + pow(lambda,-3.));
    return res;
}

/** Computes the stretching forces (Dupin2007)
    @return 0
*/
inline int calc_oif_stretching_pair_force(Particle *p1, Particle *p2,
				 Bonded_ia_parameters *iaparams, double force[3])// first-fold-then-the-same approach
{
	int i, img[3];
	double fp1[3],fp2[3];
	double AA[3];
    double dx[3],fac,dr,len2,len,lambda;

	// first find out which particle out of p1, p2 (possibly p3, p4) is not a ghost particle. In almost all cases it is p2, however, it might be other one. we call this particle reference particle.
	if (p2->l.ghost != 1) {
		//unfold non-ghost particle using image, because for physical particles, the structure p->l.i is correctly set
		memmove(fp2, p2->r.p, 3*sizeof(double));
		memmove(img, p2->l.i, 3*sizeof(int));
		unfold_position(fp2,img);
		// other coordinates are obtained from its relative positions to the reference particle
		get_mi_vector(AA, p1->r.p, fp2);
		for (i=0; i<3; i++) {
            fp1[i] = fp2[i] + AA[i];
        }
	} else {
		// in case  particle p2 is a ghost particle
		if (p1->l.ghost != 1) {
			memmove(fp1, p1->r.p, 3*sizeof(double));
			memmove(img, p1->l.i, 3*sizeof(int));
			unfold_position(fp1,img);
			get_mi_vector(AA, p2->r.p, fp1);
			for (i=0; i<3; i++) {
                fp2[i] = fp1[i] + AA[i];
            }
        } else {
			printf("Something wrong in oif_stretching.hpp: Both particles in a bond are ghost particles, impossible to unfold the positions...\n");
			return 0;
		}
	}

    // non-linear stretching
    if (iaparams->p.oif_stretching.ks > TINY_OIF_ELASTICITY_COEFFICIENT) {
        vecsub(fp2,fp1,dx);
        len2 = sqrlen(dx);
        len = sqrt(len2);
        dr = len - iaparams->p.oif_stretching.r0;
        lambda = 1.0*len/iaparams->p.oif_stretching.r0;
        fac = iaparams->p.oif_stretching.ks * ks_nonlinearity(lambda) * dr; // no normalization
        for(i=0; i<3; i++) {
            force[i] = fac*dx[i]/len;
        }
    }
    
    // linear stretching
    if (iaparams->p.oif_stretching.kslin > TINY_OIF_ELASTICITY_COEFFICIENT) {
        vecsub(fp2,fp1,dx);
        len2 = sqrlen(dx);
        len = sqrt(len2);
        dr = len - iaparams->p.oif_stretching.r0;
        fac = iaparams->p.oif_stretching.kslin * dr; // no normalization
        for(i=0; i<3; i++) {
            force[i] = fac*dx[i]/len;
        }
    }

  return 0;
}

#endif

