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
#ifndef _OBJECT_IN_FLUID_OIF_LOCAL_FORCES_H
#define _OBJECT_IN_FLUID_OIF_LOCAL_FORCES_H

/** \file oif_local_forces.hpp
 *  Routines to calculate the OIF_LOCAL_FORCES
 *  for a particle quadruple (two neighboring triangles with common edge). (Dupin2007)
 *  \ref forces.cpp
 */

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"
#include "config.hpp"
#include "integrate.hpp"

// set parameters for local forces
int oif_local_forces_set_params(int bond_type, double r0, double ks, double kslin, double phi0, double kb, double A01, double A02, double kal, double kvisc);

inline double KS(double lambda){ // Defined by (19) from Dupin2007
    double res;
    res = (pow(lambda,0.5) + pow(lambda,-2.5))/(lambda + pow(lambda,-3.));
    return res;
}

/** Computes the local forces (Dupin2007) and adds them
    to the particle forces (see \ref tclcommand_inter). 
    @param p1,p2,p3     Pointers to particles of triangle 1.
    @param p2,p3,p4     Pointers to particles of triangle 2.
    (triangles have particles p2 and p3 in common)
    @return 0
*/
inline int calc_oif_local(Particle *p2, Particle *p1, Particle *p3, Particle *p4,
				 Bonded_ia_parameters *iaparams, double force[3],
				 double force2[3], double force3[3], double force4[3])// first-fold-then-the-same approach
{
	int i, img[3];
	double fp1[3],fp2[3],fp3[3],fp4[3];
	double AA[3],BB[3],CC[3];
    double n1[3],n2[3],dn1,dn2,phi,aa;
    double dx[3],fac,dr,len2,len,lambda;
    double A,A1,A2,h[3],rh[3],hn,h1[3],h2[3],h3[3],h4[3];
    double v1,v2,footC[3],footD[3];
    double m1[3],m2[3],m3[3];
    double v[3], len_new, len_old, def_vel;
    double m1_length,m2_length,m3_length,t;

    // variables for sanity check
    double Aforce[3],Bforce[3],Cforce[3],Dforce[3];
    double Atorque[3],Btorque[3],Ctorque[3],Dtorque[3];
    double centroid[3],AtoCentroid[3],BtoCentroid[3],CtoCentroid[3],DtoCentroid[3];
    double check_force[3],check_torque[3];

	#ifdef GHOST_FLAG
	// first find out which particle out of p1, p2 (possibly p3, p4) is not a ghost particle. In almost all cases it is p2, however, it might be other one. we call this particle reference particle.
	if (p2->l.ghost != 1) {
		//unfold non-ghost particle using image, because for physical particles, the structure p->l.i is correctly set
		memmove(fp2, p2->r.p, 3*sizeof(double));
		memmove(img, p2->l.i, 3*sizeof(int));
		unfold_position(fp2,img);
		// other coordinates are obtained from its relative positions to the reference particle
		get_mi_vector(AA, p1->r.p, fp2);
		get_mi_vector(BB, p3->r.p, fp2);
		get_mi_vector(CC, p4->r.p, fp2);
		for (i=0; i<3; i++) {
            fp1[i] = fp2[i] + AA[i];
            fp3[i] = fp2[i] + BB[i];
            fp4[i] = fp2[i] + CC[i];
        }
	} else {
		// in case  particle p2 is a ghost particle
		if (p1->l.ghost != 1) {
			memmove(fp1, p1->r.p, 3*sizeof(double));
			memmove(img, p1->l.i, 3*sizeof(int));
			unfold_position(fp1,img);
			get_mi_vector(AA, p2->r.p, fp1);
			get_mi_vector(BB, p3->r.p, fp1);
			get_mi_vector(CC, p4->r.p, fp1);
			for (i=0; i<3; i++) {
                fp2[i] = fp1[i] + AA[i];
                fp3[i] = fp1[i] + BB[i];
                fp4[i] = fp1[i] + CC[i];
            }
		} else {
			// in case the first and the second particle are ghost particles
			if (p3->l.ghost != 1) {
				memmove(fp3, p3->r.p, 3*sizeof(double));
				memmove(img, p3->l.i, 3*sizeof(int));
				unfold_position(fp3,img);
				get_mi_vector(AA, p1->r.p, fp3);
				get_mi_vector(BB, p2->r.p, fp3);
				get_mi_vector(CC, p4->r.p, fp3);
				for (i=0; i<3; i++) {
                    fp1[i] = fp3[i] + AA[i];
                    fp2[i] = fp3[i] + BB[i];
                    fp4[i] = fp3[i] + CC[i];
                }
			} else {
				// in case the first and the second particle are ghost particles
				if (p4->l.ghost != 1) {
					memmove(fp4, p4->r.p, 3*sizeof(double));
					memmove(img, p4->l.i, 3*sizeof(int));
					unfold_position(fp4,img);
					get_mi_vector(AA, p1->r.p, fp4);
					get_mi_vector(BB, p2->r.p, fp4);
					get_mi_vector(CC, p3->r.p, fp4);
					for (i=0; i<3; i++) {
                        fp1[i] = fp4[i] + AA[i];
                        fp2[i] = fp4[i] + BB[i];
                        fp3[i] = fp4[i] + CC[i];
                    }
				} else {
					printf("Something wrong in oif_local_forces.hpp: All particles in a bond are ghost particles, impossible to unfold the positions...\n");
					return 0;
				}
			}
		}
	}
	#endif
	#ifndef GHOST_FLAG
		// if ghost flag was not defined we have no other option than to assume the first particle (p2) is a physical one.
		memmove(fp2, p2->r.p, 3*sizeof(double));
		memmove(img, p2->l.i, 3*sizeof(int));
		unfold_position(fp2,img);
		// other coordinates are obtained from its relative positions to the reference particle
		get_mi_vector(AA, p1->r.p, fp2);
		get_mi_vector(BB, p3->r.p, fp2);
		get_mi_vector(CC, p4->r.p, fp2);
		for (i=0; i<3; i++) {
            fp1[i] = fp2[i] + AA[i];
            fp3[i] = fp2[i] + BB[i];
            fp4[i] = fp2[i] + CC[i];
        }
	#endif
	
    for(i=0; i<3; i++) {
        force[i] = 0;
        force2[i] = 0;
        force3[i] = 0;
        force4[i] = 0;
    }

    // non-linear stretching
    if (iaparams->p.oif_local_forces.ks > TINY_OIF_ELASTICITY_COEFFICIENT) {
        vecsub(fp2,fp3,dx);
        len2 = sqrlen(dx);
        len = sqrt(len2);
        dr = len - iaparams->p.oif_local_forces.r0;
        lambda = 1.0*len/iaparams->p.oif_local_forces.r0;
        fac = -iaparams->p.oif_local_forces.ks * KS(lambda) * dr; // no normalization
        for(i=0; i<3; i++) {
            force2[i] += fac*dx[i]/len;
            force3[i] += -fac*dx[i]/len;
        }
    }
    
    // linear stretching
    if (iaparams->p.oif_local_forces.kslin > TINY_OIF_ELASTICITY_COEFFICIENT) {
        vecsub(fp2,fp3,dx);
        len2 = sqrlen(dx);
        len = sqrt(len2);
        dr = len - iaparams->p.oif_local_forces.r0;
        fac = -iaparams->p.oif_local_forces.kslin * dr; // no normalization
        for(i=0; i<3; i++) {
            force2[i] += fac*dx[i]/len;
            force3[i] += -fac*dx[i]/len;
        }
    }
    
     // viscous force
    if (iaparams->p.oif_local_forces.kvisc > TINY_OIF_ELASTICITY_COEFFICIENT) {	// to be implemented....
        //vecsub(fp2,fp3,dx);
        //len2 = sqrlen(dx);
        //len = sqrt(len2);
        //dr = len - iaparams->p.oif_local_forces.r0;
        //fac = -iaparams->p.oif_local_forces.kslin * dr; // no normalization
        //for(i=0; i<3; i++) {
            //force2[i] += fac*dx[i]/len;
            //force3[i] += -fac*dx[i]/len;
        //}
        
        // KOD od Mata:
        vecsub(fp2,fp3,dx);
        len2 = sqrlen(dx);
        len = sqrt(len2);
          
		len_new = sqrt(pow((fp2[0] +  p2->m.v[0]) - (fp3[0] +  p3->m.v[0]),2) + pow((fp2[1] +  p2->m.v[1]) - (fp3[1] +  p3->m.v[1]),2) + pow((fp2[2] +  p2->m.v[2]) - (fp3[2] +  p3->m.v[2]),2));
		len_old = sqrt(pow((fp2[0] -  p2->m.v[0]) - (fp3[0] -  p3->m.v[0]),2) + pow((fp2[1] -  p2->m.v[1]) - (fp3[1] -  p3->m.v[1]),2) + pow((fp2[2] -  p2->m.v[2]) - (fp3[2] -  p3->m.v[2]),2));

		//printf("future length %lf\n", len_new);
		//printf("old length %lf\n", len_old);    
		//printf("actual length %lf\n", len);

		if ( - len_old + len < 0 ) {
			def_vel = -sqrt(pow(p2->m.v[0] - p3->m.v[0],2) + pow(p2->m.v[1] - p3->m.v[1],2) + pow(p2->m.v[2] - p3->m.v[2],2));}
		else if ( - len_old + len > 0 ) {
			def_vel = +sqrt(pow(p2->m.v[0] - p3->m.v[0],2) + pow(p2->m.v[1] - p3->m.v[1],2) + pow(p2->m.v[2] - p3->m.v[2],2));}
		else { def_vel = 0; }
 	 
		fac = -iaparams->p.oif_local_forces.kvisc*def_vel;
  
		/* unscale velocities ! */
		v[0] = (p3->m.v[0] - p2->m.v[0])/time_step;
		v[1] = (p3->m.v[1] - p2->m.v[1])/time_step;
		v[2] = (p3->m.v[2] - p2->m.v[2])/time_step;
 
		// Variant A
		for(i=0;i<3;i++) {
			force2[i] += iaparams->p.oif_local_forces.kvisc*v[i];
			force3[i] -= iaparams->p.oif_local_forces.kvisc*v[i];
		}
		// Variant B
		//for(i=0;i<3;i++) {
			//force2[i] += fac*dx[i]/time_step;
			//force3[i] -= fac*dx[i]/time_step;
		//}
    }
    
    /* bending
       force-free and torque free bending
       m = 1/(B-A), where B = fp3, A = fp2, C = fp1, D = fp4
     */
    if (iaparams->p.oif_local_forces.kb > TINY_OIF_ELASTICITY_COEFFICIENT) {
        if (iaparams->p.oif_local_forces.phi0 < 0.001 || iaparams->p.oif_local_forces.phi0 > 2*M_PI - 0.001)
            printf("oif_local_forces.hpp, calc_oif_local: Resting angle is close to zero!\n");
        get_n_triangle(fp2,fp1,fp3,n1);
        dn1=normr(n1);
        get_n_triangle(fp2,fp3,fp4,n2);
        dn2=normr(n2);
        phi = angle_btw_triangles(fp1,fp2,fp3,fp4);
        if (phi < 0.001 || phi > 2*M_PI - 0.001) printf("oif_local_forces.hpp, calc_oif_local: Angle approaches 0 or 2*Pi\n");

        // common edge
        vecsub(fp2,fp3,h);  // A - B

        // altitude (v1) of the first triangle onto the common edge
        vecsub(fp3,fp1,h1);    // B - C
        vecsub(fp2,fp1,h2);    // A - C
        vector_product(h1,h2,h3);
        v1 = normr(h3)/normr(h);

        // foot of this altitude on the common edge
        t = - scalar(h,h2)/(normr(h)*normr(h));
        for(i=0; i<3; i++) {
            footC[i] = fp2[i] + t * h[i];
        }

        // altitude (v2) of the second triangle onto the common edge
        vecsub(fp3,fp4,h1);   // B - D
        vecsub(fp2,fp4,h2);   // A - D
        vector_product(h1,h2,h3);
        v2 = normr(h3)/normr(h);

        // foot of this altitude on the common edge
        t = - scalar(h,h2)/(normr(h)*normr(h));
        for(i=0; i<3; i++) {
            footD[i] = fp2[i] + t * h[i];
        }

        // beta
        fac = iaparams->p.oif_local_forces.kb * phi;

        // areas of the two triangles
        A1=area_triangle(fp1,fp2,fp3);
        A2=area_triangle(fp2,fp3,fp4);

        // feet of the altitudes to vertices on the common edge
        vecsub(footC,fp3,h1);
        vecsub(footD,fp3,h2);
        vecsub(footC,fp2,h3);
        vecsub(footD,fp2,h4);

        for(i=0; i<3; i++) {
            force[i] += fac * n1[i]/(dn1*v1);
            force2[i] -= fac * (normr(h1)*n1[i]/(2*A1*dn1) + normr(h2)*n2[i]/(2*A2*dn2));
            force3[i] -= fac * (normr(h3)*n1[i]/(2*A1*dn1) + normr(h4)*n2[i]/(2*A2*dn2));
            force4[i] += fac * n2[i]/(dn2*v2);
        }

        // verify force-free
        for(i=0; i<3; i++) {
            Aforce[i] = fac * n1[i]/(dn1*v1);
            Bforce[i] = - fac * (normr(h1)*n1[i]/(2*A1*dn1) + normr(h2)*n2[i]/(2*A2*dn2));
            Cforce[i] = - fac * (normr(h3)*n1[i]/(2*A1*dn1) + normr(h4)*n2[i]/(2*A2*dn2));
            Dforce[i] = fac * n2[i]/(dn2*v2);
        }
        for(i=0; i<3; i++) {
            check_force[i] = Aforce[i] + Bforce[i] + Cforce[i] + Dforce[i];
        }
        printf("oif_local_forces.hpp, bending: Total force = [%e, %e, %e] \n",check_force[0], check_force[1], check_force[2]);

        // verify torque-free
        for(i=0; i<3; i++) {
            centroid[i] = (fp1[i] + fp2[i] + fp3[i] + fp4[i])/4.0;
        }
        vecsub(fp2,centroid,AtoCentroid);
        vecsub(fp3,centroid,BtoCentroid);
        vecsub(fp1,centroid,CtoCentroid);
        vecsub(fp4,centroid,DtoCentroid);
        vector_product(Aforce,AtoCentroid,Atorque);
        vector_product(Bforce,BtoCentroid,Btorque);
        vector_product(Cforce,CtoCentroid,Ctorque);
        vector_product(Dforce,DtoCentroid,Dtorque);
        for(i=0; i<3; i++) {
            check_torque[i] = Atorque[i] + Btorque[i] + Ctorque[i] + Dtorque[i];
        }
        printf("oif_local_forces.hpp, bending: Total torque = [%e, %e, %e] \n",check_torque[0], check_torque[1], check_torque[2]);

    }

    /* local area
       for both triangles
       only 1/3 of calculated forces are added, because each triangle will enter this calculation 3 times (one time per edge)

		Proportional distribution of forces, implemented according to the article
		I.Jancigova, I.Cimrak, 
		Non-uniform force allocation for area preservation in spring network models,  
		International Journal for Numerical Methods in Biomedical Engineering, DOI: 10.1002/cnm.2757
     
    */
    if (iaparams->p.oif_local_forces.kal > TINY_OIF_ELASTICITY_COEFFICIENT) {
        
        // triangle p1,p2,p3
        for(i=0; i<3; i++){            // centroid of triangle p1,p2,p3
            h[i]=1.0/3.0 *(fp1[i]+fp2[i]+fp3[i]);
        }
        A=area_triangle(fp1,fp2,fp3); 	
        t = sqrt(A/iaparams->p.oif_local_forces.A01) - 1.0; 
        vecsub(h,fp1,m1);
        vecsub(h,fp2,m2);
        vecsub(h,fp3,m3);					
        
        m1_length = normr(m1);
		m2_length = normr(m2);
		m3_length = normr(m3);
		
        fac = iaparams->p.oif_local_forces.kal*iaparams->p.oif_local_forces.A01*(2*t+t*t)/(m1_length*m1_length + m2_length*m2_length + m3_length*m3_length);
		
		for(i=0; i<3; i++) {          // local area force for p1
			force[i] += fac*m1[i]/3.0;
		}    
		for(i=0; i<3; i++) {          // local area force for p2
			force2[i] += fac*m2[i]/3.0;
		}
		for(i=0; i<3; i++) {          // local area force for p3
			force3[i] += fac*m3[i]/3.0;
		}
            
        // triangle p2,p3,p4
        for(i=0; i<3; i++) {         // centroid of triangle p2,p3,p4
            h[i]=1.0/3.0 *(fp2[i]+fp3[i]+fp4[i]);
        }
        A=area_triangle(fp2,fp3,fp4);
        t = sqrt(A/iaparams->p.oif_local_forces.A02) - 1.0;    				////
        vecsub(h,fp2,m1);
        vecsub(h,fp3,m2);
        vecsub(h,fp4,m3);					
        
        m1_length = normr(m1);
		m2_length = normr(m2);
		m3_length = normr(m3);
		
		fac = iaparams->p.oif_local_forces.kal*iaparams->p.oif_local_forces.A02*(2*t+t*t)/(m1_length*m1_length + m2_length*m2_length + m3_length*m3_length);
		
		for(i=0; i<3; i++) {          // local area force for p2
			force2[i] += fac*m1[i]/3.0;
		}    
		for(i=0; i<3; i++) {          // local area force for p3
			force3[i] += fac*m2[i]/3.0;
		}
		for(i=0; i<3; i++) {          // local area force for p4
			force4[i] += fac*m3[i]/3.0;
		}

    }
  return 0;
}

#endif

