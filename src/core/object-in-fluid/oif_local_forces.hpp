/*
 * Copyright (C) 2012-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _OBJECT_IN_FLUID_OIF_LOCAL_FORCES_H
#define _OBJECT_IN_FLUID_OIF_LOCAL_FORCES_H

/** \file
 *  Routines to calculate the OIF local forces for a particle quadruple
 *  (two neighboring triangles with common edge).
 *  See @cite dupin07a.
 *
 *  Implementation in \ref oif_local_forces.cpp
 */

#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "grid.hpp"
//<<<<<<< HEAD
//#include "particle_data.hpp"
//=======
#include <utils/Vector.hpp>
//>>>>>>> espressomd/python
#include <utils/math/triangle_functions.hpp>

/** Set parameters for OIF local forces */
int oif_local_forces_set_params(int bond_type, double r0, double ks,
                                double kslin, double phi0, double kb,
                                double A01, double A02, double kal,
                                double kvisc);

/** @details see eq. (19) in @cite dupin07a */
inline double KS(double lambda) {
  double res;
  res = (pow(lambda, 0.5) + pow(lambda, -2.5)) / (lambda + pow(lambda, -3.));
  return res;
}

/** Compute the OIF local forces.
 *  See @cite dupin07a, @cite jancigova16a.
 *  @param p2           %Particle of triangle 1.
 *  @param p1 , p3      Particles common to triangle 1 and triangle 2.
 *  @param p4           %Particle of triangle 2.
 *  @param iaparams     Bonded parameters for the OIF interaction.
 *  @return forces on @p p1, @p p2, @p p3, @p p4
 */
//<<<<<<< HEAD
//inline int calc_oif_local(Particle *p2, Particle *p1, Particle *p3,
//                          Particle *p4, Bonded_ia_parameters *iaparams,
//                          double force[3], double force2[3], double force3[3],
//                          double force4[3]) // first-fold-then-the-same approach
//{
//  auto const fp2 = unfolded_position(*p2);
//  auto const fp1 = fp2 + get_mi_vector(p1->r.p, fp2);
//  auto const fp3 = fp2 + get_mi_vector(p3->r.p, fp2);
//  auto const fp4 = fp2 + get_mi_vector(p4->r.p, fp2);

//  for (int i = 0; i < 3; i++) {
//    force[i] = 0;
//    force2[i] = 0;
//    force3[i] = 0;
//    force4[i] = 0;
//  }
//=======
inline std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d,
                  Utils::Vector3d>
calc_oif_local(Particle const &p2, Particle const &p1, Particle const &p3,
               Particle const &p4, Bonded_ia_parameters const &iaparams) {

  // first-fold-then-the-same approach
  auto const fp2 = unfolded_position(p2.r.p, p2.l.i, box_geo.length());
  auto const fp1 = fp2 + get_mi_vector(p1.r.p, fp2, box_geo);
  auto const fp3 = fp2 + get_mi_vector(p3.r.p, fp2, box_geo);
  auto const fp4 = fp2 + get_mi_vector(p4.r.p, fp2, box_geo);

  auto const dx = fp2 - fp3;
  auto const len = dx.norm();

  if ((len > 5.0*iaparams.p.oif_local_forces.r0)) {
    return {};
  }


  Utils::Vector3d force1{}, force2{}, force3{}, force4{};
//>>>>>>> espressomd/python

  // non-linear stretching
  if (iaparams.p.oif_local_forces.ks > TINY_OIF_ELASTICITY_COEFFICIENT) {
    auto const dx = fp2 - fp3;
    auto const len = dx.norm();
    auto const dr = len - iaparams.p.oif_local_forces.r0;
    auto const lambda = 1.0 * len / iaparams.p.oif_local_forces.r0;
    auto fac =
        -iaparams.p.oif_local_forces.ks * KS(lambda) * dr; // no normalization
    double val = 0.0;
    auto const r0 = iaparams.p.oif_local_forces.r0;
    if ((len > 2.0*r0) && (len < 2.5*r0)) {
        val = 0.01*r0*r0/((2.0*r0 - len)*(2.0*r0 - len));
        if (val > 150.0) {
            val = 150.0;
        }
        fac = fac*(1.0 + val);
        //printf("%lf %lf %lf",val,r0,len);
        //printf("x");
    }        
    if ((len > 0.25*r0) && (len < 0.5*r0)) {
        val = 0.01*r0*r0/((len - 0.5*r0)*(len - 0.5*r0));
        if (val > 150.0) {
            val = 150.0;
        }
        fac = fac*(1.0 + val);
    }
    auto const f = (fac / len) * dx;
    force2 += f;
    force3 -= f;
  }

  // linear stretching
  if (iaparams.p.oif_local_forces.kslin > TINY_OIF_ELASTICITY_COEFFICIENT) {
    auto const dx = fp2 - fp3;
    auto const len = dx.norm();
    auto const dr = len - iaparams.p.oif_local_forces.r0;
    auto const fac =
        -iaparams.p.oif_local_forces.kslin * dr; // no normalization
    auto const f = (fac / len) * dx;
    force2 += f;
    force3 -= f;
  }

  // viscous force
  if (iaparams.p.oif_local_forces.kvisc >
      TINY_OIF_ELASTICITY_COEFFICIENT) { // to be implemented....
    auto const dx = fp2 - fp3;
    auto const len2 = dx.norm2();
    auto const v_ij = p3.m.v - p2.m.v;

    // Variant A
    // Here the force is in the direction of relative velocity btw points

    // Code:
    // force2 += iaparams.p.oif_local_forces.kvisc * v;
    // force3 -= iaparams.p.oif_local_forces.kvisc * v;

    // Variant B
    // Here the force is the projection of relative velocity btw points onto
    // line btw the points

    // denote p vector between p2 and p3
    // denote v the velocity difference between the points p2 and p3
    // denote alpha the angle between p and v
    // denote x the projected v onto p
    // cos alpha = |x|/|v|
    // cos alpha = (v,p)/(|v||p|)
    // together we get |x|=(v,p)/|p|
    // also, x is along p, so x = |x|.p/|p|
    // so x = p/|p| . (v,p)/|p|
    // altogether x = p . (v,p)/(|p|)^2
    // |p|^2 is stored in len2

    // Code:
    auto const fac = iaparams.p.oif_local_forces.kvisc * (dx * v_ij) / len2;
    auto const f = fac * dx;
    force2 += f;
    force3 -= f;
  }

//<<<<<<< HEAD
  /* bending
     implemented according to updated expressions in book Computational Blood Cell Mechanics, by I.Cimrak and I.Jancigova, see errata...*/
//  if (iaparams->p.oif_local_forces.kb > TINY_OIF_ELASTICITY_COEFFICIENT) {
    // how fp1 - fp4 correspond to points A,B,C,D from the book, Figure A.1:
     //    fp1 -> C
     //    fp2 -> A
     //    fp3 -> B
     //    fp4 -> D
    
//    auto const Nc = Utils::get_n_triangle(fp1, fp2, fp3);    // returns (fp2 - fp1)x(fp3 - fp1), thus Nc = (A - C)x(B - C)  
//    auto const Nd = Utils::get_n_triangle(fp4, fp3, fp2);    // returns (fp3 - fp4)x(fp2 - fp4), thus Nd = (B - D)x(A - D)  
//=======
  /* bending
     forceT1 is restoring force for triangle p1,p2,p3 and force2T restoring
     force for triangle p2,p3,p4 p1 += forceT1; p2 -= 0.5*forceT1+0.5*forceT2;
     p3 -= 0.5*forceT1+0.5*forceT2; p4 += forceT2; */
  if (iaparams.p.oif_local_forces.kb > TINY_OIF_ELASTICITY_COEFFICIENT) {
    // how fp1 - fp4 correspond to points A,B,C,D from the book, Figure A.1:
     //    fp1 -> C
     //    fp2 -> A
     //    fp3 -> B
     //    fp4 -> D
    auto const Nc = Utils::get_n_triangle(fp1, fp2, fp3);    // returns (fp2 - fp1)x(fp3 - fp1), thus Nc = (A - C)x(B - C)  
    auto const Nd = Utils::get_n_triangle(fp4, fp3, fp2);    // returns (fp3 - fp4)x(fp2 - fp4), thus Nd = (B - D)x(A - D)  
//>>>>>>> espressomd/python

    auto const phi = Utils::angle_btw_triangles(fp1, fp2, fp3, fp4);
    auto const aa = (phi - iaparams.p.oif_local_forces
                               .phi0); // no renormalization by phi0, to be
                                       // consistent with Krueger and Fedosov
//<<<<<<< HEAD
    auto const BminA = fp3 - fp2;
    auto const fac = iaparams.p.oif_local_forces.kb * aa;
    auto const factorFaNc = (fp2 - fp3) * (fp1 - fp3) / BminA.norm() / Nc.norm2();
    auto const factorFaNd = (fp2 - fp3) * (fp4 - fp3) / BminA.norm() / Nd.norm2();
    auto const factorFbNc = (fp2 - fp3) * (fp2 - fp1) / BminA.norm() / Nc.norm2();
    auto const factorFbNd = (fp2 - fp3) * (fp2 - fp4) / BminA.norm() / Nd.norm2();
    
//    for (int i = 0; i < 3; i++) {
//      force[i] -= fac * BminA.norm()/Nc.norm2() * Nc[i];                // Fc
//      force2[i] +=  fac * (factorFaNc * Nc[i] + factorFaNd * Nd[i]);    // Fa
//      force3[i] +=  fac * (factorFbNc * Nc[i] + factorFbNd * Nd[i]);    // Fb
//      force4[i] -= fac * BminA.norm()/Nd.norm2() * Nd[i];               // Fd
//    }
//=======
//    auto const fac = iaparams.p.oif_local_forces.kb * aa;
//    auto const f = 0.5 * fac * n1 + 0.5 * fac * n2;

    force1 -= fac * BminA.norm()/Nc.norm2() * Nc;
    force2 += fac * (factorFaNc * Nc + factorFaNd * Nd);
    force3 += fac * (factorFbNc * Nc + factorFbNd * Nd);
    force4 -= fac * BminA.norm()/Nd.norm2() * Nd;
//>>>>>>> espressomd/python
  }

  /* local area
   * for both triangles, only 1/3 of calculated forces are added, because each
   * triangle will enter this calculation 3 times (one time per edge).
   * Proportional distribution of forces, implemented according to
   * @cite jancigova16a.
   */
  if (iaparams.p.oif_local_forces.kal > TINY_OIF_ELASTICITY_COEFFICIENT) {

    auto handle_triangle = [](double kal, double A0, Utils::Vector3d const &fp1,
                              Utils::Vector3d const &fp2,
                              Utils::Vector3d const &fp3,
                              Utils::Vector3d &force1, Utils::Vector3d &force2,
                              Utils::Vector3d &force3) {
      auto const h = (1. / 3.) * (fp1 + fp2 + fp3);
      auto const A = Utils::area_triangle(fp1, fp2, fp3);
      auto const t = sqrt(A / A0) - 1.0;

      auto const m1 = h - fp1;
      auto const m2 = h - fp2;
      auto const m3 = h - fp3;

      auto const m1_length2 = m1.norm2();
      auto const m2_length2 = m2.norm2();
      auto const m3_length2 = m3.norm2();

      auto const fac =
          kal * A0 * (2 * t + t * t) / (m1_length2 + m2_length2 + m3_length2);

      // local area force for p1
      force1 += (fac / 3.0) * m1;
      force2 += (fac / 3.0) * m2;
      force3 += (fac / 3.0) * m3;
    };

    handle_triangle(iaparams.p.oif_local_forces.kal,
                    iaparams.p.oif_local_forces.A01, fp1, fp2, fp3, force1,
                    force2, force3);
    handle_triangle(iaparams.p.oif_local_forces.kal,
                    iaparams.p.oif_local_forces.A02, fp2, fp3, fp4, force2,
                    force3, force4);
  }
  return std::make_tuple(force2, force1, force3, force4);
}

#endif
