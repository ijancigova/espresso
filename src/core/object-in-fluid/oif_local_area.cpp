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

/** \file oif_local_area.hpp
 *  Routines to calculate the OIF_LOCAL_AREA
 *  for three particles that form a triangle. (Dupin2007)
 *  \ref forces.cpp
 */

#include "communication.hpp"
#include "oif_local_area.hpp"

/** set parameters for the OIF_LOCAL_AREA potential.*/

int oif_local_forces_set_params(int bond_type, double A0, double kal)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.oif_local_area.A0 = A0;
  bonded_ia_params[bond_type].p.oif_local_area.kal = kal;
    
  bonded_ia_params[bond_type].type = BONDED_IA_OIF_LOCAL_AREA;
  bonded_ia_params[bond_type].num  = 2;

  mpi_bcast_ia_params(bond_type, -1); 

  return ES_OK;
}

