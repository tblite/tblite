/* This file is part of tblite.
 * SPDX-Identifier: LGPL-3.0-or-later
 *
 * tblite is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tblite is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tblite.  If not, see <https://www.gnu.org/licenses/>.
**/

/**@file tblite/container.h
 * @brief
 * Provides an interaction container which can be added to a #tblite_calculator.
 */

#pragma once

#include "tblite/macros.h"
#include "tblite/structure.h"
#include "tblite/calculator.h"
#include "tblite/context.h"

/// Interaction container
typedef struct _tblite_container* tblite_container;

/// Create new electric field container
///
/// @param efield: Electric field in atomic units (Hartree/(Bohr*e)), shape: [3]
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_electric_field(double* efield);

/// Create a new external point-charge embedding container.
/// Coordinates are Cartesian in Bohr, charges in elementary charge, and
/// hardness in Hartree. The hardness is combined with each xTB shell value
/// using harmonic averaging for GFN1/IPEA1 and arithmetic averaging for GFN2.
/// Only QM-to-point-charge coupling is included; external-charge interactions
/// are handled by the embedding or MM driver.
/// @param ctx Context handle
/// @param mol Molecular structure
/// @param calc GFN1-xTB, GFN2-xTB, or IPEA1-xTB calculator
/// @param npc Number of point charges
/// @param charge Point charges, shape [npc]
/// @param xyz Point-charge coordinates, Fortran order shape [3, npc]
/// @param hardness Point-charge hardness; zero selects the method's element hardness
/// @param numbers Atomic numbers used when hardness is zero
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_point_charges(tblite_context ctx,
                         tblite_structure mol,
                         tblite_calculator calc,
                         int npc,
                         double* charge,
                         double* xyz,
                         double* hardness,
                         int* numbers);

/// Create new spin polarization container using internal parameters
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @param calc: Calculator instance
/// @param wscale: Scaling factor for spin polarization (default: 1)
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_spin_polarization(tblite_context ctx,
                             tblite_structure mol,
                             tblite_calculator calc,
                             double wscale);


/// Add container to calculator object.
///
/// Note: Ownership is transferred and container handle is destroyed after function call
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param cont: Interaction container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_calculator_push_back(tblite_context ctx,
                            tblite_calculator calc,
                            tblite_container* cont);

/// Delete container handle
///
/// @param cont: Container handle
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_container(tblite_container* cont);
