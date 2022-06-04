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

/// Interaction container
typedef struct _tblite_container* tblite_container;

/// Create new electric field container
///
/// @param efield: Electric field in atomic units (Hartree/(Bohr*e)), shape: [3]
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_electric_field(double* efield);

/// Delete container handle
///
/// @param cont: Container handle
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_container(tblite_container* cont);
