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

#pragma once

#include "tblite/macros.h"
#include "tblite/error.h"

/// Container to for storing and handling calculation results
typedef struct _tblite_result* tblite_result;

/// Create new result contains
TBLITE_API_ENTRY tblite_result TBLITE_API_CALL
tblite_new_result(void);

/// Delete a calculation environment object
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_result(tblite_result* /* res */);

/// Retrieve energy from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_energy(tblite_error /* err */,
                         tblite_result /* res */,
                         double* /* energy */);

/// Retrieve gradient from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_gradient(tblite_error /* err */,
                           tblite_result /* res */,
                           double* /* gradient[n][3] */);

/// Retrieve virial from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_virial(tblite_error /* err */,
                         tblite_result /* res */,
                         double* /* sigma[3][3] */);
