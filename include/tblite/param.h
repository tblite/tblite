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
#include "tblite/table.h"

/// Parametrization records
typedef struct _tblite_param* tblite_param;

/// Create new parametrization records object
TBLITE_API_ENTRY tblite_param TBLITE_API_CALL
tblite_new_param(void);

/// Delete a parametrization records object
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_param(tblite_param* /* param */);

/// Load parametrization records from data table
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_load_param(tblite_error /* error */,
                  tblite_param /* param */,
                  tblite_table /* table */);

/// Dump parametrization records to data table
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_dump_param(tblite_error /* error */,
                  tblite_param /* param */,
                  tblite_table /* table */);

/// Export GFN2-xTB parametrization records
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_export_gfn2_param(tblite_error /* error */,
                         tblite_param /* param */);

/// Export GFN1-xTB parametrization records
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_export_gfn1_param(tblite_error /* error */,
                         tblite_param /* param */);

/// Export IPEA1-xTB parametrization records
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_export_ipea1_param(tblite_error /* error */,
                          tblite_param /* param */);
