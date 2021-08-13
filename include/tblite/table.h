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

/// Data table object
typedef struct _tblite_table* tblite_table;

/// Create new data table object
TBLITE_API_ENTRY tblite_table TBLITE_API_CALL
tblite_new_table(void);

/// Delete a data table object
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_table(tblite_table* /* table */);

/// Set floating point number to data table
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_double(tblite_error /* error */,
                        tblite_table /* table */,
                        char[] /* key */,
                        double* /* value */,
                        int /* n */);

/// Set integer number to data table
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_long(tblite_error /* error */,
                      tblite_table /* table */,
                      char[] /* key */,
                      long* /* value */,
                      int /* n */);

/// Set boolean value to data table
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_bool(tblite_error /* error */,
                      tblite_table /* table */,
                      char[] /* key */,
                      bool* /* value */,
                      int /* n */);

/// Set character string to data table
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_char(tblite_error /* error */,
                      tblite_table /* table */,
                      char[] /* key */,
                      char(*)[] /* value */,
                      int /* n */);

/// Set data table to data table, value table reference is destroyed on success
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_table(tblite_error /* error */,
                       tblite_table /* table */,
                       char[] /* key */,
                       tblite_table* /* value */);
