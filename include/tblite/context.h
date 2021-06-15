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

/// Context manager for the library usage
typedef struct _tblite_context* tblite_context;

/// Create new calculation environment object
TBLITE_API_ENTRY tblite_context TBLITE_API_CALL
tblite_new_context(void);

/// Delete a calculation environment object
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_context(tblite_context* /* ctx */);

/// Check calculation environment status
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_check_context(tblite_context /* ctx */);

/// Get error message from calculation environment
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_context_error(tblite_context /* ctx */,
                         char* /* buffer */,
                         const int* /* buffersize */);
