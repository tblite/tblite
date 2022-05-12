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

/// Define callback function for use in custom logger
typedef void (*tblite_logger_callback)(char*, int, void*);

#ifdef TBLITE_CFFI
extern "Python" void TBLITE_API_CALL
logger_callback(char* msg, int len, void* user_data);
#endif

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

/// Set custom logger function
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_context_logger(tblite_context /* ctx */,
                          tblite_logger_callback /* callback */,
                          void* /* userdata */);

/// Enable colorful output
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_context_color(tblite_context /* ctx */,
                         int /* color */);
