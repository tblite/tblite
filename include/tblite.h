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

#include "tblite/error.h"
#include "tblite/container.h"
#include "tblite/context.h"
#include "tblite/structure.h"
#include "tblite/calculator.h"
#include "tblite/result.h"
#include "tblite/table.h"
#include "tblite/param.h"
#include "tblite/version.h"

#define tblite_delete(ptr) _Generic((ptr), \
                       tblite_error: tblite_delete_error, \
                   tblite_container: tblite_delete_container, \
                     tblite_context: tblite_delete_context, \
                   tblite_structure: tblite_delete_structure, \
                  tblite_calculator: tblite_delete_calculator, \
                      tblite_result: tblite_delete_result, \
                       tblite_table: tblite_delete_table, \
                       tblite_param: tblite_delete_param \
                                   )(&ptr)

