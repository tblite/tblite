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

#ifdef __cplusplus
#define TBLITE_API_ENTRY extern "C"
#ifndef TBLITE_CFFI
#include <cstdint>
#endif
#else
#define TBLITE_API_ENTRY extern
#ifndef TBLITE_CFFI
#include <stdbool.h>
#include <stdint.h>
#endif
#endif
#define TBLITE_API_CALL
