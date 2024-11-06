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

/**@file tblite/macros.h
 * @brief
 * General macro definitions for the tblite C API bindings.
 */

#pragma once

/** @def TBLITE_API_ENTRY
 * Defines an external function exported by the Fortran library.
 */

/** @def TBLITE_API_CALL
 * Macro to define calling convention.
 */

/** @def TBLITE_CFFI
 * Guard macro for CFFI preprocessing of the header files.
 *
 * Use this macro to conditionally enable or disable code snippets for the pass
 * of the CFFI generation step to obtain the Python extension module.
 * Notably, header includes should be removed if the macro is defined to avoid
 * creating bindings for the system library. Furthermore, callbacks need a
 * special extern "Python" declaration which should conditionally included.
 */

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

#ifndef TBLITE_API_CALL
#define TBLITE_API_CALL
#endif

#define vrg_cnt(vrg1,vrg2,vrg3,vrg4,vrg5,vrg6,vrg7,vrg8,vrgN, ...) vrgN
#define vrg_argn(...)  vrg_cnt(__VA_ARGS__, 8, 7, 6, 5, 4, 3, 2, 1, 0)
#define vrg_cat0(x,y)  x ## y
#define vrg_cat(x,y)   vrg_cat0(x,y)

#define vrg(vrg_f,...) vrg_cat(vrg_f, vrg_argn(__VA_ARGS__))(__VA_ARGS__)
