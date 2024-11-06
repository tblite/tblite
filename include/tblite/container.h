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


#define tblite_new_cpcm_solvation(ctx, mol, calc, x)            \
                        _Generic((x),                           \
                                char*                           \
                                : tblite_new_cpcm_solvation_solvent,\
                                double                          \
                                : tblite_new_cpcm_solvation_epsilon \
                                ) (ctx, mol, calc, x)

/// Create new CPCM implicit solvation container using internal parameters
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @param calc: Calculator instance
/// @param eps: epsilon value for solvent
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_cpcm_solvation_epsilon(tblite_context ctx,
                          tblite_structure mol,
                          tblite_calculator calc,
                          double eps);

/// Create new CPCM implicit solvation container using internal parameters
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @param calc: Calculator instance
/// @param solvent: Solvent to be modelled, can be given as name of solvent or epsilon value
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_cpcm_solvation_solvent(tblite_context ctx,
                          tblite_structure mol,
                          tblite_calculator calc,
                          char* solvent);

// standard/default argument for reference state of solvation model
char * std_refstr = "gsolv";
//Macro to call tblite_new_alpb_solvation with or without/default reference string
#define tblite_new_alpb_solvation(...) vrg(tblite_new_alpb_solvation, __VA_ARGS__)
//in case only 4 arguments are used add default argument to macro call
#define tblite_new_alpb_solvation4(ctx, mol, calc, x) tblite_new_alpb_solvation_(ctx, mol, calc, x, std_refstr)
#define tblite_new_alpb_solvation5(ctx, mol, calc, x, refstr) tblite_new_alpb_solvation_(ctx, mol, calc, x, refstr)

#define tblite_new_alpb_solvation_(ctx, mol, calc, x, refstr)            \
                        _Generic((x),                           \
                                char*                           \
                                : tblite_new_alpb_solvation_solvent,\
                                double                          \
                                : tblite_new_alpb_solvation_epsilon) (ctx, mol, calc, x, refstr)

/// Create new ALPB implicit solvation container using internal parameters
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @param calc: Calculator instance
/// @param solvent: Solvent to be modelled, can be given as name of solvent or epsilon value
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_alpb_solvation_solvent(tblite_context ctx,
                          tblite_structure mol,
                          tblite_calculator calc,
                          char* solvent,
                          char* refstr);

/// Create new ALPB implicit solvation container using internal parameters
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @param calc: Calculator instance
/// @param eps: epsilon value of solvent
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_alpb_solvation_epsilon(tblite_context ctx,
                          tblite_structure mol,
                          tblite_calculator calc,
                          double eps,
                          char* refstr);

//Macro to call tblite_new_aöllpbsolvation with or without/default reference string
#define tblite_new_gbsa_solvation(...) vrg(tblite_new_gbsa_solvation, __VA_ARGS__)
//in case only 4 arguments are used add default argument to macro call
#define tblite_new_gbsa_solvation4(ctx, mol, calc, x) tblite_new_gbsa_solvation_(ctx, mol, calc, x, std_refstr)
#define tblite_new_gbsa_solvation5(ctx, mol, calc, x, refstr) tblite_new_gbsa_solvation_(ctx, mol, calc, x, refstr)

#define tblite_new_gbsa_solvation_(ctx, mol, calc, x, refstr)            \
                        _Generic((x),                           \
                                char*                           \
                                : tblite_new_gbsa_solvation_solvent,\
                                double                          \
                                : tblite_new_gbsa_solvation_epsilon) (ctx, mol, calc, x, refstr)

/// Create new ALPB implicit solvation container using internal parameters
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @param calc: Calculator instance
/// @param solvent: Solvent to be modelled, can be given as name of solvent or epsilon value
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_gbsa_solvation_solvent(tblite_context ctx,
                          tblite_structure mol,
                          tblite_calculator calc,
                          char* solvent,
                          char* refstr);

/// Create new ALPB implicit solvation container using internal parameters
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @param calc: Calculator instance
/// @param eps: epsilon value of solvent
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_gbsa_solvation_epsilon(tblite_context ctx,
                          tblite_structure mol,
                          tblite_calculator calc,
                          double eps,
                          char* refstr);


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
