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

/// Create new result contains
TBLITE_API_ENTRY tblite_result TBLITE_API_CALL
tblite_copy_result(tblite_result /* res */);

/// Delete a calculation environment object
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_result(tblite_result* /* res */);

/// Retrieve number of atoms from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_number_of_atoms(tblite_error /* err */,
                                  tblite_result /* res */,
                                  int* /* nat */);

/// Retrieve number of shells from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_number_of_shells(tblite_error /* err */,
                                   tblite_result /* res */,
                                   int* /* nsh */);

/// Retrieve number of orbitals from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_number_of_orbitals(tblite_error /* err */,
                                     tblite_result /* res */,
                                     int* /* nao */);

/// Retrieve energy from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_energy(tblite_error /* err */,
                         tblite_result /* res */,
                         double* /* energy */);

/// Retrieve atom-resolved energies from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_energies(tblite_error /* err */,
                           tblite_result /* res */,
                           double* /* energies[nat] */);

/// Retrieve gradient from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_gradient(tblite_error /* err */,
                           tblite_result /* res */,
                           double* /* gradient[nat][3] */);

/// Retrieve virial from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_virial(tblite_error /* err */,
                         tblite_result /* res */,
                         double* /* sigma[3][3] */);

/// Retrieve atomic charges from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_charges(tblite_error /* err */,
                          tblite_result /* res */,
                          double* /* charges[nat] */);

/// Retrieve dipole moment from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_dipole(tblite_error /* err */,
                         tblite_result /* res */,
                         double* /* dipole[3] */);

/// Retrieve traceless quadrupole moment from result container (packed xx, xy, yy, xz, yz, zz)
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_quadrupole(tblite_error /* err */,
                             tblite_result /* res */,
                             double* /* quadrupole[6] */);

/// Retrieve orbital energies from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_orbital_energies(tblite_error /* err */,
                                   tblite_result /* res */,
                                   double* /* emo[nao] */);

/// Retrieve orbital occupations from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_orbital_occupations(tblite_error /* err */,
                                      tblite_result /* res */,
                                      double* /* occ[nao] */);

/// Retrieve orbital coefficients from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_orbital_coefficients(tblite_error /* err */,
                                       tblite_result /* res */,
                                       double* /* cmo[nao][nao] */);

/// Retrieve density matrix from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_density_matrix(tblite_error /* err */,
                                 tblite_result /* res */,
                                 double* /* pmat[nao][nao] */);

/// Retrieve overlap matrix from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_overlap_matrix(tblite_error /* err */,
                                 tblite_result /* res */,
                                 double* /* smat[nao][nao] */);

/// Retrieve Hamiltonian matrix from result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_hamiltonian_matrix(tblite_error /* err */,
                                     tblite_result /* res */,
                                     double* /* hmat[nao][nao] */);
