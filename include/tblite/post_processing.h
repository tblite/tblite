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

/**@file tblite/post_processing.h
 * @brief
 * Provides access to the post processing type as to add post processing conatiners using strings.
 *
 * Gives access to post processing container to add methods to the current post processing list.
 */
#pragma once
#include "tblite/macros.h"

// Post Processing Container
typedef struct _tblite_post_processing* tblite_post_processing;

/// Construct post processing container
///
/// @param charptr: String of the post processing desired
/// @return New post processing instance
TBLITE_API_ENTRY tblite_post_processing TBLITE_API_CALL
tblite_new_post_processing(char* charptr);

/// Push Back new conatiner to post processing construct
///
/// @param post_proc: Post Processing instance 
/// @param charptr: String of the post processing desired
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_push_back_post_processing(tblite_post_processing* post_proc, 
                                 char* charptr);

/// Delete calculator
///
/// @param post_proc: Post Processing instance
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_post_processing(tblite_post_processing* post_proc);