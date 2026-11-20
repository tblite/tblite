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

/**@file tblite/table.h
 * @brief
 * Provides a representation of a generic table data structure.
 *
 * Used to mirror the data available in the #tblite_param object. It aims to provide a
 * programmatic accessible representation of the parametrization records.
 */

#pragma once

#include "tblite/macros.h"
#include "tblite/error.h"

/// Handle for holding a table data structure
///
/// The table can either own its data or reference another table. Table references
/// can be created from existing table data structures or by adding new table entries
/// into an existing table, which returns a reference to the newly created table.
typedef struct _tblite_table* tblite_table;

/// Handle for holding an array data structure
typedef struct _tblite_array* tblite_array;

/// Supported value kinds for table entries
enum tblite_table_value_type {
    TBLITE_TABLE_VALUE_TYPE_NONE = 0,
    TBLITE_TABLE_VALUE_TYPE_BOOL = 1,
    TBLITE_TABLE_VALUE_TYPE_INT = 2,
    TBLITE_TABLE_VALUE_TYPE_DOUBLE = 3,
    TBLITE_TABLE_VALUE_TYPE_CHAR = 4,
    TBLITE_TABLE_VALUE_TYPE_ARRAY = 5,
    TBLITE_TABLE_VALUE_TYPE_TABLE = 6
};

/// Create new data table object
///
/// @param table: Table object to reference in new table (optional)
/// @return: New table data structure
TBLITE_API_ENTRY tblite_table TBLITE_API_CALL
tblite_new_table(tblite_table* table);

/// Delete a data table object
///
/// @param table: Table object to be deleted
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_table(tblite_table* table);

/// Set floating point number to data table
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Double value or array to add to table
/// @param n: Number of entries to add, 0 for adding scalars
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_double(tblite_error error,
                        tblite_table table,
                        char key[],
                        double* value,
                        int n);

/// Set integer number to data table (use int64_t rather than long)
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Integer value or array to add to table
/// @param n: Number of entries to add, 0 for adding scalars
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_int64_t(tblite_error error,
                         tblite_table table,
                         char key[],
                         int64_t* value,
                         int n);

/// Set boolean value to data table
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Boolean value or array to add to table
/// @param n: Number of entries to add, 0 for adding scalars
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_bool(tblite_error error,
                      tblite_table table,
                      char key[],
                      bool* value,
                      int n);

/// Set character string to data table
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Character value or array to add to table
/// @param n: Number of entries to add, 0 for adding scalars
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_char(tblite_error error,
                      tblite_table table,
                      char key[],
                      char (* value)[],
                      int n);

/// Create new subtable in existing data table
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to add new subtable at
/// @return New table data structure
TBLITE_API_ENTRY tblite_table TBLITE_API_CALL
tblite_table_add_table(tblite_error error,
                       tblite_table table,
                       char key[]);

/// Create new data array object
TBLITE_API_ENTRY tblite_array TBLITE_API_CALL
tblite_new_array(void);

/// Delete a data array object
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_array(tblite_array* array);

/// Append a double value to a data array
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_array_push_back_double(tblite_error error,
                              tblite_array array,
                              double value);

/// Append an integer value to a data array
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_array_push_back_int64_t(tblite_error error,
                               tblite_array array,
                               int64_t value);

/// Append a boolean value to a data array
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_array_push_back_bool(tblite_error error,
                            tblite_array array,
                            bool value);

/// Append a character string to a data array
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_array_push_back_char(tblite_error error,
                            tblite_array array,
                            char value[]);

/// Return the number of entries in a data array
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_array_size(tblite_error error,
                  tblite_array array);

/// Return the value kind of an array entry
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_array_get_type(tblite_error error,
                      tblite_array array,
                      int index);

/// Return an array entry as a double value
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_array_get_double(tblite_error error,
                        tblite_array array,
                        int index,
                        double* value);

/// Return an array entry as an integer value
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_array_get_int64_t(tblite_error error,
                         tblite_array array,
                         int index,
                         int64_t* value);

/// Return an array entry as a boolean value
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_array_get_bool(tblite_error error,
                       tblite_array array,
                       int index,
                       bool* value);

/// Return an array entry as a character string
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_array_get_char(tblite_error error,
                      tblite_array array,
                      int index,
                      char value[],
                      int n);

/// Set the value of a table entry to a data array
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_array(tblite_error error,
                       tblite_table table,
                       char key[],
                       tblite_array array);

/// Return the value kind of a table entry
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_table_get_type(tblite_error error,
                      tblite_table table,
                      char key[]);

/// Return a scalar bool value from a table entry
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_get_bool(tblite_error error,
                      tblite_table table,
                      char key[],
                      bool* value);

/// Return a scalar integer value from a table entry
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_get_int64_t(tblite_error error,
                         tblite_table table,
                         char key[],
                         int64_t* value);

/// Return a scalar floating point value from a table entry
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_get_double(tblite_error error,
                        tblite_table table,
                        char key[],
                        double* value);

/// Return a scalar character value from a table entry
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_get_char(tblite_error error,
                      tblite_table table,
                      char key[],
                      char value[],
                      int n);

/// Return a child table from a table entry
TBLITE_API_ENTRY tblite_table TBLITE_API_CALL
tblite_table_get_table(tblite_error error,
                       tblite_table table,
                       char key[]);

/// Return an array value from a table entry
TBLITE_API_ENTRY tblite_array TBLITE_API_CALL
tblite_table_get_array(tblite_error error,
                       tblite_table table,
                       char key[]);

/// Return the number of entries in a table
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_table_get_n_keys(tblite_error error,
                        tblite_table table);

/// Return the key of a table entry at the given index
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_get_key(tblite_error error,
                     tblite_table table,
                     int index,
                     char key[],
                     int n);

/// Serialize a data table to a TOML file
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param filename: Output file name
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_dump_table(tblite_error error,
                  tblite_table table,
                  char filename[]);

/*
 * Type generic macros
 */

/// Generic setter based on the type of the value
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Value to set at key
#define tblite_table_set_value(error, table, key, value, ...) \
    _Generic((value), \
     double*: tblite_table_set_double, \
    int64_t*: tblite_table_set_int64_t, \
       bool*: tblite_table_set_bool, \
   char(*)[]: tblite_table_set_char \
            )((error), (table), (key), (value), __VA_ARGS__)
