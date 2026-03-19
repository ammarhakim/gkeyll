#pragma once

// Needs more thought on how to get this right!!

// Defintions and interface for post-processing and analysis tool
// chain

#include <gkyl_array.h>
#include <gkyl_dynvec.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Opaque object to handle pgkyl internal state
typedef struct pgkyl_state_tag* pgkyl_state;

// All commands must have this signature
typedef bool (*pgkyl_cmd)(pgkyl_state ps, void *ctx);

// Command and command context pair
struct pgkyl_cmd_ctx_pair {
  pgkyl_cmd cmd; // command
  void *ctx; //command context
};

struct pgkyl_cmd_status {
  bool status; // status of command execution
  int outid; // ID of output dataset
};

// Data store for an array: holds the local range, the grid and the
// actual array data itself
struct pgkyl_array_data {
  struct gkyl_range range; // local range for array
  struct gkyl_rect_grid grid; // grid on which array is defined
  const struct gkyl_array *array; // array 
};

// Data store for dynvec.
struct pgkyl_dynvec_data {
  const gkyl_dynvec vec; // dynvector
};

/**
 * Initialize a new pgkyl state object
 *
 * @return New pgkyl state
 */
pgkyl_state pgkyl_state_new(void);

/**
 * Call this method before appending data. Finish using
 * pgkyl_finish_appending() call.
 *
 * @param ps pgkyl state
 * @return ID of new dataset
 */
int pgkyl_start_appending(pgkyl_state ps);

/**
 * Call this method after finishing appending all data.
 *
 * @param ps pgkyl state
 */
void pgkyl_finish_appending(pgkyl_state ps);

/**
 * Number of datasets in pgkyl state.
 *
 * @param pd pgkyl state
 * @return number of datasets
 */
int pgkyl_num_datasets(const pgkyl_state ps);

/**
 * Set the active dataset to the specified ID
 *
 * @param ps pgkyl state
 * @param dsid Dataset ID to activate
 * @return true if activation succeeded. False otherwise
 */
bool pgkyl_set_active_dataset(pgkyl_state ps, int dsid);

/**
 * Get the active dataset ID
 *
 * @param ps pgkyl state
 * @return ID of active dataset
 */
int pgkyl_get_active_dataset(const pgkyl_state ps);

/**
 * Appending array data to dataset. Only safe to use after a
 * pgkyl_start_appending has been called.
 *
 * @param ps pgkyl state
 * @param range Local range object
 * @param grid Grid on which array lives
 * @param arr Array object
 */
void pgkyl_append_array(pgkyl_state ps, const struct gkyl_range *range,
  const struct gkyl_rect_grid *grid, const struct gkyl_array *arr);

/**
 * Execute sequence of commands in the @a cmd_list sequentially. The
 * final entry in @a cmd_list should be { 0, 0 }
 *
 * @param ps pgkyl state
 * @param cmd Command to run
 * @return Status and dataset ID
 */
struct pgkyl_cmd_status pgkyl_execute_cmd(pgkyl_state ps, struct pgkyl_cmd_ctx_pair cmd);

/**
 * Execute sequence of commands in the @a cmd_list sequentially. The
 * final entry in @a cmd_list should be { 0, 0 }. Each command uses the dataset ID of the previous command as it input.
 *
 * @param ps pgkyl state
 * @param cmd_list List of commands to execute
 * @return Status and dataset ID of last command run
 */
struct pgkyl_cmd_status pgkyl_execute_cmd_chain(pgkyl_state ps, const struct pgkyl_cmd_ctx_pair *cmd_list);

/**
 * Release pgkyl state object
 *
 * @param ps pgkyl state to release
 */
void pgkyl_state_release(pgkyl_state ps);
