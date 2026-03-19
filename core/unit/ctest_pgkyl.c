#include <acutest.h>
#include <stdio.h>

#include <pgkyl_ifc.h>

void
test_basic(void)
{
  pgkyl_state ps = pgkyl_state_new();

  int shape[] = { 10, 10 };
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  double lower[] = { 0.0, 0.0 };
  double upper[] = { 1.0, 1.0 };
  int cells[] = { shape[0], shape[1] };
  
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  
  int aid = pgkyl_start_appending(ps);
  pgkyl_append_array(ps, &range, &grid, arr);
  pgkyl_finish_appending(ps);

  TEST_CHECK( aid == pgkyl_get_active_dataset(ps) );

  int bid = pgkyl_start_appending(ps);
  pgkyl_append_array(ps, &range, &grid, arr);
  pgkyl_append_array(ps, &range, &grid, arr);  
  pgkyl_finish_appending(ps);

  TEST_CHECK( bid == pgkyl_get_active_dataset(ps) );
  
  TEST_CHECK( 2 == pgkyl_num_datasets(ps) );
  
  gkyl_array_release(arr);
  pgkyl_state_release(ps);
}

TEST_LIST = {
  { "test_basic", test_basic },
  { NULL, NULL },
};
