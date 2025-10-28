#include <gkyl_alloc.h>
#include <gkyl_basis.h>

void
gkyl_cart_modal_basis_release(struct gkyl_basis *basis)
{
  gkyl_free(basis);
}

void
gkyl_cart_modal_basis_release_cu(struct gkyl_basis *basis)
{
  gkyl_cu_free(basis);
}

