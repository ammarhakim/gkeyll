#pragma once
// Per-surface metadata for the generated nodal surface-evaluation tables
// (gkyl_vlasov_surf_tables_*.h) and their per-node kernels. One instance per
// surface (vst_<tag>_ph_v<d>_meta / vst_<tag>_ph_x<d>_meta) is emitted next
// to the tables. Loop bounds and array extents are compiled into the node
// functions as literals; this record is for the host-side launch logic
// (block sizing, shared-memory sizing, output offsets).
struct gkyl_vlasov_surf_meta {
  int nb; // phase-basis coefficients (f/hamil length; stage-1 loop bound)
  int na; // distinct outer (chi) shapes: Cm columns, G-column length
  int nmv; // distinct inner (psi) shapes: shared V columns
  int no; // outer nodes (Cm rows; i extent)
  int ni; // inner nodes (V rows; j extent; flat node n = i*ni + j)
  int ns; // sparse-Hamiltonian support size (vel_sparse_idx length)
  int fhat_off; // base offset of this surface's block in Fhat_nodal/alpha_quad
};

// Per-surface metadata for the from-flux consumer projection tables
// (vst_<tag>_prj_{v,x}<d>_*): the nodal flux projected onto the p=polyOrder
// surface basis and applied to the volume increment. One instance per surface
// (vst_<tag>_prj_v<d>_meta / vst_<tag>_prj_x<d>_meta) is emitted next to the
// tables; the per-mode launch uses nk threads per surface (distinct surface
// modes write disjoint out coefficients, so no atomics are needed).
struct gkyl_vlasov_prj_meta {
  int nk; // surface-basis modes (mode-function index extent)
  int np; // volume-basis coefficients (out length; out_mode/out_cl/out_cr length)
  int na; // distinct outer (chi) shapes: Cw columns
  int nmv; // distinct inner (psi) shapes: Vw columns
  int no; // outer nodes (Cw rows; i extent)
  int ni; // inner nodes (Vw rows; j extent; flat node n = i*ni + j)
  int fhat_off; // base offset of this surface's block in Fhat_nodal
};
