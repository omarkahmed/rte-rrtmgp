
#pragma once

#include "const.h"


// increment 2stream by 2stream
extern "C" void inc_2stream_by_2stream_bybnd(int ncol, int nlay, int ngpt,
                                             real *tau1_p, real *ssa1_p, real *g1_p,
                                             real *tau2_p, real *ssa2_p, real *g2_p,
                                             int nbnd, int *gpt_lims_p);


// Incrementing when the second set of optical properties is defined at lower spectral resolution
//   (e.g. by band instead of by gpoint)
extern "C" void inc_1scalar_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real *tau1_p, real *tau2_p,             
                                             int nbnd, int *gpt_lims_p);


// Delta-scale
//   f = g*g
extern "C" void delta_scale_2str_k(int ncol, int nlay, int ngpt, real *tau_p, real *ssa_p, real *g_p);


// Delta-scaling, provided only for two-stream properties at present
// -------------------------------------------------------------------------------------------------
// Delta-scale
//   user-provided value of f (forward scattering)
extern "C" void delta_scale_2str_f_k(int ncol, int nlay, int ngpt, real *tau_p, real *ssa_p, real *g_p, real *f_p);


// Addition of optical properties: the first set are incremented by the second set.
//
//   There are three possible representations of optical properties (scalar = optical depth only;
//   two-stream = tau, single-scattering albedo, and asymmetry factor g, and
//   n-stream = tau, ssa, and phase function moments p.) Thus we need nine routines, three for
//   each choice of representation on the left hand side times three representations of the
//   optical properties to be added.
//
//   There are two sets of these nine routines. In the first the two sets of optical
//   properties are defined at the same spectral resolution. There is also a set of routines
//   to add properties defined at lower spectral resolution to a set defined at higher spectral
//   resolution (adding properties defined by band to those defined by g-point)
extern "C" void increment_1scalar_by_1scalar(int ncol, int nlay, int ngpt, real *tau1_p, real *tau2_p);


// increment 1scalar by 2stream
extern "C" void increment_1scalar_by_2stream(int ncol, int nlay, int ngpt, real *tau1_p, real *tau2_p, real *ssa2_p);


// increment 1scalar by nstream
extern "C" void increment_1scalar_by_nstream(int ncol, int nlay, int ngpt, real *tau1_p, real *tau2_p, real *ssa2_p);


// increment 2stream by 1scalar
extern "C" void increment_2stream_by_1scalar(int ncol, int nlay, int ngpt, real *tau1_p, real *ssa1_p, real *tau2_p);


// increment 2stream by 2stream
extern "C" void increment_2stream_by_2stream(int ncol, int nlay, int ngpt, real *tau1_p, real *ssa1_p, real *g1_p, real *tau2_p, real *ssa2_p, real *g2_p);


// increment 2stream by nstream
extern "C" void increment_2stream_by_nstream(int ncol, int nlay, int ngpt, int nmom2, real *tau1_p, real *ssa1_p, real *g1_p, real *tau2_p, real *ssa2_p, real *p2_p);


// increment nstream by 1scalar
extern "C" void increment_nstream_by_1scalar(int ncol, int nlay, int ngpt, real *tau1_p, real *ssa1_p, real *tau2_p);


// increment nstream by 2stream
extern "C" void increment_nstream_by_2stream(int ncol, int nlay, int ngpt, int nmom1, real *tau1_p, real *ssa1_p, real *p1_p, real *tau2_p, real *ssa2_p, real *g2_p);


// increment nstream by nstream
extern "C" void increment_nstream_by_nstream(int ncol, int nlay, int ngpt, int nmom1, int nmom2, real *tau1_p, real *ssa1_p, real *p1_p, real *tau2_p, real *ssa2_p, real *p2_p);


// increment 1scalar by 2stream
extern "C" void inc_1scalar_by_2stream_bybnd(int ncol, int nlay, int ngpt, real *tau1_p, real *tau2_p, real *ssa2_p, int nbnd, int *gpt_lims_p);


// increment 1scalar by nstream
extern "C" void inc_1scalar_by_nstream_bybnd(int ncol, int nlay, int ngpt, real *tau1_p, real *tau2_p, real *ssa2_p, int nbnd, int *gpt_lims_p);


// increment 2stream by 1scalar
extern "C" void inc_2stream_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real *tau1_p, real *ssa1_p, real *tau2_p, int nbnd, int *gpt_lims_p);


// increment 2stream by nstream
extern "C" void inc_2stream_by_nstream_bybnd(int ncol, int nlay, int ngpt, int nmom2, real *tau1_p, real *ssa1_p, real *g1_p, real *tau2_p, real *ssa2_p, real *p2_p, int nbnd, int *gpt_lims_p);


// increment nstream by 1scalar
extern "C" void inc_nstream_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real *tau1_p, real *ssa1_p, real *tau2_p, int nbnd, int *gpt_lims_p);


// increment nstream by 2stream
extern "C" void inc_nstream_by_2stream_bybnd(int ncol, int nlay, int ngpt, int nmom1, real *tau1_p, real *ssa1_p, real *p1_p, real *tau2_p, real *ssa2_p, real *g2_p, int nbnd, int *gpt_lims_p);


// increment nstream by nstream
extern "C" void inc_nstream_by_nstream_bybnd(int ncol, int nlay, int ngpt, int nmom1, int nmom2, real *tau1_p, real *ssa1_p, real *p1_p, real *tau2_p, real *ssa2_p, real *p2_p, int nbnd, int *gpt_lims_p);


// Subsetting, meaning extracting some portion of the 3D domain
extern "C" void extract_subset_dim1_3d(int ncol, int nlay, int ngpt, real *array_in_p, int colS, int colE, real *array_out_p);


extern "C" void extract_subset_dim2_4d(int nmom, int ncol, int nlay, int ngpt, real *array_in_p, int colS, int colE, real *array_out_p);


// Extract the absorption optical thickness which requires mulitplying by 1 - ssa
extern "C" void extract_subset_absorption_tau(int ncol, int nlay, int ngpt, real *tau_in_p, real *ssa_in_p, int colS, int colE, real *tau_out_p);


