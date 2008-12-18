/*! \file buccaneer-prep.h buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for merging overlapped Ca chains and grouping by symmetry
class Ca_prep {
 public:
  struct Rama_flt { double phi, psi, rad; };

  Ca_prep( double main_tgt_rad, double side_tgt_rad, Rama_flt rama_flt, bool correl, bool seqnc, bool debug=false ) : rama_flt_(rama_flt), main_tgt_rad_(main_tgt_rad), side_tgt_rad_(side_tgt_rad), correl_(correl), seqnc_(seqnc), debug_(debug) {}
  bool operator() ( LLK_map_target& llktgt, std::vector<LLK_map_target>& llkcls, const clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const;

  // ramachandran filter data
  static const Rama_flt rama_flt_all, rama_flt_helix, rama_flt_strand, rama_flt_nonhelix;

  static void set_cpus( int cpus ) { ncpu = cpus; }
 private:
  Rama_flt rama_flt_;
  double main_tgt_rad_, side_tgt_rad_;
  bool correl_, seqnc_, debug_;
  static int ncpu;
};
