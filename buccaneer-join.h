/*! \file buccaneer-join.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for merging overlapped Ca chains and grouping by symmetry
class Ca_join {
 public:
  Ca_join( double rad_merge = 2.0, double rad_join = 2.0 ) : rmerg(rad_merge), rjoin(rad_join) {}
  bool operator() ( clipper::MiniMol& mol ) const;

  static bool join( clipper::MiniMol& mol, const double& rmerg, const double& rjoin, const clipper::Coord_orth& com );
 private:
  static std::vector<int> longest_chain( std::vector<std::vector<int> >& fwd_ptrs );  
  double rmerg, rjoin;
};
