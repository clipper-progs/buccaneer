/*! \file buccaneer-join.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for merging overlapped Ca chains and grouping by symmetry
class Ca_join {
 public:
  Ca_join( clipper::ftype rad_merge = 2.0, clipper::ftype rad_join = 2.0 ) : rmerg(rad_merge), rjoin(rad_join) {}
  bool operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1 ) const;
 private:
  static std::vector<int> longest_chain( std::vector<std::vector<int> >& fwd_ptrs );  
  clipper::ftype rmerg, rjoin;
};
