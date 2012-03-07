/*! \file buccaneer-prune.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for pruning clashing Ca chains using density
class Ca_prune {
 public:
  Ca_prune( double rad = 3.0 ) : rad_(rad) {}
  static bool prune( clipper::MiniMol& mol, double rad = 3.0 );
  bool operator() ( clipper::MiniMol& mol ) const;
 private:
  static std::vector<int> score_positions( const clipper::MPolymer& mp );
  clipper::ftype rad_;
};
