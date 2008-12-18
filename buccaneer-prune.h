/*! \file buccaneer-prune.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for pruning clashing Ca chains using density
class Ca_prune {
 public:
  Ca_prune( clipper::ftype rad = 2.0 ) : rad_(rad) {}
  bool operator() ( clipper::MiniMol& mol ) const;
 private:
  clipper::ftype rad_;
};
