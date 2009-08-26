/*! \file buccaneer-mr.h buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for augmenting model with MR model
class Ca_mr {
 public:
  Ca_mr( clipper::ftype rad = 2.0 ) : rad_(rad) {}
  bool operator() ( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr ) const;
 private:
  clipper::ftype rad_;
};
