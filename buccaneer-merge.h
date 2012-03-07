/*! \file buccaneer-merge.h buccaneer library */
/* (C) 2009 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for augmenting model with MERGE model
class Ca_merge {
 public:
  Ca_merge( double reliability = 0.5 ) : reliability_(reliability) {}
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq ) const;
 private:
  double reliability_;
};
