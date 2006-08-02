/*! \file buccaneer-build.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for building Ca chains using density
class Ca_build {
 public:
  Ca_build( clipper::String type = "ALA" ) : newrestype( type ) {}
  bool operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap ) const;
 private:
  struct Clash { int p1, m1, p2, m2; };
  std::vector<Clash> find_clashes( clipper::MiniMol& mol ) const;
  void fix_clash  ( clipper::MMonomer& m1, clipper::MMonomer& m2,
		    const clipper::Xmap<float>& xmap ) const;
  void fix_clashes( clipper::MiniMol& mol,
		    const clipper::Xmap<float>& xmap ) const;
  clipper::String newrestype;
};
