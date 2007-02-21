/*! \file buccaneer-find.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"
#include "simplex-lib.h"


//! Class for finding Ca's from density
class Ca_find {
 public:
  Ca_find( int n_find = 500 ) : nfind( n_find ) {}
  bool operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget );
 private:
  int nfind;
  clipper::Xmap<float> resultscr;
  clipper::Xmap<int>   resultrot;
  clipper::Xmap<int>   resulttrn;
};


//! class for refining Ca groups
class Target_fn_refine_llk_map_target : Target_fn_order_zero
{
 public:
  Target_fn_refine_llk_map_target() {}
  Target_fn_refine_llk_map_target( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const double& rot_step, const double& trn_step );
  ~Target_fn_refine_llk_map_target() {}
  int num_params() const { return 6; }
  //! evaluate target function for given rotation
  double operator() ( const clipper::RTop_orth& rtop ) const;
  //! \internal evaluate target function for EulerXYZr offset from rot_
  double operator() ( const std::vector<double>& args ) const;
  //! \internal convert params to rotation
  clipper::RTop_orth rtop_orth( const std::vector<double>& args ) const;  
  //! refine rotation
  clipper::RTop_orth refine( const clipper::RTop_orth& rtop );
 private:
  const clipper::Xmap<float>* xmap_;
  const LLK_map_target* llktarget_;
  double rot_step_, trn_step_;
  clipper::RTop_orth rtop_;
};
