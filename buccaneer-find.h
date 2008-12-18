/*! \file buccaneer-find.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"
#include "simplex-lib.h"

#include <clipper/clipper-contrib.h>


// Result class
class SearchResult {
 public:
  float score; int rot; int trn;
  bool operator <( SearchResult other ) const { return score < other.score; }
};


//! Class for finding Ca's from density
class Ca_find {
 public:
  enum TYPE { LIKELIHOOD, SECSTRUC };
  Ca_find( int n_find = 500, double resol = 1.0 ) : nfind( n_find ), resol_( resol ) {}
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const TYPE type = LIKELIHOOD );
  static void set_cpus( int cpus ) { ncpu = cpus; }

 private:
  //friend class Search_threaded;
  // convenient analytial approximate distance funtion
  static double prob_dist( double x ) { return 0.999*exp(-75.0*pow(x-3.50,2.0)*pow(x,-2.5))+0.001; }
  // perform a single search from a list
  static void search_op( std::vector<SearchResult>& results, clipper::Xmap<float> xmap1, const clipper::Xmap<int>& xlookp1, const clipper::FFFear_fft<float>& srch, const LLK_map_target& llktarget, const std::vector<clipper::RTop_orth>& ops, int op );
  // FFFear map search function
  std::vector<SearchResult> search_llk( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget ) const;
  // SSfind map search function
  std::vector<SearchResult> search_sec( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget ) const;
  // prepare lateral growing prior
  void prep_prior( clipper::Xmap<float>& prior, const clipper::MiniMol& mol, const double radius=9.0 ) const;

  int nfind;
  double resol_;
  std::vector<clipper::RTop_orth> ops;
  std::vector<SearchResult> results;
  static int ncpu;
};


//! class for fast secondary structure finding (alternative to fffear)
class SSfind
{
 public:
  enum SSTYPE { ALPHA2, ALPHA3, ALPHA4, BETA2, BETA3, BETA4 };
  typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;

  void prep_xmap( const clipper::Xmap<float>& xmap, const double radius );
  void prep_target( SSfind::SSTYPE type, int num_residues );
  void set_target_score( const double score );
  const std::vector<SearchResult>& search( const std::vector<clipper::RTop_orth>& ops, const double rhocut );

  const std::vector<Pair_coord>          target_coords() { return target_cs; }
  const std::vector<clipper::Coord_orth> calpha_coords() { return calpha_cs; }
  const std::vector<SearchResult>& results() const { return rslts; }
 private:
  std::vector<Pair_coord>          target_cs;
  std::vector<clipper::Coord_orth> calpha_cs;
  std::vector<float> mapbox;
  clipper::Grid grid;
  clipper::Grid_range mxgr;
  clipper::Mat33<> grrot;
  std::vector<SearchResult> rslts;
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
