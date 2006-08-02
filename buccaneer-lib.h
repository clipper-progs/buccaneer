/*! \file buccaneer-lib.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#ifndef BUCCANEER_LIB
#define BUCCANEER_LIB

#include <clipper/clipper.h>


//! Score list class.
/*! The class keeps a list of the best fits to some target function,
  scoring both the score and an object of a user defined type
  representing the parameters whuich gave that score. On construction,
  a maximum size is specified. Only the best N matches will be
  stored.

  The implementation could be more efficient. */
template<class T> class Score_list {
 public:
  Score_list() {}  //!< null constructor
  //! constructor: takes maximum size of list
  Score_list( const int& n ) { init( n ); }
  //! initialiser: takes maximum size of list
  void init( const int& n ) { max = n; list.clear(); list.reserve( max ); }
  //! is a given score good enough to be added to the list?
  bool addable( const clipper::ftype& score ) const
    { return ( list.size() < max || score < list.back().first ); }
  //! add a score to the list, if it is good enough
  void add( const clipper::ftype& scr, const T& data ) {
    if ( addable( scr ) ) {
      if ( list.size() >= max )	list.pop_back();
      int i; for ( i = list.size()-1; i >= 0; i-- ) if ( scr > score(i) ) break;
      list.insert( list.begin()+i+1, std::pair<clipper::ftype,T>(scr,data) );
    }
  }
  //! delete a score from the list
  void del( const int& i ) { list.erase( list.begin()+i ); }
  //! access score
  const clipper::ftype& score( const int& i ) const { return( list[i].first ); }
  //! access list
  const T& operator[] ( const int& i ) const { return( list[i].second ); }
  //! list size
  int size() const { return list.size(); }
 private:
  int max;
  std::vector<std::pair<clipper::ftype,T> > list;
};


//! Log-likelihood map matching target
/*! This class is used in determining the log-likelihood fit of some
  desired density features from some region of a noisy electron
  density map. It contains methods to accumulate the log likelihood
  target from a number of sample density regions from a map with
  similar noise levels to the target map; methods for a FFFear
  6-dimensional (rotation/orientation) search of the target map; and
  methods for testing indivdual sample positions and orientations.

  Note the results from the 6-d search and the fast and full LLK
  calculations are on different scales and so cannot be compared
  directly. */
class LLK_map_target {
 public:
  //! constructor: provide radius and sampling in A for LLK target
  LLK_map_target( const clipper::ftype& rad, const clipper::ftype& sampling );
  //! prepare LLK target after accumulating density or loading map
  void prep_llk();
  //! accumulate density statistics from a sample orientation in a sample map
  void accumulate( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop );
  //! access to LLK target for load/save
  clipper::NXmap<float>& llk_target() { return target; }
  //! access to LLK weight for load/save
  clipper::NXmap<float>& llk_weight() { return weight; }
  //! perform 6-d (RT) search for LLK target in given map
  Score_list<clipper::RTop_orth> search( const clipper::Xmap<float>& xmap, const clipper::ftype& step, const int& nres, const clipper::ftype& dres ) const;
  //! calculate fast approx to LLK for given orientation
  clipper::ftype llk_approx( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop ) const;
  //! calculate full LLK for given orientation
  clipper::ftype llk       ( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop ) const;
  //! output formatted representation
  clipper::String format() const;
 private:
  clipper::ftype radius;  //!< density sphere radius
  int naccum;             //!< number of maps accumulated
  clipper::NXmap<float> target;  //!< target map
  clipper::NXmap<float> weight;  //!< weight map
  std::vector<clipper::Coord_orth> repxyz, fullxyz;  //!< fast target lists
  std::vector<clipper::ftype>      reptgt, fulltgt;
  std::vector<clipper::ftype>      repwgt, fullwgt;
};


//! Log-likelihood map classifying target
/*! This class is used in determining classifying some region of a
  noisy electron density map according to the best log-likelihood fit
  to one of a set of targets, e.g. for identification of side chains.
  It contains methods for constructing a list of targets and for
  identifying which target a particular region of density resembles.

  Targets are identified by an ordinal number, so the nature and the
  order of the targets is irrelevent.

  This llk_raw returns raw scores, which may be on different scales
  for each of the targets and so need post-processing, e.g. using a
  z-score.

  \code
  LLK_map_classifier llkmcl( 4.0, 0.5, 20 );
  // accumulate scores for map features in reference map
  for ( int i = 0; i < ref_rtops.size(); i++ )
    llkmcl.accumulate( ref_map, ref_rtops[i], res_types[i] );
  llkmcl.prep_llk()
  // get the scores for a map feature in the target map 
  std::vector<double> scores = llkmcl.llk_raw( tgt_map, tgt_rtops[i] );
  \endcode
*/
class LLK_map_classifier {
 public:
  //! constructor: provide radius and sampling and number of targets
  LLK_map_classifier( const clipper::ftype& rad, const clipper::ftype& sampling, const int& num_targets );
  //! accumulate density statistics from a sample orientation in a sample map
  void accumulate( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop, const int& target_type );
  //! prepare LLK target after accumulating density or loading map
  void prep_llk();
  //! calculate raw LLKs for given map location for each target
  std::vector<clipper::ftype> llk_raw( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop ) const;
  //! direct access to target function for an individual target
  LLK_map_target& llk_map_target( const int& target_type )
    { return llktgts[target_type]; }
  //! number of targets into which features are classified
  int num_targets() const { return llktgts.size(); }
 private:
  std::vector<LLK_map_target> llktgts;  //!< the targets
};


#endif
