/*! \file buccaneer-sequence.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for sequence Ca chains using density
class Ca_sequence {
 public:
  Ca_sequence( double reliability = 0.5 ) : reliability_(reliability) {}
  bool operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const LLK_map_classifier& llktarget, const clipper::MMoleculeSequence& seq );
  int num_sequenced() const;
  clipper::String format() const;
 private:
  static Score_list<clipper::String> sequence_matches( const std::vector<std::vector<double> >& scores, const clipper::MMoleculeSequence& seq );
  double reliability_;
  int num_seq;
  std::vector<Score_list<clipper::String> > history;
};
