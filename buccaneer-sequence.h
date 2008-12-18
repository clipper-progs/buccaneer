/*! \file buccaneer-sequence.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for sequence Ca chains using density
class Ca_sequence {
 public:
  Ca_sequence( double reliability = 0.5 ) : reliability_(reliability) {}
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq );
  int num_sequenced() const;
  clipper::String format() const;

  static double phi_approx( double z );
  static std::vector<double> get_cached_scores( clipper::MMonomer& mm, const Ca_group& ca );
  static void set_cached_scores( clipper::MMonomer& mm, const Ca_group& ca, const std::vector<double>& scr );
  static double sequence_overlap( const clipper::String& seq1, const clipper::String& seq2 );
  static double sequence_similarity( const clipper::String& seq1, const clipper::String& seq2 );
  static Score_list<clipper::String> sequence_combine( const Score_list<clipper::String>& seq, const double& reliability );
  static std::pair<double,std::pair<int,int> > sequence_score( const std::vector<std::vector<double> >& scores, const clipper::String& subseq );
  static std::vector<clipper::String> sequence_align( const std::vector<std::vector<double> >& scores, const clipper::String& seq );
  static Score_list<clipper::String> sequence_match( const std::vector<std::vector<double> >& scores, const clipper::MMoleculeSequence& seq );
  static Score_list<clipper::String> sequence_chain( clipper::MChain& chain, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llktarget, const clipper::MMoleculeSequence& seq );
  static void sequence_apply( clipper::MChain& chain, const clipper::String& seq );
  static Score_list<clipper::String> sequence( clipper::MChain& chain, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llksample, const clipper::MMoleculeSequence& seq, const double& reliability );

  static void set_semet( bool semet ) { semet_ = semet; }
  static void set_prior_model( const clipper::MiniMol& mol );

  static void set_cpus( int cpus ) { ncpu = cpus; }

  class Sequence_data {
   public:
    Sequence_data() {}
    Sequence_data( const Ca_group& c, const std::vector<double>& d ) :
      ca(c), data(d) {}
    Ca_group ca;
    std::vector<double> data;
  };
 private:
  double reliability_;
  int num_seq;
  std::vector<Score_list<clipper::String> > history;
  static int ncpu;
  static bool semet_;
  static clipper::MiniMol molprior;
};
