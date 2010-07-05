/*! \file buccaneer-tidy.h buccaneer library */
/* (C) 2002-2010 Kevin Cowtan & University of York all rights reserved */


#ifndef BUCCANEER_TIDY
#define BUCCANEER_TIDY

#include "buccaneer-prot.h"


//! Class for tidying model
class ModelTidy {
 public:
 ModelTidy( double rmsd = 1.0, double nmin = 12, bool verbose = false ) : rmsd_(rmsd), nmin_(nmin), verbose_(verbose) {}
  bool tidy( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const clipper::MMoleculeSequence& seq ) const;
 private:
  bool update_model( clipper::MiniMol& mol, std::vector<int> sources, std::vector<clipper::String> srcids ) const;
  static std::vector<int> assign_chains( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr );
  static std::vector<int> assign_chains( clipper::MiniMol& mol, const double rmsd, const int nmin );
  static bool chain_renumber( clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq );
  static bool label_map( clipper::Xmap<int>& xmap, const std::vector<clipper::Coord_orth>& coords, const double radius, const int label, const int clash, const int null );
  static std::vector<int> count_map( const clipper::Xmap<int>& xmap, const clipper::MPolymer& mp, const int nsrc );
  static bool move_chain( clipper::MPolymer& mp1, const clipper::MPolymer& mp0, const clipper::Spacegroup& spgr, const clipper::Cell& cell );
  static std::vector<clipper::Coord_orth> protein_atoms( const clipper::MPolymer& mp );
  static std::vector<int>      sequence_count( const clipper::MiniMol& mol );
  static clipper::Array2d<int> sequence_flags( const clipper::MiniMol& mol );
  static void best_closed_ncs_group( const clipper::Array2d<int>& super, const std::vector<int>& num_seq, std::vector<int>& used, std::vector<int>& used_best );
  double rmsd_;
  int nmin_;
  bool verbose_;
};


#endif
