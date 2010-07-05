/*! \file buccaneer-known.h buccaneer library */
/* (C) 2002-2010 Kevin Cowtan & University of York all rights reserved */


#ifndef BUCCANEER_KNOWN
#define BUCCANEER_KNOWN

#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


//! Class for augmenting model with known model
class KnownStructure {
 public:
  //! Constructor: from model and arguments
  KnownStructure( const clipper::MiniMol& mol,
		  const std::vector<std::pair<clipper::String,double> >& ids );
  //! Return chain ids used by known structure elements
  clipper::String chain_ids() const;
  //! Add known structure to existing structure
  void copy_to( clipper::MiniMol& mol ) const;
  //! check for clashes against known model
  bool clash( const clipper::Coord_orth& coord ) const;
  //! prune model where it clashes with known model
  bool prune( clipper::MiniMol& mol ) const;
  //! Parse and store an input argument
  static std::pair<clipper::String,double> parse( clipper::String arg );
  void debug() const;
 private:
  typedef clipper::Property<double> PROP;
  clipper::MiniMol known;
  clipper::MAtomNonBond knownnb;
  double radius_max;
};


#endif