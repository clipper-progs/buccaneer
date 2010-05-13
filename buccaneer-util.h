/*! \file buccaneer-util.h buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */


#ifndef BUCCANEER_UTIL
#define BUCCANEER_UTIL

#include <clipper/clipper-minimol.h>


class BuccaneerUtil {
 public:
  static void set_reference( clipper::String& mtz, clipper::String& pdb );
  // coordinate utilities
  static void read_model( clipper::MiniMol& mol, clipper::String file, bool verbose );
};


class BuccaneerLog {
 public:
  BuccaneerLog() : currentcpu(0.0) { log(""); }
  void log( const clipper::String& id );
  void log( const clipper::String& id, const clipper::MiniMol& mol, bool view );
  void xml( const clipper::String& file, const clipper::MiniMol& mol );
  void profile();
 private:
  std::vector<std::pair<std::string,double> > prof;
  double currentcpu;
};


#endif
