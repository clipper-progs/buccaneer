/*! \file buccaneer-util.h buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */


#ifndef BUCCANEER_UTIL
#define BUCCANEER_UTIL

#include <clipper/clipper-minimol.h>


class Buccaneer_util {
 public:
  static void set_reference( clipper::String& mtz, clipper::String& pdb );
};


class Buccaneer_log {
 public:
  Buccaneer_log() : currentcpu(0.0) { log(""); }
  void log( const clipper::String& id );
  void log( const clipper::String& id, const clipper::MiniMol& mol, bool view );
  void profile();
 private:
  std::vector<std::pair<std::string,double> > prof;
  double currentcpu;
};


#endif
