/*! \file buccaneer-util.cpp buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-util.h"

#include <fstream>
extern "C" {
#include <stdlib.h>
}


void Buccaneer_util::set_reference( clipper::String& mtz, clipper::String& pdb )
{
  const char* clibdptr = getenv( "CLIBD" );
  if ( clibdptr != NULL ) {
    clipper::String clibd( clibdptr );
    clipper::String path;
    std::ifstream file;
    if ( pdb == "NONE" ) {
      path = clibd+"/reference_structures/reference-1tqw.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( pdb == "NONE" ) {
      path = clibd+"\\reference_structures\\reference-1tqw.pdb";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) pdb = path;
    }
    if ( mtz == "NONE" ) {
      path = clibd+"/reference_structures/reference-1tqw.mtz";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) mtz = path;
    }
    if ( mtz == "NONE" ) {
      path = clibd+"\\reference_structures\\reference-1tqw.mtz";
      file.open( path.c_str(), std::ifstream::in ); file.close();
      if ( !file.fail() ) mtz = path;
    }
    if ( pdb == "NONE" || mtz == "NONE" ) 
      clipper::Message::message( clipper::Message_fatal( "No reference data specified and not in $CLIBD" ) );
  } else {
    clipper::Message::message( clipper::Message_fatal( "No reference data specified and $CLIBD not found" ) );
  }
}


#ifdef BUCCANEER_PROFILE
#include <sys/times.h>
void Buccaneer_log::log( const clipper::String& id )
{
  int i;
  tms tmst;
  times( &tmst );
  long ticks = sysconf(_SC_CLK_TCK);
  double cpu = double( tmst.tms_utime ) / double( ticks );
  double elapsed = cpu - currentcpu;
  if ( id != "" ) {
    for ( i = 0; i < prof.size(); i++ )
      if ( id == prof[i].first ) break;
    if ( i < prof.size() )
      prof[i].second += elapsed;
    else
      prof.push_back( std::pair<std::string,double>( id, elapsed ) );
  }
  currentcpu = cpu;
}
#else
void Buccaneer_log::log( const clipper::String& id ) {}
#endif


void Buccaneer_log::log( const clipper::String& id, const clipper::MiniMol& mol, bool view )
{
  if ( view ) {
    for ( int c = 0; c < mol.size(); c++ )
      for ( int r = 0; r < mol[c].size(); r++ ) {
	clipper::Coord_orth co(0.0,0.0,0.0);
	for ( int a = 0; a < mol[c][r].size(); a++ )
	  co += mol[c][r][a].coord_orth();
	co = (1.0/mol[c][r].size()) * co;
	std::cout << id << " " << c << "\t" << r << "\t" << co.format() << "\n";
	std::cout << id << " " << mol[c][r].type() << " ";
	for ( int a = 0; a < mol[c][r].size(); a++ )
	  std::cout << mol[c][r][a].id() << " ";
	std::cout << std::endl;
	int cn = mol[c][r].lookup( " N  ", clipper::MM::ANY );
	int ca = mol[c][r].lookup( " CA ", clipper::MM::ANY );
	int cc = mol[c][r].lookup( " C  ", clipper::MM::ANY );
	if ( ca >= 0 && cn >= 0 ) {
	  double d2 = ( mol[c][r][ca].coord_orth() -
			mol[c][r][cn].coord_orth() ).lengthsq();
	  if ( d2 > 6.25 )
	    std::cout << "BOND N-CA " << d2 << std::endl;
	}
	if ( ca >= 0 && cc >= 0 ) {
	  double d2 = ( mol[c][r][ca].coord_orth() -
			mol[c][r][cc].coord_orth() ).lengthsq();
	  if ( d2 > 6.25 )
	    std::cout << "BOND CA-C " << d2 << std::endl;
	}
      }
  }
  log( id );
}


void Buccaneer_log::profile()
{
  if ( prof.size() > 0 ) {
    std::cout << std::endl << "Profile:" << std::endl;
    for ( int i = 0; i < prof.size(); i++ )
      std::cout << prof[i].first << ": " << clipper::String( prof[i].second, 8 ) << " s" << std::endl;
  }
}
