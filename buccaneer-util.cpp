/*! \file buccaneer-util.cpp buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-util.h"
#include "buccaneer-prot.h"

#include <fstream>
extern "C" {
#include <stdlib.h>
}


void BuccaneerUtil::set_reference( clipper::String& mtz, clipper::String& pdb )
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


void BuccaneerUtil::read_model( clipper::MiniMol& mol, clipper::String file, bool verbose )
{
  const int mmdbflags = ( MMDBF_IgnoreBlankLines |
			  MMDBF_IgnoreDuplSeqNum |
			  MMDBF_IgnoreNonCoorPDBErrors |
			  MMDBF_IgnoreRemarks );
  clipper::MMDBfile mmdb;
  mmdb.SetFlag( mmdbflags );
  if ( file != "NONE" ) {
    std::vector<clipper::String> files = file.split(",");
    for ( int f = 0; f < files.size(); f++ ) {
      try {
	clipper::MiniMol moltmp;
	mmdb.read_file( files[f] );
	mmdb.import_minimol( moltmp );
	std::cout << "Read PDB file: " << files[f] << std::endl;
	for ( int c = 0; c < moltmp.size(); c++ )
	  if ( moltmp[c].id() != "!" ) mol.insert( moltmp[c] );	
	if ( verbose ) {
	  clipper::Atom_list atoms = moltmp.atom_list();
	  std::cout << "Number of atoms read: " << atoms.size() << std::endl;
	  for ( int i = 0; i < atoms.size(); i += atoms.size()-1 ) printf("%i6  %4s  %8.3f %8.3f %8.3f\n", i, atoms[i].element().c_str(), atoms[i].coord_orth().x(), atoms[i].coord_orth().y(), atoms[i].coord_orth().z() );
	}
      } catch ( clipper::Message_fatal ) {
	std::cout << "FAILED TO READ PDB FILE: " << file << std::endl;
      }
    }
  }
}


#ifdef BUCCANEER_PROFILE
#include <sys/times.h>
void BuccaneerLog::log( const clipper::String& id )
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
void BuccaneerLog::log( const clipper::String& id ) {}
#endif


void BuccaneerLog::log( const clipper::String& id, const clipper::MiniMol& mol, bool view )
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


void BuccaneerLog::xml( const clipper::String& file, const clipper::MiniMol& mol )
{
  int nres, nseq, nchn, nmax;
  nchn = mol.size();
  nres = nseq = nmax = 0;
  for ( int c = 0; c < mol.size(); c++ ) {
    if ( mol[c].size() > nmax ) nmax = mol[c].size();
    for ( int r = 0; r < mol[c].size(); r++ ) {
      if ( mol[c][r].lookup( " CA ", clipper::MM::ANY ) >= 0 ) nres++;
      if ( ProteinTools::residue_index_3( mol[c][r].type() ) >= 0 ) nseq++;
    }
  }

  std::ofstream f;
  f.open( file.c_str(), std::ios::out );
  f << "<BuccaneerResult>" << std::endl;
  f << "<ChainsBuilt>" << nchn << "</ChainsBuilt>" << std::endl;
  f << "<ResiduesBuilt>" << nres << "</ResiduesBuilt>" << std::endl;
  f << "<ResiduesSequenced>" << nseq << "</ResiduesSequenced>" << std::endl;
  f << "<ResiduesLongestChain>" << nmax << "</ResiduesLongestChain>" << std::endl;
  f << "</BuccaneerResult>" << std::endl;
  f.close();

}


void BuccaneerLog::profile()
{
  if ( prof.size() > 0 ) {
    std::cout << std::endl << "Profile:" << std::endl;
    for ( int i = 0; i < prof.size(); i++ )
      std::cout << prof[i].first << ": " << clipper::String( prof[i].second, 8 ) << " s" << std::endl;
  }
}
