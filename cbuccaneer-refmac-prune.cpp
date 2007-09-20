// Clipper buccaneer
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>


int main( int argc, char** argv )
{
  CCP4Program prog( "cbuccaneer-refmac-prune", "0.1.0", "$Date: 2006/12/13" );

  std::cout << "\nCopyright 2002-2006 Kevin Cowtan and University of York. All rights reserved.\n\n";
  std::cout << "Please reference:\n Cowtan K. (2006) Acta Cryst. D62, 1002-1011.\n\n";

  // defaults
  clipper::String title;
  clipper::String ipmap = "NONE";
  clipper::String ippdb = "NONE";
  clipper::String oppdb = "NONE";
  double sig_d = 3.0;
  double sig_u = 3.0;
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-mapin" ) {
      if ( ++arg < args.size() ) ipmap = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-pdbout" ) {
      if ( ++arg < args.size() ) oppdb = args[arg];
    } else if ( args[arg] == "-sigma-density" ) {
      if ( ++arg < args.size() ) sig_d = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-sigma-tempfac" ) {
      if ( ++arg < args.size() ) sig_u = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cbuccaneer-refmac-prune\n\t-mapin <filename>\t\tCOMPULSORY\n\t-pdbin <filename>\t\tCOMPULSORY\n\t-pdbout <filename>\n\t-sigma-density <value>\n\t-sigma-tempfac <value>\nAn input pdb and map are required. Residues with high temperature factors or very negative difference density are removed, and the result written to the output pdb file.\n";
    exit(1);
  }

  // get map
  clipper::CCP4MAPfile file;
  clipper::Xmap<float> xmap;
  file.open_read( ipmap );
  file.import_xmap( xmap );
  file.close_read();

  // get model
  clipper::MiniMol mol;
  clipper::MMDBfile mmdb;
  mmdb.read_file( ippdb );
  mmdb.import_minimol( mol );

  // loop over residues and accumulate stats
  clipper::Cell cell = mol.cell();
  double d, d0, d1, d2, u, u0, u1, u2;
  d0 = d1 = d2 = u0 = u1 = u2 = 0.0;
  for ( int p = 0; p < mol.size(); p++ )
    for ( int m = 0; m < mol[p].size(); m++ ) {
      int cn = mol[p][m].lookup( " N  ", clipper::MM::ANY );
      int ca = mol[p][m].lookup( " CA ", clipper::MM::ANY );
      int cc = mol[p][m].lookup( " C  ", clipper::MM::ANY );
      if ( cn >= 0 && ca >= 0 && cc >= 0 ) {
	d = xmap.interp<clipper::Interp_cubic>( mol[p][m][cn].coord_orth().coord_frac(cell) ) + xmap.interp<clipper::Interp_cubic>( mol[p][m][ca].coord_orth().coord_frac(cell) ) + xmap.interp<clipper::Interp_cubic>( mol[p][m][cc].coord_orth().coord_frac(cell) );
	u = mol[p][m][cn].u_iso() + mol[p][m][ca].u_iso() + mol[p][m][cc].u_iso();
	d0 += 1.0;
	d1 += d;
	d2 += d*d;
	u0 += 1.0;
	u1 += u;
	u2 += u*u;
      }
    }
  d1 /= d0;
  d2 /= d0;
  d2 = sqrt( clipper::Util::max( d2 - d1*d1, 0.0 ) );
  u1 /= u0;
  u2 /= u0;
  u2 = sqrt( clipper::Util::max( u2 - u1*u1, 0.0 ) );
  double cut_d = d1 - sig_d * d2;
  double cut_u = u1 + sig_u * u2;

  // loop over residues and filter residues
  for ( int p = 0; p < mol.size(); p++ )
    for ( int m = 0; m < mol[p].size(); m++ ) {
      int cn = mol[p][m].lookup( " N  ", clipper::MM::ANY );
      int ca = mol[p][m].lookup( " CA ", clipper::MM::ANY );
      int cc = mol[p][m].lookup( " C  ", clipper::MM::ANY );
      if ( cn >= 0 && ca >= 0 && cc >= 0 ) {
	d = xmap.interp<clipper::Interp_cubic>( mol[p][m][cn].coord_orth().coord_frac(cell) ) + xmap.interp<clipper::Interp_cubic>( mol[p][m][ca].coord_orth().coord_frac(cell) ) + xmap.interp<clipper::Interp_cubic>( mol[p][m][cc].coord_orth().coord_frac(cell) );
	u = mol[p][m][cn].u_iso() + mol[p][m][ca].u_iso() + mol[p][m][cc].u_iso();
	if ( d < cut_d || u > cut_u )
	  if ( mol[p][m].type() == "UNK" )
	    mol[p][m].set_type( "~~~" );
      }
    }

  // make new molecule
  clipper::MiniMol mol2( mol.spacegroup(), mol.cell() );
  for ( int p = 0; p < mol.size(); p++ ) {
    clipper::MPolymer mp;
    mp.set_id( mol[p].id() );
    for ( int m = 0; m < mol[p].size(); m++ )
      if ( mol[p][m].type() != "~~~" )
	mp.insert( mol[p][m] );
    mol2.insert( mp );
  }

  // write answers
  clipper::MMDBfile mmdb_out;
  mmdb_out.export_minimol( mol2 );
  mmdb_out.write_file( oppdb );
}
