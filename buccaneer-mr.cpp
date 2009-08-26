/*! \file buccaneer-mr.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-mr.h"

#include <clipper/clipper-contrib.h>



bool Ca_mr::operator() ( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr ) const
{
  typedef clipper::MMonomer Mm;
  clipper::ftype r2cut = rad_*rad_;

  clipper::Cell       cell = mol.cell();
  clipper::Spacegroup spgr = mol.spacegroup();

  // copy across chains from MR model
  clipper::MiniMol moltmp = mol;
  clipper::Property<bool> ptrue( true );
  for ( int c = 0; c < mol_mr.size(); c++ ) moltmp.insert( mol_mr[c] );
  for ( int c = mol.size(); c < moltmp.size(); c++ )
    for ( int r = 0; r < moltmp[c].size(); r++ )
      moltmp[c][r].set_property( "MR", ptrue );
  std::cout << "DEBUG " << mol.size() << " " << moltmp.size() << std::endl;

  // now loop over chains
  clipper::Coord_frac cf1, cf2;
  for ( int chn1 = 0; chn1 < moltmp.size()-1; chn1++ ) {
    for ( int chn2 = chn1; chn2 < moltmp.size(); chn2++ ) {
      // find any clashing residues between these chains
      for ( int res1 = 0; res1 < moltmp[chn1].size(); res1++ ) {
	for ( int res2 = 0; res2 < moltmp[chn2].size(); res2++ ) {
	  if ( ( chn1 != chn2 || res1 != res2 ) && 
	       moltmp[chn1][res1].type() != "~~~" &&
	       moltmp[chn2][res2].type() != "~~~" ) {
	    int a1 = moltmp[chn1][res1].lookup( " CA ", clipper::MM::ANY );
	    int a2 = moltmp[chn2][res2].lookup( " CA ", clipper::MM::ANY );
	    if ( a1 >= 0 && a2 >= 0 ) {
	      cf1 = moltmp[chn1][res1][a1].coord_orth().coord_frac(cell);
	      cf2 = moltmp[chn2][res2][a2].coord_orth().coord_frac(cell);
	      cf2 = cf2.symmetry_copy_near( spgr, cell, cf1 ) - cf1;
	      if ( cf2.lengthsq(cell) < r2cut ) {
		// clash found: if only one is sequenced, keep that,
		//              otherwise keep the longer.
		int scr1 = moltmp[chn1].size();
		int scr2 = moltmp[chn2].size();
		if ( moltmp[chn1][res1].exists_property( "MR") ||
		     moltmp[chn1][res1].type() == "UNK" ) scr1 += -1000000;
		if ( moltmp[chn2][res2].exists_property( "MR") ||
		     moltmp[chn2][res2].type() == "UNK" ) scr2 += -1000000;
		if ( scr1 > scr2 ) moltmp[chn2][res2].set_type( "~~~" );
		else               moltmp[chn1][res1].set_type( "~~~" );
	      }
	    }
	  }
	}
      }
    }
  }

  // eliminate any sequences of less than 6 residues
  mol = clipper::MiniMol( spgr, cell );
  clipper::MPolymer mp, mpnull;
  for ( int chn = 0; chn < moltmp.size(); chn++ ) {
    mp = mpnull;
    for ( int res = 0; res < moltmp[chn].size(); res++ ) {
      if ( moltmp[chn][res].type() != "~~~" ) {
	mp.insert( moltmp[chn][res] );
      } else {
	if ( mp.size() > 5 ) mol.insert( mp );
	mp = mpnull;
      }
    }
    if ( mp.size() > 5 ) mol.insert( mp );
  }

  // relabel MR chains and residues
  ProteinTools::chain_tidy( mol );
  int rc = 1;
  for ( int c = 0; c < mol.size(); c++ ) {
    std::cout << "DEBUG " << mol[c].id() << std::endl;
    if ( mol[c][0].exists_property( "MR") ) {
      mol[c].set_id( "!" );
      for ( int r = 0; r < mol[c].size(); r++ ) mol[c][r].set_seqnum( rc++ );
    }
    rc++;
  }

  return true;
}
