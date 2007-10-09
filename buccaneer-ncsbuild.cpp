/*! \file buccaneer-ncsbuild.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-ncsbuild.h"


#include "buccaneer-join.h"
#include "buccaneer-sequence.h"
#include "buccaneer-filter.h"


clipper::RTop_orth Ca_ncsbuild::superpose( const clipper::MPolymer& mp1, const clipper::MPolymer& mp2 ) const
{
  clipper::RTop_orth result = clipper::RTop_orth::null();
  clipper::String seq1 = ProteinTools::chain_sequence( mp1 );
  clipper::String seq2 = ProteinTools::chain_sequence( mp2 );

  // see if sequences can be aligned (crude version)
  int bestscr = -1;
  int bestoff = 0;
  const int l1 = seq1.length();
  const int l2 = seq2.length();
  for ( int off = -l2+nmin_; off <= l1-nmin_; off++ ) {
    int scr = 0;
    for ( int i1 = 0; i1 < l1; i1++ ) {
      int i2 = i1 - off;
      if ( i2 >= 0 && i2 < l2 ) {
	if ( isupper( seq1[i1] ) && isupper( seq2[i2] ) ) {
	  if ( seq1[i1] == seq2[i2] ) scr += 1;
	  else                        scr -= 2;
	}
      }
    }
    if ( scr > bestscr ) {
      bestscr = scr;
      bestoff = off;
    }
  }
  // if no good alignment, return
  if ( bestscr < nmin_ ) return result;

  // now get the coordinates
  std::vector<clipper::Coord_orth> c1, c2;
  for ( int i1 = 0; i1 < l1; i1++ ) {
    int i2 = i1 - bestoff;
    if ( i2 >= 0 && i2 < l2 ) {
      if ( isupper( seq1[i1] ) && isupper( seq2[i2] ) ) {
	if ( seq1[i1] == seq2[i2] ) {
	  int a1 = mp1[i1].lookup( " CA ", clipper::MM::ANY );
	  int a2 = mp2[i2].lookup( " CA ", clipper::MM::ANY );
	  if ( a1 >= 0 && a2 >= 0 ) {
	    c1.push_back( mp1[i1][a1].coord_orth() );
	    c2.push_back( mp2[i2][a2].coord_orth() );
	  }
	}
      }
    }
  }

  // refine the alignment
  clipper::RTop_orth rtop_tmp;
  double r2;
  for ( int c = 0; c < 5; c++ ) {
    int nc = c1.size();
    // get transformation
    rtop_tmp = clipper::RTop_orth( c1, c2 );
    // get rmsd
    std::vector<std::pair<double,int> > r2index( nc );
    r2 = 0.0;
    for ( int i = 0; i < nc; i++ ) {
      double d2 = ( rtop_tmp * c1[i] - c2[i] ).lengthsq();
      r2 += d2;
      r2index[i] = std::pair<double,int>( d2, i );
    }
    r2 /= double( nc );
    // prune the list to improve it
    std::sort( r2index.begin(), r2index.end() );
    std::vector<clipper::Coord_orth> t1, t2;
    for ( int i = 0; i < (9*r2index.size())/10; i++ ) {
      t1.push_back( c1[r2index[i].second] );
      t2.push_back( c2[r2index[i].second] );
    }
    c1 = t1;
    c2 = t2;
  }

  /*
  clipper::String s1 = seq1;
  clipper::String s2 = seq2;
  while ( bestoff < 0 ) { s1 = " " + s1; bestoff++; }
  while ( bestoff > 0 ) { s2 = " " + s2; bestoff--; }
  std::cout << s1 << std::endl;
  std::cout << s2 << std::endl;
  */

  // if a close match has been found, return it
  if ( r2 < rmsd_*rmsd_ ) result = rtop_tmp;
  return result;
}


bool Ca_ncsbuild::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq ) const
{
  typedef clipper::MMonomer Mm;

  clipper::Cell       cell = xmap.cell();
  clipper::Spacegroup spgr = xmap.spacegroup();

  // split into separate chains
  ProteinTools::chain_tidy( mol2, mol1 );

  // now loop over chains
  for ( int chn1 = 0; chn1 < mol2.size(); chn1++ ) {
    // assemble combined model for this chain
    clipper::MPolymer mp1, mp2;
    mp1 = mol2[chn1];
    clipper::MiniMol mol_wrk, mol_tmp;
    mol_wrk.init( mol2.spacegroup(), mol2.cell() );
    mol_wrk.insert( mp1 );
    // add any other matching chains
    for ( int chn2 = 0; chn2 < mol2.size(); chn2++ ) {
      if ( chn2 != chn1 ) {
	clipper::RTop_orth rtop = superpose( mol2[chn2], mol2[chn1] );
	if ( !rtop.is_null() ) {
	  clipper::MPolymer mp = mol2[chn2];
	  mp.transform( rtop );
	  mol_wrk.insert( mp );
	}
      }
    }
    // were any matches found?
    if ( mol_wrk.size() > 1 ) {
      // remove sequence
      for ( int c = 0; c < mol_wrk.size(); c++ )
	for ( int r = 0; r < mol_wrk[c].size(); r++ )
	  mol_wrk[c][r].set_type( "UNK" );
      // join
      Ca_join cajoin( 2.0 );
      cajoin( mol_tmp, mol_wrk );
      mol_wrk = clipper::MiniMol( xmap.spacegroup(), xmap.cell() );
      mol_wrk = mol_tmp;
      // sequence
      Ca_sequence caseq( reliability_ );
      caseq( mol_tmp, mol_wrk, xmap, llktarget, seq );
      mol_wrk = mol_tmp;
      // filter
      Ca_filter cafiltr( 1.0 );
      cafiltr( mol_tmp, mol_wrk, xmap );
      mol_wrk = mol_tmp;
      // tidy
      ProteinTools::chain_tidy( mol_tmp, mol_wrk );
      // make new chain
      if ( mol_tmp.size() > 0 ) mp2 = mol_tmp[0];
      mp2.copy( mp1, clipper::MM::COPY_MP );
      // test if the chain is improved
      int l0, l1, s0, s1;
      l0 = mp1.size();
      l1 = mp2.size();
      s0 = s1 = 0;
      for ( int r = 0; r < mp1.size(); r++ )
	if ( ProteinTools::residue_index( mp1[r].type() ) >= 0 ) s0++;
      for ( int r = 0; r < mp2.size(); r++ )
	if ( ProteinTools::residue_index( mp2[r].type() ) >= 0 ) s1++;
      // if new chain is better, keep it
      if ( l1 > l0 && s1 > s0 )
	mol2[chn1] = mp2;
    }
  }

  return true;
}
