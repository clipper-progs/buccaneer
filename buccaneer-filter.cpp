/*! \file buccaneer-filter.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-filter.h"

#include <clipper/clipper-contrib.h>

#include <algorithm>


std::vector<double> smooth( const std::vector<double>& x )
{
  int n = x.size();
  std::vector<double> result( n );
  result[ 0 ] = 0.25*( 3.0*x[0] + x[1] );
  for ( int i = 1; i < n-1; i++ )
    result[i] = 0.25*( x[i-1] + 2.0*x[i] + x[i+1] );
  result[n-1] = 0.25*( x[n-2] + 3.0*x[n-1] );
  return result;
}


bool Ca_filter::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap )
{
  clipper::Cell       cell = mol1.cell();
  clipper::Spacegroup spgr = mol1.spacegroup();

  // determine sigma cutoff based on map
  clipper::Map_stats stats( xmap );
  double s1 = stats.mean();
  double s2 = stats.std_dev();

  // now apply it a chain at a time
  std::vector<int> atms;
  clipper::MiniMol moltmp = mol1;
  for ( int chn = 0; chn < moltmp.size(); chn++ ) {
    // score the residues
    std::vector<double> scores( moltmp[chn].size(), 0.0 );
    for ( int res = 0; res < moltmp[chn].size(); res++ ) {
      atms.clear();
      atms.push_back( moltmp[chn][res].lookup( " N  ", clipper::MM::ANY ) );
      atms.push_back( moltmp[chn][res].lookup( " CA ", clipper::MM::ANY ) );
      atms.push_back( moltmp[chn][res].lookup( " C  ", clipper::MM::ANY ) );
      double t0, t1;
      t0 = t1 = 0.0;
      for ( int i = 0; i < atms.size(); i++ ) {
	int atm = atms[i];
	clipper::Coord_orth co = moltmp[chn][res][atm].coord_orth();
	clipper::Coord_frac cf = co.coord_frac( xmap.cell() );
	double r = xmap.interp<clipper::Interp_cubic>( cf );
	t0 += 1.0;
	t1 += r;
      }
      if ( t0 > 0.5 ) t1 /= t0;
      scores[res] = ( t1 - s1 ) / s2;
    }
    // smooth
    for ( int i = 0; i < 5; i++ ) scores = smooth( scores );

    // mark residues in sequence breaks
    int res = 0;
    while( res < moltmp[chn].size() ) {
      if ( moltmp[chn][res].type() != "UNK" ) break;
      res++;
    }
    while( res < moltmp[chn].size() ) {
      while( res < moltmp[chn].size() ) {
	if ( moltmp[chn][res].type() == "UNK" ) break;
	res++;
      }
      int res0 = res;
      while( res < moltmp[chn].size() ) {
	if ( moltmp[chn][res].type() != "UNK" ) break;
	res++;
      }
      int res1 = res;
      if ( res < moltmp[chn].size() && res1 - res0 >= 3 ) {
	int resm = res0+1;
	for ( int r = res0+1; r < res1-1; r++ )
	  if ( scores[r] < scores[resm] ) resm = r;
        moltmp[chn][resm-1].set_type( "~~~" );
	moltmp[chn][resm  ].set_type( "~~~" );
	moltmp[chn][resm+1].set_type( "~~~" );
      }
    }

    // mark residues in poor density
    for ( int res = 0; res < moltmp[chn].size(); res++ )
      if ( moltmp[chn][res].type() == "UNK" && scores[res] <= sigcut )
	moltmp[chn][res].set_type( "~~~" );

  }

  // eliminate any sequences of less than 6 residues
  mol2 = clipper::MiniMol( spgr, cell );
  clipper::MPolymer mp, mpnull;
  for ( int chn = 0; chn < moltmp.size(); chn++ ) {
    mp = mpnull;
    for ( int res = 0; res < moltmp[chn].size(); res++ ) {
      if ( moltmp[chn][res].type() != "~~~" ) {
	mp.insert( moltmp[chn][res] );
      } else {
	if ( mp.size() > 5 ) mol2.insert( mp );
	mp = mpnull;
      }
    }
    if ( mp.size() > 5 ) mol2.insert( mp );
  }

  return true;
}


int Ca_filter::num_filtered() const
{
  return num_filter;
}
