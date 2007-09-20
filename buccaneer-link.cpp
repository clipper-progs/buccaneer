/*! \file buccaneer-link.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-link.h"

#include <clipper/clipper-contrib.h>

#include <algorithm>


bool Ca_link::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget )
{
  mol2 = mol1;

  // establish map statistics to determine stopping value for building
  std::vector<double> scr;
  for ( int i = 0; i < 5000; i++ ) {
    double r = double(i);
    clipper::Euler_ccp4 rot( 2.0*r, 3.0*r, 5.0*r );
    clipper::Rotation ro( rot );
      clipper::Coord_grid cg( 23*i, 29*i, 31*i );
      clipper::Coord_orth trn =
	cg.coord_frac( xmap.grid_sampling() ).coord_orth( xmap.cell() );
      clipper::RTop_orth rtop( ro.matrix(), trn );
      double s = llktarget.llk( xmap, rtop );
      scr.push_back( s );
  }
  std::sort( scr.begin(), scr.end() );
  double cutoff = scr[25];  // select score for top 0.5% of cases

  // now do some rebuilding
  num_link = 0;
  bool cont;
  do {
    // search over possible chains and find short links:
    std::vector<std::pair<double, std::pair<int,int> > > links;
    for ( int chnn = 0; chnn < mol2.size(); chnn++ )
      for ( int chnc = 0; chnc < mol2.size(); chnc++ ) 
	if ( chnn != chnc && mol2[chnn].size() > 3 && mol2[chnc].size() > 3 ) {
	  int m = mol2[chnc].size()-1;
	  int in = mol2[chnn][0].lookup( " CA ", clipper::MM::ANY );
	  int ic = mol2[chnc][m].lookup( " CA ", clipper::MM::ANY );
	  double d2 = ( mol2[chnc][m][ic].coord_orth() -
			mol2[chnn][0][in].coord_orth() ).lengthsq();
	  if ( d2 < rlink*rlink ) links.push_back( std::pair<double,std::pair<int,int> >( d2, std::pair<int,int>( chnn, chnc ) ) );
	}

    // sort the links by length
    std::sort( links.begin(), links.end() );

    // now try each link in turn and see if it can be rebuilt
    cont = false;
    for ( int i = 0; i < links.size(); i++ ) {
      //std::cout << "Trying to link " << i << " of " << links.size() << "\t" << links[i].first << "\t" << links[i].second.first << "\t" << links[i].second.second << "\n";
      int chnn = links[i].second.first;
      int chnc = links[i].second.second;
      const clipper::MPolymer& mpn = mol2[chnn];
      const clipper::MPolymer& mpc = mol2[chnc];
      int resnbest = -1;
      int rescbest = -1;
      double scrbest = 1.0e12;
      ProteinLoop::CoordList<8> r8best;
      for ( int resc = mpc.size()-2; resc < mpc.size(); resc++ )
	for ( int resn = 0; resn < 2; resn++ ) 
	  if ( resc > 1 && resn < mol2[resn].size() - 1 ) {
	    int index_cc0 = mpc[resc-1].lookup( " C  ", clipper::MM::ANY );
	    int index_cn1 = mpc[resc  ].lookup( " N  ", clipper::MM::ANY );
	    int index_ca1 = mpc[resc  ].lookup( " CA ", clipper::MM::ANY );
	    int index_ca4 = mpn[resn  ].lookup( " CA ", clipper::MM::ANY );
	    int index_cc4 = mpn[resn  ].lookup( " C  ", clipper::MM::ANY );
	    int index_cn5 = mpn[resn+1].lookup( " N  ", clipper::MM::ANY );
	    if ( index_cc0 >= 0 && index_cn1 >= 0 && index_ca1 >= 0 &&
		 index_ca4 >= 0 && index_cc4 >= 0 && index_cn5 >= 0 ) {
	      // rebuild loop
	      ProteinLoop pl( torsion_sampling_ );
	      std::vector<ProteinLoop::CoordList<8> > r8;
	      r8 = pl.rebuild8atoms( mpc[resc-1][index_cc0].coord_orth(),
				     mpc[resc  ][index_cn1].coord_orth(),
				     mpc[resc  ][index_ca1].coord_orth(),
				     mpn[resn  ][index_ca4].coord_orth(),
				     mpn[resn  ][index_cc4].coord_orth(),
				     mpn[resn+1][index_cn5].coord_orth() );
	      for ( int i = 0; i < r8.size(); i++ ) {
		Ca_group ca1( r8[i][1], r8[i][2], r8[i][3] );
		Ca_group ca2( r8[i][4], r8[i][5], r8[i][6] );
		double scr = 0.5 * (
		  llktarget.llk( xmap, ca1.rtop_from_std_ori() ) +
		  llktarget.llk( xmap, ca2.rtop_from_std_ori() ) );
		if ( scr < scrbest ) {
		  scrbest = scr;
		  r8best = r8[i];
		  resnbest = resn;
		  rescbest = resc;
		}
	      }
	    }
	  }
      // now test whether the link is good enough
      if ( scrbest < cutoff ) {
	clipper::MPolymer mp;
	clipper::MAtom ca( clipper::MAtom::null() ), cn, cc;
	ca.set_occupancy(1.0); ca.set_u_iso(1.0); cn = cc = ca;
	cn.set_element( "N" ); cn.set_id( "N"  );
	ca.set_element( "C" ); ca.set_id( "CA" );
	cc.set_element( "C" ); cc.set_id( "C"  );
	// add first chain
	for ( int r = 0         ; r < rescbest  ; r++ ) mp.insert( mpc[r] );
	// and final residue
	clipper::MMonomer mm0 = mpc[rescbest];
	int i0 = mm0.lookup( " C  ", clipper::MM::ANY );
	if ( i0 >= 0 ) mm0[i0].set_coord_orth( r8best[0] );
	// first interpolated residue
	cn.set_coord_orth( r8best[1] );
	ca.set_coord_orth( r8best[2] );
	cc.set_coord_orth( r8best[3] );
	clipper::MMonomer mm1; mm1.set_type( "UNK" );
	mm1.insert( cn ); mm1.insert( ca ); mm1.insert( cc );
	// second interpolated residue
	cn.set_coord_orth( r8best[4] );
	ca.set_coord_orth( r8best[5] );
	cc.set_coord_orth( r8best[6] );
	clipper::MMonomer mm2; mm2.set_type( "UNK" );
	mm2.insert( cn ); mm2.insert( ca ); mm2.insert( cc );
	// first residue of next chain
	clipper::MMonomer mm3 = mpn[resnbest];
	int i3 = mm3.lookup( " N  ", clipper::MM::ANY );
	if ( i3 >= 0 ) mm3[i3].set_coord_orth( r8best[7] );
	// add the joining residues
	mp.insert( mm0 );
	mp.insert( mm1 );
	mp.insert( mm2 );
	mp.insert( mm3 );
	// and the remaining residues
	for ( int r = resnbest+1; r < mpn.size(); r++ ) mp.insert( mpn[r] );
	// store the first chain and delete the second
	int chnlo = ( chnn < chnc ) ? chnn : chnc ;
	int chnhi = ( chnn < chnc ) ? chnc : chnn ;
	clipper::MiniMol moltmp( mol2.spacegroup(), mol2.cell() );
	for ( int c = 0; c < mol2.size(); c++ )
	  if      ( c == chnlo ) moltmp.insert( mp );
	  else if ( c != chnhi ) moltmp.insert( mol2[c] );
	mol2 = moltmp;
	// and go back for another chain
	num_link++;
	cont = true;
	break;
      } // if we build a new join
    } // loop over possible joins
  } while ( cont );


  return true;
}


int Ca_link::num_linked() const
{
  return num_link;
}
