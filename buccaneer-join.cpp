/*! \file buccaneer-join.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-join.h"

#include <clipper/clipper-contrib.h>



struct Tri_residue {
  int flag;
  clipper::ftype score;
  clipper::String type[3];
  clipper::Coord_orth res[3][3];
};


// find longest single chain ignoring loops
std::vector<int> Ca_join::longest_chain( std::vector<std::vector<int> >& fwd_ptrs ) {
  // first make a list of loop starts
  std::vector<int> node_count( fwd_ptrs.size(), 0 );
  // loop starts are now marks as 0
  // declare a list of 'dirty' nodes
  std::set<int> dirty;
  for ( int i = 0; i < node_count.size(); i++ )
    if ( node_count[i] == 0 ) dirty.insert( i );
  // now propogate the values
  std::vector<int> bck_ptrs( fwd_ptrs.size(), -1 );
  while ( dirty.size() > 0 ) {
    // get a node fromt he dirty list and remove it from the list
    std::set<int>::iterator iter = dirty.begin();
    int node = *iter;
    dirty.erase( iter );
    // now check the children
    for ( int j = 0; j < fwd_ptrs[node].size(); j++ ) {
      int next_node = fwd_ptrs[node][j];
      // test whether we've found a longer route
      if ( node_count[next_node] < node_count[node]+1 ) {
	// check here for loops and broken paths
	int back_node = node;
	int back_node_next;
	while ( bck_ptrs[back_node] >= 0 ) {
	  back_node_next = bck_ptrs[back_node];
	  if ( back_node_next == next_node ) break;
	  if ( node_count[back_node_next] >= node_count[back_node] ) break;
	  back_node = back_node_next;
	}
	// if the path to this node is clean, we can update the node
	if ( bck_ptrs[back_node] < 0 ) {
	  // if this is a longer non-looped route, store it
	  node_count[next_node] = node_count[node]+1;
	  bck_ptrs[next_node] = node;
	  dirty.insert( next_node );
	}
      }
    }
  }
  // we've found all the long routes, now find the longest and back-trace it
  int node_max = 0;
  for ( int i = 1; i < node_count.size(); i++ )
    if ( node_count[i] > node_count[node_max] )
      node_max = i;
  // and back-trace
  std::vector<int> result;
  int node = node_max;
  result.push_back( node );
  while ( bck_ptrs[node] >= 0 ) {
    int next = bck_ptrs[node];
    bck_ptrs[node] = -1;
    node = next;
    result.push_back( node );
  }
  // reverse the list
  std::reverse( result.begin(), result.end() );
  return result;
}


// build chains by merging and joining tri-residue fragments
bool Ca_join::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1 ) const
{
  typedef clipper::MMonomer Mm;
  clipper::ftype r2merg = rmerg*rmerg;
  clipper::ftype r2join = rjoin*rjoin;

  clipper::Cell       cell = mol1.cell();
  clipper::Spacegroup spgr = mol1.spacegroup();

  // first calculate a convenient ASU centre for output of results
  // (cosmetic only)
  clipper::Coord_frac com;
  {
    // calc mask
    clipper::Resolution reso( 2.0 );
    clipper::Grid_sampling grid( spgr, cell, reso );
    clipper::Xmap<float> xmap( spgr, cell, grid ), xflt( spgr, cell, grid );
    clipper::EDcalc_mask<float> maskcalc( 2.0 );
    maskcalc( xmap, mol1.atom_list() );
    // calc smoothing radius
    clipper::ftype rad = 0.5 * pow( cell.volume()/spgr.num_symops(), 0.333 );
    clipper::MapFilterFn_linear fn( rad );
    clipper::MapFilter_fft<float>
      fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
    fltr( xflt, xmap ); 
    // find peak
    typedef clipper::Xmap<float>::Map_reference_index MRI;
    MRI iy = xflt.first();
    for ( MRI ix = xflt.first(); !ix.last(); ix.next() )
      if ( xflt[ix] > xflt[iy] ) iy = ix;
    com = iy.coord().coord_frac( grid );
  }

  // create 3-residue segments
  std::vector<Tri_residue> fragments;
  Tri_residue fragment;
  for ( int chn = 0; chn < mol1.size(); chn++ ) {
    for ( int res = 1; res < mol1[chn].size()-1; res++ ) {
      if ( Mm::protein_peptide_bond( mol1[chn][res-1], mol1[chn][res] ) &&
	   Mm::protein_peptide_bond( mol1[chn][res], mol1[chn][res+1] ) ) {
	int n0 = mol1[chn][res-1].lookup( " N  ", clipper::MM::ANY );
	int a0 = mol1[chn][res-1].lookup( " CA ", clipper::MM::ANY );
	int c0 = mol1[chn][res-1].lookup( " C  ", clipper::MM::ANY );
	int n1 = mol1[chn][res  ].lookup( " N  ", clipper::MM::ANY );
	int a1 = mol1[chn][res  ].lookup( " CA ", clipper::MM::ANY );
	int c1 = mol1[chn][res  ].lookup( " C  ", clipper::MM::ANY );
	int n2 = mol1[chn][res+1].lookup( " N  ", clipper::MM::ANY );
	int a2 = mol1[chn][res+1].lookup( " CA ", clipper::MM::ANY );
	int c2 = mol1[chn][res+1].lookup( " C  ", clipper::MM::ANY );
	if ( n0 >= 0 && a0 >= 0 && c0 >= 0 &&
	     n1 >= 0 && a1 >= 0 && c1 >= 0 &&
	     n2 >= 0 && a2 >= 0 && c2 >= 0 ) {
	  fragment.type[0] = mol1[chn][res-1].type();
	  fragment.type[1] = mol1[chn][res  ].type();
	  fragment.type[2] = mol1[chn][res+1].type();
	  fragment.res[0][0] = mol1[chn][res-1][n0].coord_orth();
	  fragment.res[0][1] = mol1[chn][res-1][a0].coord_orth();
	  fragment.res[0][2] = mol1[chn][res-1][c0].coord_orth();
	  fragment.res[1][0] = mol1[chn][res  ][n1].coord_orth();
	  fragment.res[1][1] = mol1[chn][res  ][a1].coord_orth();
	  fragment.res[1][2] = mol1[chn][res  ][c1].coord_orth();
	  fragment.res[2][0] = mol1[chn][res+1][n2].coord_orth();
	  fragment.res[2][1] = mol1[chn][res+1][a2].coord_orth();
	  fragment.res[2][2] = mol1[chn][res+1][c2].coord_orth();
	  fragment.flag = 1;
	  fragment.score = 1.0;
	  // upweight sequenced fragments, and flag core seqeunced regions.
	  if ( mol2[chn][res].type() != "UNK" ) {
	    int r1, r2;
	    for ( r1 = res-1; r1 >= 0; r1-- )
	      if ( mol2[chn][r1].type() == "UNK" ) break;
	    for ( r2 = res+1; r2 < mol2[chn].size(); r2++ )
	      if ( mol2[chn][r2].type() == "UNK" ) break;
	    int d = clipper::Util::min( res-r1, r2-res );
	    fragment.score += 0.1 * double(d);  // upweight sequenced
	    if ( d > 15 ) fragment.flag = 2;    // flag core sequenced
	  }
	  fragments.push_back( fragment );
	} // if mainchain atoms present
      } // if connected residues
    } // loop over residues
  } // loop over chains

  // now merge equivalent fragments
  for ( int f1 = 0; f1 < fragments.size()-1; f1++ ) {
    clipper::Coord_frac cx0 = fragments[f1].res[0][1].coord_frac(cell);
    clipper::Coord_frac cx1 = fragments[f1].res[1][1].coord_frac(cell);
    clipper::Coord_frac cx2 = fragments[f1].res[2][1].coord_frac(cell);
    for ( int f2 = f1+1; f2 < fragments.size(); f2++ ) {
      if ( fragments[f1].flag == 1 && fragments[f2].flag == 1 ) {
	clipper::Coord_frac cy0 = fragments[f2].res[0][1].coord_frac(cell);
	clipper::Coord_frac cy1 = fragments[f2].res[1][1].coord_frac(cell);
	clipper::Coord_frac cy2 = fragments[f2].res[2][1].coord_frac(cell);
	cy0 = cy0.symmetry_copy_near( spgr, cell, cx1 );
	cy1 = cy1.symmetry_copy_near( spgr, cell, cx1 );
	cy2 = cy2.symmetry_copy_near( spgr, cell, cx1 );
	if ( ( cy0 - cx0 ).lengthsq(cell) < r2merg &&
	     ( cy1 - cx1 ).lengthsq(cell) < r2merg &&
	     ( cy2 - cx2 ).lengthsq(cell) < r2merg ) {
	  clipper::ftype s1 = fragments[f1].score / (fragments[f1].score+1.0);
	  clipper::ftype s2 =                 1.0 / (fragments[f1].score+1.0);
	  for ( int r = 0; r < 3; r++ )
	    for ( int a = 0; a < 3; a++ ) {
	      cy1 = fragments[f2].res[r][a].coord_frac(cell);
	      cy1 = cy1.symmetry_copy_near( spgr, cell, cx1 );
	      fragments[f1].res[r][a] =
		s1 * fragments[f1].res[r][a] + s2 * ( cy1.coord_orth(cell) );
	    }
	  fragments[f1].score += fragments[f2].score;
	  fragments[f2].score = 0.0;
	  fragments[f2].flag = 0;
	}
      }
    }
  }

  // make a list of joins
  std::vector<std::vector<int> > joins( fragments.size() );
  for ( int f1 = 0; f1 < fragments.size(); f1++ )
    if ( fragments[f1].flag != 0 ) {
      clipper::Coord_frac cx0 = fragments[f1].res[0][1].coord_frac(cell);
      clipper::Coord_frac cx1 = fragments[f1].res[1][1].coord_frac(cell);
      clipper::Coord_frac cx2 = fragments[f1].res[2][1].coord_frac(cell);
      for ( int f2 = 0; f2 < fragments.size(); f2++ )
	if ( fragments[f2].flag != 0 ) {
	  if ( f1 != f2 ) {
	    clipper::Coord_frac cy0 = fragments[f2].res[0][1].coord_frac(cell);
	    clipper::Coord_frac cy1 = fragments[f2].res[1][1].coord_frac(cell);
	    clipper::Coord_frac cy2 = fragments[f2].res[2][1].coord_frac(cell);
	    cy0 = cy0.symmetry_copy_near( spgr, cell, cx1 );
	    cy1 = cy1.symmetry_copy_near( spgr, cell, cx1 );
	    cy2 = cy2.symmetry_copy_near( spgr, cell, cx1 );
	    if ( (cx1-cy0).lengthsq(cell) < r2join &&
		 (cx2-cy1).lengthsq(cell) < r2join ) {
	      if ( fragments[f1].flag == 1 && fragments[f2].flag == 1 )
		joins[f1].push_back( f2 );
	      else
		if ( f2 == f1+1 )
		  joins[f1].push_back( f2 );
	    }
	  }
	}
    }

  /*
  for ( int j1 = 0; j1 < joins.size(); j1++ ) {
    std::cout << j1 << "(" << fragments[j1].flag << "):\t";
    for ( int j2 = 0; j2 < joins[j1].size(); j2++ ) {
      std::cout << joins[j1][j2] << "(" << fragments[joins[j1][j2]].flag << ")\t";
    }
    std::cout << std::endl;
  }
  */

  // use threading to extract successive longest chains
  std::vector<std::vector<int> > chns;
  { // code block to avoid windows compiler bugs over use of label "chn"
    std::vector<int> chn = longest_chain( joins );
    while ( chn.size() > 5 ) {
      // add longest chain to list
      chns.push_back( chn );
      // remove used fragments
      for ( int r = 0; r < chn.size(); r++ )
	fragments[chn[r]].flag = 0;
      // remove links from used fragments
      for ( int f = 0; f < joins.size(); f++ )
	if ( fragments[f].flag == 0 )
	  joins[f].clear();
      // and links to used fragments
      for ( int f = 0; f < joins.size(); f++ )
	for ( int j = joins[f].size()-1; j >= 0; j-- )
	  if ( fragments[joins[f][j]].flag == 0 )
	    joins[f].erase( joins[f].begin() + j );
      // get longest remaining chain
      chn = longest_chain( joins );
    }
  }

  /*
  for ( int c = 0; c < chns.size(); c++ ) {
    std::cout << c << ":\t";
    for ( int r = 0; r < chns[c].size(); r++ ) {
      std::cout << chns[c][r] << " ";
    }
    std::cout << std::endl;
  }
  */

  // now join the fragments
  char atomid[3][5] = { " N  ", " CA ", " C  " };
  char atomel[3][2] = { "N", "C", "C" };

  mol2 = clipper::MiniMol( spgr, cell );
  for ( int c = 0; c < chns.size(); c++ ) {
    // chain and atom info
    clipper::MPolymer chain;
    int ires = 1;
    clipper::MAtom atom = clipper::Atom::null();
    atom.set_occupancy(1.0);
    atom.set_u_iso( 1.0 );

    const std::vector<int>& chn = chns[c];

    clipper::Coord_frac cx, cy;
    cx = com;  // reference coord

    for ( int f = -1; f < int(chn.size())+1; f++ ) {
      // residue info
      clipper::MMonomer residue;
      residue.set_type("UNK");

      // add this residue
      for ( int a = 0; a < 3; a++ ) {
	clipper::Coord_orth co( 0.0, 0.0, 0.0 );
	clipper::ftype s = 0.0;
	for ( int r = -1; r <= 1; r++ )
	  if ( f+r >= 0 && f+r < chn.size() ) {
	    s += 1.0;
	    cy = fragments[chn[f+r]].res[1-r][a].coord_frac(cell);
	    cy = cy.symmetry_copy_near( spgr, cell, cx );
	    co += cy.coord_orth(cell);
	  }
	co = (1.0/s) * co;
	atom.set_element( atomel[a] );
	atom.set_id( atomid[a] );
	atom.set_coord_orth( co );
	residue.insert( atom );
      }
      residue.set_seqnum( ires++ );
      chain.insert( residue );

      cx = residue[1].coord_orth().coord_frac(cell);  // update reference coord
    }
    mol2.insert( chain );
  }

  // tidy up the peptide bonds
  const double dmax = 1.45;
  for ( int chn = 0; chn < mol2.size(); chn++ ) {
    for ( int res = 0; res < mol2[chn].size()-1; res++ ) {
      int a1 = mol2[chn][res  ].lookup( " CA ", clipper::MM::ANY );
      int c1 = mol2[chn][res  ].lookup( " C  ", clipper::MM::ANY );
      int n2 = mol2[chn][res+1].lookup( " N  ", clipper::MM::ANY );
      int a2 = mol2[chn][res+1].lookup( " CA ", clipper::MM::ANY );
      if ( a1 >= 0 && c1 >= 0 && n2 >= 0 && a2 >= 0 ) {
	// rebuild peptide units
	clipper::Coord_orth ca1 = mol2[chn][res  ][a1].coord_orth();
	clipper::Coord_orth cc1 = mol2[chn][res  ][c1].coord_orth();
	clipper::Coord_orth cn2 = mol2[chn][res+1][n2].coord_orth();
	clipper::Coord_orth ca2 = mol2[chn][res+1][a2].coord_orth();
	// check and rebuild peptide units if necessary
	if ( (cc1-cn2).lengthsq() > dmax*dmax ) {
	  clipper::Vec3<> v = clipper::Vec3<>::cross( cn2-cc1, ca2-ca1 );
	  v = clipper::Vec3<>::cross( v, ca2-ca1 ).unit();
	  cc1 = clipper::Coord_orth( 0.63*ca1 + 0.37*ca2 + 0.57*v );
	  cn2 = clipper::Coord_orth( 0.37*ca1 + 0.63*ca2 - 0.43*v );
	}
	// check and restore C-N connectivity if necessary
	if ( (cc1-cn2).lengthsq() > dmax*dmax ) {
	  double d = sqrt( ( cc1 - cn2 ).lengthsq() );
	  double f = 0.5 * ( 1.0 - dmax / d );
	  cc1 = (1.0-f)*cc1 + f*cn2;
	  cn2 = (1.0-f)*cn2 + f*cc1;
	}
	// store
	mol2[chn][res  ][c1].set_coord_orth( cc1 );
	mol2[chn][res+1][n2].set_coord_orth( cn2 );
      }
    }
  }

  // globularise
  ProteinTools::globularise( mol2 );

  // restore the residue types, if any
  ProteinTools::copy_residue_types( mol2, mol1 );
  return true;
}
