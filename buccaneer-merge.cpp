/*! \file buccaneer-merge.cpp buccaneer library */
/* (C) 2009 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-merge.h"


#include "buccaneer-join.h"
#include "buccaneer-sequence.h"
#include "buccaneer-prune.h"
#include "buccaneer-build.h"


bool Ca_merge::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq ) const
{
  // extract the necessary bits of the likelihood targets
  std::vector<LLK_map_target::Sampled> llksample( llktarget.size() );
  for ( int t = 0; t < llktarget.size(); t++ )
    llksample[t] = llktarget[t].sampled();

  // constants
  const clipper::Spacegroup& spgr = xmap.spacegroup();
  const clipper::Cell&       cell = xmap.cell();
  const clipper::Coord_orth com( 0.0, 0.0, 0.0 );
  const double rad = 2.0;

  // update sequencing
  for ( int c1 = 0; c1 < mol.size(); c1++ ) {
    Ca_sequence::prepare_scores( mol[c1], xmap, llksample );
    Ca_sequence::sequence( mol[c1], seq, reliability_ );
  }
  ProteinTools::break_chains( mol, xmap );
  ProteinTools::chain_tidy( mol );

  // now try to augment each chain in turn to increasing the sequencing
  clipper::Coord_frac f1, f2;
  for ( int c1 = 0; c1 < mol.size()-1; c1++ ) {
    std::vector<clipper::Coord_frac> atoms;
    for ( int r1 = 0; r1 < mol[c1].size(); r1++ ) {
	int a1 = mol[c1][r1].lookup( "CA", clipper::MM::ANY );
	if ( a1 >= 0 )
	  atoms.push_back( mol[c1][r1][a1].coord_orth().coord_frac(cell) );
    }
    for ( int c2 = c1+1; c2 < mol.size(); c2++ ) {
      std::cout << c1 << " " << c2 << std::endl;

      // test if this chain overlaps the target chain
      bool overlap = false;
      for ( int r2 = 0; r2 < mol[c2].size(); r2++ ) {
	int a2 = mol[c2][r2].lookup( "CA", clipper::MM::ANY );
	if ( a2 >= 0 ) {
	  f2 = mol[c2][r2][a2].coord_orth().coord_frac(cell);
	  for ( int i1 = 0; i1 < atoms.size(); i1++ ) {
	    f1 = atoms[i1];
	    f2 = f2.symmetry_copy_near( spgr, cell, f1 );
	    double d2 = ( f2 - f1 ).lengthsq( cell );
	    if ( d2 < rad*rad ) { overlap = true; break; }
	  }
	}
	if ( overlap ) break;
      }

      // chains overlap - try combining them
      if ( overlap ) {
	clipper::MiniMol mol_wrk( spgr, cell );
	mol_wrk.insert( mol[c1] );
	mol_wrk.insert( mol[c2] );
	Ca_join::join( mol_wrk, 2.0, 2.0, com );
	if ( mol_wrk.size() > 0 ) {
	  Ca_sequence::prepare_scores( mol_wrk[0], xmap, llksample );
	  Ca_sequence::sequence( mol_wrk[0], seq, reliability_ );
	  // test if the chain is improved
	  const clipper::MPolymer& mp1 = mol[c1];
	  const clipper::MPolymer& mp2 = mol_wrk[0];
	  int l0, l1, s0, s1;
	  l0 = mp1.size();
	  l1 = mp2.size();
	  s0 = s1 = 0;
	  for ( int r = 0; r < mp1.size(); r++ )
	    if ( ProteinTools::residue_index_3( mp1[r].type() ) >= 0 ) s0++;
	  for ( int r = 0; r < mp2.size(); r++ )
	    if ( ProteinTools::residue_index_3( mp2[r].type() ) >= 0 ) s1++;
	  // if new chain is better, keep it
	  if ( s1 > s0 ) mol[c1] = mol_wrk[0];
	}
      }
    }
  }

  ProteinTools::chain_tidy( mol );  // split chains
  Ca_prune::prune( mol );
  Ca_build::build( mol, xmap );
  ProteinTools::chain_tidy( mol );  // rename chains
  return true;

  /*
    clipper::MiniMol mol_nb( spgr, cell );
    mol_nb.insert( mol[c1].select("* CA") );  !!!!Correct
    clipper::MAtomNonBond nb( mol_nb, 2.0*rad );
	  f1 = mol[c2][r2][a2].coord_orth().coord_frac(cell);
	  atoms = nb.atoms_near( mol[c2][r2][a2].coord_orth(), rad );
	  for ( int i = 0; i < atoms.size(); i++ ) {
	    const clipper::MAtom& atom = mol_nb.atom(atoms[i]);
	    f2 = atom.coord_orth().coord_frac(cell);
	    f2 = spgr.symop(atoms[i].symmetry()) * f2;
	    f2 = f2.lattice_copy_near( f1 );
	    double d2 = ( f2 - f1 ).lengthsq( cell );
	    if ( d2 < rad*rad ) overlap = true;
	  }
  */
}
