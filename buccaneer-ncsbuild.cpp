/*! \file buccaneer-ncsbuild.cpp buccaneer library */
/* (C) 2006-2020 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-ncsbuild.h"


#include "buccaneer-join.h"
#include "buccaneer-sequence.h"
#include "buccaneer-filter.h"


bool Ca_ncsbuild::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq ) const
{
  clipper::Cell       cell = xmap.cell();
  clipper::Spacegroup spgr = xmap.spacegroup();
  clipper::MiniMol mol_old = mol;

  // get center of mass
  clipper::Coord_orth com( 0.0, 0.0, 0.0 );
  clipper::Atom_list atoms = mol.atom_list();
  for ( int i = 0; i < atoms.size(); i++ ) com = com + atoms[i].coord_orth();
  com = ( 1.0 / double( atoms.size() ) ) * com;

  // extract the necessary bits of the likelihood targets
  std::vector<LLK_map_target::Sampled> llksample( llktarget.size() );
  for ( int t = 0; t < llktarget.size(); t++ )
    llksample[t] = llktarget[t].sampled();

  // split into separate chains
  ProteinTools::split_chains_at_gap( mol );

  // now loop over chains
  for ( int chn1 = 0; chn1 < mol.size(); chn1++ ) {
    // assemble combined model for this chain
    clipper::MPolymer mp1, mp2;
    mp1 = mol[chn1];
    clipper::MiniMol mol_wrk;
    mol_wrk.init( mol.spacegroup(), mol.cell() );
    mol_wrk.insert( mp1 );
    // add any other matching chains
    for ( int chn2 = 0; chn2 < mol.size(); chn2++ ) {
      if ( chn2 != chn1 ) {
        clipper::RTop_orth rtop =
          ProteinTools::superpose( mol[chn2], mol[chn1], rmsd_, nmin_, nmin_ );
        if ( !rtop.is_null() ) {
          clipper::MPolymer mp = mol[chn2];
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
      Ca_join::join( mol_wrk, 2.0, 2.0, com );
      if ( mol_wrk.size() > 0 ) {  // trap empty models
        if ( mol_wrk[0].size() > mp1.size() ) {  // optimisation
          // sequence
          Ca_sequence::prepare_scores( mol_wrk[0], xmap, llksample );
          Ca_sequence::sequence( mol_wrk[0], seq, reliability_ );
          // filter
          Ca_filter::filter( mol_wrk, xmap, 1.0 );
          // tidy
          ProteinTools::split_chains_at_gap( mol_wrk );
          // make new chain
          if ( mol_wrk.size() > 0 ) mp2 = mol_wrk[0];
          mp2.copy( mp1, clipper::MM::COPY_MP );
          // test if the chain is improved
          int l0, l1, s0, s1;
          l0 = mp1.size();
          l1 = mp2.size();
          s0 = s1 = 0;
          for ( int r = 0; r < mp1.size(); r++ )
            if ( ProteinTools::residue_index_3( mp1[r].type() ) >= 0 ) s0++;
          for ( int r = 0; r < mp2.size(); r++ )
            if ( ProteinTools::residue_index_3( mp2[r].type() ) >= 0 ) s1++;
          // if new chain is better, keep it
          if ( l1 > l0 && s1 > s0 ) {
            mol[chn1] = mp2;
          }
        }
      }
    }
  }

  // now we need to check if any later chains are duplicates of earlier, more reliable chains
  // get non-bond list for Ca's to eliminate trivial superpositions
  const double dmin = 1.5;
  clipper::MiniMol camodel( spgr, cell );
  camodel.model() = mol.select( "*/*/ CA ", clipper::MM::ANY );
  if ( camodel.size() != mol.size() ) clipper::Message::message( clipper::Message_fatal( "NCSbuild: internal error - length mismatch" ) );
  clipper::MAtomNonBond nb( camodel, 3.0 );
  std::vector<clipper::MAtomIndexSymmetry> atomnb;
  clipper::Matrix<int> contacts( mol.size(), mol.size(), 0 );
  for ( int c1 = 0; c1 < camodel.size(); c1++ ) {
    for ( int r1 = 0; r1 < camodel[c1].size(); r1++ ) {
      if ( camodel[c1][r1].size() >= 1 ) {
        atomnb = nb( camodel[c1][r1][0].coord_orth(), dmin );
        for ( int i = 0; i < atomnb.size(); i++ ) {
          const int c2 = atomnb[i].polymer();
          contacts(c1,c2) += 1;
        }
      }
    }
  }
  // for any heavily contacted chain, check if they superpose
  for ( int c2 = 0; c2 < camodel.size(); c2++ ) {
    for ( int c1 = 0; c1 < c2; c1++ ) {
      if ( contacts(c1,c2) > contacts(c2,c2)/2 ) {
        mol[c2] = mol_old[c2];
      }
    }
  }

  return true;
}
