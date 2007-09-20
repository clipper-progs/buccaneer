/*! \file buccaneer-correct.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-correct.h"

#include <clipper/clipper-contrib.h>


double Ca_correct::score_chain_sequence( const clipper::MPolymer& mp, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget )
{
  double z = 0.0;
  for ( int r = 0; r < mp.size(); r++ ) {
    // get residue type
    int t = ProteinTools::residue_index( mp[r].type() );
    // and score if recongnized
    if ( t >= 0 && t < llktarget.size() ) {
      // find ca, c, n
      int index_n  = mp[r].lookup( " N  ", clipper::MM::ANY );
      int index_ca = mp[r].lookup( " CA ", clipper::MM::ANY );
      int index_c  = mp[r].lookup( " C  ", clipper::MM::ANY );
      // if we have all three atoms, then score it
      if ( index_ca >= 0 && index_c >= 0 && index_n >= 0 ) {
	Ca_group ca( mp[r][index_n ].coord_orth(), mp[r][index_ca].coord_orth(),
		     mp[r][index_c ].coord_orth() );
        z += llktarget[t].llk( xmap, ca.rtop_beta_carbon() );
      }
    }
  }
  return z;
}


/* position is the index of the first Ca to move - the preceding C moves also */
std::pair<double,clipper::MPolymer> Ca_correct::best_rebuild_sequence( const int& position, const clipper::MPolymer& mp, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget )
{
  std::pair<double,clipper::MPolymer> best;
  best.first = 1.0e12;
  // find the atoms to rebuild
  int index_cc0 = mp[position-2].lookup( " C  ", clipper::MM::ANY );
  int index_cn1 = mp[position-1].lookup( " N  ", clipper::MM::ANY );
  int index_ca1 = mp[position-1].lookup( " CA ", clipper::MM::ANY );
  int index_cc1 = mp[position-1].lookup( " C  ", clipper::MM::ANY );
  int index_cn2 = mp[position  ].lookup( " N  ", clipper::MM::ANY );
  int index_ca2 = mp[position  ].lookup( " CA ", clipper::MM::ANY );
  int index_cc2 = mp[position  ].lookup( " C  ", clipper::MM::ANY );
  int index_cn3 = mp[position+1].lookup( " N  ", clipper::MM::ANY );
  int index_ca3 = mp[position+1].lookup( " CA ", clipper::MM::ANY );
  int index_cc3 = mp[position+1].lookup( " C  ", clipper::MM::ANY );
  int index_cn4 = mp[position+2].lookup( " N  ", clipper::MM::ANY );
  int index_ca4 = mp[position+2].lookup( " CA ", clipper::MM::ANY );
  int index_cc4 = mp[position+2].lookup( " C  ", clipper::MM::ANY );
  int index_cn5 = mp[position+3].lookup( " N  ", clipper::MM::ANY );
  if ( index_cc0 >= 0 &&
       index_cn1 >= 0 && index_ca1 >= 0 && index_cc1 >= 0 && 
       index_cn2 >= 0 && index_ca2 >= 0 && index_cc2 >= 0 && 
       index_cn3 >= 0 && index_ca3 >= 0 && index_cc3 >= 0 && 
       index_cn4 >= 0 && index_ca4 >= 0 && index_cc4 >= 0 && 
       index_cn5 >= 0 ) {
    // rebuild loop
    ProteinLoop pl( torsion_sampling_ );
    std::vector<ProteinLoop::CoordList<8> > r8;
    r8 = pl.rebuild8atoms( mp[position-2][index_cc0].coord_orth(),
			   mp[position-1][index_cn1].coord_orth(),
			   mp[position-1][index_ca1].coord_orth(),
			   mp[position+2][index_ca4].coord_orth(),
			   mp[position+2][index_cc4].coord_orth(),
			   mp[position+3][index_cn5].coord_orth() );
    // loop over results
    for ( int i = 0; i < r8.size(); i++ ) {
      // modify chain
      clipper::MPolymer mp1 = mp;
      mp1[position-1][index_cc1].set_coord_orth( r8[i][0] );
      mp1[position  ][index_cn2].set_coord_orth( r8[i][1] );
      mp1[position  ][index_ca2].set_coord_orth( r8[i][2] );
      mp1[position  ][index_cc2].set_coord_orth( r8[i][3] );
      mp1[position+1][index_cn3].set_coord_orth( r8[i][4] );
      mp1[position+1][index_ca3].set_coord_orth( r8[i][5] );
      mp1[position+1][index_cc3].set_coord_orth( r8[i][6] );
      mp1[position+2][index_cn4].set_coord_orth( r8[i][7] );
      // score the chain
      double z = score_chain_sequence( mp1, xmap, llktarget );
      // if this chain is best, then update
      if ( z < best.first ) {
	best.first  = z;
	best.second = mp1;
      }
    }  // done loop over fitted loops
  }  // if
  return best;
}


bool Ca_correct::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq )
{
  typedef clipper::MMonomer Mm;

  // split into separate chains
  ProteinTools::chain_tidy( mol2, mol1 );

  num_cor = 0;

  // now loop over chains
  int last_pos;
  for ( int chn = 0; chn < mol2.size(); chn++ ) {
    // extract the sequence of this chain
    clipper::MPolymer mpwrk = mol2[chn];
    clipper::String seq0 = ProteinTools::chain_sequence( mpwrk );

    // deal with insertions
    last_pos = 0;
    while ( seq0.find( "+", last_pos ) != std::string::npos ) {
      const int offset  = 8;    // how far to search for insertion/deletion
      const int buffer  = 6;    // number of residues padding at end
      const double zoff = 0.1;  // fudge factor to favour rebuilding
      int ins1 = seq0.find_first_of( "+", last_pos );
      int ins2 = seq0.find_first_not_of( "+", ins1 );
      int nins = ins2 - ins1;
      int len = mpwrk.size();
      if ( ins2 == std::string::npos ) break;
      // pick a range of positions to modify
      int r1, r2;
      int o1 = clipper::Util::max(ins1-offset,0);
      int o2 = clipper::Util::min(ins2+offset,int(seq0.length())-1);
      for ( r1 = ins1-1; r1 >  o1; r1-- ) if ( seq0[r1] == '?' ) break;
      for ( r2 = ins2;   r2 <= o2; r2++ ) if ( seq0[r2] == '?' ) break;
      r1 = clipper::Util::max( r1, buffer );
      r2 = clipper::Util::min( r2, int(seq0.length())-buffer );
      // OPTIMISATION: mask the sequence
      for ( int r = 0; r < mpwrk.size(); r++ )
	if ( r < r1-2 || r > r2+2 ) mpwrk[r].set_type( mpwrk[r].type() + "*" );
      // create null sequence
      clipper::MPolymer mpbest = mpwrk;
      double zbest = 0.0;
      double znull = score_chain_sequence( mpbest, xmap, llktarget ) + zoff;
      for ( int rl = r1; rl < r2; rl++ ) {
	// make a new chain with one less residue
	clipper::MPolymer mp;
	mp.set_id( mpwrk.id() );
	for ( int r = 0;       r < rl; r++ )   mp.insert( mpwrk[r] );
	for ( int r = rl+nins; r < len;  r++ ) mp.insert( mpwrk[r] );
	// set sequence
	for ( int r = 0; r < ins1; r++ )
	  mp[r].set_type( mpwrk[r].type() );
	for ( int r = ins2; r < len; r++ )
	  mp[r-nins].set_type( mpwrk[r].type() );
	// rebuild and score
	std::pair<double,clipper::MPolymer> result =
	  best_rebuild_sequence( rl, mp, xmap, llktarget );
	double z = result.first - znull;
	if ( z < zbest ) {
	  zbest = z;
	  mpbest = result.second;
	}
      }  // done loop over fitting positions
      // store the sequence
      if ( zbest < 0.0 ) {
	mpwrk = mpbest;
	num_cor++;
      }
      // OPTIMISATION: unmask the sequence
      for ( int r = 0; r < mpwrk.size(); r++ )
	mpwrk[r].set_type( mpwrk[r].type().substr(0,3) );
      // and update
      clipper::String seq0 = ProteinTools::chain_sequence( mpwrk );
      last_pos = ins2 + 1;
    }  // done insertion correction

    // deal with deletions
    last_pos = 0;
    while ( seq0.find( "-", last_pos ) != std::string::npos ) {
      const int offset  = 8;    // how far to search for insertion/deletion
      const int buffer  = 6;    // number of residues padding at end
      const double zoff = 0.1;  // fudge factor to favour rebuilding
      int del1 = seq0.find_first_of( "-", last_pos );
      int del2 = seq0.find_first_not_of( "-", del1 );
      int ndel = del2 - del1;
      int len = mpwrk.size();
      if ( del2 == std::string::npos ) break;
      // look for the missing sequence elements in the sequence data
      int i1, i2, f1, f2, e1, e2;
      for ( i1 = del1-1;   i1 > 0; i1-- ) if ( seq0[i1-1] == '?' ) break;
      for ( i2 = del2; i2 < len-1; i2++ ) if ( seq0[i2+1] == '?' ) break;
      clipper::String s1 = seq0.substr( i1  , del1-i1 );
      clipper::String s2 = seq0.substr( del2, i2-del2 );
      clipper::String seqx;
      for ( int c = 0; c < seq.size(); c++ ) {
	f1 = seq[c].sequence().find(s1);
	f2 = seq[c].sequence().find(s2);
	e1 = f1 + s1.length(); e2 = f2 + s2.length();
	if ( f1 != std::string::npos && f2 != std::string::npos )
	  if ( f2 > e1 && f2 <= f1 + s1.length() + 4 ) {
	    seqx = seq[c].sequence().substr( e1, f2-e1 );
	    break;
	  }
      }
      seqx += "????????";
      // pick a range of positions to modify
      int r1, r2;
      int o1 = clipper::Util::max(del1-offset,0);
      int o2 = clipper::Util::min(del2+offset,int(seq0.length())-1);
      for ( r1 = del1-1; r1 >  o1; r1-- ) if ( seq0[r1] == '?' ) break;
      for ( r2 = del2;   r2 <= o2; r2++ ) if ( seq0[r2] == '?' ) break;
      r1 = clipper::Util::max( r1, buffer );
      r2 = clipper::Util::min( r2, int(seq0.length())-buffer );
      // OPTIMISATION: mask the sequence
      for ( int r = 0; r < mpwrk.size(); r++ )
	if ( r < r1-2 || r > r2+2 ) mpwrk[r].set_type( mpwrk[r].type() + "*" );
      // create null sequence
      clipper::MPolymer mpbest = mpwrk;
      double zbest = 0.0;
      double znull = score_chain_sequence( mpbest, xmap, llktarget ) + zoff;
      for ( int rl = r1; rl < r2; rl++ ) {
	// make a new chain with one more residue
	clipper::MPolymer mp;
	mp.set_id( mpwrk.id() );
	for ( int r = 0;       r < rl; r++ )   mp.insert( mpwrk[r] );
	for ( int r = rl-ndel; r < len;  r++ ) mp.insert( mpwrk[r] );
	// set sequence
	for ( int r = 0; r < del1; r++ )
	  mp[r].set_type( mpwrk[r].type() );
	for ( int r = del1; r < del2+ndel; r++ ) {
	  int t = ProteinTools::residue_index( seqx.substr(r-del1,1) );
	  if ( t >= 0 ) mp[r].set_type( ProteinTools::residue_code_3( t ) );
	  else          mp[r].set_type( "UNK" );
	}
	for ( int r = del2; r < len; r++ )
	  mp[r+ndel].set_type( mpwrk[r].type() );
	// rebuild and score
	std::pair<double,clipper::MPolymer> result =
	  best_rebuild_sequence( rl, mp, xmap, llktarget );
	double z = result.first - znull;
	if ( z < zbest ) {
	  zbest = z;
	  mpbest = result.second;
	}
      }  // done loop over fitting positions
      // store the sequence
      if ( zbest < 0.0 ) {
	mpwrk = mpbest;
	num_cor++;
      }
      // OPTIMISATION: unmask the sequence
      for ( int r = 0; r < mpwrk.size(); r++ )
	mpwrk[r].set_type( mpwrk[r].type().substr(0,3) );
      // and update
      clipper::String seq0 = ProteinTools::chain_sequence( mpwrk );
      last_pos = del2 + ndel + 1;
    }  // done deletion correction

    // store the chain again
    mol2[chn] = mpwrk;
  }  // next chain

  return true;
}


int Ca_correct::num_corrected() const
{
  return num_cor;
}
