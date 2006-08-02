/*! \file buccaneer-sequence.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-sequence.h"

#include <clipper/clipper-contrib.h>


Score_list<clipper::String> Ca_sequence::sequence_matches( const std::vector<std::vector<double> >& scores, const clipper::MMoleculeSequence& seq )
{
  // prepare data for partial sequence match
  clipper::String nulseq = std::string( scores.size(), '?' );
  // find the best match of this chain against the sequence
  std::vector<double> score( scores.size() );
  std::vector<double> sccum( scores.size()+1 );
  Score_list<clipper::String> matches_tmp(100);
  for ( int seqchn = 0; seqchn < seq.size(); seqchn++ ) {
    clipper::String chnseq = nulseq + seq[seqchn].sequence() + nulseq;
    int maxoff = chnseq.length()-scores.size();
    for ( int seqoff = 0; seqoff <= maxoff; seqoff++ ) {
      // score the whole sequence
      clipper::String subseq = chnseq.substr( seqoff, scores.size() );
      for ( int r = 0; r < scores.size(); r++ ) {
	int t = ProteinTools::residue_index( subseq.substr( r, 1 ) );
	if ( t >= 0 ) score[r] = scores[r][t];  // add residue z-score
	else          score[r] = 0.10;          // or a penalty of +0.10
      }
      // calculate cumulative score
      sccum[0] = 0.0;
      for ( int r = 0; r < score.size(); r++ ) sccum[r+1] = sccum[r]+score[r];
      // now find grestest subsequence score
      int minseq = 1;
      int r1min = 0;
      int r2min = sccum.size()-1;
      double scmin = 0.0;
      for ( int r1 = 0; r1 < sccum.size()-minseq; r1++ )
	for ( int r2 = r1+minseq; r2 < sccum.size(); r2++ ) {
	  double l1 = ( r2 - r1 )/50.0;  // downweight very long sequences
	  double sc = ( sccum[r2]-sccum[r1] ) / pow( 1.0+l1*l1 , 0.25 );
	  if ( sc < scmin ) { r1min = r1; r2min = r2; scmin = sc; }
	}
      subseq = ( nulseq.substr(0,r1min) +
		 subseq.substr(r1min,r2min-r1min) +
		 nulseq.substr(r2min) );
      matches_tmp.add( scmin, subseq );
    }
  }

  // eliminate any duplicates
  Score_list<clipper::String> matches(10);
  for ( int i = 0; i < matches_tmp.size(); i++ ) {
    bool clash = false;
    for ( int j = 0; j < matches.size(); j++ )
      if ( matches[j] == matches_tmp[i] ) clash = true;
    if ( !clash )
      matches.add( matches_tmp.score(i), matches_tmp[i] );
  }

  return matches;
}


bool Ca_sequence::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const LLK_map_classifier& llktarget, const clipper::MMoleculeSequence& seq )
{
  typedef clipper::MMonomer Mm;

  clipper::Cell       cell = mol1.cell();
  clipper::Spacegroup spgr = mol1.spacegroup();

  // split into separate chains
  ProteinTools::chain_tidy( mol2, mol1 );

  // set relibility criteria and docking count
  double score_cutoff = -20.0 * reliability_;  // min score
  double screl_cutoff =   0.3 * reliability_;  // min diff between 1st & 2nd
  num_seq = 0;

  // now loop over chains
  clipper::Coord_orth coord_n, coord_ca, coord_c;
  for ( int chn = 0; chn < mol2.size(); chn++ ) {
    // for each chain, classify each residue against the density
    std::vector<double> scores_null( llktarget.num_targets(), 0.0 );
    std::vector<std::vector<double> > scores;
    for ( int res = 0; res < mol2[chn].size(); res++ ) {
      // find ca, c, n
      int index_n  = mol2[chn][res].lookup( " N  ", clipper::MM::ANY );
      int index_ca = mol2[chn][res].lookup( " CA ", clipper::MM::ANY );
      int index_c  = mol2[chn][res].lookup( " C  ", clipper::MM::ANY );
      // if we have all three atoms, then add residue
      if ( index_ca >= 0 && index_c >= 0 && index_n >= 0 ) {
	coord_n  = mol2[chn][res][index_n].coord_orth();
	coord_ca = mol2[chn][res][index_ca].coord_orth();
	coord_c  = mol2[chn][res][index_c].coord_orth();
	Ca_group ca( coord_n, coord_ca, coord_c );
	scores.push_back( llktarget.llk_raw( xmap, ca.rtop_beta_carbon() ) );
      } else {
	scores.push_back( scores_null );
      }
    }

    // normalise down columns by mean and variance to get z-scores
    double s0, s1, s2;
    s0 = double( scores.size() );
    for ( int t = 0; t < llktarget.num_targets(); t++ ) {
      s1 = s2 = 0.0;
      for ( int r = 0; r < scores.size(); r++ ) {
	s1 += scores[r][t];
	s2 += scores[r][t]*scores[r][t];
      }
      s1 /= s0;
      s2 /= s0;
      s2 = sqrt( s2 - s1*s1 );
      for ( int r = 0; r < scores.size(); r++ )
	scores[r][t] = ( scores[r][t] - s1 ) / s2;
    }

    // do sequence match
    Score_list<clipper::String> matches = sequence_matches( scores, seq );

    // test the reliability of the sequence match
    if ( matches.size() >= 2 ) {
      if ( matches.score(0) < score_cutoff &&
	   (1.0-screl_cutoff)*matches.score(0) < matches.score(1) ) {
	// copy the residue types from the best match into the model
	clipper::String bestseq = matches[0];
	for ( int res = 0; res < mol2[chn].size(); res++ ) {
	  int t = ProteinTools::residue_index( bestseq.substr(res,1) );
	  if ( t >= 0 ) {
	    mol2[chn][res].set_type( ProteinTools::residue_code_3( t ) );
	    num_seq++;
	  }
	}
      }
    }

    // save some info for future output
    history.push_back( matches );
  }

  return true;
}


int Ca_sequence::num_sequenced() const
{
  return num_seq;
}


clipper::String Ca_sequence::format() const
{
  clipper::String result = "";
  for ( int chn = 0; chn < history.size(); chn++ ) {
    result += "Chain number: " + clipper::String( chn, 4 ) + "    length: " + clipper::String( int(history[chn][0].length()) ) + "\n";
    for ( int res = 0; res < history[chn].size(); res++ ) {
      result += history[chn][res] + " \t" +
	clipper::String( history[chn].score(res), 10, 6 ) + "\n";
    }
  }
  return result;
}
