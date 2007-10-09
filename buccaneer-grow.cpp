/*! \file buccaneer-grow.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-grow.h"

#include <clipper/clipper-contrib.h>


Ca_grow::Ca_grow( int n_grow )
{
  ngrow = n_grow;
  max_conf1 = 50;
  max_conf2 = 30;
  rama1 = clipper::Ramachandran( clipper::Ramachandran::All );
  rama2 = clipper::Ramachandran( clipper::Ramachandran::NonGly );
}


bool Ca_grow::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget ) const
{
  // loop over chains, finding starting chains to expand
  std::vector<Ca_chain> chains;
  Ca_chain chain;
  for ( int chn = 0; chn < mol1.size(); chn++ )
    for ( int res = 0; res < mol1[chn].size(); res++ ) {
      // find ca, c, n
      int index_n  = mol1[chn][res].lookup( " N  ", clipper::MM::ANY );
      int index_ca = mol1[chn][res].lookup( " CA ", clipper::MM::ANY );
      int index_c  = mol1[chn][res].lookup( " C  ", clipper::MM::ANY );
      // if we have all three atoms, then add residue
      if ( index_ca >= 0 && index_c >= 0 && index_n >= 0 ) {
	clipper::Coord_orth coord_n  = mol1[chn][res][index_n].coord_orth();
	clipper::Coord_orth coord_ca = mol1[chn][res][index_ca].coord_orth();
	clipper::Coord_orth coord_c  = mol1[chn][res][index_c].coord_orth();
	// check if we must start a new chain
	if ( chain.size() > 0 )
	  if ( (chain.back().coord_c()-coord_n).lengthsq() > 2.25 ) {
	    chains.push_back( chain );
	    chain.clear();
	  }
	// add this atom to chain
	chain.push_back( Ca_group( coord_n, coord_ca, coord_c ) );
      }
    }
  if ( chain.size() > 0 ) chains.push_back( chain );

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
  double cutoff = scr[50];  // select score for top 1% of cases

  // now build up to 20 residues each way
  Ca_group ca;
  for ( int chn = 0; chn < chains.size(); chn++ ) {
    for ( int i = 0; i < ngrow; i++ ) {
      ca = next_ca_group( chains[chn], xmap, llktarget );
      if ( llktarget.llk( xmap, ca.rtop_from_std_ori() ) > cutoff ) break;
      chains[chn].push_back( ca );
    }
    for ( int i = 0; i < ngrow; i++ ) {
      ca = prev_ca_group( chains[chn], xmap, llktarget );
      if ( llktarget.llk( xmap, ca.rtop_from_std_ori() ) > cutoff ) break;
      chains[chn].push_front( ca );
    }
  }

  // make a new mmdb
  mol2 = clipper::MiniMol( mol1.spacegroup(), mol1.cell() );
  for ( int chn = 0; chn < chains.size(); chn++ ) {
    clipper::MPolymer chain;
    int ires = 1;
    for ( int res = 0; res < chains[chn].size(); res++ ) {
      clipper::MMonomer residue;
      residue.set_type("UNK");
      clipper::MAtom atom = clipper::Atom::null();
      atom.set_occupancy(1.0);
      atom.set_u_iso( 0.5 );

      atom.set_element( "N" );
      atom.set_id( "N" );
      atom.set_coord_orth( chains[chn][res].coord_n() );
      residue.insert( atom );

      atom.set_element( "C" );
      atom.set_id( "CA" );
      atom.set_coord_orth( chains[chn][res].coord_ca() );
      residue.insert( atom );

      atom.set_id( "C" );
      atom.set_coord_orth( chains[chn][res].coord_c() );
      residue.insert( atom );

      residue.set_seqnum( ires++ );
      chain.insert( residue );
    }
    mol2.insert( chain );
  }

  // restore the residue types, if any
  ProteinTools::copy_residue_types( mol2, mol1 );
  return true;
}


Ca_group Ca_grow::next_ca_group( const Ca_chain& chain, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget ) const
{
  Rama_ang1 conf1; Rama_ang2 conf2;
  Score_list<Rama_ang1> scores_l1( max_conf1 );
  Score_list<Rama_ang2> scores_l2( max_conf2 );
  Ca_group ca0, ca1, ca2;
  double r1, r2;
  double phi0 = -9.999; // prev Ramachandran angle
  if ( chain.size() > 1 ) phi0 = chain.ramachandran_phi( chain.size()-1 );
  ca0 = chain.back();  // start residue
  const double deg360 = clipper::Util::d2rad(359.0);
  const double deg20  = clipper::Util::d2rad( 20.0);
  const double deg30  = clipper::Util::d2rad( 30.0);
  // search all conformations of first residue
  for ( conf1.r1.psi = 0.0; conf1.r1.psi < deg360; conf1.r1.psi += deg20 )
    for ( conf1.r1.phi = 0.0; conf1.r1.phi < deg360; conf1.r1.phi += deg20 )
      if ( phi0 < -6.283 || rama1.allowed( phi0, conf1.r1.psi ) ) {
	ca1 = ca0.next_ca_group( conf1.r1.psi, conf1.r1.phi );
	r1 = llktarget.llk_approx( xmap, ca1.rtop_from_std_ori() );
	scores_l1.add( r1, conf1 );
      }
  // seach all conformations of second residue using best confirmations of first
  for ( int l1 = 0; l1 < scores_l1.size(); l1++ ) {
    r1 = scores_l1.score(l1);
    conf2.r1 = scores_l1[l1].r1;
    ca1 = ca0.next_ca_group( conf2.r1.psi, conf2.r1.phi );
    for ( conf2.r2.psi = 0.0; conf2.r2.psi < deg360; conf2.r2.psi += deg20 )
      for ( conf2.r2.phi = 0.0; conf2.r2.phi < deg360; conf2.r2.phi += deg30 )
	if ( rama2.favored( conf2.r1.phi, conf2.r2.psi ) ) {
	  ca2 = ca1.next_ca_group( conf2.r2.psi, conf2.r2.phi );
	  r2 = llktarget.llk_approx( xmap, ca2.rtop_from_std_ori() );
	  scores_l2.add( r1+r2, conf2 );
	}
  }
  //return ca0.next_ca_group( scores_l2[0].r1.psi, scores_l2[0].r1.phi );
  // now calculate full likelihood scores and pick best
  double ll_best = 1.0e6;
  Rama_ang2 ra_best = scores_l2[0];
  for ( int l2 = 0; l2 < scores_l2.size(); l2++ ) {
    ca1 = ca0.next_ca_group( scores_l2[l2].r1.psi, scores_l2[l2].r1.phi );
    ca2 = ca1.next_ca_group( scores_l2[l2].r2.psi, scores_l2[l2].r2.phi );
    r1 = ( llktarget.llk( xmap, ca1.rtop_from_std_ori() ) +
	   llktarget.llk( xmap, ca2.rtop_from_std_ori() ) );
    if ( r1 < ll_best ) {
      ll_best = r1;
      ra_best = scores_l2[l2];
    }
  }
  // refine the result
  Target_fn_refine_c_terminal_build tgt( xmap, llktarget, rama1, rama2, 0.1 );
  std::vector<double> args( tgt.num_params() );
  args[0] = ra_best.r1.psi;
  args[1] = ra_best.r1.phi;
  args[2] = ra_best.r2.psi;
  args[3] = ra_best.r2.phi;
  args = tgt.refine( chain, args );
      ca1 = ca0.next_ca_group( args[0], args[1] );
      ca2 = ca1.next_ca_group( args[2], args[3] );
      r2 = ( llktarget.llk( xmap, ca1.rtop_from_std_ori() ) +
	     llktarget.llk( xmap, ca2.rtop_from_std_ori() ) );
  // and return it
  return ca1;
}


Ca_group Ca_grow::prev_ca_group( const Ca_chain& chain, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget ) const
{
  Rama_ang1 conf1; Rama_ang2 conf2;
  Score_list<Rama_ang1> scores_l1( max_conf1 );
  Score_list<Rama_ang2> scores_l2( max_conf2 );
  Ca_group ca0, ca1, ca2;
  double r1, r2;
  double psi0 = -9.999; // prev Ramachandran angle
  if ( chain.size() > 1 ) psi0 = chain.ramachandran_psi( 0 );
  ca0 = chain.front();  // start residue
  const double deg360 = clipper::Util::d2rad(359.0);
  const double deg20  = clipper::Util::d2rad( 20.0);
  const double deg30  = clipper::Util::d2rad( 30.0);
  // search all conformations of first residue
  for ( conf1.r1.phi = 0.0; conf1.r1.phi < deg360; conf1.r1.phi += deg20 ) 
    for ( conf1.r1.psi = 0.0; conf1.r1.psi < deg360; conf1.r1.psi += deg20 )
      if ( psi0 < -6.283 || rama1.allowed( conf1.r1.phi, psi0 ) ) {
      ca1 = ca0.prev_ca_group( conf1.r1.phi, conf1.r1.psi );
      r1 = llktarget.llk_approx( xmap, ca1.rtop_from_std_ori() );
      scores_l1.add( r1, conf1 );
    }
  // seach all conformations of second residue using best confirmations of first
  for ( int l1 = 0; l1 < scores_l1.size(); l1++ ) {
    r1 = scores_l1.score(l1);
    conf2.r1 = scores_l1[l1].r1;
    ca1 = ca0.prev_ca_group( conf2.r1.phi, conf2.r1.psi );
    for ( conf2.r2.phi = 0.0; conf2.r2.phi < deg360; conf2.r2.phi += deg20 )
      for ( conf2.r2.psi = 0.0; conf2.r2.psi < deg360; conf2.r2.psi += deg30 )
	if ( rama2.favored( conf2.r2.phi, conf2.r1.psi) ) {
	  ca2 = ca1.prev_ca_group( conf2.r2.phi, conf2.r2.psi );
	  r2 = llktarget.llk_approx( xmap, ca2.rtop_from_std_ori() );
	  scores_l2.add( r1+r2, conf2 );
	}
  }
  //return ca0.prev_ca_group( scores_l2[0].r1.phi, scores_l2[0].r1.psi );
  // now calculate full likelihood scores and pick best
  double ll_best = 1.0e6;
  Rama_ang2 ra_best = scores_l2[0];
  for ( int l2 = 0; l2 < scores_l2.size(); l2++ ) {
    ca1 = ca0.prev_ca_group( scores_l2[l2].r1.phi, scores_l2[l2].r1.psi );
    ca2 = ca1.prev_ca_group( scores_l2[l2].r2.phi, scores_l2[l2].r2.psi );
    r1 = ( llktarget.llk( xmap, ca1.rtop_from_std_ori() ) +
	   llktarget.llk( xmap, ca2.rtop_from_std_ori() ) );
    if ( r1 < ll_best ) {
      ll_best = r1;
      ra_best = scores_l2[l2];
    }
  }
  // refine the result
  Target_fn_refine_n_terminal_build tgt( xmap, llktarget, rama1, rama2, 0.1 );
  std::vector<double> args( tgt.num_params() );
  args[0] = ra_best.r1.phi;
  args[1] = ra_best.r1.psi;
  args[2] = ra_best.r2.phi;
  args[3] = ra_best.r2.psi;
  args = tgt.refine( chain, args );
      ca1 = ca0.prev_ca_group( args[0], args[1] );
      ca2 = ca1.prev_ca_group( args[2], args[3] );
      r2 = ( llktarget.llk( xmap, ca1.rtop_from_std_ori() ) +
	     llktarget.llk( xmap, ca2.rtop_from_std_ori() ) );
  // and return it
  return ca1;
}


Target_fn_refine_n_terminal_build::Target_fn_refine_n_terminal_build( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const double& rot_step )
{
  xmap_ = &xmap;
  llktarget_ = &llktarget;
  rama1_ = &rama1;
  rama2_ = &rama2;
  rot_step_ = rot_step;
}

double Target_fn_refine_n_terminal_build::operator() ( const std::vector<double>& args ) const
{
  const Ca_chain& chain = *chain_;
  const clipper::Ramachandran& rama1 = *rama1_;
  const clipper::Ramachandran& rama2 = *rama2_;
  double psi0 = -9.999; // prev Ramachandran angle
  if ( chain.size() > 1 ) psi0 = chain.ramachandran_psi( 0 );
  const Ca_group ca0 = chain.front();
  const Ca_group ca1 = ca0.prev_ca_group( args[0], args[1] );
  const Ca_group ca2 = ca1.prev_ca_group( args[2], args[3] );
  double r = ( (*llktarget_).llk( *xmap_, ca1.rtop_from_std_ori() ) +
	       (*llktarget_).llk( *xmap_, ca2.rtop_from_std_ori() ) );
  if ( !rama1.allowed( args[0], psi0 ) ) 
    if ( psi0 > -6.283 ) r += 10.0;
  if ( !rama2.favored( args[2], args[1] ) )
    r += 10.0;
  return r;
}

std::vector<double> Target_fn_refine_n_terminal_build::refine( const Ca_chain& chain, const std::vector<double>& args )
{
  // store initial chain
  chain_ = &chain;
  // calculate initial params
  std::vector<double> arg_init;
  std::vector<std::vector<double> > args_init;
  args_init.push_back( args );
  for ( int i = 0; i < num_params(); i++ ) {
    arg_init = args;
    arg_init[i] += rot_step_;
    args_init.push_back( arg_init );
  }
  // simple refinement
  double tol = 0.005 * (*this)( args_init[0] );
  Optimiser_simplex os( tol, 50, Optimiser_simplex::GRADIENT );
  return os( *this, args_init );
}


Target_fn_refine_c_terminal_build::Target_fn_refine_c_terminal_build( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const double& rot_step )
{
  xmap_ = &xmap;
  llktarget_ = &llktarget;
  rama1_ = &rama1;
  rama2_ = &rama2;
  rot_step_ = rot_step;
}

double Target_fn_refine_c_terminal_build::operator() ( const std::vector<double>& args ) const
{
  const Ca_chain& chain = *chain_;
  const clipper::Ramachandran& rama1 = *rama1_;
  const clipper::Ramachandran& rama2 = *rama2_;
  double phi0 = -9.999; // prev Ramachandran angle
  if ( chain.size() > 1 ) phi0 = chain.ramachandran_phi( chain.size()-1 );
  const Ca_group ca0 = chain.back();
  const Ca_group ca1 = ca0.next_ca_group( args[0], args[1] );
  const Ca_group ca2 = ca1.next_ca_group( args[2], args[3] );
  double r = ( (*llktarget_).llk( *xmap_, ca1.rtop_from_std_ori() ) +
	       (*llktarget_).llk( *xmap_, ca2.rtop_from_std_ori() ) );
  if ( !rama1.allowed( phi0,    args[0] ) ) 
    if ( phi0 > -6.283 ) r += 10.0;
  if ( !rama2.favored( args[1], args[2] ) )
    r += 10.0;
  return r;
}

std::vector<double> Target_fn_refine_c_terminal_build::refine( const Ca_chain& chain, const std::vector<double>& args )
{
  // store initial chain
  chain_ = &chain;
  // calculate initial params
  std::vector<double> arg_init;
  std::vector<std::vector<double> > args_init;
  args_init.push_back( args );
  for ( int i = 0; i < num_params(); i++ ) {
    arg_init = args;
    arg_init[i] += rot_step_;
    args_init.push_back( arg_init );
  }
  // simple refinement
  double tol = 0.005 * (*this)( args_init[0] );
  Optimiser_simplex os( tol, 50, Optimiser_simplex::GRADIENT );
  return os( *this, args_init );
}
