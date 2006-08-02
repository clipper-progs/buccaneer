/*! \file buccaneer-find.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-find.h"


bool Ca_find::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget ) const
{
  // do the fffear search
  Score_list<clipper::RTop_orth> score_temp =
    llktarget.search( xmap, 24.0, nfind, 5.0 );

  // now refine the best matches
  Score_list<clipper::RTop_orth> score_list( score_temp.size() );
  for ( int i = 0; i < score_temp.size(); i++ ) {
    Target_fn_refine_llk_map_target tgt( xmap, llktarget, 0.2, 0.2 );
    clipper::RTop_orth rtop = tgt.refine( score_temp[i] );
    double score = llktarget.llk( xmap, rtop );
    score_list.add( score, rtop );
  }

  // Now we build a model for output
  // Make an mmdb
  mol = clipper::MiniMol( xmap.spacegroup(), xmap.cell() );
  clipper::MPolymer chain;
  chain.set_id("A");
  mol.insert( chain );

  int ires = 0;
  clipper::Coord_orth ca( 0.0, 0.0, 0.0 );
  clipper::Coord_orth c( 0.85, 0.0, 1.2 );
  clipper::Coord_orth n( 0.85, 0.0, -1.2 );
  for ( int i = 0; i < score_list.size(); i++ ) {
    clipper::MMonomer residue;
    residue.set_type("UNK");
    clipper::MAtom atom = clipper::Atom::null();
    atom.set_occupancy(1.0);
    atom.set_u_iso( exp(-10.0*(score_list.score(i)-score_list.score(0))) );

    atom.set_element( "N" );
    atom.set_id( "N" );
    atom.set_coord_orth( score_list[i] * n );
    residue.insert( atom );

    atom.set_element( "C" );
    atom.set_id( "CA" );
    atom.set_coord_orth( score_list[i] * ca );
    residue.insert( atom );

    atom.set_id( "C" );
    atom.set_coord_orth( score_list[i] * c );
    residue.insert( atom );

    residue.set_seqnum( ires += 2 );
    mol[0].insert( residue );
  }

  return true;
}


// methods for refinement target


Target_fn_refine_llk_map_target::Target_fn_refine_llk_map_target( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const double& rot_step, const double& trn_step )
{
  xmap_ = &xmap;
  llktarget_ = &llktarget;
  rot_step_ = rot_step;
  trn_step_ = trn_step;
}

double Target_fn_refine_llk_map_target::operator() ( const clipper::RTop_orth& rtop ) const 
{
  return (*llktarget_).llk( *xmap_, rtop );
}

double Target_fn_refine_llk_map_target::operator() ( const std::vector<double>& args ) const
{
  return (*this)( rtop_orth( args ) );
}

clipper::RTop_orth Target_fn_refine_llk_map_target::rtop_orth( const std::vector<double>& args ) const
{
  return clipper::RTop_orth( clipper::Euler<clipper::Rotation::EulerXYZs>(args[0],args[1],args[2]).rotation().matrix() * rtop_.rot(), clipper::Coord_orth(args[3],args[4],args[5]) + rtop_.trn() );
}

clipper::RTop_orth Target_fn_refine_llk_map_target::refine( const clipper::RTop_orth& rtop )
{
  // store initial rtop
  rtop_ = rtop;

  // calculate initial params
  std::vector<std::vector<double> > args_init;
  std::vector<double> arg(6,0.0);
  // identity
  clipper::Euler<clipper::Rotation::EulerXYZs> euler( 0.0, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args_init.push_back( arg );
  // rotation steps
  double step = 0.5 * rot_step_;
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( step, 0.0, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args_init.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, step, 0.0 );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args_init.push_back( arg );
  euler = clipper::Euler<clipper::Rotation::EulerXYZs>( 0.0, 0.0, step );
  arg[0] = euler.alpha();
  arg[1] = euler.beta();
  arg[2] = euler.gamma();
  args_init.push_back( arg );
  // translation steps
  step = 0.5 * trn_step_;
  arg = args_init[0];
  arg[3] = step;
  arg[4] = 0.0;
  arg[5] = 0.0;
  args_init.push_back( arg );
  arg[3] = 0.0;
  arg[4] = step;
  arg[5] = 0.0;
  args_init.push_back( arg );
  arg[3] = 0.0;
  arg[4] = 0.0;
  arg[5] = step;
  args_init.push_back( arg );
  // simple refinement
  double tol = 0.005 * (*this)( args_init[0] );
  Optimiser_simplex os( tol, 50, Optimiser_simplex::GRADIENT );
  return rtop_orth( os( *this, args_init ) );
}
