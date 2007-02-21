/*! \file buccaneer-find.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-find.h"


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


// find Ca groups


// convenient analytial approximate distance funtion
double prob_dist( double x ) { return 0.999*exp(-75.0*pow(x-3.50,2.0)*pow(x,-2.5))+0.001; }


bool Ca_find::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget )
{
  const clipper::Spacegroup&    spgr = xmap.spacegroup();
  const clipper::Cell&          cell = xmap.cell();
  const clipper::Grid_sampling& grid = xmap.grid_sampling();

  // make prior map
  clipper::Xmap<float> prior( spgr, cell, grid );

  // do non-bond search
  clipper::MAtomNonBond nb( mol1, 12.0 );
  std::vector<clipper::MAtomIndexSymmetry> atoms;
  clipper::Coord_orth o1, o2;
  clipper::Coord_frac f1, f2;
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  for ( MRI ix = prior.first(); !ix.last(); ix.next() ) {
    double d2min = pow( 12.0, 2.0 );
    f1 = ix.coord().coord_frac( grid );
    o1 = f1.coord_orth( cell );
    atoms = nb( o1, 12.0 );
    for ( int i = 0; i < atoms.size(); i++ ) {
      const clipper::MAtom& atom =
	mol1[atoms[i].polymer()][atoms[i].monomer()][atoms[i].atom()];
      o2 = atom.coord_orth();
      f2 = o2.coord_frac( cell );
      f2 = spgr.symop(atoms[i].symmetry()) * f2;
      f2 = f2.lattice_copy_near( f1 );
      double d2 = ( f2 - f1 ).lengthsq( cell );
      if ( d2 < d2min ) d2min = d2;
    }
    prior[ix] = prob_dist( sqrt(d2min) );
  }

  // turn the prior into a z-score
  for ( MRI ix = prior.first(); !ix.last(); ix.next() )
    prior[ix] = -log( prior[ix] );
  clipper::Map_stats zstats( prior );
  double std_dev = clipper::Util::max( zstats.std_dev(), 1.0e-12 );
  for ( MRI ix = prior.first(); !ix.last(); ix.next() )
    prior[ix] = ( prior[ix] - zstats.mean() ) / std_dev;

  // do the fffear search:
  const double step = 24.0;
  const double dres = 5.0;
  // make a list of rotation ops to try
  std::vector<clipper::RTop_orth> ops = llktarget.rtop_list( spgr, step );

  // now search for ML target in each orientation in turn
  if ( resultscr.is_null() ) {
    resultscr.init( spgr, cell, grid );
    resultrot.init( spgr, cell, grid );
    resulttrn.init( spgr, cell, grid );
    llktarget.search( resultscr, resultrot, resulttrn, xmap, ops );
  }

  /*
  clipper::Xmap<float> resultscr( spgr, cell, grid );
  clipper::Xmap<int>   resultrot( spgr, cell, grid );
  clipper::Xmap<int>   resulttrn( spgr, cell, grid );
  clipper::Xmap<float> resultp1 ( clipper::Spacegroup::p1(), cell, grid );
  clipper::Xmap<float>::Map_reference_index i1(resultp1);
  clipper::Xmap<float>::Map_reference_coord ix(resultscr);
  resultscr = 1.0e20;

  // set up z scoring
  clipper::FFFear_fft<float> srch( xmap );
  LLK_map_target llktgt = llktarget;
  clipper::NX_operator nxop( xmap, llktgt.llk_target(), ops[0] );
  srch( resultp1, llktgt.llk_target(), llktgt.llk_weight(), nxop );
  zstats = clipper::Map_stats( resultp1 );

  // loop over orienatations
  for ( int op = 0; op < ops.size(); op++ ) {
    // do the fffear search
    clipper::NX_operator nxop( xmap, llktgt.llk_target(), ops[op].inverse() );
    srch( resultp1, llktgt.llk_target(), llktgt.llk_weight(), nxop );

    // store best scores
    for ( i1 = resultp1.first(); !i1.last(); i1.next() ) {
      ix.set_coord( i1.coord() );
      float score = ( resultp1[i1] - zstats.mean() ) / zstats.std_dev();
      score += prior[ix];
      if ( score < resultscr[ix] ) {
	resultscr[ix] = score;
	resultrot[ix] = op;
	resulttrn[ix] = grid.index( i1.coord() );
      }
    }
  }

  // now create a long scores list from the maps of results
  Score_list<clipper::RTop_orth> score_long( 5*nfind );
  for ( MRI ix = resultscr.first(); !ix.last(); ix.next() ) {
    float score = resultscr[ix];
    clipper::RTop_orth  rtop = ops[ resultrot[ix] ].inverse();
    clipper::Coord_grid cg = grid.deindex( resulttrn[ix] );
    rtop.trn() = xmap.coord_orth( cg.coord_map() );
    score_long.add( score, rtop );
  }
  */

  // now create a long scores list from the maps of results
  Score_list<clipper::RTop_orth> score_long( 5*nfind );
  for ( MRI ix = resultscr.first(); !ix.last(); ix.next() ) {
    float score = resultscr[ix] + prior[ix];
    clipper::RTop_orth  rtop = ops[ resultrot[ix] ].inverse();
    clipper::Coord_grid cg = grid.deindex( resulttrn[ix] );
    rtop.trn() = xmap.coord_orth( cg.coord_map() );
    score_long.add( score, rtop );
  }

  // create a pruned scores list omitting near-clashes
  Score_list<clipper::RTop_orth> score_trim( nfind );
  for ( int i = 0; i < score_long.size(); i++ ) {
    bool clash = false;
    for ( int j = 0; j < score_trim.size(); j++ )
      if ( clipper::Coord_orth( (score_trim[j].trn()-score_long[i].trn()) ).lengthsq() < dres*dres ) clash = true;
    if ( !clash ) score_trim.add( score_long.score(i), score_long[i] );
  }

  // now refine the best matches
  Score_list<clipper::RTop_orth> score_list( score_trim.size() );
  for ( int i = 0; i < score_trim.size(); i++ ) {
    Target_fn_refine_llk_map_target tgt( xmap, llktarget, 0.2, 0.2 );
    clipper::RTop_orth rtop = tgt.refine( score_trim[i] );
    double score = llktarget.llk( xmap, rtop );
    score_list.add( score, rtop );
  }

  // Now we build a model for output
  mol2 = clipper::MiniMol( xmap.spacegroup(), xmap.cell() );
  for ( int p = 0; p < mol1.size(); p++ ) mol2.insert( mol1[p] );

  clipper::MPolymer chain;
  chain.set_id(" ");
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
    chain.insert( residue );
  }
  mol2.insert( chain );

  return true;
}
