/*! \file buccaneer-lib.cpp buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-lib.h"

#include <clipper/clipper-contrib.h>


/*! The target is constructed and initialised to zero. It can then be
  filled by accumulation or by loading the target and weight.
  \param rad Radius in which LLK function will be used, in Angstroms.
  \param sampling Sampling spacing for the LLK function, in Angstroms. */
LLK_map_target::LLK_map_target( const clipper::ftype& rad, const clipper::ftype& sampling )
{
  radius = rad;

  clipper::ftype extent = radius + 2.0 * sampling;  // model box size
  
  clipper::ftype ng = rint( extent / sampling );
  clipper::Grid grid( 2*int(ng)+1, 2*int(ng)+1, 2*int(ng)+1 );
  clipper::RTop<> rtop( clipper::Mat33<>( 1.0/sampling, 0.0, 0.0,
					  0.0, 1.0/sampling, 0.0,
					  0.0, 0.0, 1.0/sampling ),
			clipper::Vec3<>( ng, ng, ng ) );

  target.init( grid, rtop );
  weight.init( grid, rtop );

  naccum = 0;  // number of samples
}

/*! This most be called after loading and before saving or using the
  LLK targets in any way. i.e. if you accumulate the log likelihood
  target in one program and use it in another, then you must call
  prep_llk() both before saving the targets in the first program and
  after loading in the second program. */
void LLK_map_target::prep_llk()
{
  clipper::NXmap_base::Map_reference_index ix;

  // first make a llk target, if necessary
  if ( naccum != 0 ) {
    // aliases for maps
    clipper::NXmap<float>& mrho = target;
    clipper::NXmap<float>& mrho2 = weight;
    // calculate density moments for whole map
    clipper::ftype sn = 0.0, sr = 0.0, sr2 = 0.0;
    for ( ix = mrho.first(); !ix.last(); ix.next() )
      if ( mrho2[ix] > 0.0 ) {
	sn  += naccum;        // overall desnity stats
	sr  += mrho[ix];
	sr2 += mrho2[ix];
      }
    float rmap = sr/sn;
    float smap = sqrt( sn*sr2 - sr*sr )/sn;
    // and for individual positions within map
    for ( ix = mrho.first(); !ix.last(); ix.next() )
      if ( mrho2[ix] > 0.0 ) {
	mrho[ix]  /= naccum;  // local density stats
	mrho2[ix] /= naccum;  // limit std to no less than 3% of map std
	mrho2[ix] = sqrt( clipper::Util::max( mrho2[ix] - mrho[ix]*mrho[ix],
					      0.001f*smap*smap ) );
      }
    // convert density moments to llk target
    float rho, std, v1, v2, w1, w2;
    for ( ix = mrho.first(); !ix.last(); ix.next() )
      if ( mrho2[ix] > 0.0 ) {
	rho = mrho[ix];
	std = mrho2[ix];
	v1 = std*std;
	v2 = smap*smap;
	w1 = clipper::Util::max( v2/v1-1.0, 0.001 );  // weight must be +ve
	w2 = clipper::Util::min( 1.0/w1, 2.0 );  // limit density upweighting
	target[ix] = rho + w2*(rho-rmap);
	weight[ix] = w1 * 0.5/v2;
      }
    // done
    naccum = 0;
  }

  // now truncate the functions at the limit radii
  for ( ix = target.first(); !ix.last(); ix.next() )
    if ( ix.coord_orth().lengthsq() > radius * radius )
      weight[ix] = target[ix] = 0.0;

  // now assemble optimised llk lists
  clipper::Coord_grid cg;
  clipper::Coord_orth co;
  clipper::Coord_map cm;
  const clipper::ftype r0 = radius * 3.0 / 8.0;

  // Make list of sites for the approximate LLK fn
  for ( cg.u() = -1; cg.u() <= 1; cg.u()++ )
    for ( cg.v() = -1; cg.v() <= 1; cg.v()++ )
      for ( cg.w() = -1; cg.w() <= 1; cg.w()++ )
	if ( (cg.u()+cg.v()+cg.w())%2 == 0 ) {
	  co = clipper::Coord_orth( r0 * cg.coord_map() );
	  cm = target.coord_map( co );
	  cm = clipper::Coord_map( target.operator_orth_grid() * co );
	  repxyz.push_back( co );
	  reptgt.push_back( target.interp<clipper::Interp_cubic>( cm ) );
	  repwgt.push_back( weight.interp<clipper::Interp_cubic>( cm ) );
	}
  // Make list of sites for the accurate LLK fn
  for ( ix = target.first(); !ix.last(); ix.next() ) {
    cg = ix.coord();
    if ( (cg.u()+cg.v()+cg.w())%2 == 0 ) {  // only use alternate grids
      co = target.coord_orth( cg.coord_map() );
      if ( co.lengthsq() <= radius*radius ) {  // within a sphere
	fullxyz.push_back( co );
	fulltgt.push_back( target[ix] );
	fullwgt.push_back( weight[ix] );
	if ( clipper::Util::isnan(fullwgt.back()) ) std::cout << "Error2 " << ix.coord().format() << "\n";
      }
    }
  }
}

/*! Accumulate the statistics for a log-likelihood target using a known
  map and operator maping into that map. After accumulation is
  completed, you must call prep_llk().
  \param xmap The known map from which to accumulate density statistics.
  \param rtop The operator from the target map at the origin into the xmap. */
void LLK_map_target::accumulate( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop )
{
  // aliases for maps
  clipper::NXmap<float>& mrho = target;
  clipper::NXmap<float>& mrho2 = weight;

  // zero maps if necessary
  if ( naccum == 0.0 ) target = weight = float(0.0);
  naccum++;

  // accumulate stats
  clipper::ftype extentsq =
    pow( mrho.operator_grid_orth().rot()(0,0) * (mrho.grid().nu()-1)/2, 2 );
  clipper::NXmap_base::Map_reference_index ix;
  float rho;
  for ( ix = mrho.first(); !ix.last(); ix.next() )
    if ( ix.coord_orth().lengthsq() <= extentsq ) {
      rho = xmap.interp<clipper::Interp_cubic>( (rtop*ix.coord_orth()).coord_frac(xmap.cell()) );
      mrho[ix]  += rho;
      mrho2[ix] += rho*rho;
    }
}

/*! A log-likelihood FFFear search is performed for the target in the given map.
  \param xmap The map to search.
  \param step The search step in degrees.
  \param nres The number of results to return.
  \return A list of the best orientations and their scores. */
Score_list<clipper::RTop_orth> LLK_map_target::search( const clipper::Xmap<float>& xmap, const clipper::ftype& step, const int& nres, const clipper::ftype& dres ) const
{
  std::vector<clipper::RTop_orth>  ops;
  const clipper::Spacegroup&    spgr = xmap.spacegroup();
  const clipper::Cell&          cell = xmap.cell();
  const clipper::Grid_sampling& grid = xmap.grid_sampling();
  // make a list of rotation ops to try
  float glim = 360.0;  // gamma
  float blim = 180.0;  // beta
  float alim = 360.0;  // alpha
  // reduce search angles by symmetry rotations
  alim /= float( spgr.order_of_symmetry_about_axis( clipper::Spacegroup::C ) );
  if ( spgr.order_of_symmetry_about_axis( clipper::Spacegroup::A ) % 2 == 0 ||
       spgr.order_of_symmetry_about_axis( clipper::Spacegroup::B ) % 2 == 0 )
    blim /= 2.0;
  // do a uniformly sampled search of orientation space
  float anglim = clipper::Util::min( alim, glim );
  for ( float bdeg=step/2; bdeg < 180.0; bdeg += step ) {
    float beta = clipper::Util::d2rad(bdeg);
    float spl = anglim/clipper::Util::intf(cos(0.5*beta)*anglim/step+1);
    float smi = anglim/clipper::Util::intf(sin(0.5*beta)*anglim/step+1);
    for ( float thpl=spl/2; thpl < 720.0; thpl += spl )
      for ( float thmi=smi/2; thmi < 360.0; thmi += smi ) {
	float adeg = clipper::Util::mod(0.5*(thpl+thmi),360.0);
	float gdeg = clipper::Util::mod(0.5*(thpl-thmi),360.0);
	if ( adeg <= alim && bdeg <= blim && gdeg <= glim ) {
	  float alpha = clipper::Util::d2rad(adeg);
	  float gamma = clipper::Util::d2rad(gdeg);
	  clipper::Euler_ccp4 euler( alpha, beta, gamma );
	  ops.push_back(clipper::RTop_orth(clipper::Rotation(euler).matrix()));
	}
      }
  }

  // now search for ML target in each orientation in turn
  clipper::Grid_map tg( cell, grid, radius );
  clipper::NXmap<float> target_rot( cell, grid, tg );
  clipper::NXmap<float> weight_rot( cell, grid, tg );
  clipper::Xmap<float> resultscr( spgr, cell, grid );
  clipper::Xmap<int>   resultrot( spgr, cell, grid );
  clipper::Xmap<int>   resulttrn( spgr, cell, grid );
  clipper::Xmap<float> resultp1 ( clipper::Spacegroup::p1(), cell, grid );
  clipper::Xmap<float>::Map_reference_index i1(resultp1);
  clipper::Xmap<float>::Map_reference_coord ix(resultscr);
  resultscr = 1.0e20;

  // loop over orienatations
  clipper::FFFear_fft<float> srch( xmap );
  for ( int op = 0; op < ops.size(); op++ ) {
    // do the fffear search
    clipper::NX_operator nxop( xmap, target, ops[op].inverse() );
    srch( resultp1, target, weight, nxop );

    // store best scores
    for ( i1 = resultp1.first(); !i1.last(); i1.next() ) {
      ix.set_coord( i1.coord() );
      float score = resultp1[i1];
      if ( score < resultscr[ix] ) {
	resultscr[ix] = score;
	resultrot[ix] = op;
	resulttrn[ix] = grid.index( i1.coord() );
      }
    }
  }

  // now create a long scores list from the maps of results
  Score_list<clipper::RTop_orth> score_list( 5*nres );
  for ( i1 = resultscr.first(); !i1.last(); i1.next() ) {
    float score = resultscr[i1];
    clipper::RTop_orth  rtop = ops[ resultrot[i1] ].inverse();
    clipper::Coord_grid cg = grid.deindex( resulttrn[i1] );
    rtop.trn() = xmap.coord_orth( cg.coord_map() );
    score_list.add( score, rtop );
  }

  // create a pruned scores list omitting near-clashes
  Score_list<clipper::RTop_orth> score_trim( nres );
  for ( int i = 0; i < score_list.size(); i++ ) {
    bool clash = false;
    for ( int j = 0; j < score_trim.size(); j++ )
      if ( clipper::Coord_orth( (score_trim[j].trn()-score_list[i].trn()) ).lengthsq() < dres*dres ) clash = true;
    if ( !clash ) score_trim.add( score_list.score(i), score_list[i] );
  }

  return score_trim;
}

clipper::ftype LLK_map_target::llk_approx( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop ) const
{
  clipper::ftype r( 0.0 ), s( 0.0 );
  for ( int i = 0; i < repxyz.size(); i++ ) {
    r += repwgt[i] * pow( xmap.interp<clipper::Interp_linear>( (rtop*repxyz[i]).coord_frac(xmap.cell()) ) - reptgt[i], 2 );
    s += repwgt[i];
  }
  return r/s;
}

clipper::ftype LLK_map_target::llk       ( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop ) const
{
  clipper::ftype r( 0.0 ), s( 0.0 );
  for ( int i = 0; i < fullxyz.size(); i++ ) {
    r += fullwgt[i] * pow( xmap.interp<clipper::Interp_linear>( (rtop*fullxyz[i]).coord_frac(xmap.cell()) ) - fulltgt[i], 2 );
    s += fullwgt[i];
  }
  return r/s;
}

clipper::String LLK_map_target::format() const
{
  clipper::String result = "";
  clipper::NXmap_base::Map_reference_index ix;
  clipper::ftype y = 0.0;
  for ( ix = target.first(); !ix.last(); ix.next() ) y += pow(target[ix],2);
  y = sqrt( y / double( target.grid().nu() * target.grid().nv() * 
		        target.grid().nw() / 2 ) );
  for ( int w = 0; w < target.grid().nw(); w++ ) {
    result += "\n w=" + clipper::String(w,3);
    for ( int v = 0; v < target.grid().nv(); v++ ) {
      result += "\n  v=" + clipper::String(v,3) + " ";
      for ( int u = 0; u < target.grid().nu(); u++ ) {
	clipper::Coord_grid g( u, v, w );
	float x = target.get_data(g) / y;
	if ( x < 0.0 ) result += "-";
	else if ( x > 1.0 ) result += "#";
	else if ( x > 0.3 ) result += "+";
	else if ( x > 0.1 ) result += ".";
	else result += " ";
      }
    }
  }
  result += "\n";
  return result;
}

/*! The classifier is constructed to classify between a range of targets
  \param rad Radius in which LLK function will be used, in Angstroms.
  \param sampling Sampling spacing for the LLK function, in Angstroms.
  \param num_targets The number of types to be classified. */
LLK_map_classifier::LLK_map_classifier( const clipper::ftype& rad, const clipper::ftype& sampling, const int& num_targets )
{
  LLK_map_target blank_tgt( rad, sampling );
  llktgts.clear();
  llktgts.resize( num_targets, blank_tgt );
}

/*! Accumulate the statistics for a log-likelihood target using a known
  map and operator maping into that map. After accumulation is
  completed, you must call prep_llk().
  \param xmap The known map from which to accumulate density statistics.
  \param rtop The operator from the target map at the origin into the xmap.
  \param target_type The type code of given density example. */
void LLK_map_classifier::accumulate( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop, const int& target_type )
{ llktgts[target_type].accumulate( xmap, rtop ); }

void LLK_map_classifier::prep_llk()
{ for ( int i = 0; i < llktgts.size(); i++ ) llktgts[i].prep_llk(); }

std::vector<clipper::ftype> LLK_map_classifier::llk_raw( const clipper::Xmap<float>& xmap, const clipper::RTop_orth rtop ) const
{
  std::vector<clipper::ftype> result( llktgts.size() );
  for ( int i = 0; i < llktgts.size(); i++ )
    result[i] = llktgts[i].llk( xmap, rtop );
  return result;
}
