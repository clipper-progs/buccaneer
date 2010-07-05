/*! \file buccaneer-tidy.cpp buccaneer library */
/* (C) 2010 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-tidy.h"
#include "buccaneer-prot.h"

#include <algorithm>


bool ModelTidy::tidy( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const clipper::MMoleculeSequence& seq ) const
{
  chain_renumber( mol, seq );
  std::vector<int> sources;
  std::vector<clipper::String> srcids;
  if ( mol_mr.size() > 0 ) {
    sources = assign_chains( mol, mol_mr );
    for ( int c = 0; c < mol_mr.size(); c++ )
      srcids.push_back( mol_mr[c].id() );
  } else {
    sources = assign_chains( mol, rmsd_, nmin_ );
  }
  return update_model( mol, sources, srcids );
}


bool ModelTidy::update_model( clipper::MiniMol& mol, std::vector<int> sources, std::vector<clipper::String> srcids ) const
{
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();

  // current number of chains
  const int nchn = mol.size();
  // ideal number  of chains (sources)
  int nsrc = 0;
  for ( int i = 0; i < sources.size(); i++ )
    nsrc = std::max( nsrc, sources[i]+1 );
  for ( int i = srcids.size(); i < nsrc; i++ )
    srcids.push_back( "" );

  // assemble new chains
  clipper::MiniMol molnew( spgr, cell );
  for ( int c2 = 0; c2 < nsrc; c2++ ) {
    clipper::MPolymer mp;
    mp.set_id( srcids[c2] );

    // globularise around the longest chain
    clipper::MPolymer mptmp;
    for ( int c1 = 0; c1 < nchn; c1++ )
      if ( sources[c1] == c2 ) {
	if ( mptmp.size() > 0 )
	  move_chain( mol[c1], mptmp, spgr, cell );
	for ( int r1 = 0; r1 < mol[c1].size(); r1++ )
	  mptmp.insert( mol[c1][r1] );
      }

    // sort fragments by sequence number
    std::vector<std::pair<int,int> > srcsort;
    for ( int c1 = 0; c1 < nchn; c1++ ) {
      if ( sources[c1] == c2 ) {
	srcsort.push_back(std::pair<int,int>(mol[c1][0].seqnum(),c1));
      }
    }
    std::sort( srcsort.begin(), srcsort.end() );
    for ( int i = 0; i < srcsort.size(); i++ ) {
      int c1 = srcsort[i].second;
      for ( int r1 = 0; r1 < mol[c1].size(); r1++ )
	mp.insert( mol[c1][r1] );
      if ( verbose_ ) std::cout << "Adding chain " << mol[c1].id() << " to chain " << c1 << std::endl;
    }

    // record chain gaps
    std::vector<bool> gaps( mp.size(), false );
    for ( int i = 1; i < gaps.size(); i++ )
      gaps[i] = ! clipper::MMonomer::protein_peptide_bond( mp[i-1], mp[i] );

    // renumber the model, with inscodes if necessary
    const clipper::String inscodes = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const int l = inscodes.length();
    std::vector<int> nums( mp.size() );
    for ( int i = 0; i < nums.size(); i++ ) nums[i] = mp[i].seqnum();
    for ( int i = 0; i < nums.size(); i++ ) nums[i] = nums[i] * l;
    for ( int i = 1; i < nums.size(); i++ )
      if ( gaps[i] ) nums[i] = std::max( nums[i], nums[i-1]-nums[i-1]%l+2*l );
      else           nums[i] = std::max( nums[i], nums[i-1]+1 );
    for ( int i = 0; i < mp.size(); i++ ) {
      const int nres = (nums[i]/l) + 1;
      const int nins = clipper::Util::mod( nums[i], l );
      if ( nins == 0 ) mp[i].set_seqnum( nres );
      else             mp[i].set_seqnum( nres, inscodes.substr( nins, 1 ) );
    }

    if ( mp.size() > 0 ) molnew.insert( mp );
  }

  // assemble unassigned chains
  clipper::MiniMol molunk( spgr, cell );
  for ( int c = 0; c < nchn; c++ )
    if ( sources[c] < 0 ) molunk.insert( mol[c] );
  if ( verbose_ ) std::cout << "Chains -  Old: " << nchn << "  Assigned: " << molnew.size() << "  Unassigned: " << molunk.size() << std::endl;
  // copy chains across
  for ( int c = 0; c < molunk.size(); c++ ) molunk[c].set_id( "" );
  for ( int c = 0; c < molunk.size(); c++ ) molnew.insert( molunk[c] );

  mol = molnew;
  return true;
}


std::vector<int> ModelTidy::assign_chains( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr )
{
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();

  const int nchn = mol.size();
  const int nsrc = mol_mr.size();

  if ( nchn == 0 ) return std::vector<int>();

  // first move the chains onto the MR model
  ProteinTools::symm_match( mol, mol_mr );

  // set up sources
  // source chain map
  const int fnull(-1), fclsh(-2);
  clipper::Resolution reso( 4.0 );
  clipper::Grid_sampling grid( spgr, cell, reso );
  clipper::Xmap<int> xflg( spgr, cell, grid );
  xflg = fnull;
  for ( int c = 0; c < nsrc; c++ )
    label_map( xflg, protein_atoms(mol_mr[c]), 3.0, c, fclsh, fnull );
  //for (clipper::Xmap<int>::Map_reference_index ix=xflg.first();!ix.last();ix.next())std::cout<<xflg[ix];std::cout<<std::endl;

  // count references
  std::vector<int> sources( nchn, -1 );
  for ( int c1 = 0; c1 < nchn; c1++ ) {
    std::vector<int> refcount = count_map( xflg, mol[c1], nsrc );
    int c2 = 0;
    for ( int c = 0; c < nsrc; c++ )
      if ( refcount[c] > refcount[c2] ) c2 = c;
    if ( refcount[c2] > 5 ) sources[c1] = c2;
  }

  // rearrange the chains to match the sources
  std::vector<clipper::String> srcids( nsrc );
  for ( int c2 = 0; c2 < nsrc; c2++ ) srcids[c2] = mol_mr[c2].id();
  return sources;
}


std::vector<int> ModelTidy::assign_chains( clipper::MiniMol& mol, const double rmsd, const int nmin )
{
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();

  const int nchn = mol.size();

  // source chain clash data
  const std::vector<int> num_seq = sequence_count( mol );
  const clipper::Array2d<int> seqflg = sequence_flags( mol );
  const int nseq = seqflg.cols();

  // build matrix of which chains superpose
  clipper::Array2d<int> super( nchn, nchn, true );
  for ( int c1 = 0; c1 < nchn-1; c1++ )
    for ( int c2 = c1 + 1; c2 < nchn; c2++ ) {
      const clipper::RTop_orth rtop =
	ProteinTools::superpose( mol[c2], mol[c1], rmsd, nmin, nmin );
      const bool s = ! rtop.is_null();
      super(c1,c2) = s; super(c2,c1) = s;
    }

  /*
  for ( int c1 = 0; c1 < nchn; c1++ ) {
    for ( int c2 = 0; c2 < nchn; c2++ ) std::cout << super(c1,c2);
    std::cout << std::endl;
  }
  */

  // find the best closed group of ncs related chains
  std::vector<int> used, used_best( nchn, false );
  best_closed_ncs_group( super, num_seq, used, used_best );
  used = used_best;

  // make initial sources list
  std::vector<int> sources( nchn, -1 );
  int nsrc = 0;
  for ( int c = 0; c < nchn; c++ )
    if ( used[c] ) sources[c] = nsrc++;

  // source chain map
  const int fnull(-1), fclsh(-2);
  clipper::Resolution reso( 4.0 );
  clipper::Grid_sampling grid( spgr, cell, reso );
  clipper::Xmap<int> xflg( spgr, cell, grid );
  xflg = fnull;
  for ( int c = 0; c < nchn; c++ )
    if ( sources[c] >= 0 )
      label_map( xflg, protein_atoms(mol[c]), 4.0, sources[c], fclsh, fnull );

  // source chain clash map
  clipper::Array2d<int> srcflg( nsrc, nseq, false );
  for ( int c = 0; c < nchn; c++ )
    if ( sources[c] >= 0 )
      for ( int r = 0; r < nseq; r++ )
	srcflg(sources[c],r) = seqflg(c,r);

  // now expand sources list
  for ( int cyc = 0; cyc < nchn; cyc++ ) {
    std::vector<double> scr1(nchn,0.0), scr2(nchn,0.0), scr(nchn,0.0);
    std::vector<int> src(nchn,-1);

    for ( int c = 0; c < nchn; c++ ) 
      if ( sources[c] < 0 && num_seq[c] > 0 ) {
	std::vector<double> sscr1(nsrc,0.0), sscr2(nsrc,0.0), sscr(nsrc,0.0);
	// score by intimacy to known chains
	std::vector<int> refcount = count_map( xflg, mol[c], nsrc );
	for ( int s = 0; s < nsrc; s++ )
	  sscr1[s] = double(refcount[s]) / double(num_seq[c]);

	// score by overlap with known chains
	std::vector<int> clscount( nsrc, 0 );
	for ( int c1 = 0; c1 < nchn; c1++ )
	  if ( sources[c1] >= 0 )
	    for ( int r = 0; r < nseq; r++ )
	      if ( seqflg(c,r) && seqflg(c1,r) ) clscount[sources[c1]]++;
	for ( int s = 0; s < nsrc; s++ )
	  sscr2[s] = double(clscount[s]) / double(num_seq[c]);

	// total score
	for ( int s = 0; s < nsrc; s++ ) sscr[s] = sscr1[s] - 2.0*sscr2[s];

	// get best score for this chain
	int schn = 0;
	for ( int s = 0; s < nsrc; s++ )
	  if ( sscr[s] > sscr[schn] ) schn = s;
	src[c] = schn;
	scr[c] = sscr[schn];
	scr1[c] = sscr1[schn];
	scr2[c] = sscr2[schn];
      }

    // find best chain to assign
    int c = 0;
    for ( int c1 = 0; c1 < nchn; c1++ ) if ( scr[c1] > scr[c] ) c = c1;

    if ( scr[c] <= 0.0 ) break;

    // add new chain to source
    sources[c] = src[c];
    // update label map
    label_map( xflg, protein_atoms(mol[c]), 4.0, sources[c], fclsh, fnull );
    // update clash map
    for ( int r = 0; r < nseq; r++ )
      srcflg(sources[c],r) = srcflg(sources[c],r) || seqflg(c,r);
  }

  // rearrange the chains to match the sources
  std::vector<clipper::String> srcids( nsrc, "" );
  return sources;
}


bool ModelTidy::chain_renumber( clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq )
{
  for ( int c = 0; c < mol.size(); c++ ) {
    clipper::MPolymer& mp = mol[c];

    // convert sequences to unique strings
    clipper::String chnseq = ProteinTools::chain_sequence( mp );
    std::vector<clipper::String> seqs( seq.size() );
    for ( int chn = 0; chn < seq.size(); chn++ ) {
      clipper::String s = "";
      for ( int res = 0; res < seq[chn].sequence().length(); res++ )
        s += ProteinTools::residue_code_1(ProteinTools::residue_index_1(seq[chn].sequence().substr(res,1)));
      seqs[chn] = s;
    }
  
    // now find best match
    int bestchn = -1;
    int bestscr = 0;
    clipper::MSequenceAlign align( clipper::MSequenceAlign::LOCAL,
  				 1.0, -1000.0, -4.0 );
    std::pair<std::vector<int>,std::vector<int> > result;
    for ( int chn = 0; chn < seqs.size(); chn++ ) {
      const clipper::String& seqseq = seqs[chn];
      result = align( chnseq, seqseq );
      int scr = 0;
      for ( int i = 0; i < result.first.size(); i++ ) {
        if ( result.first[i] >= 0 ) {
  	if ( chnseq[i] == seqseq[result.first[i]] )
  	  scr++;
  	else
  	  scr--;
        }
      }
      if ( scr > bestscr ) {
        bestchn = chn;
        bestscr = scr;
      }
    }
    if ( bestchn < 0 ) return false;
  
    // now number residues
    clipper::String truseq = seqs[bestchn];
    result = align( chnseq, truseq );
    std::vector<int> nums = result.first;
  
    /*
    std::cout << bestchn << " " << bestscr << std::endl;
    for ( int i = 0; i < chnseq.size(); i++ )
      std::cout << chnseq[i];
    std::cout << std::endl;
    for ( int i = 0; i < chnseq.size(); i++ )
      std::cout << (nums[i]>=0&&nums[i]<truseq.length() ? truseq[nums[i]] : '-');
    std::cout << std::endl;
    */
  
    // find bounds of sequenced region
    int i0, i1, i;
    for ( i0 = 0;    i0 < nums.size(); i0++ ) if ( nums[i0] >= 0 ) break;
    for ( i1 = nums.size()-1; i1 >= 0; i1-- ) if ( nums[i1] >= 0 ) break;
    if ( i0 < nums.size() )
      for ( i = i0 - 1; i >= 0;          i-- ) nums[i] = nums[i+1] - 1;
    if ( i1 >= 0 )
      for ( i = i1 + 1; i < nums.size(); i++ ) nums[i] = nums[i-1] + 1;
  
    // renumber the model, with inscodes if necessary
    const clipper::String inscodes = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const int l = inscodes.length();
    for ( i = 0; i < nums.size(); i++ ) nums[i] = nums[i] * l;
    for ( i = 1; i < nums.size(); i++ )
      if ( nums[i] <= nums[i-1] ) nums[i] = nums[i-1]+1;
    for ( i = 0; i < mp.size(); i++ ) {
      const int nres = (nums[i]/l) + 1;
      const int nins = clipper::Util::mod( nums[i], l );
      if ( nins == 0 ) mp[i].set_seqnum( nres );
      else             mp[i].set_seqnum( nres, inscodes.substr( nins, 1 ) );
    }
  
    /*
    for ( int i = 0; i < chnseq.size()-1; i++ )
      if ( nums[i+1] != nums[i]+1 ) std::cout << "! " << mp.id() << " " << nums[i] << " " << nums[i+1] << std::endl;
    */
  }
  return true;
}


bool ModelTidy::label_map( clipper::Xmap<int>& xmap, const std::vector<clipper::Coord_orth>& coords, const double radius, const int label, const int clash, const int null )
{
  const clipper::Cell& cell          = xmap.cell();
  const clipper::Grid_sampling& grid = xmap.grid_sampling();

  clipper::Coord_orth xyz;
  clipper::Coord_grid g0, g1;
  clipper::Grid_range gd( cell, grid, radius );
  clipper::Xmap<int>::Map_reference_coord i0, iu, iv, iw;
  for ( int i = 0; i < coords.size(); i++ ) {
    const clipper::Coord_orth& xyz = coords[i];
    g0 = xmap.coord_map( xyz ).coord_grid() + gd.min();
    g1 = xmap.coord_map( xyz ).coord_grid() + gd.max();
    i0 = clipper::Xmap<int>::Map_reference_coord( xmap, g0 );
    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
	  if ( ( xyz - iw.coord_orth() ).lengthsq() < radius*radius ) {
	    if        ( xmap[iw] == null  ) {
	      xmap[iw] = label;
	    } else if ( xmap[iw] != label ) {
	      xmap[iw] = clash;
	    }
	  }
  }

  return true;
}


std::vector<int> ModelTidy::count_map( const clipper::Xmap<int>& xmap, const clipper::MPolymer& mp, const int nsrc )
{
  const clipper::Cell& cell          = xmap.cell();
  const clipper::Grid_sampling& grid = xmap.grid_sampling();
  
  std::vector<int> refcount( nsrc, 0 );
  for ( int r1 = 0; r1 < mp.size(); r1++ ) {
    if ( ProteinTools::residue_index_3( mp[r1].type() ) >= 0 ) {
      int a1 = mp[r1].lookup( " CA ", clipper::MM::ANY );
      if ( a1 >= 0 ) {
	const clipper::Coord_orth& co = mp[r1][a1].coord_orth();
	const clipper::Coord_grid cg = co.coord_frac(cell).coord_grid(grid);
	int c = xmap.get_data( cg );
	if ( c >= 0 && c < nsrc ) refcount[c] += 1;  // count chain refs
      }
    }
  }
  return refcount;
}


bool ModelTidy::move_chain( clipper::MPolymer& mp1, const clipper::MPolymer& mp0, const clipper::Spacegroup& spgr, const clipper::Cell& cell )
{
  clipper::Coord_orth co;
  double no;

  // get COM of target chain
  co = clipper::Coord_orth( 0.0, 0.0, 0.0 );
  no = 0.0;
  for ( int r1 = 0; r1 < mp0.size(); r1++ )
    for ( int a1 = 0; a1 < mp0[r1].size(); a1++ ) {
      co += mp0[r1][a1].coord_orth();
      no += 1.0;
    }
  clipper::Coord_orth com0 = (1.0/no) * co;

  // get COM of moving chain
  co = clipper::Coord_orth( 0.0, 0.0, 0.0 );
  no = 0.0;
  for ( int r1 = 0; r1 < mp1.size(); r1++ )
    for ( int a1 = 0; a1 < mp1[r1].size(); a1++ ) {
      co += mp1[r1][a1].coord_orth();
      no += 1.0;
    }
  clipper::Coord_orth com1 = (1.0/no) * co;

  // find symop and lattice shift which brings it close to com0
  clipper::Coord_frac cf0 = com0.coord_frac( cell );
  clipper::Coord_frac cf1 = com1.coord_frac( cell );
  double r2min = 1.0e9;
  int symin = 0;
  clipper::Coord_frac cfmin( 0.0, 0.0, 0.0 );
  for ( int s = 0; s < spgr.num_symops(); s++ ) {
    clipper::Coord_frac cf2 = spgr.symop(s) * cf1;
    clipper::Coord_frac cf3 = cf2.lattice_copy_near( cf0 );
    double r2 = (cf3-cf0).lengthsq(cell);
    if ( r2 < r2min ) {
      r2min = r2;
      symin = s;
      cfmin = cf3 - cf2;
    }
  }

  // apply the shifts
  for ( int r1 = 0; r1 < mp1.size(); r1++ )
    for ( int a1 = 0; a1 < mp1[r1].size(); a1++ ) {
      cf1 = mp1[r1][a1].coord_orth().coord_frac( cell );
      cf1 = spgr.symop(symin) * cf1 + cfmin;
      mp1[r1][a1].set_coord_orth( cf1.coord_orth( cell ) );
    }

  return true;
}


std::vector<clipper::Coord_orth> ModelTidy::protein_atoms( const clipper::MPolymer& mp )
{
  std::vector<clipper::Coord_orth> coords;
  for ( int r = 0; r < mp.size(); r++ )
    if ( ProteinTools::residue_index_3( mp[r].type() ) >= 0 )
      for ( int a = 0; a < mp[r].size(); a++ )
	coords.push_back( mp[r][a].coord_orth() );
  return coords;
}


std::vector<int> ModelTidy::sequence_count( const clipper::MiniMol& mol )
{
  std::vector<int> result( mol.size(), 0 );
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ )
      if ( ProteinTools::residue_index_3( mol[c][r].type() ) >= 0 )
	result[c]++;
  return result;
}


clipper::Array2d<int> ModelTidy::sequence_flags( const clipper::MiniMol& mol )
{
  const int nchn = mol.size();

  // make list of used residue numbers
  int smin(1000000), smax(-1000000);
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ )
      if ( ProteinTools::residue_index_3( mol[c][r].type() ) >= 0 ) {
	int s = mol[c][r].seqnum();
	smin = std::min( smin, s );
	smax = std::max( smax, s );
      }
  int nseq = smax - smin + 1;
  if ( nseq <= 0 ) return clipper::Array2d<int>( 0, 0 );

  // make flags of used residue numbers
  clipper::Array2d<int> seqflg( nchn, nseq, false );
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ )
      if ( ProteinTools::residue_index_3( mol[c][r].type() ) >= 0 ) {
	int s = mol[c][r].seqnum();
	seqflg(c,s-smin) = true;
      }

  return seqflg;
}


void ModelTidy::best_closed_ncs_group( const clipper::Array2d<int>& super, const std::vector<int>& num_seq, std::vector<int>& used, std::vector<int>& used_best )
{
  const int nchn = used_best.size();
  const int chn  = used.size();
  if ( chn == nchn ) {
    // finished recursion - test if better solution found
    int n1(0),n2(0),r1(0),r2(0);
    for ( int c = 0; c < nchn; c++ ) {
      if ( used[c] )      { n1++; r1+=num_seq[c]; }
      if ( used_best[c] ) { n2++; r2+=num_seq[c]; }
    }
    if ( n1 > n2 || ( n1 == n2 && r1 > r2 ) ) used_best = used;
    //for ( int c = 0; c < nchn; c++ ) std::cout << used[c];
    //std::cout << " " << n1 << " " << r1 << std::endl;
  } else {
    // check if we can use this chain
    bool can_use = true;
    for ( int c = 0; c < chn; c++ )
      if ( used[c] ) can_use = can_use && super(c,chn);
    if ( can_use ) {
      // if so, try using it
      used.push_back( true );
      best_closed_ncs_group( super, num_seq, used, used_best );
      used.pop_back();
    }
    // try not using it
    used.push_back( false );
    best_closed_ncs_group( super, num_seq, used, used_best );
    used.pop_back();
  }
}
