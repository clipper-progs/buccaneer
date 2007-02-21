/*! \file buccaneer-build.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-build.h"

#include <clipper/clipper-contrib.h>


bool Ca_build::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap ) const
{
  typedef clipper::MMonomer Mm;

  // copy molecule
  mol2 = mol1;

  // grow side chains
  for ( int chn = 0; chn < mol2.size(); chn++ ) {
    // build sidechains
    for ( int res = 0; res < mol2[chn].size(); res++ ) {
      // preserve residue name, but build "UNK" as newrestype
      clipper::String oldrestype = mol2[chn][res].type();
      if ( oldrestype == "UNK" ) mol2[chn][res].set_type( newrestype );
      // search for best rotomer
      clipper::MMonomer mm;
      int minr = 0;
      double mins = 0.0;
      for ( int r = 0; r < mm.protein_sidechain_number_of_rotomers(); r++ ) {
        mm = mol2[chn][res];
        mm.protein_sidechain_build_rotomer( r );
        double s = 0.0;
        for ( int atm = 0; atm < mm.size(); atm++ )
          s += xmap.interp<clipper::Interp_cubic>( mm[atm].coord_orth().coord_frac( xmap.cell() ) );
        if ( s > mins ) {
          mins = s;
          minr = r;
        }
      }
      // build best
      mol2[chn][res].protein_sidechain_build_rotomer( minr );
      // restore name
      mol2[chn][res].set_type( oldrestype );
    }

    // build oxygens
    for ( int res = 0; res < mol2[chn].size() - 1; res++ )
      if ( Mm::protein_peptide_bond( mol2[chn][res], mol2[chn][res+1] ) )
        mol2[chn][res].protein_mainchain_build_carbonyl_oxygen( mol2[chn][res+1] );
  }

  // fix clashes
  fix_clashes( mol2, xmap );

  return true;
}




std::vector<Ca_build::Clash> Ca_build::find_clashes( clipper::MiniMol& mol ) const
{
  // now search for any clashes
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();
  Clash clash;
  std::vector<Clash> clashes;
  clipper::MAtomNonBond nb( mol, 8.0 );
  std::vector<clipper::MAtomIndexSymmetry> atoms;
  std::vector<clipper::Coord_orth> catoms;
  clipper::Coord_orth o1, o2;
  clipper::Coord_frac f1, f2;
  for ( int p = 0; p < mol.size(); p++ ) {
    for ( int m = 0; m < mol[p].size(); m++ ) {
      int a = mol[p][m].lookup( " CA ", clipper::MM::ANY );
      if ( a >= 0 ) {
	o1 = mol[p][m][a].coord_orth();
	f1 = o1.coord_frac( cell );
	atoms = nb( o1, 8.0 );
	catoms.resize( atoms.size() );
	for ( int i = 0; i < atoms.size(); i++ ) {
	  o2 = mol[atoms[i].polymer()][atoms[i].monomer()][atoms[i].atom()]
	    .coord_orth();
	  f2 = o2.coord_frac( cell );
	  f2 = spgr.symop(atoms[i].symmetry()) * f2;
	  f2 = f2.lattice_copy_near( f1 );
	  catoms[i] = f2.coord_orth(cell);
	}
	// search over atoms in this residue
	double d2min = 1.0e9;
	int i2min = 0;
	for ( int i1 = 0; i1 < atoms.size(); i1++ )
	  if ( atoms[i1].polymer() == p && atoms[i1].monomer() == m ) {
	    // is native atom is movable?
	    clipper::String id = mol[atoms[i1].polymer()]
	      [atoms[i1].monomer()][atoms[i1].atom()].name();
	    if ( id != " CA " && id != " N  " &&
		 id != " C  " && id != " O  " && id != " CB " ) {
	      // and atoms from elsewhere
	      for ( int i2 = 0; i2 < atoms.size(); i2++ )
		if ( atoms[i2].polymer() != p || atoms[i2].monomer() != m ) {
		  // check for a clash
		  double d2 = (catoms[i1]-catoms[i2]).lengthsq();
		  if ( d2 < d2min ) {
		    d2min = d2;
		    i2min = i2;
		  }
		}
	    }
	  }
	// now check the closest clash
	if ( sqrt(d2min) < 1.6 ) {
	  clash.p1 = p;
	  clash.m1 = m;
	  clash.p2 = atoms[i2min].polymer();
	  clash.m2 = atoms[i2min].monomer();
	  clashes.push_back( clash );
	}
      }
    }
  }
  return clashes;
}


void Ca_build::fix_clash  ( clipper::MMonomer& m1, clipper::MMonomer& m2, const clipper::Xmap<float>& xmap ) const
{
  // useful data
  const clipper::Spacegroup& spgr = xmap.spacegroup();
  const clipper::Cell&       cell = xmap.cell();

  // deal with unknow residues
  clipper::String oldrestype1 = m1.type();
  clipper::String oldrestype2 = m2.type();
  if ( oldrestype1 == "UNK" ) m1.set_type( newrestype );
  if ( oldrestype2 == "UNK" ) m2.set_type( newrestype );

  // set up two monomers
  clipper::MMonomer mm1(m1), mm2(m2);
  // move second to be close to first
  clipper::Coord_orth co( 0.0, 0.0, 0.0 );
  for ( int i1 = 0; i1 < mm1.size(); i1++ ) co += mm1[i1].coord_orth();
  co = (1.0/double(mm1.size())) * co;
  clipper::Coord_frac cf1 = co.coord_frac( cell );
  for ( int i2 = 0; i2 < mm2.size(); i2++ ) {
    clipper::Coord_frac cf2 = mm2[i2].coord_orth().coord_frac( cell );
    cf2 = cf2.symmetry_copy_near( spgr, cell, cf1 );
    mm2[i2].set_coord_orth( cf2.coord_orth( cell ) );
  }
  // now search over orientations
  double smax = -1.0e9;
  int r1max(0), r2max(0);
  int nr1 = mm1.protein_sidechain_number_of_rotomers();
  int nr2 = mm2.protein_sidechain_number_of_rotomers();
  for ( int r1 = 0; r1 < nr1; r1++ ) {
    // build and score first rotomer
    mm1.protein_sidechain_build_rotomer( r1 );
    double s1 = 0.0;
    for ( int a1 = 0; a1 < mm1.size(); a1++ )
      s1 += xmap.interp<clipper::Interp_cubic>( mm1[a1].coord_orth().coord_frac( cell ) );
    s1 /= double( mm1.size() );
    for ( int r2 = 0; r2 < nr2; r2++ ) {
      // build and score second rotomer
      mm2.protein_sidechain_build_rotomer( r2 );
      double s2 = 0.0;
      for ( int a2 = 0; a2 < mm2.size(); a2++ )
	s2 += xmap.interp<clipper::Interp_cubic>( mm2[a2].coord_orth().coord_frac( cell ) );
      s2 /= double( mm2.size() );
      // score for clashes
      double d2min = 1.0e9;
      for ( int a1 = 0; a1 < mm1.size(); a1++ )
	for ( int a2 = 0; a2 < mm2.size(); a2++ ) {
	  double d2 = (mm1[a1].coord_orth()-mm2[a2].coord_orth()).lengthsq();
	  if ( d2 < d2min ) d2min = d2;
	}
      double s = s1 + s2;
      if ( d2min < 1.6 ) s = s - 10.0;
      // keep the best
      if ( s > smax ) {
	smax = s;
	r1max = r1;
	r2max = r2;
      }
    }
  }

  // rebuild
  m1.protein_sidechain_build_rotomer( r1max );
  m2.protein_sidechain_build_rotomer( r2max );
  m1.set_type( oldrestype1 );
  m2.set_type( oldrestype2 );
}


void Ca_build::fix_clashes( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const
{
  // check for clashes
  std::vector<Clash> clashes = find_clashes( mol );

  // and try and fix them
  for ( int i = 0; i < clashes.size(); i++ )
    fix_clash( mol[clashes[i].p1][clashes[i].m1], mol[clashes[i].p2][clashes[i].m2], xmap );
}
