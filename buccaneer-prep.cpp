/*! \file buccaneer-prep.cpp buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-prep.h"


const Ca_prep::Rama_flt Ca_prep::rama_flt_all      = {  0.0,  0.0, 10.0 };
const Ca_prep::Rama_flt Ca_prep::rama_flt_helix    = { -1.5, -1.0,  1.5 };
const Ca_prep::Rama_flt Ca_prep::rama_flt_strand   = { -2.0,  2.5,  1.5 };
const Ca_prep::Rama_flt Ca_prep::rama_flt_nonhelix = { -1.5, -1.0, -1.5 };
int Ca_prep::ncpu = 0;


bool Ca_prep::operator() ( LLK_map_target& llktgt, std::vector<LLK_map_target>& llkcls, const clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const
{
  typedef clipper::NXmap<float>::Map_reference_index MRI;
  const int nmain = 1;
  const int nside = 20;
  std::vector<std::vector<clipper::RTop_orth> > rtops(nmain+nside);
  std::vector<LLK_map_target> targets(nmain+nside);
  const int type_gly = ProteinTools::residue_index_3( "GLY" );
  const int type_pro = ProteinTools::residue_index_3( "PRO" );

  LLK_map_target::TYPE tgttyp =
    correl_ ? LLK_map_target::CORREL : LLK_map_target::NORMAL;
  for ( int t = 0; t < nmain; t++ )
    targets[t      ].init( main_tgt_rad_, 0.5, tgttyp );
  for ( int t = 0; t < nside; t++ )
    targets[t+nmain].init( side_tgt_rad_, 0.5, tgttyp );

  // prepare llk targets
  typedef clipper::MMonomer MM;
  for ( int chn = 0; chn < mol.size(); chn++ )
    for ( int res = 1; res < mol[chn].size()-1; res++ ) {
      const clipper::MMonomer& mm1 = mol[chn][res-1];
      const clipper::MMonomer& mm2 = mol[chn][res  ];
      const clipper::MMonomer& mm3 = mol[chn][res+1];
      bool b1 = MM::protein_peptide_bond( mm1, mm2 );
      bool b2 = MM::protein_peptide_bond( mm2, mm3 );
      int type = ProteinTools::residue_index_3( mm2.type() );
      if ( type >= 0 && type < nside ) {
        Ca_group ca( mm2 );
        if ( b1 && b2 && !ca.is_null() ) {
          // main chain target: check residue types and bonding
          if ( type != type_gly && type != type_pro ) {
            double phi = MM::protein_ramachandran_phi( mm1, mm2 );
            double psi = MM::protein_ramachandran_psi( mm2, mm3 );
            if ( !clipper::Util::is_nan(phi) && !clipper::Util::is_nan(psi) ) {
              double r2 = ( pow(acos(cos(phi-rama_flt_.phi)),2.0) +
                            pow(acos(cos(psi-rama_flt_.psi)),2.0) );
              // negative rad allows condition to be reversed
              if ( r2/rama_flt_.rad < rama_flt_.rad )
                rtops[res%nmain].push_back( ca.rtop_from_std_ori() );
            }
          }
          // side chain target:
          if ( seqnc_ )
            rtops[type+nmain].push_back( ca.rtop_beta_carbon() );
        }
      }
    }

  for ( int t = 0; t < rtops.size(); t++ )
    for ( int r = 0; r < rtops[t].size(); r++ )
      targets[t].accumulate( xmap, rtops[t][r] );
  /*
  Prep_threaded prep( targets, xmap, rtops );
  prep( ncpu );
  targets = prep.result();
  */

  llktgt.init( main_tgt_rad_, 0.5, tgttyp );
  llkcls.resize(nside);
  for ( int t = 0; t < nmain; t++ ) {
    for ( MRI ix = llktgt.llk_target().first(); !ix.last(); ix.next() )
      llktgt.llk_target()[ix] += targets[t].llk_target()[ix];
    for ( MRI ix = llktgt.llk_weight().first(); !ix.last(); ix.next() )
      llktgt.llk_weight()[ix] += targets[t].llk_weight()[ix];
    llktgt.num_samples() += targets[t].num_samples();
  }
  for ( int t = 0; t < nside; t++ ) {
    if ( seqnc_ ) llkcls[t] = targets[t+nmain];
  }

  if ( debug_ )
    std::cout << "Target:" << std::endl << llktgt.format() << std::endl;

  // convert to llk
  llktgt.prep_llk();
  if ( seqnc_ )
    for ( int t = 0; t < nside; t++ ) llkcls[t].prep_llk();

  return true;
}
