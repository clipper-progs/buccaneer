// Clipper buccaneer
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "simulate-lib.h"
#include "buccaneer-find.h"
#include "buccaneer-grow.h"
#include "buccaneer-join.h"
#include "buccaneer-prune.h"
#include "buccaneer-sequence.h"
#include "buccaneer-build.h"


// ramachandran filter data
struct Rama_flt { double phi, psi, rad; };
const Rama_flt rama_flt_all      = {  0.0,  0.0, 10.0 };
const Rama_flt rama_flt_helix    = { -1.5, -1.0,  1.5 };
const Rama_flt rama_flt_strand   = { -2.0, -2.5,  1.5 };
const Rama_flt rama_flt_nonhelix = { -1.5, -1.0, -1.5 };


int main( int argc, char** argv )
{
  CCP4Program prog( "cbuccaneer", "0.4.2", "$Date: 2006/06/16" );

  std::cout << "\nCopyright 2002-2006 Kevin Cowtan and University of York\n";
  std::cout << "All rights reserved. Please reference:\n";
  std::cout << " Cowtan K. (2001) Acta Cryst. D57, 1435-1444.\n";

  // defaults
  clipper::String title;
  clipper::String ipmtz_ref = "NONE";
  clipper::String ipmtz_wrk = "NONE";
  clipper::String ippdb_ref = "NONE";
  clipper::String ippdb_wrk = "NONE";
  clipper::String ipseq_wrk = "NONE";
  clipper::String ipcol_ref_fo = "/*/*/FP";
  clipper::String ipcol_ref_hl = "/*/*/FC";
  clipper::String ipcol_wrk_fo = "NONE";
  clipper::String ipcol_wrk_hl = "NONE";
  clipper::String ipcol_wrk_fc = "NONE";
  clipper::String ipcol_wrk_fr = "NONE";
  clipper::String oppdb = "buccaneer.pdb";
  clipper::String opmap = "NONE";
  clipper::String newresname = "UNK";
  clipper::String newrestype = "ALA";
  double res_in = 1.0;         // Resolution limit
  int nfrag  = 500;
  int nfragr = 20;
  bool find  = false;
  bool grow  = false;
  bool join  = false;
  bool seqnc = false;
  bool prune = false;
  bool build = false;
  double main_tgt_rad = 4.0;
  double side_tgt_rad = 5.5;
  double seq_rel = 0.5;
  double moffset=0.0;
  std::vector<Rama_flt> rama_flt;
  int verbose = 0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-mtzin-ref" ) {
      if ( ++arg < args.size() ) ipmtz_ref = args[arg];
    } else if ( args[arg] == "-mtzin-wrk" ) {
      if ( ++arg < args.size() ) ipmtz_wrk = args[arg];
    } else if ( args[arg] == "-pdbin-ref" ) {
      if ( ++arg < args.size() ) ippdb_ref = args[arg];
    } else if ( args[arg] == "-pdbin-wrk" ) {
      if ( ++arg < args.size() ) ippdb_wrk = args[arg];
    } else if ( args[arg] == "-seqin-wrk" ) {
      if ( ++arg < args.size() ) ipseq_wrk = args[arg];
    } else if ( args[arg] == "-pdbout-wrk" ) {
      if ( ++arg < args.size() ) oppdb = args[arg];
    } else if ( args[arg] == "-mapout" ) {
      if ( ++arg < args.size() ) opmap  = args[arg];
    } else if ( args[arg] == "-colin-ref-fo" ) {
      if ( ++arg < args.size() ) ipcol_ref_fo = args[arg];
    } else if ( args[arg] == "-colin-ref-hl" ) {
      if ( ++arg < args.size() ) ipcol_ref_hl = args[arg];
    } else if ( args[arg] == "-colin-wrk-fo" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fo = args[arg];
    } else if ( args[arg] == "-colin-wrk-hl" ) {
      if ( ++arg < args.size() ) ipcol_wrk_hl = args[arg];
    } else if ( args[arg] == "-colin-wrk-fc" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fc = args[arg];
    } else if ( args[arg] == "-colin-wrk-free" ) {
      if ( ++arg < args.size() ) ipcol_wrk_fr = args[arg];
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) res_in = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-find" ) {
      find = true;
    } else if ( args[arg] == "-grow" ) {
      grow = true;
    } else if ( args[arg] == "-join" ) {
      join = true;
    } else if ( args[arg] == "-sequence" ) {
      seqnc = true;
    } else if ( args[arg] == "-prune" ) {
      prune = true;
    } else if ( args[arg] == "-rebuild" ) {
      build = true;
    } else if ( args[arg] == "-fragments" ) {
      if ( ++arg < args.size() ) nfrag  = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-fragments-per-100-residues" ) {
      if ( ++arg < args.size() ) nfragr = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-ramachandran-filter" ) {
      if ( ++arg < args.size() ) {
	std::vector<clipper::String> words =
	  clipper::String(args[arg]).split( "," );
	for ( int i = 0; i < words.size(); i++ ) {
	  if ( words[i] == "all"      ) rama_flt.push_back( rama_flt_all );
	  if ( words[i] == "helix"    ) rama_flt.push_back( rama_flt_helix );
	  if ( words[i] == "strand"   ) rama_flt.push_back( rama_flt_strand );
	  if ( words[i] == "nonhelix" ) rama_flt.push_back( rama_flt_nonhelix );
	}
      }
    } else if ( args[arg] == "-main-chain-likelihood-radius" ) {
      if ( ++arg < args.size() ) main_tgt_rad = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-side-chain-likelihood-radius" ) {
      if ( ++arg < args.size() ) side_tgt_rad = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-sequence-reliability" ) {
      if ( ++arg < args.size() ) seq_rel = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-new-residue-name" ) {
      if ( ++arg < args.size() ) newresname = args[arg];
    } else if ( args[arg] == "-new-residue-type" ) {
      if ( ++arg < args.size() ) newrestype = args[arg];
    } else if ( args[arg] == "-offset" ) {
      if ( ++arg < args.size() ) moffset = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cbuccaneer\n\t-mtzin-ref <filename>\t\tCOMPULSORY\n\t-pdbin-ref <filename>\t\tCOMPULSORY\n\t-mtzin-wrk <filename>\t\tCOMPULSORY\n\t-pdbin-wrk <filename>\n\t-seqin-wrk <filename>\n\t-pdbout-wrk <filename>\n\t-colin-ref-fo <colpath>\n\t-colin-ref-hl <colpath>\n\t-colin-wrk-fo <colpath>\t\tCOMPULSORY\n\t-colin-wrk-hl <colpath>\t\tCOMPULSORY\n\t-colin-wrk-fc <colpath>\n\t-colin-wrk-free <colpath>\n\t-resolution <resolution/A>\n\t-find\n\t-grow\n\t-join\n\t-sequence\n\t-prune\n\t-rebuild\n\t-fragments <max_fragments>\n\t-fragments-per-100-residues <num_fragments>\n\t-ramachandran-filter <type>\n\t-main-chain-likelihood-radius <radius/A>\n\t-side-chain-likelihood-radius <radius/A>\n\t-sequence-reliability <value>\n\t-new-residue-name <type>\n\t-new-residue-type <type>\nAn input pdb and mtz are required for the reference structure, and \nan input mtz file for the work strructure. Chains will be located and \ngrown for the work structure and written to the output pdb file. \nThis involves 4 steps: finding, growing, joining, and pruning. \nIf the optional input pdb file is provided for the work structure, \nthen finding is skipped and the input chains are grown. \nThe other steps may be disabled using the appropriate options.\n";
    exit(1);
  }

  // other initialisations
  if ( rama_flt.size() == 0 ) rama_flt.push_back( rama_flt_all );
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
 
  // Get resolution for calculation
  mtzfile.open_read( ipmtz_ref );
  double res_ref = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  mtzfile.open_read( ipmtz_wrk );
  double res_wrk = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  mtzfile.close_read();
  resol = clipper::Resolution( clipper::Util::max( res_ref, res_wrk ) );
  if ( res_ref > res_wrk ) std::cout << "\nWARNING: resolution of work structure truncated to reference:\n Ref: " << res_ref << " Wrk: " << res_wrk << std::endl;

  // Get reference reflection data
  clipper::HKL_info hkls_ref;
  mtzfile.open_read( ipmtz_ref );
  hkls_ref.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<clipper::data32::F_sigF> ref_f( hkls_ref );
  clipper::HKL_data<clipper::data32::ABCD> ref_hl( hkls_ref );
  mtzfile.import_hkl_data( ref_f, ipcol_ref_fo );
  mtzfile.import_hkl_data( ref_hl, ipcol_ref_hl );
  mtzfile.close_read();

  // Get work reflection data
  clipper::HKL_info hkls_wrk;
  mtzfile.open_read( ipmtz_wrk );
  hkls_wrk.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<clipper::data32::F_sigF> wrk_f ( hkls_wrk );
  clipper::HKL_data<clipper::data32::ABCD>   wrk_hl( hkls_wrk );
  clipper::HKL_data<clipper::data32::F_phi>  fphi( hkls_wrk );
  clipper::HKL_data<clipper::data32::Flag>   flag( hkls_wrk );
  mtzfile.import_hkl_data( wrk_f , ipcol_wrk_fo );
  mtzfile.import_hkl_data( wrk_hl, ipcol_wrk_hl );
  if ( ipcol_wrk_fc != "NONE" ) mtzfile.import_hkl_data( fphi, ipcol_wrk_fc );
  if ( ipcol_wrk_fr != "NONE" ) mtzfile.import_hkl_data( flag, ipcol_wrk_fr );
  mtzfile.close_read();

  // apply free flag
  clipper::HKL_data<clipper::data32::F_sigF> wrk_f1 = wrk_f;
  wrk_f1.mask( flag != 0 );

  // Get reference model
  clipper::MiniMol mol_ref;
  clipper::MMDBfile mmdb_ref;
  mmdb_ref.read_file( ippdb_ref );
  mmdb_ref.import_minimol( mol_ref );

  // Get work model (optional)
  clipper::MiniMol mol_wrk;
  if ( ippdb_wrk != "NONE" ) {
    clipper::MMDBfile mmdb_wrk;
    mmdb_wrk.read_file( ippdb_wrk );
    mmdb_wrk.import_minimol( mol_wrk );
  }

  // Get work sequence (optional)
  clipper::MMoleculeSequence seq_wrk;
  if ( ipseq_wrk != "NONE" ) {
    clipper::SEQfile seqf_wrk;
    seqf_wrk.read_file( ipseq_wrk );
    seqf_wrk.import_molecule_sequence( seq_wrk );
  }

  // check input files match mode
  if ( !find && mol_wrk.is_null() )
    clipper::Message::message(clipper::Message_fatal("Missing work model."));
  if ( seqnc && seq_wrk.is_null() )
    clipper::Message::message(clipper::Message_fatal("Missing work sequence."));

  // DO INITIAL MAP SIMULATION
  clipper::HKL_data<clipper::data32::F_sigF> sim_f( hkls_ref );
  clipper::HKL_data<clipper::data32::ABCD> sim_hl( hkls_ref );
  MapSimulate mapsim( 100, 20 );
  mapsim( sim_f, sim_hl, ref_f, ref_hl, wrk_f1, wrk_hl );

  // make llk target objects
  LLK_map_classifier llktgt( main_tgt_rad, 0.5, rama_flt.size() );
  LLK_map_classifier llkcls( side_tgt_rad, 0.5, 20 );

  // STAGE 1: Calculate target from reference data

  {
    // reference map
    clipper::HKL_data<clipper::data32::F_phi> fphi( hkls_ref );
    clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls_ref );
    phiw.compute( sim_hl, clipper::data32::Compute_phifom_from_abcd() );
    fphi.compute( sim_f, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( hkls_ref.spacegroup(), hkls_ref.cell(), hkls_ref.resolution() );
    clipper::Xmap<float> xref( hkls_ref.spacegroup(), hkls_ref.cell(), grid );
    xref.fft_from( fphi );

    // prepare llk targets
    typedef clipper::MMonomer MM;
    for ( int chn = 0; chn < mol_ref.size(); chn++ )
      for ( int res = 1; res < mol_ref[chn].size()-1; res++ ) {
	const clipper::MMonomer& mm1 = mol_ref[chn][res-1];
	const clipper::MMonomer& mm2 = mol_ref[chn][res  ];
	const clipper::MMonomer& mm3 = mol_ref[chn][res+1];
	bool b1 = MM::protein_peptide_bond( mm1, mm2 );
	bool b2 = MM::protein_peptide_bond( mm2, mm3 );
	int index_n  = mm2.lookup( " N  ", clipper::MM::ANY );
	int index_ca = mm2.lookup( " CA ", clipper::MM::ANY );
	int index_c  = mm2.lookup( " C  ", clipper::MM::ANY );
	if ( b1 && b2 && index_ca >= 0 && index_c >= 0 && index_n >= 0 ) {
	  Ca_group ca( mm2[index_n ].coord_orth(),
		       mm2[index_ca].coord_orth(),
		       mm2[index_c ].coord_orth() );
	  // main chain target: check residue types and bonding
	  if ( mm2.type() != "GLY" && mm2.type() != "PRO" ) {
	    double phi = MM::protein_ramachandran_phi( mm1, mm2 );
	    double psi = MM::protein_ramachandran_psi( mm2, mm3 );
	    if ( !clipper::Util::is_nan(phi) && !clipper::Util::is_nan(psi) ) {
	      for ( int t = 0; t < rama_flt.size(); t++ ) {
		double r2 = pow(acos(cos(phi-rama_flt[t].phi)),2.0) +
		            pow(acos(cos(psi-rama_flt[t].psi)),2.0);
		if ( r2/rama_flt[t].rad < rama_flt[t].rad )
		  // accumulate target
		  llktgt.accumulate( xref, ca.rtop_from_std_ori(), t );
	      }
	    }
	  }
	  // side chain target:
	  if ( seqnc ) {
	    int type = ProteinTools::residue_index( mm2.type() );
	    if ( type >= 0 )
	      llkcls.accumulate( xref, ca.rtop_beta_carbon(), type );
	  }
	}
      }

    if ( verbose >= 3 )
      for ( int t = 0; t < llktgt.num_targets(); t++ )
	std::cout << "Target:" << t << std::endl
		  << llktgt.llk_map_target(t).format() << std::endl;

    // convert to llk
    llktgt.prep_llk();
    if ( seqnc ) llkcls.prep_llk();
  }


  // STAGE 2: Apply target to work data

  {
    // work map
    if ( ipcol_wrk_fc == "NONE" ) {
      clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls_wrk );
      phiw.compute( wrk_hl, clipper::data32::Compute_phifom_from_abcd() );
      fphi.compute( wrk_f1, phiw, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    }
    clipper::Grid_sampling grid( hkls_wrk.spacegroup(), hkls_wrk.cell(), hkls_wrk.resolution() );
    clipper::Xmap<float> xwrk( hkls_wrk.spacegroup(), hkls_wrk.cell(), grid );
    xwrk.fft_from( fphi );

    // optionlly write work map
    if ( opmap != "NONE" ) {
      clipper::CCP4MAPfile mapfile;
      mapfile.open_write( opmap );
      mapfile.export_xmap( xwrk );
      mapfile.close_write();
    }

    // number of residues to find
    double vol = xwrk.cell().volume() / double(xwrk.spacegroup().num_symops());
    int nres   = int( vol / 320.0 );  // 320A^3/residue on average (inc solvent)
    nfrag = clipper::Util::min( nfrag, (nfragr*nres)/100 );
    clipper::MiniMol mol_all( hkls_wrk.spacegroup(), hkls_wrk.cell() );
    clipper::MiniMol mol_tgt, mol_tmp;

    // offset the map density
    clipper::Map_stats stats( xwrk );
    clipper::Xmap_base::Map_reference_index ix;
    for ( ix = xwrk.first(); !ix.last(); ix.next() )
      xwrk[ix] += moffset * stats.std_dev();

    // tidy input model
    ProteinTools::chain_tidy( mol_tmp, mol_wrk );
    mol_wrk = mol_tmp;

    // now loop over all the chosen target types to find, grow, join
    // and sequence chain fragments
    for ( int t = 0; t < llktgt.num_targets(); t++ ) {
      mol_tgt = mol_wrk;

      // find C-alphas
      if ( find ) {
	Ca_find cafind( nfrag );
	cafind( mol_tgt, xwrk, llktgt.llk_map_target(t) );
	std::cout << "\nC-alphas after finding:    " << mol_tgt.select("*/*/CA").atom_list().size() << std::endl;
      }

      // grow C-alphas
      if ( grow ) {
	Ca_grow cagrow( 25 );
	cagrow( mol_tmp, mol_tgt, xwrk, llktgt.llk_map_target(t) );
	mol_tgt = mol_tmp;
	std::cout << "\nC-alphas after growing:    " << mol_tgt.select("*/*/CA").atom_list().size() << std::endl;
      }
    
      // join C-alphas
      if ( join ) {
	Ca_join cajoin( 2.0 );
	cajoin( mol_tmp, mol_tgt );
	mol_tgt = mol_tmp;
	std::cout << "\nC-alphas after joining:    " << mol_tgt.select("*/*/CA").atom_list().size() << std::endl;
      }
    
      // assign sequences
      if ( seqnc ) {
	Ca_sequence caseq( seq_rel );
	caseq( mol_tmp, mol_tgt, xwrk, llkcls, seq_wrk );
	mol_tgt = mol_tmp;
	if ( verbose >= 1 ) std::cout << caseq.format() << std::endl;
	std::cout << "\nC-alphas sequenced:        " << caseq.num_sequenced() << std::endl;
      }

      // accumulate sequenced fragments
      for ( int chn = 0; chn < mol_tgt.size(); chn++ )
	mol_all.insert( mol_tgt[chn] );
    } // end loop over targets

    // prune C-alphas
    if ( prune ) {
      Ca_prune caprune( 3.0 );
      caprune( mol_tmp, mol_all );
      mol_all = mol_tmp;
      std::cout << "\nC-alphas after pruning:    " << mol_all.select("*/*/CA").atom_list().size() << std::endl;
    }

    // build side chains/atoms
    if ( build ) {
      Ca_build cabuild( newrestype );
      cabuild( mol_tmp, mol_all, xwrk );
      mol_all = mol_tmp;      
      std::cout << "\nC-alphas after rebuilding: " << mol_all.select("*/*/CA").atom_list().size() << std::endl;
    }

    // tidy output model
    ProteinTools::chain_tidy( mol_tmp, mol_all );
    mol_all = mol_tmp;
    
    // adjust residue names
    for ( int c = 0; c < mol_all.size(); c++ )
      for ( int r = 0; r < mol_all[c].size(); r++ )
	if ( mol_all[c][r].type() == "UNK" )
	  mol_all[c][r].set_type( newresname );

    // write answers
    clipper::MMDBfile mmdb;
    mmdb.export_minimol( mol_all );
    mmdb.write_file( oppdb );
  }
}
