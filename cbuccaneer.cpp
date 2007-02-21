// Clipper buccaneer
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "simulate-lib.h"
#include "buccaneer-find.h"
#include "buccaneer-grow.h"
#include "buccaneer-join.h"
#include "buccaneer-link.h"
#include "buccaneer-sequence.h"
#include "buccaneer-correct.h"
#include "buccaneer-prune.h"
#include "buccaneer-filter.h"
#include "buccaneer-build.h"


// ramachandran filter data
struct Rama_flt { double phi, psi, rad; };
const Rama_flt rama_flt_all      = {  0.0,  0.0, 10.0 };
const Rama_flt rama_flt_helix    = { -1.5, -1.0,  1.5 };
const Rama_flt rama_flt_strand   = { -2.0, -2.5,  1.5 };
const Rama_flt rama_flt_nonhelix = { -1.5, -1.0, -1.5 };


int main( int argc, char** argv )
{
  CCP4Program prog( "cbuccaneer", "0.7.1", "$Date: 2006/12/14" );

  std::cout << "\nCopyright 2002-2006 Kevin Cowtan and University of York. All rights reserved.\n\n";
  std::cout << "Please reference:\n Cowtan K. (2006) Acta Cryst. D62, 1002-1011.\n\n";

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
  int ncyc = 1;
  bool find  = false;
  bool grow  = false;
  bool join  = false;
  bool link  = false;
  bool seqnc = false;
  bool corct = false;
  bool prune = false;
  bool filtr = false;
  bool build = false;
  double main_tgt_rad = 4.0;
  double side_tgt_rad = 5.5;
  double seq_rel = 0.5;
  double moffset = 0.0;
  Rama_flt rama_flt = rama_flt_all;
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
    } else if ( args[arg] == "-cycles" ) {
      if ( ++arg < args.size() ) ncyc  = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-find" ) {
      find = true;
    } else if ( args[arg] == "-grow" ) {
      grow = true;
    } else if ( args[arg] == "-join" ) {
      join = true;
    } else if ( args[arg] == "-link" ) {
      link = true;
    } else if ( args[arg] == "-sequence" ) {
      seqnc = true;
    } else if ( args[arg] == "-correct" ) {
      corct = true;
    } else if ( args[arg] == "-filter" ) {
      filtr = true;
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
	if ( args[arg] == "all"      ) rama_flt = rama_flt_all;
	if ( args[arg] == "helix"    ) rama_flt = rama_flt_helix;
	if ( args[arg] == "strand"   ) rama_flt = rama_flt_strand;
	if ( args[arg] == "nonhelix" ) rama_flt = rama_flt_nonhelix;
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
    std::cout << "\nUsage: cbuccaneer\n\t-mtzin-ref <filename>\t\tCOMPULSORY\n\t-pdbin-ref <filename>\t\tCOMPULSORY\n\t-mtzin-wrk <filename>\t\tCOMPULSORY\n\t-pdbin-wrk <filename>\n\t-seqin-wrk <filename>\n\t-pdbout-wrk <filename>\n\t-colin-ref-fo <colpath>\n\t-colin-ref-hl <colpath>\n\t-colin-wrk-fo <colpath>\t\tCOMPULSORY\n\t-colin-wrk-hl <colpath>\t\tCOMPULSORY\n\t-colin-wrk-fc <colpath>\n\t-colin-wrk-free <colpath>\n\t-resolution <resolution/A>\n\t-find\n\t-grow\n\t-join\n\t-link\n\t-sequence\n\t-correct\n\t-filter\n\t-prune\n\t-rebuild\n\t-cycles <num_cycles>\n\t-fragments <max_fragments>\n\t-fragments-per-100-residues <num_fragments>\n\t-ramachandran-filter <type>\n\t-main-chain-likelihood-radius <radius/A>\n\t-side-chain-likelihood-radius <radius/A>\n\t-sequence-reliability <value>\n\t-new-residue-name <type>\n\t-new-residue-type <type>\nAn input pdb and mtz are required for the reference structure, and \nan input mtz file for the work strructure. Chains will be located and \ngrown for the work structure and written to the output pdb file. \nThis involves 6 main steps:\n finding, growing, joining, sequencing, pruning, and rebuilding. \nIf the optional input pdb file is provided for the work structure, \nthen the input chains are grown.\n";
    exit(1);
  }

  // other initialisations
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
  clipper::MiniMol mol_wrk( hkls_wrk.spacegroup(), hkls_wrk.cell() );
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
  LLK_map_target llktgt;
  std::vector<LLK_map_target> llkcls( 20 );
  llktgt.init( main_tgt_rad, 0.5 );
  if (seqnc) for ( int t = 0; t < 20; t++ ) llkcls[t].init( side_tgt_rad, 0.5 );

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
	      double r2 = pow(acos(cos(phi-rama_flt.phi)),2.0) +
		          pow(acos(cos(psi-rama_flt.psi)),2.0);
	      if ( r2/rama_flt.rad < rama_flt.rad )
		// accumulate target
		llktgt.accumulate( xref, ca.rtop_from_std_ori() );
	    }
	  }
	  // side chain target:
	  if ( seqnc ) {
	    int type = ProteinTools::residue_index( mm2.type() );
	    if ( type >= 0 )
	      llkcls[type].accumulate( xref, ca.rtop_beta_carbon() );
	  }
	}
      }

    if ( verbose >= 3 )
      std::cout << "Target:" << std::endl << llktgt.format() << std::endl;

    // convert to llk
    llktgt.prep_llk();
    if ( seqnc )
      for ( int t = 0; t < llkcls.size(); t++ ) llkcls[t].prep_llk();
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

    // and sequence chain fragments
    Ca_find cafind( nfrag );

    // model building loop
    for ( int c = 0; c < ncyc; c++ ) {
      std::cout << "Cycle: " << c+1 << std::endl;
      Ca_sequence::History seq_history;

      // find C-alphas
      if ( find ) {
	cafind( mol_tmp, mol_wrk, xwrk, llktgt );
	mol_wrk = mol_tmp;
	std::cout << "\nC-alphas after finding:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      }

      // grow C-alphas
      if ( grow ) {
	Ca_grow cagrow( 25 );
	cagrow( mol_tmp, mol_wrk, xwrk, llktgt );
	mol_wrk = mol_tmp;
	std::cout << "\nC-alphas after growing:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      }
    
      // join C-alphas
      if ( join ) {
	Ca_join cajoin( 2.0 );
	cajoin( mol_tmp, mol_wrk );
	mol_wrk = mol_tmp;
	std::cout << "\nC-alphas after joining:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      }

      // link C-alphas
      if ( link ) {
	Ca_link calnk( 10.0, 24 );
	calnk( mol_tmp, mol_wrk, xwrk, llktgt );
	mol_wrk = mol_tmp;
	std::cout << "\nC-alphas linked:           " << calnk.num_linked() << std::endl;
      }
    
      // assign sequences
      if ( seqnc ) {
	Ca_sequence caseq( seq_rel );
	caseq( mol_tmp, mol_wrk, xwrk, llkcls, seq_wrk );
	mol_wrk = mol_tmp;
	std::cout << "\nC-alphas sequenced:        " << caseq.num_sequenced() << std::endl;
	seq_history.append( caseq );
      }

      // correct insertions/deletions
      if ( corct ) {
	Ca_correct cacor( 12 );
	cacor( mol_tmp, mol_wrk, xwrk, llkcls, seq_wrk );
	mol_wrk = mol_tmp;
	std::cout << "\nC-alphas corrected:        " << cacor.num_corrected() << std::endl;
      }

      // filter poor density
      if ( filtr ) {
	Ca_filter cafiltr( 1.0 );
	cafiltr( mol_tmp, mol_wrk, xwrk );
	mol_wrk = mol_tmp;
	std::cout << "\nC-alphas after filtering:  " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      }

      // prune C-alphas
      if ( prune ) {
	Ca_prune caprune( 3.0 );
	caprune( mol_tmp, mol_wrk );
	mol_wrk = mol_tmp;
	std::cout << "\nC-alphas after pruning:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      }

      // build side chains/atoms
      if ( build ) {
	Ca_build cabuild( newrestype );
	cabuild( mol_tmp, mol_wrk, xwrk );
	mol_wrk = mol_tmp;
	std::cout << "\nC-alphas after rebuilding: " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
      }

      // tidy output model
      ProteinTools::chain_tidy( mol_tmp, mol_wrk );
      mol_wrk = mol_tmp;
    
      // user output
      if ( verbose >= 0 ) std::cout << std::endl << seq_history.format( mol_wrk );

    } // next cycle

    // adjust residue names
    for ( int c = 0; c < mol_wrk.size(); c++ )
      for ( int r = 0; r < mol_wrk[c].size(); r++ )
	if ( mol_wrk[c][r].type() == "UNK" || mol_wrk[c][r].type() == "???" ||
	     mol_wrk[c][r].type() == "+++" || mol_wrk[c][r].type() == "---" )
	  mol_wrk[c][r].set_type( newresname );

    // adjust residue numbers
    for ( int c = 0; c < mol_wrk.size(); c++ )
      if ( c < 26 ) ProteinTools::chain_renumber( mol_wrk[c], seq_wrk );

    // write answers
    clipper::MMDBfile mmdb;
    mmdb.export_minimol( mol_wrk );
    mmdb.write_file( oppdb );

  }
}
