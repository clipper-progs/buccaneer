// Clipper buccaneer
/* Copyright 2002-2008 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "simulate-lib.h"
#include "buccaneer-prep.h"
#include "buccaneer-find.h"
#include "buccaneer-grow.h"
#include "buccaneer-join.h"
#include "buccaneer-link.h"
#include "buccaneer-sequence.h"
#include "buccaneer-correct.h"
#include "buccaneer-filter.h"
#include "buccaneer-ncsbuild.h"
#include "buccaneer-prune.h"
#include "buccaneer-build.h"
#include "buccaneer-util.h"


int main( int argc, char** argv )
{
  CCP4Program prog( "cbuccaneer", "1.1.9", "$Date: 2008/10/01" );
  prog.set_termination_message( "Failed" );

  std::cout << std::endl << "Copyright 2002-2008 Kevin Cowtan and University of York." << std::endl << std::endl;
  prog.summary_beg();
  std::cout << "$TEXT:Reference: $$ Please reference $$" << std::endl << std::endl << " 'The Buccaneer software for automated model building'" << std::endl << " Cowtan K. (2006) Acta Cryst. D62, 1002-1011." << std::endl << std::endl << "$$";
  prog.summary_end();

  // defaults
  clipper::String title;
  clipper::String ipmtz_ref = "NONE";
  clipper::String ipmtz_wrk = "NONE";
  clipper::String ippdb_ref = "NONE";
  clipper::String ippdb_wrk = "NONE";
  clipper::String ipseq_wrk = "NONE";
  clipper::String ippdb_seq = "NONE";
  clipper::String ipcol_ref_fo = "/*/*/FP";
  clipper::String ipcol_ref_hl = "/*/*/FC";
  clipper::String ipcol_wrk_fo = "NONE";
  clipper::String ipcol_wrk_hl = "NONE";
  clipper::String ipcol_wrk_pw = "NONE";
  clipper::String ipcol_wrk_fc = "NONE";
  clipper::String ipcol_wrk_fr = "NONE";
  clipper::String oppdb = "buccaneer.pdb";
  clipper::String opmap = "NONE";
  clipper::String newresname = "UNK";
  clipper::String newrestype = "ALA";
  double res_in = 1.0;         // Resolution limit
  int nfrag  = 500;
  int nfragr = 20;
  int ncyc = 3;
  int ncpu = 0;
  bool find  = false;  // calculation steps
  bool grow  = false;
  bool join  = false;
  bool link  = false;
  bool seqnc = false;
  bool corct = false;
  bool filtr = false;
  bool ncsbd = false;
  bool prune = false;
  bool build = false;
  bool fast   = false;  // further options
  bool semet  = false;
  bool doanis = false;
  bool optemp = false;
  bool fixpos = false;
  double main_tgt_rad = 4.0;
  double side_tgt_rad = 5.5;
  double seq_rel = 0.95;
  double moffset = 0.0;
  bool correl = false;
  Ca_prep::Rama_flt rama_flt = Ca_prep::rama_flt_all;
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
    } else if ( args[arg] == "-seqin-wrk" ) {
      if ( ++arg < args.size() ) ipseq_wrk = args[arg];
    } else if ( args[arg] == "-pdbin-wrk" ) {
      if ( ++arg < args.size() ) ippdb_wrk = args[arg];
    } else if ( args[arg] == "-pdbin-wrk-sequence-prior" ) {
      if ( ++arg < args.size() ) ippdb_seq = args[arg];
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
    } else if ( args[arg] == "-colin-wrk-phifom" ) {
      if ( ++arg < args.size() ) ipcol_wrk_pw = args[arg];
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
    } else if ( args[arg] == "-ncsbuild" ) {
      ncsbd = true;
    } else if ( args[arg] == "-prune" ) {
      prune = true;
    } else if ( args[arg] == "-rebuild" ) {
      build = true;
    } else if ( args[arg] == "-fast" ) {
      fast = true;
    } else if ( args[arg] == "-build-semet" ) {
      semet = true;
    } else if ( args[arg] == "-anisotropy-correction" ) {
      doanis = true;
    } else if ( args[arg] == "-fix-position" ) {
      fixpos = true;
    } else if ( args[arg] == "-output-intermediate-models" ) {
      optemp = true;
    } else if ( args[arg] == "-fragments" ) {
      if ( ++arg < args.size() ) nfrag  = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-fragments-per-100-residues" ) {
      if ( ++arg < args.size() ) nfragr = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-ramachandran-filter" ) {
      if ( ++arg < args.size() ) {
	if ( args[arg] == "all"      ) rama_flt = Ca_prep::rama_flt_all;
	if ( args[arg] == "helix"    ) rama_flt = Ca_prep::rama_flt_helix;
	if ( args[arg] == "strand"   ) rama_flt = Ca_prep::rama_flt_strand;
	if ( args[arg] == "nonhelix" ) rama_flt = Ca_prep::rama_flt_nonhelix;
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
    } else if ( args[arg] == "-correlation-mode" ) {
      correl = true;
    } else if ( args[arg] == "-jobs" || args[arg] == "-j" ) {
      if ( ++arg < args.size() ) ncpu = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "\nUsage: cbuccaneer\n\t-mtzin-ref <filename>\n\t-pdbin-ref <filename>\n\t-mtzin-wrk <filename>\t\tCOMPULSORY\n\t-seqin-wrk <filename>\n\t-pdbin-wrk <filename>\n\t-pdbin-wrk-sequence-prior <filename>\n\t-pdbout-wrk <filename>\n\t-colin-ref-fo <colpath>\n\t-colin-ref-hl <colpath>\n\t-colin-wrk-fo <colpath>\t\tCOMPULSORY\n\t-colin-wrk-hl <colpath> or -colin-wrk-phifom <colpath>\tCOMPULSORY\n\t-colin-wrk-fc <colpath>\n\t-colin-wrk-free <colpath>\n\t-resolution <resolution/A>\n\t-find\n\t-grow\n\t-join\n\t-link\n\t-sequence\n\t-correct\n\t-filter\n\t-ncsbuild\n\t-prune\n\t-rebuild\n\t-fast\n\t-anisotropy-correction\n\t-build-semet\n\t-fix-position\n\t-cycles <num_cycles>\n\t-fragments <max_fragments>\n\t-fragments-per-100-residues <num_fragments>\n\t-ramachandran-filter <type>\n\t-main-chain-likelihood-radius <radius/A>\n\t-side-chain-likelihood-radius <radius/A>\n\t-sequence-reliability <value>\n\t-new-residue-name <type>\n\t-new-residue-type <type>\n\t-correlation-mode\n\t-jobs <CPUs>\nAn input pdb and mtz are required for the reference structure, and \nan input mtz file for the work structure. Chains will be located and \ngrown for the work structure and written to the output pdb file. \nThis involves 6 main steps:\n finding, growing, joining, sequencing, pruning, and rebuilding. \nIf the optional input pdb file is provided for the work structure, \nthen the input chains are grown.\n";
    exit(1);
  }

  // other initialisations
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  typedef clipper::Xmap_base::Map_reference_index MRI;
  typedef clipper::NXmap_base::Map_reference_index NRI;
  using clipper::data32::Compute_abcd_from_phifom;
  using clipper::data32::Compute_phifom_from_abcd;
  using clipper::data32::Compute_fphi_from_fsigf_phifom;
  using clipper::data32::Compute_scale_u_iso_fsigf;
  using clipper::data32::Compute_scale_u_aniso_fphi;
  using clipper::data32::F_sigF;
  using clipper::data32::F_phi;
  using clipper::data32::Phi_fom;
  using clipper::data32::ABCD;
  using clipper::data32::Flag;
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
  std::string msg;
  const int mmdbflags = MMDBF_IgnoreBlankLines | MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreNonCoorPDBErrors | MMDBF_IgnoreRemarks;
  ProteinTools proteintools = ProteinTools();
  Ca_prep::set_cpus( ncpu );
  Ca_find::set_cpus( ncpu );
  Ca_grow::set_cpus( ncpu );
  Ca_sequence::set_cpus( ncpu );
  Ca_sequence::set_semet( semet );
  Ca_find::TYPE findtype = fast ? Ca_find::SECSTRUC : Ca_find::LIKELIHOOD;
  if ( !(find||grow||join||link||seqnc||corct||filtr||ncsbd||prune||build) )
    find=grow=join=link=seqnc=corct=filtr=ncsbd=prune=build=true;
  if ( ipmtz_ref == "NONE" || ippdb_ref == "NONE" )
    Buccaneer_util::set_reference( ipmtz_ref, ippdb_ref );
  Buccaneer_log log;

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
  clipper::HKL_data<F_sigF> ref_f( hkls_ref );
  clipper::HKL_data<ABCD>   ref_hl( hkls_ref );
  mtzfile.import_hkl_data( ref_f,  ipcol_ref_fo );
  mtzfile.import_hkl_data( ref_hl, ipcol_ref_hl );
  mtzfile.close_read();

  // Get work reflection data
  clipper::MTZcrystal cxtl;
  clipper::HKL_info hkls_wrk;
  mtzfile.open_read( ipmtz_wrk );
  hkls_wrk.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  mtzfile.import_crystal( cxtl, ipcol_wrk_fo+".F_sigF.F" );
  clipper::HKL_data<F_sigF>  wrk_f ( hkls_wrk, cxtl );
  clipper::HKL_data<ABCD>    wrk_hl( hkls_wrk, cxtl );
  clipper::HKL_data<Phi_fom> wrk_pw( hkls_wrk, cxtl );
  clipper::HKL_data<F_phi>   wrk_fp( hkls_wrk, cxtl );
  clipper::HKL_data<Flag>    flag( hkls_wrk, cxtl );
  mtzfile.import_hkl_data( wrk_f , ipcol_wrk_fo );
  if ( ipcol_wrk_hl != "NONE" ) mtzfile.import_hkl_data( wrk_hl,ipcol_wrk_hl );
  if ( ipcol_wrk_pw != "NONE" ) mtzfile.import_hkl_data( wrk_pw,ipcol_wrk_pw );
  if ( ipcol_wrk_fc != "NONE" ) mtzfile.import_hkl_data( wrk_fp,ipcol_wrk_fc );
  if ( ipcol_wrk_fr != "NONE" ) mtzfile.import_hkl_data( flag,  ipcol_wrk_fr );
  mtzfile.close_read();

  // do anisotropy correction
  clipper::U_aniso_orth uaniso( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
  if ( doanis ) {
    // scale obs data
    typedef clipper::SFscale_aniso<float> SFscale;
    SFscale sfscl( 3.0, SFscale::SHARPEN );
    sfscl( wrk_f );
    uaniso = sfscl.u_aniso_orth( SFscale::F );
    // scale map coeffs
    Compute_scale_u_aniso_fphi compute_aniso( 1.0, -uaniso );
    if ( ipcol_wrk_fc != "NONE" ) wrk_fp.compute( wrk_fp, compute_aniso );    
    // output
    std::cout << std::endl << "Applying anisotropy correction:"
	      << std::endl << uaniso.format() << std::endl;
  }

  // apply free flag
  clipper::HKL_data<F_sigF> wrk_fwrk = wrk_f;
  //wrk_fwrk.mask( flag != 0 );
  for ( clipper::HKL_data_base::HKL_reference_index ih = hkls_wrk.first(); !ih.last(); ih.next() ) if ( flag[ih].flag() == 0 ) wrk_fwrk[ih] = F_sigF();  //ugly hack for broken SGI compilers
  // and fill in hl
  if ( ipcol_wrk_hl == "NONE" )
    wrk_hl.compute( wrk_pw, Compute_abcd_from_phifom() );

  // Get reference model
  clipper::MiniMol mol_ref, mol_wrk, mol_seq, mol_wrk_in, mol_tmp;
  clipper::MMDBfile mmdb_ref;
  mmdb_ref.SetFlag( mmdbflags );
  mmdb_ref.read_file( ippdb_ref );
  mmdb_ref.import_minimol( mol_ref );

  // Get work model (optional)
  mol_wrk.init( hkls_wrk.spacegroup(), cxtl );
  if ( ippdb_wrk != "NONE" ) {
    clipper::MMDBfile mmdb_wrk;
    mmdb_wrk.SetFlag( mmdbflags );
    mmdb_wrk.read_file( ippdb_wrk );
    mmdb_wrk.import_minimol( mol_tmp );
    mol_wrk.copy( mol_tmp, clipper::MM::COPY_MPC );
  }
  mol_wrk_in = mol_wrk;
  // Get sequencing model - heavy atoms or MR (optional)
  mol_seq.init( hkls_wrk.spacegroup(), cxtl );
  if ( ippdb_seq != "NONE" ) {
    clipper::MMDBfile mmdb_wrk;
    mmdb_wrk.SetFlag( mmdbflags );
    mmdb_wrk.read_file( ippdb_seq );
    mmdb_wrk.import_minimol( mol_tmp );
    mol_seq.copy( mol_tmp, clipper::MM::COPY_MPC );
    Ca_sequence::set_prior_model( mol_seq );
    log.log( "", mol_seq, verbose>9 );
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
  clipper::HKL_data<F_sigF> sim_f( hkls_ref );
  clipper::HKL_data<ABCD> sim_hl( hkls_ref );
  MapSimulate mapsim( 100, 20 );
  mapsim( sim_f, sim_hl, ref_f, ref_hl, wrk_f, wrk_hl );

  // make llk target objects
  LLK_map_target llktgt;
  std::vector<LLK_map_target> llkcls( 20 );

  // STAGE 1: Calculate target from reference data

  {

    // reference map
    clipper::HKL_data<F_phi>   ref_fp( hkls_ref );
    clipper::HKL_data<Phi_fom> ref_pw( hkls_ref );
    ref_pw.compute( sim_hl, Compute_phifom_from_abcd() );
    ref_fp.compute( sim_f, ref_pw, Compute_fphi_from_fsigf_phifom() );
    clipper::Grid_sampling grid( hkls_ref.spacegroup(), hkls_ref.cell(), hkls_ref.resolution() );
    clipper::Xmap<float> xref( hkls_ref.spacegroup(), hkls_ref.cell(), grid );
    xref.fft_from( ref_fp );

    // prepare llk targets
    Ca_prep caprep( main_tgt_rad, side_tgt_rad, rama_flt, correl, seqnc,
		    verbose>3 );
    caprep( llktgt, llkcls, mol_ref, xref );

    log.log( "PREP" );
  }


  // STAGE 2: Apply target to work data

  {
    // work map
    wrk_pw.compute( wrk_hl, Compute_phifom_from_abcd() );
    if ( ipcol_wrk_fc == "NONE" )
      wrk_fp.compute( wrk_fwrk, wrk_pw, Compute_fphi_from_fsigf_phifom() );
    clipper::Spacegroup cspg = hkls_wrk.spacegroup();
    clipper::Grid_sampling grid( cspg, cxtl, hkls_wrk.resolution() );
    clipper::Xmap<float>   xwrk( cspg, cxtl, grid );
    xwrk.fft_from( wrk_fp );

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

    // offset the map density
    clipper::Map_stats stats( xwrk );
    clipper::Xmap_base::Map_reference_index ix;
    for ( ix = xwrk.first(); !ix.last(); ix.next() )
      xwrk[ix] += moffset * stats.std_dev();

    // generate distribution of llk target values for cutoffs
    llktgt.prep_llk_distribution( xwrk );

    // tidy input model
    ProteinTools::chain_tidy( mol_wrk );

    // and sequence chain fragments
    Ca_find cafind( nfrag, resol.limit() );

    // model building loop
    for ( int cyc = 0; cyc < ncyc; cyc++ ) {
      std::cout << std::endl << "Cycle: " << cyc+1 << std::endl << std::endl;
      clipper::String history = "";

      // find C-alphas by slow likelihood search
      if ( find ) {
	cafind( mol_wrk, xwrk, llktgt, findtype );
	std::cout << " C-alphas after finding:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
	log.log( "FIND", mol_wrk, verbose>9 );
      }

      // grow C-alphas
      if ( grow ) {
	Ca_grow cagrow( 25 );
	cagrow( mol_wrk, xwrk, llktgt );
	std::cout << " C-alphas after growing:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
	log.log( "GROW", mol_wrk, verbose>9 );
      }
    
      // join C-alphas
      if ( join ) {
	Ca_join cajoin( 2.0 );
	cajoin( mol_wrk );
	std::cout << " C-alphas after joining:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
	log.log( "JOIN", mol_wrk, verbose>9 );
      }

      // link C-alphas
      if ( link ) {
	Ca_link calnk( 10.0, 24 );
	calnk( mol_wrk, xwrk, llktgt );
	std::cout << " C-alphas linked:           " << calnk.num_linked() << std::endl;
	log.log( "LINK", mol_wrk, verbose>9 );
      }
    
      // assign sequences
      if ( seqnc ) {
	Ca_sequence caseq( seq_rel );
	caseq( mol_wrk, xwrk, llkcls, seq_wrk );
	std::cout << " C-alphas sequenced:        " << caseq.num_sequenced() << std::endl;
	history = caseq.format();
	log.log( "SEQU", mol_wrk, verbose>9 );
      }

      // correct insertions/deletions
      if ( corct ) {
	Ca_correct cacor( 12 );
	cacor( mol_wrk, xwrk, llkcls, seq_wrk );
	std::cout << " C-alphas corrected:        " << cacor.num_corrected() << std::endl;
	log.log( "CORR", mol_wrk, verbose>9 );
      }

      // filter poor density
      if ( filtr ) {
	Ca_filter cafiltr( 1.0 );
	cafiltr( mol_wrk, xwrk );
	std::cout << " C-alphas after filtering:  " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
	log.log( "FILT", mol_wrk, verbose>9 );
      }

      // ncsbuild C-alphas
      if ( ncsbd ) {
	Ca_ncsbuild cancsbuild( seq_rel, 1.0, 12 );
	cancsbuild( mol_wrk, xwrk, llkcls, seq_wrk );
	std::cout << " C-alphas after NCS build:  " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
	log.log( "NCSB", mol_wrk, verbose>9 );
      }

      // prune C-alphas
      if ( prune ) {
	Ca_prune caprune( 3.0 );
	caprune( mol_wrk );
	std::cout << " C-alphas after pruning:    " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
	log.log( "PRUN", mol_wrk, verbose>9 );
      }

      // build side chains/atoms
      if ( build ) {
	Ca_build cabuild( newrestype );
	cabuild( mol_wrk, xwrk );
	std::cout << " C-alphas after rebuilding: " << mol_wrk.select("*/*/CA").atom_list().size() << std::endl;
	log.log( "REBU", mol_wrk, verbose>9 );
      }

      // tidy output model
      ProteinTools::chain_tidy( mol_wrk );

      // user output
      std::cout << std::endl;
      if ( verbose >= 1 ) std::cout << history << std::endl;
      int nres, nseq, nchn, nmax;
      nchn = mol_wrk.size();
      nres = nseq = nmax = 0;
      for ( int c = 0; c < mol_wrk.size(); c++ ) {
	if ( mol_wrk[c].size() > nmax ) nmax = mol_wrk[c].size();
	for ( int r = 0; r < mol_wrk[c].size(); r++ ) {
	  if ( mol_wrk[c][r].lookup( " CA ", clipper::MM::ANY ) >= 0 ) nres++;
	  if ( ProteinTools::residue_index_3( mol_wrk[c][r].type() ) >= 0 ) nseq++;
	}
      }
      prog.summary_beg();
      std::cout << "Internal cycle " << clipper::String( cyc+1, 3 ) << std::endl;
      msg = ( clipper::String( nres, 5 ) + " residues were built in " +
	      clipper::String( nchn, 2 ) + " chains, the longest having " +
	      clipper::String( nmax, 4 ) + " residues." + "\n" +
	      clipper::String( nseq, 5 ) + " residues were sequenced, " +
	      "after pruning.\n" );
      std::cout << msg << std::endl;
      prog.summary_end();

      // temporary file output
      if ( optemp ) {
	clipper::String c( cyc, 1 );
	clipper::MMDBfile mmdb;
	mmdb.export_minimol( mol_wrk );
	mmdb.write_file( oppdb.substr(0,oppdb.rfind(".")+1) + c + ".pdb" );
      }
    } // next cycle

    // move model
    if ( fixpos )
      ProteinTools::symm_match( mol_wrk, mol_wrk_in );

    // adjust residue names
    if ( newresname != "NONE" )
      for ( int c = 0; c < mol_wrk.size(); c++ )
	for ( int r = 0; r < mol_wrk[c].size(); r++ )
	  if ( ProteinTools::residue_index_3( mol_wrk[c][r].type() ) < 0 )
	    mol_wrk[c][r].set_type( newresname );

    // adjust residue numbers
    for ( int c = 0; c < mol_wrk.size(); c++ )
      if ( c < 26 ) ProteinTools::chain_renumber( mol_wrk[c], seq_wrk );

    // write answers
    clipper::MMDBfile mmdb;
    mmdb.export_minimol( mol_wrk );
    mmdb.write_file( oppdb );
  }

  std::cout << "$TEXT:Result: $$ $$" << std::endl << msg << "$$" << std::endl;
  log.profile();
  prog.set_termination_message( "Normal termination" );
}
