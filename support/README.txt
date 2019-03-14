
This GISTIC2 release has the following directory structure.

<gistic_install_directory>                    where the archive was unpacked
|-- README.txt                                this file
|-- GISTICDocumentation_standalone.htm        HTML documentation
|-- GISTICDocumentation_standalone_files      documentation images
|   `-- {16 image files}
|-- INSTALL.txt                               how install the compiled program
|-- LICENSE.txt                               licensing information
|-- MATLAB_Compiler_Runtime {empty}           suggested location of MCR installation
|-- MCR_Installer 
|   `-- MCR_R2014a_glnxa64_installer.zip      Matlab MCR v8.3 installer code
|-- example_results {empty}                   output directory for example             
|-- gp_gistic2_from_seg                       compiled GISTIC 2.0 executable
|-- gistic2                                   shell script wrapper that 
|-- run_gistic_example                        shell script to run the example
|-- examplefiles                              example input files
|   |-- arraylistfile.txt
|   |-- cnvfile.txt
|   |-- markersfile.txt
|   `-- segmentationfile.txt
|-- refgenefiles                              various reference genomes
|   |-- hg16.mat
|   |-- hg17.mat
|   |-- hg18.mat
|   |-- hg19.mat                              (hg19 from previous releases)
|   |-- hg19.UCSC.add_miR.140312.refgene.mat  (updated hg19)
|   `-- hg38.UCSC.add_miR.160920.refgene.mat  new! hg38
`-- source                                    MATLAB source code
    |-- {~200 *.m GISTIC 2 source files}      source code for GISTIC 2.0
    `-- @SegArray                             source code for the SegArray class
        `-- {~100 *.m SegArray class files}

This release contains both the compiled executable and the source code from 
which it was compiled, as well as documentation and a runnable example. 

Please read LICENSE.txt to understand the terms under which GISTIC 2.0
is licensed to end users. (The license agreement was amended as of release 
2.0.22.)

The executable is the gene pattern module 'gp_gistic2_from_seg'. This
module depends on the libraries in the MATLAB Compiler Runtime (MCR)
version 8.3 (R2014a). If your Linux system does not already have
the MCR installed, you may do so by running the provided
MCRInstaller.bin from The MathWorks.  The example script,
run_gistic_example,  See INSTALL.txt for details.

The wrapper script 'gistic2' is responsible for setting up the MCR
environment and calling 'gp_gistic2_from_seg'. The script assumes that 
the MCR will be installed in the MATLAB_Compiler_Runtime subdirectory. 
If it is not, you must edit the MCR_ROOT definition in 'gistic2'.

The source code is in the 'source' subdirectory. The
<gistic_install_directory>/source directory must be added to the
matlab path before running or compiling the GISTIC source code. The top-
level module for the executable program provided is gp_gistic2_from_seg.m.

The minimum requirement for MATLAB is R2010b (7.14).  The GISTIC
software was developed on x86-64 Linux systems and has not been tested
on other MATLAB platforms.

HTML documentation for the input parameters and important output files
is available by browsing GISTICDocumentation_standalone.htm. 
The source code is extensively documented and most functions will return 
information in response to a 'help <function>' at the matlab prompt.

This release contains input files to support the run_gistic_example in
the directories examplefiles and refgenefiles.


KNOWN ISSUES WITH THIS RELEASE
------------------------------
Centromeres are marked more than once when using (unsupported) Mus musculus
refebrence genomes.


REVISION HISTORY
----------------
2.00.18546 (2011-03-24)  - initial source code release for Genome Biology publication
2.01.19298 (2011-04-20)  - removed dependencies on external files on the Broad filesystem
2.02.20156 (2011-05-03)  - (internal) added arm level peeloff option for peak resolution
2.029.20718 (2011-05-18) - (internal) fixed data-dependent gene gistic bug
2.03.21133 (2011-05-23)  - (internal) modified SegArray so that large "across the grain" 
                           operations do not take as much memory
2.0.4 (2011-06-08)       - (internal) fixed a SegArray issue with logical indexing, 
      			   SegArray version 1.02, new <maj>.<min>.<bugfix> version numbering
2.0.5 (2011-07-01)       - (internal) fixed remove_cnv function so that Mac-formatted CNV 
      			   files work properly
2.0.6 (2011-07-05)       - (internal) additional patch to remove_cnv to work with the mixed 
                           line delimiter conventions used in the "combined verified" CNV lists
2.0.7 (2011-07-19)       - (internal) output changes: added prefixed "base name"; shortened 
      			   and unified some output names; gene gistic deletion plot is gene 
			   gistic q-value; fully compressed gene tables (one row per gene); 
			   amp/del gene lists are sorted by residual q-values  
2.0.8 (2011-08-03)	 - (internal) smoothing (interpolation) of gene gistic deletion 
                           q-plot; fixes to scales on q-plots; some factoring of code in 
			   gistic_plots
2.0.9 (2011-08-04)       - (internal) fixed issues with plots introduced with version 2.0.8
2.0.10 (2011-09-01)      - (internal) misc. bug fixes: (1) broad length command line argument 
       			   now works in gp_gistic2_from_seg; (2) remove_cnv now works with 
			   pathological mixed-platform input files; (3) run_gistic2_from_seg 
			   throws exception when all the data is eliminated.

2.0.11gp (2011-10-05)	 - GenePattern branch created. Bug fix for detecting by-marker CNV 
	 		   lists in remove_cnv.
2.0.12gp (2011-10-14)    - scores.gistic output reflects interpolated gene gistic deletion
                           plot 
2.0.13gp (2011-12-07)	 - (1) fix dependency on alphabetized-by-symb refgene; (2) consistent 
	 		   q-values in all_lesions output; (3) fix intermittent 'no mtimes for 
			   SegArray' bug in plot_snp_score; (4) write out gistic_inputs.mat 
			   for reproducibility.
2.0.14gp (2011-12-14)	 - (1) fix issues where all q-values are zero for either amps or dels; 
	 		   (2) fix issues where some samples have no events; (3) fix capseg case 
			   with no intergenic markers; (4) less cbs chatter and gistic completion 
			   messages; (5) optimize gene_score_permutations for large number of 
			   samples; (6) fix SegArray fencepost issues; (7) '-armpeel' arm-level
			   peel-off command line option added to gp_gistic2_from_seg executable
2.0.15gp (2012-01-30)    - (1) fixed fencepost bug in gene score interpolation; (2) minor 
	 		   tweaks to gistic_plots
2.0.16gp (2012-05-07)	 - Can read Macintosh generated seg file, broad analysis safe for data 
	 		   sets with few broad events, many miscellaneous minor bug fixes. 
--------------------------
2.0.21 (2104-01-31)	 - End the gene-pattern branch of GISTIC core. This version includes the
       			   following changes incorporated in internal versions:

  2.0.16a (2012-08-17)	 - minor fixes: make_D_from_seg invalid chromosome message error; 
  	  		   identify_peaks_by_arbitration issue with peak at end of chromosome
  2.0.17 (2012-09-11)	 - memory and performance optimization of peak identification code. 
  	 		   SegArray version 1.06 with some Mex files added to improve performance. 
			   Fix gistic_plots chromosome shading for q-value 0. 
  2.0.17a (2012-10-15)	 - fix bug in 2.0.17 where output file "raw_copy_number.pdf" was being named 
  	  		   "[pathname '.pdf']"
  2.0.18 (2013-07-26)	 - Add gene_collapse_method parameter to allow control of marker-to-gene 
  	 		   compression method. Eliminate dependence on ps2pdf function which has 
			   hardwired server paths. Allow asymmetric cap internally.
  2.0.19 (2013-09-18)	 - Minor bug fixes: (1) limit absurdly high CN values to +/- 1e6 in smooth_cbs();
  	 		   (2) fix line numbers reported by make_D_from_seg() when segment shortened 
			   to 0 markers; (3) allow single- and zero-marker arms in find_med_chrarm() 
			   for broad analysis.
  2.0.20 (2013-11-20)	 - Fix hang in peak arbitration when all chromosomes (including Y) are present.
  	 		   Test combinations of res, alpha parameters and maximum copy level that could 
			   overflow memory by using too many bins in scoring.
  2.0.21 (2014-01-31)	 - Cytoband labels for peaks listed in amp_genes and del_genes outputs use 
  	 		   center of maximal segment rather than start of wide peak region (consistent 
			   with cytobands in all_lesions output). D.sdesc can now be a row or a column 
			   cell array.
---------------------------
2.0.22	(2014-07-07)	- Option -gcm (gene_collapse_method) and -scent (sample centering) added.
			- Option -armpeel, which was accidentally disabled in release 2.0.21, has been 
			  reinstated.
			- support for non-human reference genomes.
			- genes with same symbol on different chromosomes given distinct names for
			  gene table generation.
			- new sample_seg_counts.txt output lists sample segment counts and which
			  samples were excluded because they exceded the maximum.
			- input segments that are shortened to zero markers are removed with a warning
			  rather than generating an error message. 
			- updated license agreement
---------------------------
2.0.22-MCR8 (2015-05-07)  Compiled against MCR 8.0 (R1012b) instead of 7.14 (R2010b)
---------------------------
2.0.23  (2017-03-27)    - The markers file input is now optional - if omitted, pseudo-markers will be
			  generated to satisfy GISTIC's input requirements while ensuring reasonably
			  uniform coverage of the genome.
                        - The "broad analysis" of arm-level events has been revised:
 			  (1) arm-level events are now called from a single broad copy number profile 
			  instead of separate amplification and deletion profiles, which had led to 
			  arms counterintuitively called as amplified and deleted on the same sample; 
			  (2) the frequency scores used to determine z-scores and q-values, which excludes
			  arms with the opposite call from the denominator, are now in a column called
			  "frequency score". A new column called "frequncy" gives the intuitive frequency
			  with the denominator inluding arms from all the samples. The analysis results
			  for the same data will be different from that of previous GISTIC versions.
			- Error handling messages have been improved. In particular, many informative
			  error messages were masked by an "Index exceeds matrix dimensions" error 
			  in the exception handler itself.
			- An hg38 reference genome is included with this release.
			- The gp_gistic2_from_seg binary executable is now compiled for MCR 8.3 
			  (Matlab R2014a). The source code is compatible with versions of Matlab up to
			  R2016a, however, the appearance of output graphics may be altered for Matlab
			  versions R2015a and later.
			- This release adds the convenient 'gistic2' wrapper function which sets up
                          the MCR and passes its command line argument to the executable. Scripts have
			  been converted from the C-shell to the Bourne shell.
