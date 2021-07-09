package org.smogserver.scm;

import java.util.*;
import org.smogserver.util.*;
import org.smogserver.io.*;
import org.smogserver.scm.*;
import argparser.*;
import java.io.*;

public class ShadowMain {
	static GroTrajectory grofile;
	static GroTopSkinny topo;
	static PDBfile pdbfile;
	static Atom[] atoms;
	static BondedList list;
	final static int GRO = 1, PDB = 2;
 
	static StringHolder chains;
	static StringHolder top;
	static StringHolder gro;
	static StringHolder type;
	static StringHolder pdb;
	static StringHolder pdbIn;
	static StringHolder outputFile;
	static IntHolder ndec;
	static DoubleHolder atomsizeArg;
	static DoubleHolder cutoffArg;
	static DoubleHolder bondedRadiusArg;
	static IntHolder rnaDeltaArg;
	static IntHolder proteinDeltaArg;
	static BooleanHolder defaultArg;
	static BooleanHolder onlyPDB;
	static BooleanHolder reindex;
	static BooleanHolder distance;
	static String[] compare;
	static String[] differ;
	static BooleanHolder version;
	static BooleanHolder noDihedral;
	static BooleanHolder ignoreHydrogen;
	static StringHolder coarse, bifFile;
	static IntHolder multiple;
	static StringHolder gaussFile;
	static BooleanHolder smogErrors;
	static BooleanHolder smog2output;
	static BooleanHolder showProgress;
	static BooleanHolder volume;
	static BooleanHolder runCorrectedShadow; 
	static BooleanHolder freeformcoords;
	
	static double versionNum = 1.33;
	
	static String distString = ""+
	"*****************************************************************\n"+
    "This software is part of SMOG.  http://smog-server.org\n"+
    "Direct questions to: info@smog-server.org\n"+
    "Work utilizing SMOG should cite:\n\n"+
    "Noel JK, Levi M, Raghunathan M, Lammert H, Hayes R, Onuchic JN, and Whitford PC. (2016)\n"+
	"SMOG V2, A Versatile Software Package for Generating Structure-Based Models.\n"+
	"PLoS Comput Biol 12(3): e1004794. doi:10.1371/journal.pcbi.1004794\n\n"+
    "Work using the Shadow contact map should cite:\n\n"+
    "Noel JK, Whitford PC & Onuchic JN (2012)\n"+
    "The shadow map: a general contact definition for capturing the \n"+
    "dynamics of biomolecular folding and function. \n"+
    "Journal of Physical Chemistry B 116, 8692-8702.\n"+
    "*****************************************************************";
	
	
	
	public static void main(String[] args) {

		//Handle Arguments
		top = new StringHolder();
		gro = new StringHolder();
		type = new StringHolder();
		pdb = new StringHolder();
		pdbIn = new StringHolder();
		chains = new StringHolder();
		outputFile = new StringHolder();
		ndec = new IntHolder();
		atomsizeArg = new DoubleHolder();
		cutoffArg = new DoubleHolder();
		proteinDeltaArg = new IntHolder();
		rnaDeltaArg = new IntHolder();
		bondedRadiusArg = new DoubleHolder();
		defaultArg = new BooleanHolder();
		distance = new BooleanHolder();
		onlyPDB = new BooleanHolder();
		reindex = new BooleanHolder();
		compare = new String[2];
		differ = new String[2];
		version = new BooleanHolder();
		coarse = new StringHolder();
		bifFile = new StringHolder();
		multiple = new IntHolder();
        gaussFile = new StringHolder();
        noDihedral = new BooleanHolder();
        ignoreHydrogen = new BooleanHolder();
        smogErrors = new BooleanHolder();
        showProgress = new BooleanHolder();
        volume = new BooleanHolder();
        smog2output = new BooleanHolder();
        runCorrectedShadow = new BooleanHolder();
		freeformcoords = new BooleanHolder();
        
		ArgParser parser = new ArgParser("java -jar SCM.jar");		
		parser.addOption("-g,--gro %s #REQUIRED#Gromacs .gro file created by SMOG (for coordinates)", gro);
		parser.addOption("-t,--top %s #REQUIRED#Gromacs .top file created by SMOG (for connectivity info)", top);
		parser.addOption("-o,--output %s #REQUIRED#Output file for contacts",outputFile);
		parser.addOption("-m,--map %s {shadow,cutoff} #Choose contact map type", type);
		parser.addOption("--default %v #Specifies to use default parameters (-m shadow -s 1.0 -c 6.0) ", defaultArg);
		parser.addOption("-s,--size %f #Sets the atom size", atomsizeArg);
		parser.addOption("-c,--cutoff %f #Sets the contact cutoff",cutoffArg);
		parser.addOption("-rd,--rnaDelta %d {[0,2]} #Sets |i-j|>rnaDelta for RNA",rnaDeltaArg);
		parser.addOption("-pd,--proteinDelta %d {[0,10]} #Sets |i-j|>proteinDelta for protein",proteinDeltaArg);
		parser.addOption("-br,--bondedRadius %f #Changes the shadow size for bonded atoms.  Default is 0.5 A.",bondedRadiusArg);
		parser.addOption("-ch,--chain %s #Input chain file",chains);
		parser.addOption("-b,--bif,-bif %s #Input bif file",bifFile);
		parser.addOption("--reindex %v #Start each chain at residue 1 and atom 1. (only important if you have ridiculous numbers of atoms)",reindex);
		parser.addOption("-po,--PDBoutput %s #Output file for PDB",pdb);
		parser.addOption("-p,--PDBinput %s #Input pdb file",pdbIn);
		parser.addOption("--distance %v #Specifies that SCM should output contact distances as well.", distance);
		parser.addOption("-og,--outputGauss %s #Specifies that SCM should output Gaussian SMOG [ pairs ] using the input .gro to this file.", gaussFile);
		parser.addOption("--onlyPDB %v #Specifies that SCM should convert .gro to PDB and exit", onlyPDB);
		//parser.addOption("--compare %sX2 #Prints the shared contacts between the files", compare);
		//parser.addOption("--differ %sX2 #Prints the contacts not in the 1st file", differ);
		parser.addOption("--multiple %d #Sets X and loads X snapshots from the specified grofile. Assumes the format of the grofile is that of `trjconv -f traj.xtc -o traj.gro.",multiple);
		parser.addOption("--version %v #Print the version and citation info and quit", version);
		parser.addOption("--coarse %s {CA,AA,AACA,residue,atomic} #Specifies the level of coarse-graining.  "+
				"CA/residue gives residue-residue contacts.  AA/atomic gives atom-atom contacts.  "+
				"AACA give residue-residue contacts with atomic numbering. (AA default)", coarse);
		parser.addOption("--ignoreDihedrals,--noD %v #Doesn't read dihedrals from top, so no exclusion of contacts if dihedral", noDihedral);
		parser.addOption("--ignoreH %v #Removes hydrogens from contact calculations.", ignoreHydrogen);
		parser.addOption("-ndec %d #Gromacs .gro precision. Default is 3.", ndec);
		parser.addOption("--smogErrors %v #Turns on errors and checking for the SMOG server code. Not useful to most users.",smogErrors);
		parser.addOption("--smog2output %v #Customizes output as SMOG2 would like.",smog2output);
		parser.addOption("--showProgress %v #Prints incremental progress if you are running a big system and get worried nothing is happening.",showProgress);
		parser.addOption("--volume %v #Prints the atomic volume.",volume);
		parser.addOption("--correctedShadow %v #uses asin instead of atan.",runCorrectedShadow);
		parser.addOption("--freecoor,-freecoor %v #coordinates in gro file are assumed space delimited floats",freeformcoords);
		
		//some default values;
		ShadowSettings.SHADOW_RADIUS = 1.0;
		ShadowSettings.BONDED_RADIUS = 0.5;
		ShadowSettings.DISTANCE_CUTOFF = 6.0;
		ShadowSettings.PROTEIN_DELTA = 3;
		ShadowSettings.NUCLEIC_DELTA = 0;
		ShadowSettings.USE_SHADOW_MAP = true;
		//set up a way to recognize that these options were given
		rnaDeltaArg.value = -1;
		proteinDeltaArg.value = -1;		
		bondedRadiusArg.value = -1;
		atomsizeArg.value = -1;
		cutoffArg.value = -1;
		coarse.value = "atomic";
		ndec.value = 3;
		//read all the arguments (-1 will be overwritten if the arg exists)
		parser.matchAllArgs(args);

		org.smogserver.scm.ContactMap map = null;

        double millis = System.currentTimeMillis();
		if(version.value == true){ 
			printOpener();
			System.exit(0);
		}
		if(smogErrors.value){ ShadowSettings.SMOG_ERRORS_ON = true;	}
		if(showProgress.value){ ShadowSettings.SHOW_PROGRESS = true; }
		if(ignoreHydrogen.value) { ShadowSettings.IGNORE_HYDROGEN = true; }
		if(smog2output.value) { ShadowSettings.SMOG2_OUTPUT_ON = true; }
		if(runCorrectedShadow.value) { ShadowSettings.LEGACY_SHADOW = false; }
		if(freeformcoords.value) { ShadowSettings.FREE_FORM_COORDINATES = true; }
		ShadowSettings.GRO_PRECISION = ndec.value;
		if(bifFile.value != null) { 
			ShadowSettings.BIF_PARSING = true;
			Residue.parseBif(new FileIO(bifFile.value,FileIO.BUFFERED_READING));
		}
		if(gro.value == null && pdbIn.value == null){ parser.printErrorAndExit("Gromacs .gro is required! Abort."); }
		if(compare[0] != null){ 
		    //compareContacts(); 
		    print2screen("--compare is not implemented in this version. Exiting.");
		    System.exit(1);
		}
		if(differ[0] != null){ 
		    //differContacts(); 
		    print2screen("--differ is not implemented in this version. Exiting.");
    		System.exit(1);
		}
		if(onlyPDB.value == true){ //print PDB only
			if(pdb.value == null){ parser.printErrorAndExit("What is pdb output file? (-p)? Abort."); }
			else {
				grofile = new GroTrajectory(gro.value);
				atoms = grofile.getAtoms();
				print2screen("Found "+atoms.length+" atoms.");
				if(chains.value!=null){
					print2screen("Reading chains from "+chains.value+".");
					grofile.setChains(chains.value);
					print2screen("Found "+grofile.getNumChains()+" chains.");
				}
				FileIO file = new FileIO(pdb.value,FileIO.WRITING);
				if(reindex.value == true) {
					print2screen("Printing reindexed PDB to file: "+pdb.value+".");
					print2screen("Each chain will start with residue 1 and atom 1.");
					PDBfile.writePDBreindex(atoms,file);
				} else {
					print2screen("Printing PDB to file: "+pdb.value+".");
					PDBfile.writePDB(atoms,file);
				}
				file.close();
			}
			System.exit(0);
		}
		if(top.value == null){ parser.printErrorAndExit("Gromacs .top is required! (-t). Abort."); }
		if(defaultArg.value) { } //keep the defaults, -s or -c will override this
		else { 
			if(type.value==null){  parser.printErrorAndExit("Specify --default or -m [type]. Abort."); } 
			if(type.value.compareTo("shadow")==0){ ShadowSettings.USE_SHADOW_MAP = true; ShadowSettings.USE_CUTOFF_MAP = false; }
			else if (type.value.compareTo("cutoff")==0){ ShadowSettings.USE_CUTOFF_MAP = true; ShadowSettings.USE_SHADOW_MAP = false; } 
			else { parser.printErrorAndExit("'-m' Need to specify a valid map type! Abort."); }
		}
		if(proteinDeltaArg.value != -1) { ShadowSettings.PROTEIN_DELTA = proteinDeltaArg.value; } 
		if(rnaDeltaArg.value != -1) { 
			ShadowSettings.NUCLEIC_DELTA = rnaDeltaArg.value; 
			if(ShadowSettings.BIF_PARSING) { print2screen("** Be aware that the -bif option only parses 'amino' residues. "+
				"So -rd will have no effect because only types are amino or other."); }
		} 
		if(atomsizeArg.value != -1) { 
			ShadowSettings.SHADOW_RADIUS = atomsizeArg.value; 
			if(ShadowSettings.USE_CUTOFF_MAP) { print2screen("** WARNING ** You have set shadow radius and are using a cutoff map where shadow radius is unused."); }
		} 
		if(bondedRadiusArg.value != -1) { 
			ShadowSettings.BONDED_RADIUS = bondedRadiusArg.value; 
			if(ShadowSettings.USE_CUTOFF_MAP) { print2screen("** WARNING ** You have set bonded radius and are using a cutoff map where bonded radius is unused."); }
			if(ShadowSettings.BONDED_RADIUS > ShadowSettings.SHADOW_RADIUS) {
				print2screen("** WARNING ** bonded radius is greater than shadowing radius, probably not intended."); 
			}
		}
		if(cutoffArg.value != -1) { ShadowSettings.DISTANCE_CUTOFF = cutoffArg.value; } 
		if(outputFile.value==null){ parser.printErrorAndExit("Please specify output file (-o). Abort."); }
		String outputFileName = outputFile.value;
		FileIO outputFileIO = new FileIO(outputFileName,FileIO.WRITING);
		if(multiple.value == 0){
			if(gro.value != null) initialize(gro.value,ShadowMain.GRO);
			else initialize(pdbIn.value,ShadowMain.PDB);
			if(pdb.value != null){ 
				FileIO file = new FileIO(pdb.value,FileIO.WRITING);
				if(reindex.value == true) {
					print2screen("Printing reindexed PDB to file: "+pdb.value+".");
					print2screen("Each chain will start with residue 1 and atom 1.");
					PDBfile.writePDBreindex(atoms,file);
				} else {
					print2screen("Printing PDB to file: "+pdb.value+".");
					PDBfile.writePDB(atoms,file);
				}
				file.close();
			}
			if(ShadowSettings.IGNORE_HYDROGEN) {
			    atoms = GroGro.removeHydrogen(atoms);
			}
			//calculating the map
			print2screen("Calculating contacts.");
			map = org.smogserver.scm.ContactMap.createContactMap(atoms,list,ShadowSettings.SHADOW_RADIUS,ShadowSettings.BONDED_RADIUS,
				ShadowSettings.DISTANCE_CUTOFF,ShadowSettings.NUCLEIC_DELTA,ShadowSettings.PROTEIN_DELTA);
			if(distance.value) map.setPrintDistance(true);
			if(coarse.value.compareTo("AA")==0 || coarse.value.compareTo("atomic")==0){
				print2screen("Printing contacts to file: "+outputFileName+".");
				if(reindex.value == true) {
					print2screen("Reindexing to correspond with the reindexed PDB.");
					map.print(outputFileIO,grofile.getChains());
				} else {
					map.print(outputFileIO);
				}
			} else if(coarse.value.compareTo("CA")==0 || coarse.value.compareTo("residue")==0){
				print2screen("Printing residue-residue contacts to file: "+outputFileName+".");
				//org.smogserver.scm.ContactMap caMap = map.getResidueContacts();
				if(reindex.value == true) {
					map.printCoarse(ContactMap.CA,outputFileIO,grofile.getChains());
				} else { 
					map.printCoarse(ContactMap.CA,outputFileIO);
				}
			} else if(coarse.value.compareTo("AACA")==0){
				print2screen("Printing residue-residue contacts\n     "+
						"with atom numbering to file: "+outputFileName+".");
				//org.smogserver.scm.ContactMap caMap = map.getResidueContacts();
				if(reindex.value == true) {
					map.printCoarse(ContactMap.AACA,outputFileIO,grofile.getChains());
				} else {
					map.printCoarse(ContactMap.AACA,outputFileIO);
				}
			}
			if(gaussFile.value != null) { //print a Gaussian [ pairs ] section
			    FileIO gaussFileIO = new FileIO(gaussFile.value,FileIO.WRITING);
			    map.printGauss6(gaussFileIO,grofile.getChains());
			}
		} else { //multiple files to read
			initialize(gro.value,ShadowMain.GRO);
			for(int groNum = 1; groNum <= multiple.value; groNum++){
				//String groName = gro.value+"_"+groNum+".gro";
				//if(!(new File(groName)).exists()){ parser.printErrorAndExit("File: \""+groName+"\" does not exist.  Abort."); }
				print2screen("Reading snapshot "+groNum+".");
				grofile.getNext();				
				atoms = grofile.getAtoms();
				if(ShadowSettings.IGNORE_HYDROGEN) {
    			    atoms = GroGro.removeHydrogen(atoms);
    			}
				print2screen("Calculating contacts.");
				map = org.smogserver.scm.ContactMap.createContactMap(atoms,list,ShadowSettings.SHADOW_RADIUS,ShadowSettings.BONDED_RADIUS,
					ShadowSettings.DISTANCE_CUTOFF,ShadowSettings.NUCLEIC_DELTA,ShadowSettings.PROTEIN_DELTA);
				if(distance.value) map.setPrintDistance(true);
				if(pdb.value != null){ 
				    print2screen("Not willing to print multiple PDB files. Remove options \"-p,--PDB\".");
				    System.exit(1);
                    // print2screen("Printing PDB to file: "+pdb.value+".");
                    // FileIO file = new FileIO(pdb.value,FileIO.WRITING);
                    // PDBfile.writePDB(atoms,file);
                    // file.close();
				}
				if(coarse.value.compareTo("AA")==0 || coarse.value.compareTo("atomic")==0){
					print2screen("Printing contacts to file: "+outputFileName+".");
					if(reindex.value == true) {
						print2screen("Reindexing to correspond with the reindexed PDB.");
						map.print(outputFileIO,grofile.getChains());
					} else {
						map.print(outputFileIO);
					}
				} else if(coarse.value.compareTo("CA")==0 || coarse.value.compareTo("residue")==0){
					print2screen("Printing residue-residue contacts to file: "+outputFileName+".");
					//org.smogserver.scm.ContactMap caMap = map.getResidueContacts();
					if(reindex.value == true) {
    					map.printCoarse(ContactMap.CA,outputFileIO,grofile.getChains());
					} else {
    					map.printCoarse(ContactMap.CA,outputFileIO);
					}
				} else if(coarse.value.compareTo("AACA")==0){
					print2screen("Printing residue-residue contacts\n     "+
							"with atom numbering to file: "+outputFileName+".");
					//org.smogserver.scm.ContactMap caMap = map.getResidueContacts();
					if(reindex.value == true) {
    					map.printCoarse(ContactMap.AACA,outputFileIO,grofile.getChains());
					} else {
    					map.printCoarse(ContactMap.AACA,outputFileIO);
					}
				} else {
					print2screen("Printing contacts to file: "+outputFileName+".");
					map.print(outputFileIO);
				}

				outputFileIO.write("END snap "+groNum+"\n");
			}
		}
		print2screen("******* Parameters used for determination of contacts *********");
		if(ShadowSettings.USE_SHADOW_MAP) {
			print2screen("Using type: Shadow");
			print2screen("Atom shadowing size: "+ShadowSettings.SHADOW_RADIUS+" angstrom");
			print2screen("Shadowing size of bonded atom: "+ShadowSettings.BONDED_RADIUS+" angstrom");
		}
		if(ShadowSettings.USE_CUTOFF_MAP) {
			print2screen("Using type: Cutoff");
		}
		print2screen("Cutoff radius: "+ShadowSettings.DISTANCE_CUTOFF+" angstrom");
		print2screen("Minimum sequence distance between protein contacts: |i-j|>"+(ShadowSettings.PROTEIN_DELTA));
		if(!ShadowSettings.SMOG2_OUTPUT_ON) print2screen("Minimum sequence distance between RNA/DNA/LIGAND contacts: |i-j|>"+(ShadowSettings.NUCLEIC_DELTA));
		print2screen("***************************************************************");
		if(!ShadowSettings.SMOG2_OUTPUT_ON) print2screen("Total time: "+((System.currentTimeMillis()-millis)/1000)+" seconds.");
		if(!ShadowSettings.SMOG2_OUTPUT_ON) print2screen(distString);
		if(!ShadowSettings.SMOG2_OUTPUT_ON) print2screen("Finished. \n\n\n");
		else print2screen("Finished with contacts. \n");
	}


	private static void initialize(String fileName, int fileType) {
		if(fileType == ShadowMain.GRO) {
			printOpener();
			
			if(ShadowSettings.BIF_PARSING) print2screen("Using protein residue definitions from supplied .bif.");
			print2screen("Reading grofile: "+fileName+".");
			grofile = new GroTrajectory(fileName);
			atoms = grofile.getAtoms();
			if(!ShadowSettings.SMOG2_OUTPUT_ON) { print2screen("Found "+atoms.length+" atoms."); }
			if(chains.value!=null){
				print2screen("Reading chains from "+chains.value+".");
				grofile.setChains(chains.value);
				if(grofile.getNumSplices() > 0) {
				    print2screen("Found "+grofile.getNumChains()+" chains and "+grofile.getNumSplices()+" splices.");
				} else {
				    print2screen("Found "+grofile.getNumChains()+" chains.");
			    }
			}
			topo = new GroTopSkinny(top.value,grofile,noDihedral.value);
			if( topo.getNumTopAtoms() != atoms.length ) { 
				print2screen("Number of atoms in .top ("+topo.getNumTopAtoms()+") and .gro ("+atoms.length+") disagree. Quitting.");
				System.exit(1);
			}
		} else if(fileType == ShadowMain.PDB) {
			printOpener();
			
			if(ShadowSettings.BIF_PARSING) print2screen("Using protein residue definitions from supplied .bif.");
			print2screen("Reading pdbfile: "+fileName+".");
			//grofile = new GroTrajectory(fileName);
			pdbfile = new PDBfile(fileName);
			atoms = pdbfile.getAtoms();
			print2screen("Found "+atoms.length+" atoms.");
			if(chains.value!=null){
				print2screen("Reading chains from "+chains.value+".");
				grofile.setChains(chains.value);
				if(grofile.getNumSplices() > 0) {
				    print2screen("Found "+grofile.getNumChains()+" chains and "+grofile.getNumSplices()+" splices.");
				} else {
				    print2screen("Found "+grofile.getNumChains()+" chains.");
			    }
			}
			topo = new GroTopSkinny(top.value,pdbfile,noDihedral.value);
			if( topo.getNumTopAtoms() != atoms.length ) { 
				print2screen("Number of atoms in .top ("+topo.getNumTopAtoms()+") and .pdb ("+atoms.length+") disagree. Quitting.");
				System.exit(1);
			}
		}
		if(volume.value) print2screen("Total atomic volume: "+((int)topo.getVolume())+" A^3.");
		if(!ShadowSettings.SMOG2_OUTPUT_ON) { 
			print2screen("Found "+topo.atomtypes+" atom types.");
			print2screen("Found "+topo.conCount+" pairs.");
			print2screen("Found "+topo.bondCount+" bonds.");
			print2screen("Found "+topo.angCount+" angles.");
			print2screen("Found "+topo.dihCount+" dihedrals.");
		}
		list = topo.getBondedList();
	}
	
	static void printOpener() {
		print2screen("");
		print2screen("Shadow Contact Map Java Application (SCM)");
		print2screen("More information at: http://smog-server.org/Shadow.html and the SMOG2 manual");
		print2screen("Version "+versionNum);
	}

	public static void print2screen(String s) {
		if(ShadowSettings.SMOG2_OUTPUT_ON) { System.out.println("\t"+s); }
		else { System.out.println(s); }
	}
	public static void print2screenNoNewline(String s) {
		System.out.print(s);
	}
}

