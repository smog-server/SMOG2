package org.smogserver.scm;

import org.smogserver.io.*;
import org.smogserver.util.*;
import java.util.*;

class Residue {
	public static final int NUCLEIC_ACID = 1;
	public static final int PROTEIN = 2;
	public static final int LIGAND = 3;
	public static boolean useRestypeEnum = true; //if parsing a .bif make this false
	//public static thing bifResTypes; 
	
	public String id;
	public int position;
	Atom[] atoms;
	public int type;
	

	public Residue(String id, int position, Atom[] atoms) {
		this.id = id;
		this.type = Residue.getType(id); 
		this.position = position;
		this.atoms = atoms;
	}
	
	private static Hashtable residueTypeHash;
	
	/**
	* Reads a smog2 bif file and creates non-null values for the residue names that
	* have type amino. Otherwise everything defaults to type Residue.LIGAND	
	*/
	public static void parseBif(FileIO bif) {
		if(!ShadowSettings.BIF_PARSING) { 
			 System.out.println("Shadow was asked to parse a bif but you didn't set the parameter -bif."+
				 " Not sure what is going on. Exiting.\n\n");
			 System.exit(1);
		} 
		residueTypeHash = new Hashtable<String,Integer>();
		String line = bif.readLine();
		try {
			//probably should check that second line is <bif> or something, otherwise might read a totally bogus file
			while(line != null) {
				if(line.contains("residue")) { //this is a residue
					if(line.contains("amino")) { //this is a protein residue
						//format
						//   <residue name="ALA" type="amino" atomCount="5">
						// take second token space delimited and grab what is between ""
						String[] tokens = FileIO.getTokens(line," ");
						String name = FileIO.getTokens(tokens[1],"\"")[1];
						residueTypeHash.put(name,Residue.PROTEIN);
						//System.out.println(name);
						//System.out.println(Residue.getType(name));
					}
				}
				line = bif.readLine();
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
		//System.exit(0);
	}
	public static int getType(String resname) {
		if(ShadowSettings.BIF_PARSING){ //return either protein or ligand
			if(residueTypeHash.get(resname) != null) return Residue.PROTEIN; //means it was in the bif
			else return Residue.LIGAND;
		} else {
			return Restype.restype(resname);  //use the native types hard coded in Restype.java
		}
	}
	
	// private int getType(String id) {
// 		//read whatever data structure you use to store the .bif data
// 	}
	
	/* Grabs the resType associated with residue names in the xml .bif 
	* @return errorCode
	*/
/*
	public static int parseBif(String filename) {
		//parse the bif
		//do it stupidly for now
		FileIO bif = new FileIO(filename,FileIO.BUFFERED_READING);
		System.out.println(bif.exists());
		String line = bif.readLine();
		while(line!=null) {
			String[] tokens = FileIO.getTokens(line," <>=\"\t");
			//find <residue ...> tags and grab name and type
			if(tokens.length > 0) {
				if(tokens[0].equals("residue") && tokens[1].equals("name") && tokens[3].equals("type")) { 
					thing.add(tokens[2],translateBifType(tokens[4]));
				}
			}
			line=bif.readLine();
		}
	}
*/
	//assumes amino==PROTEIN nucleic==NUCLEIC_ACID anything else == LIGAND
	//this is because you can set special rules for PROTEIN and NUCLEIC_ACID only at the moment
	//private static String translateBifType(String bifResidueType) {
	//	if(bifResidueType.equals("amino")) return Residue.PROTEIN;
	//	if(bifResidueType.equals("amino")) return Residue.PROT;	}
	
}
