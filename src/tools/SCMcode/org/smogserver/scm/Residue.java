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
			//format
			//   <residue name="ALA" residueType="amino" atomCount="5">
			while(line != null) {
				if(line.contains("<residue ")) { //this is a residue
					//space delimited
					String[] tokens = FileIO.getTokens(line," ");
					String[][] tokens2 = new String[tokens.length][];
					//read in all attributes
					for(int i = 0; i< tokens.length; i++) {
						if(tokens[i].contains("=")) { //this is an attribute
							tokens2[i] = new String[2];
							tokens2[i] = FileIO.getTokens(tokens[i],"=");
						}
					}
					//check for type="amino"
					boolean proteinResidue = false;
					for(int i = 0; i< tokens.length; i++) {
						if(tokens2[i] != null) {
							if(tokens2[i][0].contains("type")) { //this is type attribute
								if(tokens2[i][1].contains("amino")) { //this is a protein residue!
									proteinResidue = true;
								}
							}
						}
					}
					String name = null;
					if(proteinResidue) { //grab the name
						for(int i = 0; i< tokens.length; i++) {
							if(tokens2[i] != null) {
								if(tokens2[i][0].contains("name")) { //this is name attribute
									name = FileIO.getTokens(tokens2[i][1],"\"")[1];
								}
							}
						}
						if(name == null) {
							System.out.println("In "+bif+" there is a residue of type amino with no name attribute");
							System.out.println("while parsing line:");
							System.out.println(line);						
							System.exit(1);
						} else {
							residueTypeHash.put(name,Residue.PROTEIN);
							System.out.println(name);
						}
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
	
}
