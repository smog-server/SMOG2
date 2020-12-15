package org.smogserver.scm;

import org.smogserver.io.*;
import org.smogserver.util.math.*;
import org.smogserver.util.*;
import java.util.*;

/**
* Parses a GROMACS .gro file.
*/
public class GroGro extends java.io.File implements Structure {

	String grofile;
	Atom[] atom;
	int numAtom = 0;
	int numChains = 1;
	int numSplices = 0;
	Vector chains;
	ConnectionList splices;

	/**
	 * Gives a list of CA atoms.
	 * @return the list of atoms
	 */
    public Atom[] getCA() {
     int numRes = atom[atom.length-1].getResNum();
     Atom[] ca = new Atom[numRes];
     int count = 0;
     for(int i=0; i<atom.length; i++){
         if(atom[i].getID().equals("CA")){
             ca[count]=atom[i];
             count++;
         }
     }
     if(count < numRes) {
         System.out.println("Number of residues is "+numRes+" but only found "+count+" Calpha atoms.\nERROR: One or more residues is missing a Calpha atom!");
         System.exit(1);
     }
     return ca;
    }
    /**
     * Gives a list of CA atoms.
     * @return the list of atoms
     */
    public static Atom[] getCA(Atom[] atom) {
     int numRes = atom[atom.length-1].getResNum();
     Atom[] ca = new Atom[numRes];
     int count = 0;
     for(int i=0; i<atom.length; i++){
         if(atom[i].getID().equals("CA")){
             ca[count]=atom[i];
             count++;
         }
     }
     if(count < numRes) {
         System.out.println("Number of residues is "+numRes+" but only found "+count+" Calpha atoms.\nERROR: One or more residues is missing a Calpha atom!");
         System.exit(1);
     }
     return ca;
    }
	/**
	 * Gives a list of CA or N1 atoms.  N1 for RNA (P does not exist for first residues in RNA), if residue has both, use CA.
	 * If neither exists, use the first atom in the residue.
	 * @return the list of atoms
	 */
	public static Atom[] getCAorN1(Atom[] atom) {
		Atom[] ca = new Atom[atom[atom.length-1].getResNum()];
		int count = -1;
		int lastRes = -1;
		int countCA = -1;
		boolean isCAadded = false;
		for(int i=0; i<atom.length; i++){
		    if(lastRes != atom[i].getResNum()) {
		        lastRes = atom[i].getResNum();
				//if(count>=0) System.out.println(ca[count]);
				count++;
				//add the first atom in this residue
				ca[count]=atom[i];
				isCAadded = false;
		    }
			if(atom[i].getID().equals("CA")) {
				isCAadded = true;
				countCA++;
				ca[count]=atom[i]; //overwrite with CA
			}
			if(!isCAadded && atom[i].getID().equals("N1")) {
				ca[count]=atom[i];
			}
		}
		if(count != countCA) {
			System.out.println("Warning: You are coarse-graining and not all residues contain a CA atom. This is ok"+
				" as long as you are aware that direct use of a coarse-grained contact map with anything other than"+
					" protein is not recommended.");
		}
		return ca;
	}
	/**
	 * Removes all hydrogens from an array of atoms. Determines if it is hydrogen with isHydrogen routine in Atom class.
	 * @return the list of atoms
	 */
	public static Atom[] removeHydrogen(Atom[] atom) {
		//check how many non-hydrogen atoms there are
		int numNonHydrogen = 0;
		for(int i = 0; i < atom.length; i++) { if(!atom[i].isHydrogen()) numNonHydrogen++; }
		Atom[] nonHydrogen = new Atom[numNonHydrogen];
		numNonHydrogen = 0;
		for(int i = 0; i < atom.length; i++) { 
		    if(!atom[i].isHydrogen()) {
		        nonHydrogen[numNonHydrogen++] = atom[i];
		    }
	    }
	    return nonHydrogen;
	}
	/**
	 * The list of atoms contained in the .gro file.
	 * @return the list of atoms
	 */
	public Atom[] getAtoms() {
		return atom;
	}
	/**
	 * Parses a GROMACS .gro file.  Assumes only one chain. 
	 * @param grofile The gro filename
	 */
	// public GroGro(String grofile) {
	// 	super(grofile);
	//     splices = new ConnectionList(0);
	//     this.grofile = grofile;
	// 	FileIO file = new FileIO(grofile,FileIO.BUFFERED_READING);
	// 	String line = file.readLine();
	// 	line = file.readLine();
	// 	numAtom = Integer.parseInt(line.trim());
	// 	atom = new Atom[numAtom];
	// 	int correction = 0; //for the wrapping of gro atomNumbers due to limited space
	// 	int lastAtomNum = 0;
	// 	try{
	// 		for (int i = 0; i < numAtom; i++) {
	// 			line = file.readLine();
	// 			int atomNum = Integer.parseInt(line.substring(15,20).trim());
	// 			if(lastAtomNum>atomNum){correction += (lastAtomNum+1);}
	// 			lastAtomNum = atomNum;
	// 			String id =line.substring(10,15).trim();
	// 			String resName = line.substring(5,10).trim();
	// 			int resNum = Integer.parseInt(line.substring(0,5).trim());
	// 			String[] coor = FileIO.getTokens(line.substring(20,line.length())," ");
	// 			//convert to angstroms
	// 			double x = 10*Double.parseDouble(coor[0]);
	// 			double y = 10*Double.parseDouble(coor[1]);
	// 			double z = 10*Double.parseDouble(coor[2]);
	// 			Coordinates coords = new Coordinates(x,y,z);
	// 			atom[i] = new Atom(id, atomNum+correction, coords, resName, resNum,i);
	// 			atom[i].setChain(1);
	// 		}
	// 		//make a dummy chain and add it to all the atoms
	// 		Chain chain = new Chain(new Vector(Arrays.asList(atom)),1);
	// 		for(int i = 0; i < atom.length; i++) { //add this chain to all the members
	// 		    atom[i].setChainObject(chain);
	// 		}
	// 	} catch (Exception e) {
	// 		System.out.println("Problem at line: \""+line+"\".\nCheck formatting.  Exiting.");
	// 		System.exit(1);
	// 	}
	// }

	/**
	 * Parses a GROMACS .gro file.  Assumes only one chain, specify precision. 
	 * @param grofile The gro filename
	 */
	public GroGro(String grofile) {
		super(grofile);
	    splices = new ConnectionList(0);
		this.grofile = grofile;
		FileIO file = new FileIO(grofile,FileIO.BUFFERED_READING);
		String line = file.readLine();
		line = file.readLine();
		numAtom = Integer.parseInt(line.trim());
		atom = new Atom[numAtom];
		int prec = ShadowSettings.GRO_PRECISION; //3 is standard
		int correction = 0; //for the wrapping of gro atomNumbers due to limited space 
		int correctionR = 0; //for the wrapping of gro resNumbers due to limited space 
		int lastAtomNum = 0;
		int lastResNum = 0;
		try{
			for (int i = 0; i < numAtom; i++) {
				line = file.readLine();
				int atomNum = Integer.parseInt(line.substring(15,20).trim());
				if(lastAtomNum>atomNum){correction += (lastAtomNum+1);}
				lastAtomNum = atomNum; 			
				String id =line.substring(10,15).trim(); 
				String resName = line.substring(5,10).trim(); 
				int resNum = Integer.parseInt(line.substring(0,5).trim()); 
				if(lastResNum>resNum){correctionR += (lastResNum+1);}
				lastResNum = resNum; 			
				//convert to angstroms
				int coordLength = 5+prec;
				double x = 10*Double.parseDouble(line.substring(20,20+coordLength).trim()); 
				double y = 10*Double.parseDouble(line.substring(20+coordLength,20+coordLength*2).trim()); 
				double z = 10*Double.parseDouble(line.substring(20+coordLength*2,20+coordLength*3).trim()); 
				Coordinates coords = new Coordinates(x,y,z);
				atom[i] = new Atom(id, atomNum+correction, coords, resName, resNum+correctionR,i);
				atom[i].setChain(1);
			}
			//make a dummy chain and add it to all the atoms
			Chain chain = new Chain(new Vector(Arrays.asList(atom)),1);
			for(int i = 0; i < atom.length; i++) { //add this chain to all the members
			    atom[i].setChainObject(chain); 
			}
		} catch (Exception e) {
			System.out.println("Problem at line: \""+line+"\".\nCheck formatting.  Exiting.");
			System.exit(1);
		}
	}
    /**
    * Reads a chain file.  Format is: <code>
    *   <br>[ 1 ] (chain number)
    *   <br>1
    *   <br>2 (atom number)
    *   <br>3
    *   <br>.
    *   <br>.
    *   <br>.
    *   <br>46
    *   <br>[ 2 ]  
    *   <br>47
    *   <br>48
    *   <br>[ 3 ] ; SPLICE 
    *   <br>49
    *   <br>50
    *   <br>.
    *   <br>.
    *   <br>.
    *   </code>
    *   <br> 
    * <code>SPLICE</code> means a connection between protein and nucleic acids.
    * These are important for outputting to a PDB file, replacing <code>TER</code>.  
    * Here all residues in contact with a <code>SPLICE</code> in between are treated like
    * RNA/DNA, i.e. RNA contact delta is used.
    * @param chainfile filename containing chains
    */
	public void setChains(String chainfile) {
		chains = new Vector();
		//read in the chains
		FileIO cfile = new FileIO(chainfile,FileIO.BUFFERED_READING);
		String line = cfile.readLine();
		int chainNum = -1;
		//list of atoms that are part of the chain
		Vector<Atom> newChain = new Vector<Atom>();
		boolean splice = false;		
		try {
			while(line != null) {
				String[] tokens = FileIO.getTokens(line," ");
				if(tokens.length > 1) {
					if(tokens[0].equals("[")) {
						if(chainNum > 0) {
							numChains++;
							//add prior constructed chain
							Chain chainToAdd = new Chain(newChain,chainNum);
							if(splice) {
							    chainToAdd.splice = true;
							    chainToAdd.splicedTo = ((Chain) chains.lastElement()).getChainIndex();
							    splice = false;
							}
							chains.add(chainToAdd);
							Atom[] atoms = chainToAdd.getAtoms();
							for(int i = 0; i < atoms.length; i++) { //add this chain to all the members
							    atoms[i].setChainObject(chainToAdd); 
							}
						}
						newChain = new Vector();
						chainNum = Integer.parseInt(tokens[1]); 
					}
					if(line.indexOf("SPLICE")>0) { //this is a splice!
					    splice = true;
					    numSplices++;
					}
				}
				else if(tokens.length == 1) { 
					int num = Integer.parseInt(tokens[0]);
					if(num > atom.length) {
						System.out.println("Problem at line: \""+line+"\".\n"+num+" is too large an index.  Exiting.");
						System.exit(1);
					}
					int indexToAdd = Integer.parseInt(tokens[0])-1; 
					atom[indexToAdd].setChain(chainNum); 
					newChain.add(atom[indexToAdd]);
				} 
				line = cfile.readLine();
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Problem at line: \""+line+"\".\nCheck formatting.  Exiting.");
		}
		Chain chainToAdd = new Chain(newChain,chainNum);
		if(splice) {
		    chainToAdd.splice = true;
		    chainToAdd.splicedTo = ((Chain) chains.lastElement()).getChainIndex();
		    splice = false;
		}
		chains.add(chainToAdd);
		Atom[] atoms = chainToAdd.getAtoms();
		for(int i = 0; i < atoms.length; i++) { //add this chain to all the members
		    atoms[i].setChainObject(chainToAdd); 
		}
	}
	public Vector getChains() {
		return chains;
    }
	public int getNumChains() {
		return numChains-numSplices;
	}
	public int getNumSplices() {
		return numSplices;
	}
	public int getNumAtoms() {
		return getAtoms().length;
	}
	public static void main(String[] args) {
		/*GroGro gro = new GroGro("rop.gro");
		  Atom[] atom = gro.getAtoms();
		  for (int i = 0; i < atom.length; i++) {
		//System.out.println(atom[i]);
		}
		GroTop top = new GroTop("rop.top",gro);
		 */
	}
	
	/**
    * Parses a GROMACS .gro file. This method will be memory intensive! 
    * @param grofile The gro filename
    */
    /*	public GroGro(String grofile, PDBfile pdb) {
    		super(grofile); 
    		this.grofile = grofile;
    		FileIO file = new FileIO(grofile,FileIO.BUFFERED_READING);
    		String line = file.readLine();
    		line = file.readLine();
    		numAtom = Integer.parseInt(line.trim());
    		atom = new Atom[numAtom];
    		int correction = 0; //for the wrapping of gro atomNumbers due to limited space 
    		try {
    			int lastAtomNum = 0;
    			for (int i = 0; i < numAtom; i++) {
    				line = file.readLine();
    				int atomNum = Integer.parseInt(line.substring(15,20).trim());
    				if(lastAtomNum<atomNum){correction += (lastAtomNum+1);}
    				lastAtomNum = atomNum; 			
    				String id =line.substring(12,15).trim(); 
    				String resName = line.substring(6,12).trim(); 
    				int resNum = Integer.parseInt(line.substring(0,5).trim()); 
    				//convert to angstroms
    				double x = 10*Double.parseDouble(line.substring(20,29).trim()); 
    				double y = 10*Double.parseDouble(line.substring(29,37).trim()); 
    				double z = 10*Double.parseDouble(line.substring(37,44).trim()); 
    				Coordinates coords = new Coordinates(x,y,z);
    				atom[i] = new Atom(id, atomNum+correction, coords, resName, resNum);
    				System.out.println(atom[i]);
    			}
    		}
    		catch (Exception e) {
    			System.out.println(e);
    			file = new FileIO(grofile,FileIO.BUFFERED_READING);
    			line = file.readLine();
    			line = file.readLine();

    			for (int i = 0; i < numAtom; i++) {
    				line = file.readLine();
    				String[] tokens = FileIO.getTokens(line," ");
    				int atomNum = Integer.parseInt(tokens[3]);
    				String id =tokens[2].trim();
    				String resName = tokens[1].trim();
    				int resNum = Integer.parseInt(tokens[0].trim());
    				//convert to angstroms
    				double x = 10*Double.parseDouble(tokens[4].trim());
    				double y = 10*Double.parseDouble(tokens[5].trim());
    				double z = 10*Double.parseDouble(tokens[6].trim());
    				Coordinates coords = new Coordinates(x,y,z);
    				atom[i] = new Atom(id, atomNum, coords, resName, resNum);
    			}	
    		}
    		//assign chain numbers from pdb
    		Atom[] pdbAtom = pdb.getAtoms();
    		for (int i = 0; i < pdbAtom.length; i++) {
    			atom[i].setChain(pdbAtom[i].getChain());
    			//atom[Atom.getAtomNumber(atom,pdbAtom[i].getID(),pdbAtom[i].getResNum())-1].setChain(pdbAtom[i].getChain());
    		}
    		//atom=pdbAtom;
    	}*/
}

