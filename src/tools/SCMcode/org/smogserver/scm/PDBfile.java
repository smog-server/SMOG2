package org.smogserver.scm;

import org.smogserver.io.*;
import org.smogserver.util.math.*;
import org.smogserver.util.*;
import java.util.*;
/**
 * Parses a PDB file.
 */
public class PDBfile extends java.io.File implements Structure {

	static String pdbFormat = "ATOM%7d %4s %3s%6d    %8.3f%8.3f%8.3f\n";

	String pdbfile;
	Atom[] atom;
	int numAtom = 0;
	int numChains = 1;
	HashMap chainsH;
	Vector chains;

	/**
	 * Parses a PDB file.  Uses a Vector of atoms initially since we don't know 
	 * the total number in a PDB file.  Doesn't support SPLICES or add chain objects to the atoms.
	 * 
	 * @param pdbfile The PDB filename
	 */
	public PDBfile(String pdbfile) {
		super(pdbfile); 
		this.pdbfile = pdbfile;
		FileIO file = new FileIO(pdbfile,FileIO.BUFFERED_READING);
		String line = file.readLine();
		Vector atomTmp = new Vector(1000); //temp storage for all the atoms, will make an array at end
		this.chains = new Vector(1); //store the chains
	//	HashMap chainsH = new HashMap(1);
		Vector newChain = new Vector(1000); //stores atoms for each chain
		int correction = 0; //for the wrapping of gro atomNumbers due to limited space 
		int lastAtomNum = 0;
        while (line != null) {
            //chain code
            if (line.indexOf("TER") != -1 || line.indexOf("END") != -1) {
				//add prior constructed chain of atoms
				System.out.println(numChains);
				chains.add(new Chain(newChain,numChains));
				numChains++;
				newChain = new Vector();
			}
            //atom code
			else if (line.indexOf("ATOM") != -1 ) {
			    int atomNum = Integer.parseInt(line.substring(6,11).trim());
				if(lastAtomNum==99999){correction += (lastAtomNum+1); System.out.println(lastAtomNum+" "+atomNum+" "+correction);}
				lastAtomNum = atomNum;
				String id =line.substring(12,16).trim(); 
				String resName = line.substring(17,20).trim(); 
				int resNum = Integer.parseInt(line.substring(22,26).trim()); 
				double x = Double.parseDouble(line.substring(30,38).trim()); 
				double y = Double.parseDouble(line.substring(38,46).trim()); 
				double z = Double.parseDouble(line.substring(46,54).trim()); 
				Coordinates coords = new Coordinates(x,y,z);
				Atom newAtom = new Atom(id, atomNum+correction, coords, resName, resNum);
				newAtom.setChain(numChains);
				atomTmp.add(newAtom);
				newChain.add(newAtom);
				numAtom++;
				//System.out.println(coords);
			} else {
			    System.out.println("Unknown identifier.  Line must begin with 'ATOM', 'TER', or 'END'");
			    System.out.println(line);
			    System.exit(0);
			}
			line = file.readLine();
    		//System.exit(0);
        }
        atom = new Atom[atomTmp.size()];
        for (int i=0; i<atomTmp.size(); i++) {
            atom[i] = (Atom) atomTmp.elementAt(i);
        }
	}
	/**
	 * The list of chains contained in the PDB file.
	 * @return the list of chains
	 */
	public Vector getChains() {
		return chains;
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
		    System.out.println("in set chains");
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
	public int getNumChains() {
	    return chains.size();
	}
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
				ca[count++]=atom[i];
			}
		}
		return ca;
	}
	/**
	 * The list of atoms contained in the PDB file.
	 * @return the list of atoms
	 */

	public Atom[] getAtoms() {
		return atom;
	}

	/**
	 * Outputs in PDB format a set of atoms and reindexes
	 * @param atoms the atom list
	 * @param file the FileIO.WRITING object to send PDB file
	 */
	public static void writePDBreindex(Atom[] atoms, FileIO file) {
		int chain = atoms[0].getChain();
		int residueCorrection = 0;
		for(int i = 0; i < atoms.length; i++) {
			//chain code
			if(chain != atoms[i].getChain()) {
				chain = atoms[i].getChain();
			    //System.out.println("hello "+chain);
			    if(atoms[i].getChainObject().splice) file.write("SPLICE\n");
				else file.write("TER\n");
				residueCorrection = atoms[i].getResNum()-1;
			}
			//the '99999' is because PDB format is fixed width and
			//only allows for 5 digits for atom and residue numbers
			//System.out.println(((Chain)chains.elementAt(chain-1)).firstAtom+1);
//			System.out.println((Chain)chainsH.get(new Integer(chain)));
			file.write(String.format(pdbFormat,
					((atoms[i].getPosition()-atoms[i].getChainObject().firstAtom+1)%99999 + 99999) % 99999,
					atoms[i].getPDBname(),
					atoms[i].getResName(),
					((atoms[i].getResNum()-atoms[i].getChainObject().firstRes+1)%99999 + 99999) % 99999,
					atoms[i].getX(),
					atoms[i].getY(),
					atoms[i].getZ()));
		}
		file.write("END\n");
	}
	
	/**
	 * Outputs in PDB format a set of atoms.
	 * @param atoms the atom list
	 * @param file the FileIO.WRITING object to send PDB file
	 */
	public static void writePDB(Atom[] atoms, FileIO file) {
		int chain = atoms[0].getChain();
		int residueCorrection = 0;
		for(int i = 0; i < atoms.length; i++) {
			//chain code
			if(chain != atoms[i].getChain()) {
				chain = atoms[i].getChain();
			    //System.out.println("hello "+chain);
			    if(atoms[i].getChainObject().splice) file.write("SPLICE\n");
				else file.write("TER\n");
				residueCorrection = atoms[i].getResNum()-1;
			}
			file.write(String.format(pdbFormat,
					atoms[i].getPosition(),
					atoms[i].getPDBname(),
					atoms[i].getResName(),
					atoms[i].getResNum(),
					atoms[i].getX(),
					atoms[i].getY(),
					atoms[i].getZ()));
		}
		file.write("END\n");
	}
	public static void main(String args[]){
		//PDBfile file = new PDBfile("test.pdb");
		PDBfile file = new PDBfile("/Users/jknoel/jknoel/Ribosome/ribtest/rib.pdb");
		Atom[] atoms = file.getAtoms();
		//PDBfile.writePDB(atoms,file.getChains(),new FileIO("testWrite.pdb",FileIO.WRITING));
		//PDBfile.writePDB(atoms,new FileIO("testWrite.pdb",FileIO.WRITING));
		//GroGro gro = new GroGro("test.gro",5);
		//PDBfile.writePDB(gro.getAtoms(),new FileIO("test1.pdb",FileIO.WRITING));
	}
}
