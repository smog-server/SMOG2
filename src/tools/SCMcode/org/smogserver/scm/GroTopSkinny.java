package org.smogserver.scm;

import org.smogserver.io.*;
import org.smogserver.util.*;
import org.smogserver.util.math.*;
import org.smogserver.util.geom.*;
import java.util.*;

/**
 * Parses a GROMACS AA-Go topology file created from the Read program.  
 */
public class GroTopSkinny extends java.io.File implements Topology {
	public int numBonds, numAngles, numDih, numCon, numRep, numAtom, numResidue;
	public int conCount = 0, dihCount = 0, angCount = 0, bondCount = 0, atomtypes = 0, atomCount = 0;
	String top;
	BondedList bondedList;
	Atom[] atom;
	Hashtable atomTypeRadii;
	boolean noDihedral;
	
	/**
	 * Creates a new GROMACS topology from the .top and the .gro.
	 * @param top GROMACS .top file
	 * @param gro GROMACS .gro file object
	 */
	public GroTopSkinny(String top, Structure gro, boolean noDihedral) {
		super(top);
		numAtom = 0;
		this.noDihedral = noDihedral;
		this.top = top;
		atom = gro.getAtoms();
		atomTypeRadii = new Hashtable();
		FileIO file = new FileIO(top,FileIO.BUFFERED_READING);
		if (file == null) { //file does not exist
			System.out.println("Gromacs .top, "+file+", does not exist!  Exiting.");
			System.exit(1);
		} else {
			if(!ShadowSettings.SMOG2_OUTPUT_ON) { ShadowMain.print2screen("Reading Topology."); }
			String line = file.readLine();
			while(line != null) {
				///************ AtomTypes **********************//
			    if (line.indexOf("[ atomtypes ]") != -1) {
			        line = file.readLine();//comment
					line = file.readLine();
					boolean done = false;
					while (!done) {
						String[] tokens = FileIO.getTokens(line," \t");
						if(tokens.length != 6) { done = true; }
						else {
					        atomtypes++;
        					double c12 = Double.parseDouble(tokens[5]);
                            atomTypeRadii.put(tokens[0],new Double(Math.pow(c12,1.0/12)*10));
        					line = file.readLine();
        				}
					}
				}
				///************ Atoms **********************//
				if (line.indexOf("[ atoms ]") != -1) {
					line = file.readLine();
					line = file.readLine();
					boolean done = false;
					while (!done) {
						String[] tokens = FileIO.getTokens(line," ");
						try {
							atomCount = Integer.parseInt(tokens[0]);
							numAtom++;
						} catch (Exception e) {	done = true; }
						line = file.readLine();
					}
				}
				///************ Contacts **********************//
				if (line.indexOf("[ pairs ]") != -1) {
					line = file.readLine();
					line = file.readLine();
					boolean done = false;
					while (!done) {
						String[] tokens = FileIO.getTokens(line," \t");
						try {
							Integer.parseInt(tokens[1]);
							conCount++;
						} catch (Exception e) {	done = true; 
							//System.out.println("Found "+conCount+" pairs.");
						}
						line = file.readLine();
					}
				}
				///************ bonds **********************//
				if (line.indexOf("[ bonds ]") != -1) {
					line = file.readLine();//comment
					line = file.readLine();
					boolean done = false;
					while (!done) {
						String[] tokens = FileIO.getTokens(line," \t");
						try {
							Integer.parseInt(tokens[1]);
							bondCount++;
						} catch (Exception e) {	done = true; 
							//System.out.println("Found "+bondCount+" bonds.");
						}
						line = file.readLine();
					}
				}
				///************ Angles **********************//
				if (line.indexOf("[ angles ]") != -1) {
					line = file.readLine();//comment
					line = file.readLine();
					boolean done = false;
					while (!done) {
						String[] tokens = FileIO.getTokens(line," \t");
						try {
							Integer.parseInt(tokens[2]);
							angCount++;
						} catch (Exception e) {	done = true; 
							//System.out.println("Found "+angCount+" angles.");
						}
						line = file.readLine();
					}
				}
				///************ Dihedrals **********************//
				if (line.indexOf("[ dihedrals ]") != -1) {
					line = file.readLine();//comment
					line = file.readLine();
					boolean done = false;
					while (!done) {
						String[] tokens = FileIO.getTokens(line," \t");
						try {
							Integer.parseInt(tokens[1]);
							dihCount++;
						} catch (Exception e) {	done = true; 
							//System.out.println("Found "+dihCount+" dihedrals.");
						}
						line = file.readLine();
					}
				}
				line = file.readLine();
			}
		}
		if (bondCount==0 || dihCount==0 || angCount==0) { 
			if(!ShadowSettings.SMOG2_OUTPUT_ON) ShadowMain.print2screen("Warning: bonds, dihedrals or angles are missing.");
		}
		if(!ShadowSettings.SMOG2_OUTPUT_ON) { ShadowMain.print2screen("Initializing connected list.");	}
		bondedList = new BondedList((Topology)this,noDihedral);
	}

	public int[][] getDihedrals() {
		int[][] temp = new int[dihCount][4];
		FileIO file = new FileIO(top,FileIO.BUFFERED_READING);
		String line = file.readLine();
		while (line != null) {
			if (line.indexOf("[ dihedrals ]") != -1) {
				line = file.readLine();//comment
				line = file.readLine();
				boolean done = false;
				int i = 0;
				while (!done) {
					String[] tokens = FileIO.getTokens(line," \t");
					try {
						temp[i][0] = Integer.parseInt(tokens[0]);
						temp[i][1] = Integer.parseInt(tokens[1]);
						temp[i][2] = Integer.parseInt(tokens[2]);
						temp[i][3] = Integer.parseInt(tokens[3]);
					} catch (Exception e) { done = true; 
						if(i!=dihCount){ System.out.println("Only parsed "+i+" dihedrals.  Problem in topology at line:"+line+".  Exiting.");  System.exit(1); }
					}
					line = file.readLine();
					i++;
				}
			}
			line = file.readLine();
		}
		return temp;
	}
	public int[][] getAngles() {
		int[][] temp = new int[angCount][3];
		FileIO file = new FileIO(top,FileIO.BUFFERED_READING);
		String line = file.readLine();
		while (line != null) {
			if (line.indexOf("[ angles ]") != -1) {
				line = file.readLine();//comment
				line = file.readLine();
				boolean done = false;
				int i = 0;
				while (!done) {
					String[] tokens = FileIO.getTokens(line," \t");
					try {
						temp[i][0] = Integer.parseInt(tokens[0]);
						temp[i][1] = Integer.parseInt(tokens[1]);
						temp[i][2] = Integer.parseInt(tokens[2]);
					} catch (Exception e) { done = true; 
						if(i!=angCount){ System.out.println("Only parsed "+i+" angles.  Problem in topology.  Exiting.");  System.exit(1);}
					}
					line = file.readLine();
					i++;
				}
			}
			line = file.readLine();
		}
		return temp;
	}
    /**
     * Returns the bond information in the top file. [i][j][type]
     * @return int[][] [i][j][type]
     */
	public int[][] getBonds() {
		int[][] temp = new int[bondCount][3];
		FileIO file = new FileIO(top,FileIO.BUFFERED_READING);
		String line = file.readLine();
		while (line != null) {
			if (line.indexOf("[ bonds ]") != -1) {
				line = file.readLine();//comment
				line = file.readLine();
				boolean done = false;
				int i = 0;
				while (!done) {
					String[] tokens = FileIO.getTokens(line," \t");
					try {
						if(tokens.length == 2) { //a bond that will be filled in by bondtypes, give assumed type 1
							temp[i][0] = Integer.parseInt(tokens[0]);
							temp[i][1] = Integer.parseInt(tokens[1]);
							temp[i][2] = 1;
						} else {
							temp[i][0] = Integer.parseInt(tokens[0]);
							temp[i][1] = Integer.parseInt(tokens[1]);
							temp[i][2] = Integer.parseInt(tokens[2]);
						}							
					} catch (Exception e) { done = true; 
						if(i!=bondCount){ System.out.println("Only parsed "+i+" bonds.  Problem in topology at line:"+line+".  Exiting.");  System.exit(1);}
					}
					line = file.readLine();
					i++;
				}
			}
			line = file.readLine();
		}
		return temp;
	}
	
	/**
	 * Determines the volume of a given topology of atoms.  Uses the atomtypes information to get atom sizes.
	 * Takes each atom and adds its volume, then subtracts half of any overlapping volume from bonded neighbors.
	 * Finds the bonded neighbors fromt the bondedList.
	 * @return volume of atoms
	 */
	public double getVolume() {
	    double volume = 0;
	    double overlap = 0;
	    String id1="", id2="";
        try {
            for (int i=0; i<atom.length; i++) {
                Vector bonded = bondedList.bondedTo(atom[i].getPosition());
                overlap = 0;
                for (int j=0; j<bonded.size(); j++) {
                    int otherIndex = ((Integer) bonded.elementAt(j)).intValue() - 1;
                    id1 = atom[i].getID();
                    id2 = atom[otherIndex].getID();
                    double intersect = Sphere.intersectionVolume(
                        (Double)atomTypeRadii.get(atom[i].getID()),
                        (Double)atomTypeRadii.get(atom[otherIndex].getID()),
                        Coordinates.distance(atom[i].getCoords(),atom[otherIndex].getCoords()));
                        overlap += 0; //intersect;                
                }
                // the volume contribution of an atom is the volume of the atom minus any overlaps
                // to account for double counting we divide the overlap by 2 (each overlap is counted twice)
                volume += 4.0/3 * Math.PI * Math.pow((Double)atomTypeRadii.get(atom[i].getID()),3) - overlap/2;
            }
        } catch (Exception e) {
            System.out.println("Problem with formatting.  An atom ID does not agree between gro and top.\n"+id1+" or "+id2);
            System.exit(1);
        }
	    return volume;
	}
	
	/**
	 * Sets this GroTop with a new contact map contained in a file.
	 */
	/*public void setContactsFromFile(String conFile) {
	  FileIO file = new FileIO(conFile,FileIO.BUFFERED_READING);
	  String line = file.readLine();
	  int numContacts=0;
	  if (line.indexOf(",") >= 0) {
	  numContacts = Integer.parseInt((FileIO.getTokens(line,", "))[0]);
	  }
	  else {
	  file.close();
	  file = new FileIO(conFile,FileIO.BUFFERED_READING);
	  numContacts = file.readAllLines(0).length;
	  file.close();
	  file = new FileIO(conFile,FileIO.BUFFERED_READING);
	  }
	  Contact[] newCon = new Contact[numContacts];
	  for (int k = 0; k < numContacts; k++) {
	  line = file.readLine();
	  String[] tokens = FileIO.getTokens(line," ");
	  int i = 0,j=0;
	  try {
	  i = Integer.parseInt(tokens[1]);
	  j = Integer.parseInt(tokens[3]);
	  } catch (Exception e) { System.out.println("You forgot your chain indentifiers!"); }
	  double d = Coordinates.distanceSq(atom[i-1].getCoords(),atom[j-1].getCoords());
	  newCon[k] = new Contact(atom[i-1],atom[j-1],d,1);
	  }
	  file.close();
	  this.contacts = newCon;
	  }*/

	/**
	 * Returns the contacts with the native distances those in the GroGro.  In the future it will compute 
	 * epsilons.
	 */
	/*public Contact[] getAtomContacts() {
	  Contact[] con = new Contact[contacts.length];
	  for (int i = 0; i < con.length; i++) {
	  double dist = Coordinates.distanceSq(contacts[i].left.getCoords(),contacts[i].right.getCoords());
	  con[i] = new Contact(contacts[i].left,contacts[i].right,dist,1);
	  }
	  return con;
	  }*/
	public Atom[] getAtoms() {
		return atom;
	}
	public BondedList getBondedList() {
		return bondedList;
	}
	public int getNumTopAtoms() {
		return numAtom;
	}
	/**
	 * Gives the contacts in which <code>residueNumber</code> is a member.  First residue is number 1.
	 * @param residueNumber index of the residue in question
	 * @return array of Contact indices
	 */
	/*public int[] getContactIndicesOfResidue(int residueNumber) {
	  int[] list = Utilities.initializeWithZeros(new int[numAtom]);
	  int count = 0;
	  for (int i = 0; i < contacts.length; i++) {
	  if (contacts[i].left.getResNum() == residueNumber || contacts[i].right.getResNum() == residueNumber) {
	  list[count]=i;
	  count++;
	  }
	  }
	  int[] conList = new int[count];
	  for (int i = 0; i < count; i++) {
	  conList[i] = list[i];
	  }
	  return conList;
	  }*/
	/**
	 * Gives the contacts in which <code>residueNumber</code> is a member.
	 * @param residueNumber index of the residue in question
	 * @return array of Contacts
	 */
	/*public Contact[] getContactsOfResidue(int residueNumber) {
	  int[] list = Utilities.initializeWithZeros(new int[numAtom]);
	  int count = 0;
	  for (int i = 0; i < contacts.length; i++) {
	  if (contacts[i].left.getResNum() == residueNumber || contacts[i].right.getResNum() == residueNumber) {
	  list[count]=i;
	  count++;
	  }
	  }
	  Contact[] conList = new Contact[count];
	  for (int i = 0; i < count; i++) {
	  conList[i] = contacts[list[i]];
	  }
	  return conList;
	  }*/

}
