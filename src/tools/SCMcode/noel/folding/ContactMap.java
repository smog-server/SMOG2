package noel.folding;
import noel.util.math.*;
import noel.util.geom.*;
import noel.util.*;
import noel.io.*;
import java.util.*;
import java.text.*;
/**
 * Routines related to my way of classifying contacts.
 */
public class ContactMap {	

	int numContacts;
	AtomList list;
	Atom[] atom;
	final static int CA = 1;
	final static int AACA = 2;
	
	public ContactMap(Atom[] atom) {
		numContacts = 0;
		this.atom = atom;
		//assumes that the last atoms position is the highest possible
		//and that position is 1 indexed
		list = new AtomList(atom[atom.length-1].getPosition(),1); 
	}

    /* informs the output routines whether to output contact distances */
    static boolean PRINT_DISTANCES = false;
	    
	//scripts :)
    // public static void main(String args[]){
    //  FileIO file1 = new FileIO(args[0],FileIO.BUFFERED_READING);
    //  FileIO file2 = new FileIO(args[1],FileIO.BUFFERED_READING);
    //  GroGro gro = new GroGro(args[2]);
    //  Atom[] atoms = gro.getAtoms();
    //  ContactMap map1 = new ContactMap(atoms);
    //  ContactMap map2 = new ContactMap(atoms);
    //  String s= file1.readLine();
    //  while(s!=null){
    //      //System.out.println(s);
    //      String[] tokens = FileIO.getTokens( s," ");
    //      map1.addContact(Integer.parseInt(tokens[0]),Integer.parseInt(tokens[1]));
    //      s= file1.readLine();
    //  }
    //  s= file2.readLine();
    //  while(s!=null){
    //      //System.out.println(s);
    //      String[] tokens = FileIO.getTokens( s," ");
    //      map2.addContact(Integer.parseInt(tokens[1]),Integer.parseInt(tokens[3]));
    //      s= file2.readLine();
    // 
    //  }
    //  //System.out.println(map1.numContacts);
    //  //System.out.println(map2.numContacts);
    //  Contact[] mA1 = map1.getContacts();
    //  for (int i = 0; i < mA1.length; i++) {
    //      System.out.println(mA1[i].skinnyToString());
    //      //System.out.println(atoms[mA1[i].leftIndex-1]+" "+atoms[mA1[i].rightIndex-1]);
    //  }   
    // }



	/**
	 * Returns a set of contacts based on shadow definition.
	 * Determines which atoms are in contact using the following algorithm:
	 * <li>1. Checks if two atoms are within <code>cutoff</code> angstroms
	 * <li>2. Checks if two atoms are greater than 3 residues apart or in
	 *        separate chains.
	 * <li>3. If so then checks to see if any atoms are between the two.  This is determined by
	 *	putting a light bulb at the center of each of the atoms and seeing if any of the
	 *	other atoms cast a shadow on the other.  If so then the contact is thrown out.
	 * <li>3b. The size of the atoms is determined by <code>radius</code>.
	 * <li>3c. A potentially shadowing atom is given a size of 0.5 angstrom when it is bonded
	 * 	to the shadowee.  This is because the C-C bond length is 1.5 angstrom and
	 * 	bonded neighbors should not shadow.
	 * @param atom atom list containing contacting atoms
	 * @param list {@link noel.folding.BondedList BondedList} from the .top file
	 * @param radius the radius of all the shadowing atoms
	 * @param bondedRadius radius of a shadowing atom if it is bonded to one of the atoms being shadowed
	 * @param cutoff the maximum distance to search for contacts
	 * @param proteinDelta enforces contacts to be |i-j|>proteinDelta
	 * @param rnaDelta enforces contacts to be |i-j|>rnaDelta
	 * @param type ShadowSettings.USE_SHADOW_MAP ShadowSettings.USE_CUTOFF_MAP are defined options
	 * @return list of {@link noel.folding.ContactMap ContactMap}
	 */ 
	public static ContactMap createContactMap(Atom[] atom, BondedList list, double radius, double bondedRadius, double cutoff, int rnaDelta, int proteinDelta) {		
		if(radius <= 0 && ShadowSettings.USE_SHADOW_MAP) {
			System.out.println("Shadowing radius is <= 0, treating instead as a cutoff map.");
			ShadowSettings.USE_SHADOW_MAP = false;
			ShadowSettings.USE_CUTOFF_MAP = true;
			ShadowSettings.BONDED_RADIUS = 0;
			ShadowSettings.SHADOW_RADIUS = 0;
		}	
		ContactMap map = new ContactMap(atom);
		JGrid3D grid = new JGrid3D(atom,cutoff);
		int total = grid.getNumGrids();
		Contact[] con = new Contact[0];
		int count = 0;
		int progressCount = 1;
		for(int g = 0; g < total; g++) { //for all grids
		    if(ShadowSettings.SHOW_PROGRESS) { 
				if(((float)g)/total*10 > progressCount) {
					progressCount++;
					ShadowMain.print2screen("\t"+(int)(((float)g)/total*100)+"% done and "+map.numContacts+" contacts."); 
				}
			}
			Vector atoms = grid.getGrid(g);
			if(atoms.size() > 0) {
				Coordinates here = ((Atom)atoms.get(0)).getCoords();
				Vector others = grid.getNeighbors(here);
				//System.out.println("grid #: "+g+" numMembers: "+atoms.size()+" numNeigh: "+others.size()+
				//" Coords: "+here);
				for(int i=0;i<atoms.size();i++) {
					Atom atom1=(Atom)atoms.get(i);
					int pos1=atom1.getPosition();
					for(int j=0;j<others.size();j++) {
						Atom atom2 = (Atom)others.get(j);			
						int pos2=atom2.getPosition();
						double dist = Coordinates.distance(atom1.getCoords(),atom2.getCoords());
						if (dist < cutoff) {
							if(!list.connected(pos1,pos2)) { //not bonded or dihedraled!
								if (atom1.getChain() == atom2.getChain()) {
									if (atom1.type == Residue.NUCLEIC_ACID) { //RNA-RNA or DNA-DNA
										if (atom1.residueNumber+rnaDelta < atom2.residueNumber) { 
											if (!map.connected(atom1,atom2)) {
												if (ShadowSettings.USE_SHADOW_MAP) {
													if (!shadowed(atom1,atom2,others,list,radius,bondedRadius)) {
														map.addContact(atom1,atom2);
													}
												} else if(ShadowSettings.USE_CUTOFF_MAP) { map.addContact(atom1,atom2); }
												else { ShadowSettings.throwError(1); }
											}
										}
									} else if(atom1.type == Residue.PROTEIN) { //Protein-Protein
										if ((atom1.residueNumber+proteinDelta) < atom2.residueNumber) { 
											if (!map.connected(atom1,atom2)) {
												if (ShadowSettings.USE_SHADOW_MAP) {
													if (!shadowed(atom1,atom2,others,list,radius,bondedRadius)) {
														map.addContact(atom1,atom2);
													}
												} else if(ShadowSettings.USE_CUTOFF_MAP) { map.addContact(atom1,atom2); }
												else { ShadowSettings.throwError(1); }
											}
										}
									} else { //ligand or unknown
										if (atom1.residueNumber < atom2.residueNumber) { 
											if (!map.connected(atom1,atom2)) {
												if (ShadowSettings.USE_SHADOW_MAP) {
													if(!shadowed(atom1,atom2,others,list,radius,bondedRadius)) {
														map.addContact(atom1,atom2);
													}
												} else if(ShadowSettings.USE_CUTOFF_MAP) { map.addContact(atom1,atom2); }
												else { ShadowSettings.throwError(1); }
											}
										}
									}
								//} else if (SPLICE) //add in splice information here, for now I am just treating the
								                    //spliced chains as separate.  No contacts if connected, but otherwise
								                    //they are ok regardless of deltas.
								
								} else { //different chains
									if (!map.connected(atom1,atom2)) {
										if (ShadowSettings.USE_SHADOW_MAP) {
											if(!shadowed(atom1,atom2,others,list,radius,bondedRadius)) {
												map.addContact(atom1,atom2);
											}
										} else if(ShadowSettings.USE_CUTOFF_MAP) { map.addContact(atom1,atom2); }
										else { ShadowSettings.throwError(1); }
									}

								}
							}//not bonded

						}//within same distance						
					}
				}

			}
		} 
		return map;
	}

	/** 
	 * Determines whether any atom in <code>Vector atom</code> shadows <code>Atom j</code>
	 * from <code>Atom k</code>.	
	 * <li>1. Checks to see if any atoms are between the two.  This is determined by
	 *	putting a light bulb at the center of each of the atoms and seeing if any of the
	 *	other atoms cast a shadow on the other.  The atoms are treated as spheres of given
	 * radius.  If so then the contact is thrown out.
	 * <li>2. The size of the atoms is determined by <code>radius</code>.
	 * <li>3. A potentially shadowing atom is given a size of 0.5 angstrom when it is bonded
	 * 	to the shadowee.  This is because the C-C bond length is 1.5 angstrom and
	 * 	bonded neighbors should not shadow.
	 * 
	 * @param j left index
	 * @param k right index
	 * @param atom atom list containing contacting atoms
	 * @param list list of bonds consistent with atoms
	 * @param radius the diameter of the atoms
	 * @param bondedRadius radius of a shadowing atom if it is bonded to one of the atoms being shadowed
	 */
	private static boolean shadowed(Atom j, Atom k, Vector atom, BondedList list, double radius, double bondedRadius) {
		int posj=j.getPosition();
		int posk=k.getPosition();
		double dist = Coordinates.distance(j.getCoords(),k.getCoords());
		for (int i = 0; i < atom.size(); i++) {
			Atom test = (Atom)atom.get(i);
			int posi = test.getPosition();
			double r = radius;
			double smallR = bondedRadius;
			try{		
				if (posi != posj && posi != posk) {
					if (dist > Coordinates.distance(test.getCoords(),j.getCoords())) {
						if (dist > Coordinates.distance(test.getCoords(),k.getCoords())) {
							//if(list.connected(posi,posk) || list.connected(posi,posj)) r = smallR;
							if(list.bonded(posi,posk) || list.bonded(posi,posj)) r = smallR;
							//check if j is shadowed from k
							Coordinates p = j.getCoords();
							Sphere a = new Sphere(r,test.getCoords());
							Sphere b = new Sphere(radius,k.getCoords());
							if (ShadowSettings.LEGACY_SHADOW) { 
								if (Geometry.isShadowedLegacy(p,a,b)) {
								    return true;
							    }
							} else {
								if (Geometry.isShadowed(p,a,b)) {
								    return true;
							    }
							}

							//check if k is shadowed from j
							p = k.getCoords();
							a = new Sphere(r,test.getCoords());
							b = new Sphere(radius,j.getCoords());
							if (ShadowSettings.LEGACY_SHADOW) { 
								if (Geometry.isShadowedLegacy(p,a,b)) {
								    return true;
							    }
							} else {
								if (Geometry.isShadowed(p,a,b)) {
								    return true;
							    }
							}
						}
					}
				}
			}catch(Exception e){
				System.out.println("Exception "+posi+" "+posj+" "+posk+" ");
			}
		}
		return false;
	}



	public void addContact(Atom a, Atom b){
		list.addEdge(a,b);
		numContacts++;
	}
	
    // public void addContactFromFile(String filename) {
    //  FileIO file = new FileIO(filename,FileIO.BUFFERED_READING);
    //  String line = file.readLine();
    //  if(line.indexOf(',')>0){ line = file.readLine(); } //skip the num line
    //  int i,j;
    //  while(line != null) {
    //      String[] tokens = FileIO.getTokens(line," ");
    //      if(tokens.length < 4){ //assume only one chain
    //          i = Integer.parseInt(tokens[0]);
    //          j = Integer.parseInt(tokens[1]);
    //      } else {
    //          i = Integer.parseInt(tokens[1]);
    //          j = Integer.parseInt(tokens[3]);
    //      }
    //      if(!connected(i,j)) { addContact(i,j); }
    //      line = file.readLine();
    //  }
    //  file.close();
    // }
    
	public boolean connected(Atom a, Atom b) {
        return list.isEdge(a,b);
	}
	
	public Atom[] getConnected(Atom a){
		Vector connected = list.neighborsOf(a);
		Atom[] con = new Atom[connected.size()];
		for(int i = 0; i < connected.size(); i++){
			con[i]=(Atom)connected.elementAt(i);
		}
		Arrays.sort(con);
		return con;
	}
	
	public static void setPrintDistance(boolean print) {
	    PRINT_DISTANCES = print;
	}

	/*
	 * Prints the contacts to standard output.  Format is "chain_i  i  chain_j  j"
	 */
	public void print(){
		for(int i = 0; i < atom.length; i++) {
		    int pos1 = atom[i].getPosition();
			Atom[] con = getConnected(atom[i]);
			for(int j = 0; j < con.length; j++){
				if( pos1 < con[j].getPosition()) System.out.println(atom[i].getChain()+" "+pos1+" "+con[j].getChain()+" "+con[j].getPosition());
			}
		}
	}
	/*
     * Prints the contacts to file and reindexes.  Format is "chain_i  i  chain_j  j"
     * @param output output file
     */
    public void print(FileIO output,Vector chains) {
            for(int i = 0; i < atom.length; i++) {
                    Atom[] con = getConnected(atom[i]);
                    //System.out.println(((Chain)chains.elementAt(atom[0].getChain())).firstAtom);
                    //System.exit(0);
                    for(int j = 0; j < con.length; j++) {
                        if( atom[i].getPosition() < con[j].getPosition()) {
                                output.write(atom[i].getChain()+" "+ (atom[i].getPosition()-((Chain)chains.elementAt(atom[i].getChain()-1)).firstAtom+1) +" "+con[j].getChain()+" "+(con[j].getPosition()-((Chain)chains.elementAt(con[j].getChain()-1)).firstAtom+1)+"\n");
                        }
                    }
            }
    }
	/*
     * Prints [ pairs ] section for Gaussian function type 6 
     * @param output output file
     */
    public void printGauss6(FileIO output,Vector chains) {
            for(int i = 0; i < atom.length; i++) {
                    Atom[] con = getConnected(atom[i]);
                    //System.out.println(((Chain)chains.elementAt(atom[0].getChain())).firstAtom);
                    //System.exit(0);
                    for(int j = 0; j < con.length; j++) {
                        int left = atom[i].getPosition();
                        int right = con[j].getPosition();
                        if( left < right) {                        
                            int ftype = 6;
                            float distance = (float) (Coordinates.distance(atom[i].getCoords(),con[j].getCoords()) / 10);
                            output.write(left+" "+right+" 6 DEPTH "+distance+" "+((float)(Math.sqrt(distance*distance/50/Math.log(2))))+" EXCLUDED \n");
                        }
                    }
            }
        }        
/*
	 * Prints the residue contacts to file and reindexes.  Format is "chain_i  res_i  chain_j  res_j"
	 * @param output output file
	 */


	public void printCoarse(int type,FileIO output,Vector chains) {
	    ContactMap caMap = getResidueContacts();
		if(PRINT_DISTANCES) caMap.setPrintDistance(true);
	    if(type == ContactMap.CA) caMap.printRes(output,chains);
	    else if(type == ContactMap.AACA) caMap.print(output,chains);
	}
	public void printCoarse(int type,FileIO output) {
	    ContactMap caMap = getResidueContacts();
		if(PRINT_DISTANCES) caMap.setPrintDistance(true);
	    if(type == ContactMap.CA) caMap.printRes(output);
	    else if(type == ContactMap.AACA) caMap.print(output);
	}
	    
    public void printRes(FileIO output,Vector chains) {
		for(int i = 0; i < atom.length; i++) {
			Atom[] con = getConnected(atom[i]);
			//System.out.println(((Chain)chains.elementAt(atom[0].getChain())).firstAtom);
			//System.exit(0);
			for(int j = 0; j < con.length; j++){
				if( atom[i].getPosition() < con[j].getPosition()) { 
				    output.write(atom[i].getChain()+" "+ (atom[i].getResNum()-((Chain)chains.elementAt(atom[i].getChain()-1)).firstRes+1) +" "+
				        con[j].getChain()+" "+(con[j].getResNum()-((Chain)chains.elementAt(con[j].getChain()-1)).firstRes+1)+"\n");
			    }
			}
		}
	}
	/*
	 * Prints the residue contacts to file.  Format is "chain_i  res_i  chain_j  res_j"
	 * @param output output file
	 */
	public void printRes(FileIO output) {
		for(int i = 0; i < atom.length; i++) {
			Atom[] con = getConnected(atom[i]);
			//System.out.println(((Chain)chains.elementAt(atom[0].getChain())).firstAtom);
			//System.exit(0);
			for(int j = 0; j < con.length; j++){
				if( atom[i].getPosition() < con[j].getPosition()) { 
					if (PRINT_DISTANCES) {
						String outputDistances = String.format(atom[i].getChain()+" "+atom[i].getResNum()+" "+con[j].getChain()+" "+con[j].getResNum()+" "+
							"%.3f\n",(Coordinates.distance(atom[i].getCoords(),con[j].getCoords())/10));
						output.write(outputDistances);
					} else output.write(atom[i].getChain()+" "+atom[i].getResNum()+" "+con[j].getChain()+" "+con[j].getResNum()+"\n");
			    }
			}
		}
	}
	
	/*
	 * Prints the contacts to file.  Format is "chain_i  i  chain_j  j"
	 * @param output output file
	 */
	public void print(FileIO output) {
	    //System.out.println("PRINTING!");
		for(int i = 0; i < atom.length; i++) {
			Atom[] con = getConnected(atom[i]);
			for(int j = 0; j < con.length; j++){
				if( atom[i].getPosition() < con[j].getPosition()) {
				     if (PRINT_DISTANCES) {
						 String outputDistances = String.format(atom[i].getChain()+" "+atom[i].getPosition()+" "+con[j].getChain()+" "+con[j].getPosition()+" "+
				            "%.3f\n",(Coordinates.distance(atom[i].getCoords(),con[j].getCoords())/10));
						 output.write(outputDistances);
					 } else { 
						 output.write(atom[i].getChain()+" "+atom[i].getPosition()+" "+con[j].getChain()+" "+con[j].getPosition()+"\n");
					 }
			    }
			}
		}
	}

	/*
	 * Creates a residue-residue contact map.  First grabs all the CA atoms.
	 * Then adds a contact between all residues with at least one atom-atom 
	 * interaction.
	 * @return The CA contact map
	 */
    public ContactMap getResidueContacts(){
         Atom[] ca = GroGro.getCAorN1(atom);
         ContactMap caMap = new ContactMap(ca);
         Contact[] AAcon = getContacts();
         for(int k=0;k<AAcon.length;k++){
             int i = AAcon[k].left.getResNum()-1;
             int j = AAcon[k].right.getResNum()-1;
             if(!caMap.connected(ca[i],ca[j])){ 
                 caMap.addContact(ca[i],ca[j]); 
             }
         }
         return caMap;
    }

	/*
	 * Returns an array which gives array[AAindex]=CAindex
	 * @return conversion array
	 */
    // public static int[] getAAtoCAindexConversion(ContactMap aaMap, ContactMap caMap){
    //  Contact[] aacon = aaMap.getContacts();
    //  Contact[] cacon = caMap.getContacts();
    //  int[] conv = Utilities.initializeWithZeros(new int[aacon.length+1]);
    //  for(int k=0;k<aacon.length;k++){
    //      int i = aacon[k].left.getResNum();
    //      int j = aacon[k].right.getResNum();
    //      for(int m=0;m<cacon.length;m++){
    //          if(i==cacon[m].left.getResNum()&&j==cacon[m].right.getResNum()){ 
    //              conv[k+1]=m+1;
    //          }
    //      }
    //  }
    //  return conv;
    // }
	/**
	 * Returns the contact array wrapped as a <code>Contact[]</code>
	 */
	public Contact[] getContacts() {
		Contact[] map = new Contact[numContacts];
		int index = 0;
		for(int i = 0; i < atom.length; i++) {
			Atom[] con = getConnected(atom[i]);
			for(int j = 0; j < con.length; j++){
				if( atom[i].getPosition() < con[j].getPosition()) { 
					map[index] = new Contact(atom[i],con[j],1,1);
					index++;
				}
			}
		}
		return map;
	}
	/**
	 * Returns a contact map based on a cutoff definition.  All contacts which are less than
	 * cutoff apart and between atoms in residues greater than three apart are returned.
	 */
	public static ContactMap createCutoffContactMap(Atom[] atom, BondedList list, double cutoff, int rnaDelta, int proteinDelta) {		
		ContactMap map = new ContactMap(atom);
		JGrid3D grid = new JGrid3D(atom,cutoff);
		int total = grid.getNumGrids();
		//System.out.println(total);
		Contact[] con = new Contact[0];
		int count = 0;
		for(int g = 0; g < total; g++) { //for all grids
			Vector atoms = grid.getGrid(g);
			if(atoms.size() > 0) {
				Coordinates here = ((Atom)atoms.get(0)).getCoords();
				Vector others = grid.getNeighbors(here);
				//System.out.println("grid #: "+g+" numMembers: "+atoms.size()+" numNeigh: "+others.size()+
				//" Coords: "+here);
				for(int i=0;i<atoms.size();i++) {
					Atom atom1=(Atom)atoms.get(i);
					int pos1=atom1.getPosition();
					for(int j=0;j<others.size();j++) {
						Atom atom2 = (Atom)others.get(j);			
						int pos2=atom2.getPosition();
						double dist = Coordinates.distance(atom1.getCoords(),atom2.getCoords());
						if (dist < cutoff) {
							if(!list.connected(pos1,pos2)) { //not bonded!
								if (atom1.getChain() == atom2.getChain()) {
									if (atom1.type == Residue.NUCLEIC_ACID) { //RNA-RNA or DNA-DNA
									//if (atom1.getResName().length() == 1) { //RNA-RNA
										if (atom1.residueNumber+rnaDelta < atom2.residueNumber) { 
										    if (!map.connected(atom1,atom2)) {
											    map.addContact(atom1,atom2);
										    }
											//System.out.println(atom1.getPosition()+" "+atom2.getPosition());
										}

									} else { //Protein-Protein or ligand or unknown
										if ((atom1.residueNumber+proteinDelta) < atom2.residueNumber) { 
											if (!map.connected(atom1,atom2)) {
												map.addContact(atom1,atom2);
												//System.out.println(atom1.getPosition()+" "+atom2.getPosition());
											}
										}
									}
								} else { //different chains
									if (!map.connected(atom1,atom2)) {
										map.addContact(atom1,atom2);
										//System.out.println(atom1.getPosition()+" "+atom2.getPosition());
									}

								}
							}//not bonded

						}//within same distance
					}
				}

			}
		} 
		return map;
	}


	/**********************************************
	 * OLD CODE
	 ***********************************************/


	/**
	 * Returns a contact map based on a cutoff definition.  All contacts which are less than
	 * cutoff apart and between atoms in residues greater than three apart are returned.
	 */
	public static Contact[] createCutoffContactMap2(Atom[] atom, BondedList list, double cutoff) {	
		Contact[] con = new Contact[0];
		int count = 0;
		for (int i = 0; i < atom.length; i++) {
			for (int j = i+1; j < atom.length; j++) {
				double dist = Coordinates.distance(atom[i].getCoords(),atom[j].getCoords());
				if (dist < cutoff) {
					if ((atom[i].residueNumber+3) < atom[j].residueNumber) {
						if (!inContact(i+1,j+1,con)) {
							con=addContact(con,atom,i,j,dist,1);
						}
					}
				}
			}
		}
		return con;
	}
	/** 
	 * Removes the repeat contacts in a contact list.
	 * @param con List to trim
	 * @return the shorter list
	 */
	public static Contact[] removeDoubles(Contact[] con) {
		int count = 0;
		int[] list = new int[con.length];
		for (int i = 0; i < con.length; i++) {
			for (int j = i+1; j < con.length; j++) {
				if (con[i].leftIndex == con[j].leftIndex && con[i].rightIndex == con[j].rightIndex) {
					list[count]=i;
					count++;
					break;
				}
			}
		}
		int numUnique = con.length - count;
		count = 0;
		Contact[] conNoD = new Contact[numUnique];
		for (int i = 0; i < con.length; i++) {
			if (list[count] == i) {count++;}
			else {
				conNoD[i-count] = con[i];
			}
		}
		return conNoD;
	}


	/**
	 * The list of contacts.
	 */
	//int[][] contacts;
	//Matrix map;
	/**
	 * Reads in a contact map which has the first two columns
	 * the indices of the contacts.
	 * @param String contactFile
	 * @param int numAtoms the dimension of the map
	 */
	/*public ContactMap(String contactFile, int numAtoms) {
	  String[] file = FileIO.readFromFile(contactFile,0);
	  contacts = new int[file.length][2];
	  map = new Matrix(numAtoms,numAtoms);
	  for (int i = 0; i < file.length; i++) {
	  String[] tokens = FileIO.getTokens(file[i]," ");
	  int j = Integer.parseInt(tokens[0]);
	  int k = Integer.parseInt(tokens[1]);
	  contacts[i][0] = j;
	  contacts[i][1] = k;
	  map.set(j,k,1.0);
	  map.set(k,j,1.0);
	  }
	  }*/
	/**
	 * Returns the contact array
	 */
	/*public int[][] getContactArray() {
	  return contacts;
	  }*/
	/**
	 * Returns the contact array
	 */
	/*public Matrix getContactMatrix() {
	  return map;
	  }*/	

	/*
	 * Same as <code>createContactMap</code> except for one chain and doesnt use grids.
	 */ 
	public static Contact[] createContactMapSlow(Atom[] atom, BondedList list, double radius, double cutoff) {
		Contact[] con = new Contact[0];
		int count = 0;
		for (int i = 0; i < atom.length; i++) {
			for (int j = i+1; j < atom.length; j++) {
				double dist = Coordinates.distance(atom[i].getCoords(),atom[j].getCoords());
				if (dist < cutoff) {
					if ((atom[i].residueNumber+3) < atom[j].residueNumber ||
							atom[i].getChain() != atom[j].getChain()) {
						if (!inContact(i+1,j+1,con)) {
							if (!shadowedSlow(i,j,atom,list,radius)) {
								con=addContact(con,atom,i,j,dist,1);
							}
						}
					}
				}
			}
		}
		return con;
	}
	private static boolean shadowedSlow(int j, int k, Atom[] atom, BondedList list, double radius) {
		boolean test = false;
		double dist = Coordinates.distance(atom[k].getCoords(),atom[j].getCoords());
		for (int i = 0; i < atom.length; i++) {
			double r = radius;
			double smallR = 0.5;
			if (i != j && i != k) {
				if (dist > Coordinates.distance(atom[i].getCoords(),atom[j].getCoords())) {
					if (dist > Coordinates.distance(atom[i].getCoords(),atom[k].getCoords())) {
						if(list.bonded(i,k) || list.bonded(i,j)) r = smallR;
						//check if j is shadowed from k
						Coordinates p = atom[j].getCoords();
						Sphere a = new Sphere(r,atom[i].getCoords());
						Sphere b = new Sphere(radius,atom[k].getCoords());
						if (Geometry.isShadowed(p,a,b)) return true;

						//check if k is shadowed from j
						p = atom[k].getCoords();
						a = new Sphere(r,atom[i].getCoords());
						b = new Sphere(radius,atom[j].getCoords());
						if (Geometry.isShadowed(p,a,b)) return true;
					}
				}
			}
		}
		return test;
	}
	/**
	 * Utility method to decrease the size of a <code>Contact[]</code> array.
	 * @param con Array of contacts.
	 * @param i index of contact to remove (zero based indexing)
	 * @return The shortened contact list
	 */
	public static Contact[] removeContact(Contact[] con, int i) {
		Contact[] con2 = new Contact[con.length-1];
		for (int k = 0; k < i; k++) {con2[k] = con[k];}
		for (int k = i+1; k < con.length; k++) {con2[k-1] = con[k];}
		return con2;
	}
	/**
	 * Utility method to increase the size of a <code>Contact[]</code> array.
	 * @param con Array of contacts.
	 * @param atom Array of atoms.
	 * @param i left index
	 * @param j right index
	 * @param dist Distance of contact.
	 * @param strength Epsilon of contact.
	 * @return The lengthened contact list
	 */
	public static Contact[] addContact(Contact[] con, Atom[] atom, int i, int j, 
			double dist, double strength) {
		Contact[] con2 = new Contact[con.length+1];
		for (int k = 0; k < con.length; k++) {
			con2[k] = con[k];
		}
		con2[con.length] = new Contact(atom[i],atom[j],dist,strength);
		return con2;
	}
	/**
	 * Utility method to determine if two atoms are in contact from a list of contacts
	 */
	public static boolean inContact(int a, int b, Contact[] con) {
		for (int i = 0; i < con.length; i++) {
			if (con[i].leftIndex == a && con[i].rightIndex == b) return true;
			if (con[i].leftIndex == b && con[i].rightIndex == a) return true;
		}
		return false;
	}
}
