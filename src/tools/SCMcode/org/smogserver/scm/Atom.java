package org.smogserver.scm;
import org.smogserver.util.math.*;


/**
* Represents the information about an atom contained in a PDB file.
*/
public class Atom implements org.smogserver.util.geom.JGridable3D, Comparable {
	Coordinates coords;
	String id;
	int position;
	String residueName;
	int residueNumber;
	final static String pdbIdentifier = "ATOM";
	int chainNum = -1;
	int type; //Nucleic Acid, Protein or Ligand
	Chain chain = null;

	public Atom(String id, int position, Coordinates coords, String residueName,
						int residueNumber) {
		this.id = id;
		this.position = position;
		this.coords = coords;
		this.residueName = residueName;
		this.residueNumber = residueNumber;
		this.type = Residue.getType(residueName); 
		//System.out.println(isHydrogen()+" "+id);
	}
	public Atom(String id, int position, Coordinates coords, String residueName,
						int residueNumber, int index) {
		this.id = id;
		this.position = position;
		this.coords = coords;
		this.residueName = residueName;
		this.residueNumber = residueNumber;
		this.type = Residue.getType(residueName); 
		//System.out.println(this);
		//System.out.println("my type is: "+type);
		//System.out.println(isHydrogen()+" "+id);
	}

	public Coordinates getCoords() { return coords; }
	public double getX() {	return coords.x; }
	public double getY() {	return coords.y; }
	public double getZ() {	return coords.z; }
	public String getID() { return id; }
	public int getPosition() { return position; }
	/**
	 * PDB identifier for an atom is <code>"ATOM"</code>.
	 */
	public String getPDBIdentifier() { return pdbIdentifier; }
	public String getResName() { return residueName; }
	public int getResNum() { return residueNumber; }
	public String toString() {
		return id+" "+getChain()+" "+position+" "+residueNumber+" "+residueName+" "+getX()+" "+getY()+" "+getZ();
	}
	public String getKey() { return new String("id getChain position residueNumber residueName X Y Z");}
	public void setChain(int chain) {
		chainNum = chain;
	}
	public int getChain() {
		return chainNum;
	}
	public void setChainObject(Chain chain) {
	    this.chain = chain;
	}
	public Chain getChainObject() {
		return chain;
	}
	public void updateCoords(double x, double y, double z) {
	    coords.update(x,y,z);
	}

	public static int getAtomNumber(Atom[] atom, String id, int resNum) {
		for (int i = 0; i < atom.length; i++) {
			if (atom[i].getResNum() == resNum) {
				if (atom[i].getID().compareTo(id) == 0) {
					return atom[i].getPosition();
				}
			}
		}
		return -1;
	}
	
	//this should be consistent with SMOG namings
	static String[] pdbNames = {" CA "," CB "," C  "," CG "," CG1"," CG2"," CD "," CD1"," CD2"," CE "," CE1"," CE2"," CE3"," CZ "," CZ2"," CZ3"," CH2"," O  "," OXT"," OG "," OG1"," OD1"," OD2"," OE1"," OE2"," OH "," N  "," ND1"," ND2"," NE "," NE1"," NE2"," NZ "," NH1"," NH2"," SG "," SD "," O1P"," O2P"," O3P"," P  "," C1*"," C2*"," C3*"," C4*"," C5*"," O1*"," O2*"," O3*"," O4*"," O5*"," N1 "," N2 "," N3 "," N4 "," N5 "," N6 "," N7 "," N8 "," N9 "," O1 "," O2 "," O3 "," O4 "," O5 "," O6 "," C1 "," C2 "," C3 "," C4 "," C5 "," C6 "," C7 "," C8 "," PA "," O1A"," O2A"," O3A"," PB "," O1B"," O2B"," O3B"," PG "," O1G"," O2G"," O3G"," PD "," O1D"," O2D"," O3D"," PE "," O1E"," O2E"," O5F"," C5F"," C4F"," O4F"," C3F"," O3F"," C2F"," O2F"," C1F"," N9A"," C8A"," N7A"," C5A"," C6A"," N6A"," N1A"," C2A"," N3A"," C4A"," O5J"," C5J"," C4J"," O4J"," C3J"," O3J"," C2J"," O2J"," C1J"," N9B"," C8B"," N7B"," C5B"," C6B"," N6B"," N1B"," C2B"," N3B"," C4B"," C10"," C11"," C12"," C13"," C14"," C15"," C16"," C17"," C19"," C20"," C21"," C22"," C23"," C25"," C26"," C27"," C28"," C29"," C30"," C31"," C32"," C33"," C34"," C35"," C39"," C43"," C48"," C49"," C51"," C52"," C53"," C54"," C56"," C57"," C9 "," C98"," C99"," N18"," N24"," N37"," N50"," O36"," O40"," O41"," O44"," O45"," O46"," O47"," O7 "," S38"," S42"," S8 "," C18"," C24"," C36"," C40"," C45"," C46"," C47"," C50"," C55"," C58"," N19"," N25"," N57"," O12"," O13"," O37"," O38"," O42"," O43"," O56"," O8 "," S35"," S39"," BMG"," CY "," CAY"," OY "," C37"," C38"," C41"," C42"," C44"," N10"," N11"," N12"," O10"," O11"," O14"," S1 "," S2 "," O9 "," C7L"," C8L"," C6*"," O6*"," CC1"," CC2"," CL1"," CL2"," NC1"," NC2"," OC1"," OC2"," S10"," C1P"," C1R"," C2P"," C2R"," C3P"," C3R"," C4R"," C5M"," C5R"," C60"," C61"," C6M"," C7B"," C9B","CO  "," N21"," N22"," N23"," N29"," N33"," N40"," N45"," N52"," N59"," N62"," O28"," O34"," O39"," O51"," O58"," O63"," O6R"," O7R"," O8R","ZN  "," NG "," CH3"," C1A"," C1B"," C3A"," C3B"," C7A"," N2A"," N2B"," O4A"," O5A"," O5B"," O6A"," O6B"," O7A"," O7B"};
	
	//uses the 'id' field to get the right name
	public String getPDBname() {
	    for (String name: Atom.pdbNames) {
	        if(name.trim().equals(id)) return name;
	    }
	    if(ShadowSettings.SMOG_ERRORS_ON) { //we want to die so that SMOG doesn't get confused with "BAD" atom names
	        System.out.println("Unrecongnized atom type: "+id+".  This is only a problem in printing a PDB file (option -p) since the spacing convention needs to be known."); 
	        System.out.println("To ignore this error remove the flag --smogErrors.  The atom type will be changed to \"BAD \".");
	        System.out.println("If you are on smog-server: to have your atom type added inquire at info@smog-server.org");
	        System.exit(1);
	    }    
	    return "BAD ";
	}
	
	/**
	* Determines if an atom is hydrogen by checking if the first character in the atom id is "H"
	* @return true if hydrogen
	*/
	public boolean isHydrogen() {
	    return id.substring(0,1).equals("H");
	}
	public int compareTo(Object a) {
	    return this.getPosition() - ((Atom)a).getPosition();
	}
}
