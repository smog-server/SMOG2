package noel.folding;
import noel.util.*;
import java.util.*;
/** 
* Represents a monomer or chain.
*/
public class Chain {
	public int firstAtom;
	public int lastAtom;
	public int firstRes;
	public int lastRes;
	Vector<Atom> atom;
	int chainIndex;
	public boolean splice = false; //true if this chain is actually a splice with the previous chain
	public int splicedTo = -1; //set to the chain that this chain is spliced to
	
	/**
	* Create new chain.
	* @param atoms the atoms in the chain
	* @param index the chain number
	 */
	public Chain(Vector<Atom> atoms, int index) {
		firstAtom = ((Atom) atoms.firstElement()).getPosition();
		lastAtom = ((Atom) atoms.lastElement()).getPosition();
		firstRes = ((Atom) atoms.firstElement()).getResNum();
		lastRes = ((Atom) atoms.lastElement()).getResNum();
		//System.out.println(firstAtom+" "+lastAtom);
		atom = atoms;
		chainIndex = index;
	}
	/**
	* Create new chain.
	* @param atoms the atoms in the chain
	* @param index the chain number
	 */
	public Chain(Atom[] atoms, int index) {
		atom = new Vector();
		for (int i = 0; i < atoms.length; i++) { atom.add(atoms[i]); }
		firstAtom = ((Atom) atom.firstElement()).getPosition();
		lastAtom = ((Atom) atom.lastElement()).getPosition();
		chainIndex = index;
	}
	/**
	* Computes the center of mass assuming all atoms have equal mass.
	* @return x,y,z center of mass
	*/
	public double[] getCM() {
		Atom[] a = (Atom[]) atom.toArray();
		double[] cm = Utilities.initializeWithZeros(new double[3]);
		for (int i = 0; i < a.length; i++) {
			cm[0]+=a[i].getX();
			cm[1]+=a[i].getY();
			cm[2]+=a[i].getZ();
		}
		cm[0]/=a.length;
		cm[1]/=a.length;
		cm[2]/=a.length;
		return cm;
	}

	/**
	 * Return the chain index.
	 */
	public int getChainIndex() {
		return chainIndex;
	}

	/**
	 * Return the atoms in this chain.
	 */
	public Atom[] getAtoms() {
		return  atom.toArray(new Atom[0]);
	}
}
