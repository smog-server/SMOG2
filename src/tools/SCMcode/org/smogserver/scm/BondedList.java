package org.smogserver.scm;
import java.util.*;
import org.smogserver.util.*;

/**
* BondedList contains a list of atoms interacting locally with bonds and dihedrals.  It uses {@link ConnectionList} as a data structure.
*/ 
public class BondedList {
	Topology topo;
	ConnectionList list; //not type 6
	ConnectionList listBondsType6; //type_6 (BMG contacts)
	ConnectionList listBonds; //just type 1
	boolean noDihedral;
	
	public BondedList(Topology topo, boolean noDihedral) {
	    this.noDihedral = noDihedral;
		this.topo = topo;
		int numAtom = topo.getAtoms().length;
		//set up topo structure
		list = new ConnectionList(numAtom,1);
		listBonds = new ConnectionList(numAtom,1);
		//listBondsType6 = new ConnectionList(numAtom,1);
		int[][] bonds = topo.getBonds();
		for (int i = 0; i < bonds.length; i++) {
		    if(bonds[i][2] == 6) { //BMG bond (fake bond), ignore it?
		        //listBondsType6.addEdge(bonds[i][0],bonds[i][1]);
		    }
		    else {
    		    list.addEdge(bonds[i][0],bonds[i][1]);
    		    listBonds.addEdge(bonds[i][0],bonds[i][1]);
		    }
		}
		bonds = null;
		if(!noDihedral) {
    		int[][] dihs = topo.getDihedrals();
    		for (int i =0; i < dihs.length; i++) {
    		    list.addEdge(dihs[i][0],dihs[i][1]);
    		    list.addEdge(dihs[i][0],dihs[i][2]);
    		    list.addEdge(dihs[i][0],dihs[i][3]);
    		    list.addEdge(dihs[i][1],dihs[i][2]);
    		    list.addEdge(dihs[i][1],dihs[i][3]);
    		    list.addEdge(dihs[i][2],dihs[i][3]);
    		}
    		dihs = null;
		}
	}
	/**
	* Returns whether two atoms in a Topology file are connected by a bond.
	* @param a index of atom 1
	* @param b index of atom 2
	* @return true if connected
	*/
	public boolean bonded(int a, int b) {
        return listBonds.isEdge(a,b);
	}
	/**
	* Returns whether two atoms in a Topology file are connected by a bond or dihedral.
	* @param a index of atom 1
	* @param b index of atom 2
	* @return true if connected
	*/
	public boolean connected(int a, int b) {
        if(noDihedral) return listBonds.isEdge(a,b);
        else return list.isEdge(a,b);	
    }
	/**
	* Returns a Vector of the atoms bonded to the given atom.
	* @param a index of atom 1
	*/
	public Vector bondedTo(int a) {
	    return listBonds.neighborsOf(a);
    }
    /**
	* Returns a Vector of the atoms connected to the given atom.
	* @param a index of atom 1
	*/
	public Vector connectedTo(int a) {
	    if(noDihedral) return listBonds.neighborsOf(a);
    	else return list.neighborsOf(a);
	}
	public String toString() {
		String out = "";
		int[][] bonds=topo.getBonds();
		for (int i = 0; i < bonds.length; i++) {
			out = out+bonds[i][0]+" "+bonds[i][1]+"\n";
		}
		return out;
	}
}


