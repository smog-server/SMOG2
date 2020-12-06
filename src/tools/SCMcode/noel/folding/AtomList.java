package noel.folding;

import java.util.*;

/**
* AtomList defines a set of connections between Atoms. It is represented by a list of lists to minimize (memory use)*(query time).
*/
public class AtomList extends ArrayList {
    
    private int startIndex = 0;
    
    public AtomList(int numNodes) {
        super();
        for (int i = 0; i < numNodes; i++) {
			add(new Vector());
		}
    }
    /**
    * Assumes that the first element is <code>startIndex</code>
    */
    public AtomList(int numNodes, int startIndex) {
        this(numNodes);
        this.startIndex = startIndex;
    }
    /**
    * Creates an edge between two atoms
    * @param a atom 1
    * @param b atom 2
    * @return true if edge already exists
    */
    public boolean addEdge(Atom a, Atom b) {
        if (!isEdge(a,b)) {
			((Vector)get(a.getPosition()-startIndex)).add(b);
			((Vector)get(b.getPosition()-startIndex)).add(a);
			return false;
		} else { return true; }
    }
    
    /**
	* Returns whether two atoms are connected
	* @param a atom 1
	* @param b atom 2
	* @return true if connected
	*/
	public boolean isEdge(Atom a, Atom b) {
		//check a row
		for (int i =  0; i < ((Vector)get(a.getPosition()-startIndex)).size(); i++) {
			Atom ai = ((Atom)((Vector)get(a.getPosition()-startIndex)).elementAt(i));
			if (b.getPosition() == ai.getPosition()) return true;
		}
		return false;
	}
	/**
	* Returns a Vector of the neighbor atoms of the given node.
	* @param a atom
	*/
	public Vector neighborsOf(Atom a) {
	    return (Vector)get(a.getPosition()-startIndex);
    }
    
    
}