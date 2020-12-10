package org.smogserver.util;

import java.util.*;

/**
* ConnectionList is a set of nodes that can be connected by edges.  It is represented by a list of lists to minimize (memory use)*(query time).
*/
public class ConnectionList extends ArrayList {
    
    private int startIndex = 0;
    
    public ConnectionList(int numNodes) {
        super();
        for (int i = 0; i < numNodes; i++) {
			add(new Vector());
		}
    }
    /**
    * Assumes that the first element is <code>startIndex</code>
    */
    public ConnectionList(int numNodes, int startIndex) {
        this(numNodes);
        this.startIndex = startIndex;
    }
    /**
    * Creates an edge between two nodes
    * @param a index of node 1
    * @param b index of node 2
    * @return true if edge already exists
    */
    public boolean addEdge(int a, int b) {
        if (!isEdge(a,b)) {
			((Vector)get(a-startIndex)).add(new Integer(b));
			((Vector)get(b-startIndex)).add(new Integer(a));
			return false;
		} else { return true; }
    }
    
    /**
	* Returns whether two nodes are connected
	* @param a index of node 1
	* @param b index of node 2
	* @return true if connected
	*/
	public boolean isEdge(int a, int b) {
		//check a row
		for (int i =  0; i < ((Vector)get(a-startIndex)).size(); i++) {
			int ai = ((Integer)((Vector)get(a-startIndex)).elementAt(i)).intValue();
			if (b == ai) return true;
		}
		return false;
	}
	/**
	* Returns a Vector of the neighbors of the given node.
	* @param a index of node
	*/
	public Vector neighborsOf(int a) {
	    return (Vector)get(a-startIndex);
    }
    
    
}