package org.smogserver.util;

import java.util.*;
import org.smogserver.io.*;
/* This class is to help implement an arbitrarily length any-dimension 
 * array.  (Just making an array can be too memory intensive, and is not
 * general size).  You give this class a point contained in a 1-D int[].  All 
 * points must be of the same size and equal to depth.  Example: point (3,2,4) 
 * has depth of three.
 */
 
public class SparseIntegerArray {
	
	public NumberedList theList; //high level list
	NumberedList list; //utility list
	public int depth; //dimension
	boolean isIndexing; //should we index?
	Vector index;
	int density = 0; //number of distinct sites in the array
	
	public SparseIntegerArray(int depth, boolean indexing) {
		theList = new NumberedList(0);
		this.depth = depth;
		this.isIndexing = indexing;
		if (indexing) index = new Vector();
	}
	public SparseIntegerArray(int depth) {
		this(depth,true);
	}	
	//returns number of elements
	public int size() {
		return density;
	}
	public void increment(int[] p) {
		increment(p,1);
	}
	/**
	* Tells the index into the integerIndex of this point.
	* @return the index, -1 if p does not exist
	*/
	public int probeIndex(int[] p) {
	    int count = 0;
    	list = theList;
    	while (count < depth) {
            int found = Arrays.binarySearch(list.toArray(),new NumberedList(p[count]));
			if (found >=0 ) list = (NumberedList)list.get(found); //next layer
			else return -1;
			count++;
    	}
    	return ((Integer)list.get(1)).intValue();
	}
	/**
	* Tells the interger value of this point.
	* @return the value, -1 if p does not exist
	*/
	public int probeValue(int[] p) {
	    int count = 0;
    	list = theList;
    	while (count < depth) {
            int found = Arrays.binarySearch(list.toArray(),new NumberedList(p[count]));
			if (found >=0 ) list = (NumberedList)list.get(found); //next layer
			else return -1;
			count++;
    	}
    	return ((Integer)list.get(0)).intValue();
	}
	
	public void increment(int[] p, int amt) {
    	int count = 0;
    	list = theList;
    	while (count < depth) {
            int found = Arrays.binarySearch(list.toArray(),new NumberedList(p[count]));
			if (found >= 0) list = (NumberedList)list.get(found); //next layer
			else { //add new point
				NumberedList numL = new NumberedList(p[count]);
				list.add(numL);
				Collections.sort(list);
				list = numL;
			}
			count++;
		}
		if (list.size() < 1) {
			list.add(new Integer(0));
			density++;
			if (isIndexing) {
				int[] pcopy = new int[p.length]; //index it as well!
    			System.arraycopy(p,0,pcopy,0,p.length);
    			index.add(pcopy);
    			list.add(new Integer(index.size()-1)); //put the index into the list as the second element
			} 
		}
		int i = ((Integer)list.get(0)).intValue();
		i += amt;
		//list.clear();
		list.set(0,new Integer(i));
    }
    
    public String toString() {
    	int[][] points = getIntegerIndex();
    	StringBuffer buf = new StringBuffer();
    	for (int i = 0; i < index.size(); i++) {
    		buf.append(i+" "+points[i][depth]+" ");
	    	buf.append(Utilities.getPrint(points[i],false));
    	}	
    	return buf.toString();
    }
    
    public String printChars() {
    	int[][] points = getIntegerIndex();
    	StringBuffer buf = new StringBuffer();
    	for (int i = 0; i < index.size(); i++) {
    		buf.append(i+" "+points[i][depth]+" ");
	    	for (int j = 0; j < points[i].length-1; j++) {
	    		buf.append(CharValues.getChar((15-points[i][j]))+" ");
	    	}
	    	buf.append("\n");
    	}	
    	return buf.toString();
    }
    
    public String getStats() {
    	return "Contains - "+index.size()+" points";
    }
    
    public int[][] getIntegerIndex() {
    	int[] point;
    	int[][] newPoints = new int[index.size()][depth+1];
    	int row = 0;
    	Enumeration points = index.elements();
    	while (points.hasMoreElements()) {
    		point = (int[])points.nextElement();
    		for (int i = 0; i < point.length; i++) {
    			newPoints[row][i] = point[i];
    		}
    		int count = 0;
	    	list = theList;
	    	while (count < depth) {
                int found = Arrays.binarySearch(list.toArray(),new NumberedList(point[count]));
    			list = (NumberedList)list.get(found); //next layer
    			count++;
	    	}
	    	newPoints[row][depth] = ((Integer)list.get(0)).intValue();
	    	row++;
	    }
	    //if (removePermutations) newPoints = removePermutations(newPoints);
	    return newPoints;
    }
    
    public int[][] removePermutations(int[][] array) {
    	int depth = (array[0].length-1);
    	int numPermutations = 0;
    	//check for common between first in 0 and any place in j
    	for (int m = 0; m < array.length; m++) {//go through each
    	System.out.println("removed for "+m+" out of "+array.length);
	    	for (int j = m+1; j < array.length; j++) {//go through each other one (no doubles!)
	    		//only check first element, necessary condition for permutation is other must contain this element somewhere
				for (int k = 0; k < depth; k++) {//go through each element of others
					if (array[m][0] == array[j][k]) {//test for permutation!
	    				boolean permutation = true;
	    				for (int i = 1; i < depth; i++) {//test other places for sameness
	    					if (array[j][(k+i)%depth] != array[m][i]) permutation = false;
	    				}
	    				if (permutation) {//remove it!
	    					int[][] array2 = new int[array.length-1][];
	    					for (int n = 0; n < (array.length); n++) {
	    						if (n < j) array2[n] = array[n];
	    						if (n > j) array2[n-1] = array[n];
	    					}
	    					array2[m][depth] += array[j][depth];
	    					array = array2;
	    					//System.out.println("permutation! "+m+" and "+j);
	    					numPermutations++;
	    				}
	    			}
				}
			}
		}
		return array;
	}
	public NumberedList getList() {
		return theList;
	}
}
