package noel.util;

import java.util.*;
import noel.io.*;
import ral.*;


/**
 * This class is to help implement an arbitrary depth any-dimension
 * array.  (Just making an array can be too memory intensive, and is not
 * general size).  You give this class a point contained in a 1-D int[].  All 
 * points must be of the same size and equal to depth.  Example: point (3,2,4) 
 * has depth of three.
 */
public class SparseRealArray {
	
	public NumberedList theList; //high level list
	NumberedList list; //utility list
	public int depth; //dimension
	boolean isIndexing; //should we index?
	Vector index;
	int density = 0; //number of distinct sites in the array
	
	public SparseRealArray(int depth, boolean indexing) {
		theList = new NumberedList(0);
		this.depth = depth;
		this.isIndexing = indexing;
		if (indexing) index = new Vector();
	}
	public SparseRealArray(int depth) {
		this(depth,true);
	}	
	//returns the list corresponding to the given point p and turns
	//it into an int[].  Returns null if p.length > depth.  Get the
	//top list with getTopList
    public int[] getSparseList(int[] p) {
        if(p.length > depth) { 
            System.out.println("This sparse list has only "+depth+" dimensions, you asked for "+p.length);
            System.exit(1);
        }
        int count = 0;
        list = theList;
        while (count < p.length) {
            int found = -1;
            for (int i = 0; i < list.size(); i++) {
                NumberedList numL = (NumberedList)list.get(i);
                if (p[count] == numL.num) {
                    found = i;
                    break;
                }
            }
            if (found >=0 ) list = (NumberedList)list.get(found); //next layer
            else {
                int[] a = {-1}; 
                return a;
            }
            count++;
        }
        //turn the NumberedList into an int[]
        int [] returnList = new int[list.size()];
        for(int i=0; i < list.size(); i++) {
            returnList[i] = ((NumberedList)list.get(i)).num;
        }
        Arrays.sort(returnList);
        return returnList;
    }
    //used in conjunction with getSparseList, this returns the int[]
    //of the top list
    public int[] getTopList() {
        int [] returnList = new int[theList.size()];
        for(int i=0; i < theList.size(); i++) {
            returnList[i] = ((NumberedList)(theList).get(i)).num;
        }
        Arrays.sort(returnList);
        return returnList;
    }	
	//returns number of elements
	public int size() {
		return density;
	}
	public void increment(int[] p) {
		increment(p,new Real(1));
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
	* Tells the Real value of this point.
	* @return the value, -1 if p does not exist
	*/
	public Real probeValue(int[] p) {
	    int count = 0;
    	list = theList;
    	while (count < depth) {
            int found = Arrays.binarySearch(list.toArray(),new NumberedList(p[count]));
			if (found >=0 ) list = (NumberedList)list.get(found); //next layer
			else return new Real(-1);
			count++;
    	}
    	return (Real)list.get(0);
	}
	
    //maybe use sorts to help searches...
    public void increment(int[] p, Real amt) {
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
			list.add(new Real(0));
			density++;
			if (isIndexing) {
				int[] pcopy = new int[p.length]; //index it as well!
    			System.arraycopy(p,0,pcopy,0,p.length);
    			index.add(pcopy);
    			list.add(new Integer(index.size()-1)); //put the index into the list as the second element
			} 
		}
		((Real)list.get(0)).add(amt);
		//list.set(0,new Double(i));
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
    
    //since this is doubles, we need an integer index and a separate 
    //list of the actual double values at those indices which you obtain 
    //using getIndexValues()
    public int[][] getIntegerIndex() {
    	int[] point;
    	int[][] newPoints = new int[index.size()][depth];
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
                // int found = -1;
                // for (int i = 0; i < list.size(); i++) {
                //                  NumberedList numL = (NumberedList)list.get(i);
                //                  if (point[count] == numL.num) {
                //                      found = i;
                //                      break;
                //                  }
                //              }
                int found = Arrays.binarySearch(list.toArray(),new NumberedList(point[count]));
    			list = (NumberedList)list.get(found); //next layer
    			count++;
	    	}
	    	//newPoints[row][depth] = ((Integer)list.get(0)).intValue();
	    	row++;
	    }
	    //if (removePermutations) newPoints = removePermutations(newPoints);
	    return newPoints;
    }
    
    //since this is doubles, we need an integer index and a separate 
    //list of the actual double values at those indices which you obtain 
    //using getIndexValues()
    public Real[] getIndexValues() {
    	int[] point;
    	Real[] values = new Real[index.size()];
    	int row = 0;
    	Enumeration points = index.elements();
    	while (points.hasMoreElements()) {
    		point = (int[])points.nextElement();
    		int count = 0;
	    	list = theList;
	    	while (count < depth) {
                // int found = -1;
                // for (int i = 0; i < list.size(); i++) {
                //                  NumberedList numL = (NumberedList)list.get(i);
                //                  if (point[count] == numL.num) {
                //                      found = i;
                //                      break;
                //                  }
                //              }
                int found = Arrays.binarySearch(list.toArray(),new NumberedList(point[count]));
    			
    			list = (NumberedList)list.get(found); //next layer
    			count++;
	    	}
	    	values[row] = ((Real)list.get(0));
	    	row++;
	    }
	    //if (removePermutations) newPoints = removePermutations(newPoints);
	    return values;
    }    
    
	public NumberedList getList() {
		return theList;
	}
	
	/*****************
	* Utility methods
	******************/
	
	//This method switches the order of two indices in the array.  It is used to
	//facilitate quick integrations.  j goes to where i was and i goes to 
	//where j was.  It does this by creating a completely new array.
	public static SparseRealArray switchIndices(SparseRealArray array, int i, int j) {
	    int[][] index = array.getIntegerIndex();
	    Real[] values = array.getIndexValues();
	    SparseRealArray newArray = new SparseRealArray(array.depth);
	    for(int k = 0; k < index.length; k++) {
	        //Utilities.print(index[k],false);
	        int temp = index[k][i];
	        index[k][i] = index[k][j];
	        index[k][j] = temp;
	        //Utilities.print(index[k],false);
	        newArray.increment(index[k],values[k]);
	    }
	    return newArray;
	}
	public static SparseRealArray sumOverLast(SparseRealArray array) {
	    int[][] index = array.getIntegerIndex();
	    Real[] values = array.getIndexValues();
	    SparseRealArray newArray = new SparseRealArray(array.depth-1);
	    for(int k = 0; k < index.length; k++) {
	        int[] point = new int[array.depth-1];
	        for(int p = 0; p < point.length; p++) point[p] = index[k][p]; //grab a subarray
            newArray.increment(point,values[k]);
        }
        return newArray;
    }
	
}
