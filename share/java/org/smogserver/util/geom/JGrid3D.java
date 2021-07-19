package org.smogserver.util.geom;

import org.smogserver.util.math.*;
import java.util.*;

/*
 * This only works as a static grid for now.  No addition of elements.
 */
public class JGrid3D {

	private double GRIDSIZE;
	private double xmax,xmin,ymax,ymin,zmax,zmin;
	private int XN,YN,ZN;
	private int TOTAL_GRID;
	private ArrayList grid;

	public JGrid3D(JGridable3D[] items, double gridsize) {
		this.GRIDSIZE=gridsize; //gridsize;
		//find extremes
		for(int i=0;i<items.length;i++){
			Coordinates c=items[i].getCoords();
			if(i==0){ //initialize max/min
				xmax=c.x;
				xmin=c.x;
				ymax=c.y;
				ymin=c.y;
				zmax=c.z;
				zmin=c.z;
			}
			if(xmax<c.x)xmax=c.x;
			if(xmin>c.x)xmin=c.x;
			if(ymax<c.y)ymax=c.y;
			if(ymin>c.y)ymin=c.y;
			if(zmax<c.z)zmax=c.z;
			if(zmin>c.z)zmin=c.z;
		}
		//System.out.println(xmin+" "+xmax+" "+ymin+" "+ymax+" "+zmin+" "+zmax+" q");
		double xw=xmax-xmin;
		double yw=ymax-ymin;
		double zw=zmax-zmin;
		//Ensure that Paul doesn't try to give a ridiculous system that breaks
		//the gridding by increasing the gridsize to compensate
		long testGridNum = 1000000000;
		int timesthrough = 0;
		while(testGridNum > 100000000) {
			if(timesthrough>0) {
				System.out.println("Too many grids: doubling gridsize from "+gridsize+" to "+(gridsize*2));
				gridsize=gridsize*2;
			}
			XN=(int)(xw/gridsize)+1;
			YN=(int)(yw/gridsize)+1;
			ZN=(int)(zw/gridsize)+1;
			testGridNum = ((long)XN)*((long)YN)*((long)ZN);
			timesthrough++;
		}
		this.GRIDSIZE = gridsize;
		TOTAL_GRID=XN*YN*ZN;
		//System.out.println(XN+" "+YN+" "+ZN+" Total grids= "+TOTAL_GRID);

		grid = new ArrayList(0);
		//make grid
		for(int i=0;i<TOTAL_GRID;i++){
			grid.add(new ArrayList(0));
		}
		//add items 
		for(int i=0;i<items.length;i++){
			ArrayList list = (ArrayList)grid.get(getGridIndex(items[i].getCoords()));
			list.add(items[i]);
		}
		//for(int i=0;i<TOTAL_GRID;i++){
			//ArrayList list = (ArrayList)grid.get(i);
			//System.out.println(i+" "+((ArrayList)grid.get(i)).size());
			//for(int j=0;j<list.size();j++){
				//System.out.println(((Atom)list.get(j)).getCoords());
			//}
		//}
		//System.out.println(((Atom)(((ArrayList)grid.get(5)).get(0))).getCoords());
		//Vector friends = getGrid(((Atom)(((ArrayList)grid.get(5)).get(0))).getCoords());
		//Atom[] f = (Atom[])friends.toArray();
		//System.out.println(friends.toArray()[4]);
		//for(int i=0;i<friends.size();i++){
			//System.out.println(((Atom)friends.get(i)).getCoords());
		//}
		//Vector v = getNeighbors(((Atom)(((ArrayList)grid.get(5)).get(0))).getCoords());
		//Vector v = getNeighbors(((Atom)(((ArrayList)grid.get(5)).get(0))).getCoords());
		//for(int i=0;i<v.size();i++){
			//System.out.println(((Atom)v.get(i)));
		//}

	}//end of constructor

	/*
	 * Returns all the members of the grid containing c.
	 */
	public Vector getGrid(Coordinates c) {
		int index = getGridIndex(c);
		return getGrid(index); 
	}
	/*
	 * Returns grid with number index.
	 */
	public Vector getGrid(int index) {
		return new Vector((Collection)grid.get(index));
	}
	/*
	 * Returns grid and surrounding grids of grid containing c. 
	 */
	public Vector getNeighbors(Coordinates c) {
		int[] indices = getNeighborIndices(c);
		Vector neigh = new Vector();
		int index = 0;
		while(indices[index] >= 0) {
			neigh.addAll(getGrid(indices[index]));
//System.out.print(" |"+index+" "+indices[index]);
			index++;
		}
		return neigh;
	}
	/*
	 * Returns the number of grids
	 */
	public long getNumGrids() {
		return TOTAL_GRID;
	}
	/*
	 * Gives grid index for a coordinate
	 */
	private int getGridIndex(Coordinates c) {
		int xi,yi,zi;
		xi=(int)((c.x-xmin)/GRIDSIZE);
		yi=(int)((c.y-ymin)/GRIDSIZE);
		zi=(int)((c.z-zmin)/GRIDSIZE);
		return xi*YN*ZN+yi*ZN+zi;
	}

	private int[] getNeighborIndices(Coordinates c) {
		int[] indices = new int[28];
		for(int i=0;i<28;i++)indices[i]=-1; //initialize with negative
		int count=0;
		int xi,yi,zi;
		xi=(int)((c.x-xmin)/GRIDSIZE);
		yi=(int)((c.y-ymin)/GRIDSIZE);
		zi=(int)((c.z-zmin)/GRIDSIZE);
		for(int x=xi-1;x<=xi+1;x++){
			if(x>=0 && x<XN){
				for(int y=yi-1;y<=yi+1;y++){
					if(y>=0 && y<YN){
						for(int z=zi-1;z<=zi+1;z++){
							if(z>=0 && z<ZN){
								indices[count++]= x*YN*ZN+y*ZN+z;
							}
						}
					}
				}
			}
		}
		return indices;
	}


}//end of class
