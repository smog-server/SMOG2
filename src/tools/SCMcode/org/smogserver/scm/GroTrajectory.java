package org.smogserver.scm;

import org.smogserver.io.*;
import org.smogserver.util.math.*;
import org.smogserver.util.*;

/**
 * Parses a GROMACS .gro file.
 */
public class GroTrajectory extends GroGro {

	FileIO file;
	int bufferSize=1;
	Coordinates[][] coordinateSet;
	int frameNum = 0;
	int framesLoaded = 0;
	/**
	 * Parses a GROMACS .gro file with multiple snapshots.  Assumes only one chain. 
	 * @param grofile The gro filename
	 */
	public GroTrajectory(String grofile){
		super(grofile);
		this.file = new FileIO(grofile,FileIO.BUFFERED_READING);
	}
	public GroTrajectory(String grofile, int bufferSize){
		super(grofile);
		frameNum = 1;
		framesLoaded = 1;
		this.file = new FileIO(grofile,FileIO.BUFFERED_READING);
		this.bufferSize = bufferSize;
	}
	/**
	 * Reads the next set of atomic coordinates.  Makes new atom array.
	 */
	public void getNext(){
	    frameNum++;
	    framesLoaded++;
		String line = null; 
		try {
			line = file.readLine(); //comment
			line = file.readLine(); 
			numAtom = Integer.parseInt(line.trim());
			atom = new Atom[numAtom];
			int correction = 0; //for the wrapping of gro atomNumbers due to limited space
			int lastAtomNum = 0;
			for (int i = 0; i < numAtom; i++) {
				line = file.readLine();
				int atomNum = Integer.parseInt(line.substring(15,20).trim());
				if(lastAtomNum>atomNum){correction += (lastAtomNum+1);}
				lastAtomNum = atomNum;
				String id =line.substring(10,15).trim();
				String resName = line.substring(5,10).trim();
				int resNum = Integer.parseInt(line.substring(0,5).trim());
				String[] coor = FileIO.getTokens(line.substring(20,line.length()-1)," ");
				//convert to angstroms
				double x = 10*Double.parseDouble(coor[0]);
				double y = 10*Double.parseDouble(coor[1]);
				double z = 10*Double.parseDouble(coor[2]);
				Coordinates coords = new Coordinates(x,y,z);
				atom[i] = new Atom(id, atomNum+correction, coords, resName, resNum,i);
				atom[i].setChain(1);
			}
			line = file.readLine(); //box size
		} catch (Exception e) {
			System.out.println("Problem at line: \""+line+"\".\nCheck formatting.  Exiting.");
			System.exit(1);
		}
	}
	public Coordinates[] getNextFrameCoordinates() {
	    frameNum++;
	    if(frameNum>framesLoaded){ updateCoords(bufferSize); }
	    return coordinateSet[bufferSize-(framesLoaded-frameNum)-1];
	}
	/**
	 * Reads the next set of atomic coordinates.
	 */
	private void updateCoords(int number){
	    framesLoaded+=number;
		String line = null; 
		//System.gc();
		coordinateSet = new Coordinates[number][numAtom];
		for(int j = 0;j < number; j++) {
    	    
		try {
			line = file.readLine(); //comment
			line = file.readLine(); 
			//numAtom = Integer.parseInt(line.trim());
			//atom = new Atom[numAtom];
			int correction = 0; //for the wrapping of gro atomNumbers due to limited space
			int lastAtomNum = 0;
			for (int i = 0; i < numAtom; i++) {
				line = file.readLine();
				String[] coor = FileIO.getTokens(line.substring(20,line.length()-1)," ");
				//convert to angstroms
				double x = 10*Double.parseDouble(coor[0]);
				double y = 10*Double.parseDouble(coor[1]);
				double z = 10*Double.parseDouble(coor[2]);
				coordinateSet[j][i] = new Coordinates(x,y,z);    		    
				//atom[i] = new Atom(id, atomNum+correction, coords, resName, resNum,i);
				//atom[i].updateCoords(x,y,z);
				//atom[i].setChain(1);
			}
			line = file.readLine(); //box size
		} catch (Exception e) {
			System.out.println("Problem at line: \""+line+"\".\nCheck formatting.  Exiting.");
			System.exit(1);
		}
	    }
	}
	/**
	 * Returns the QWolynes of the current set of atoms relative to input parameter atoms.
	 * @param atoms The atom coordinates to compare.
	 */
	public double getQWolynes(Atom[] atomsNative){
		int numAtoms = atomsNative.length;
		double factor = 2.0 / ((numAtoms - 1)*(numAtoms - 2)); //normalization for sum
		double QW = 0;
		for(int i = 0; i < numAtoms; i++){
			for(int j = 0; j < i - 1; j++){
				double dist = Coordinates.distance(
						atom[i].getCoords(),
						atom[j].getCoords());
				double distNative = Coordinates.distance(
						atomsNative[i].getCoords(),
						atomsNative[j].getCoords());
				//QW += Math.exp( - (dist - distNative)*(dist - distNative) / Math.pow(i - j,0.3) );
				QW += Math.exp( - (dist - distNative)*(dist - distNative) * Math.pow(i - j,0.3) );
			}
		}
		return QW * factor;
	}
	public static void main(String args[]){
		GroTrajectory traj = new GroTrajectory("traj.gro");
		System.out.println(traj.getAtoms().length+" "+traj.getAtoms()[0].getX());
		traj.getNext();
		System.out.println(traj.getAtoms().length+" "+traj.getAtoms()[0].getX());
		traj.getNext();
		System.out.println(traj.getAtoms().length+" "+traj.getAtoms()[0].getX());

	}
}
