package org.smogserver.wham;

import java.io.*;
import java.util.*;

import argparser.*;
import org.smogserver.io.*;
import org.smogserver.util.*;
import ral.*;

public class WhamMain {
	
	static String versionNum = "1.10";
	static boolean distribution = true;
	static final double SMALL_SHIFT = 0.000001; //makes the binning not insane
	static String distString = ""+
	"*****************************************************************\n"+
    "This software is part of SMOG. Direct questions to: info@smog-server.org\n"+
	"Documentation available in the SMOG2 manual.\n"+
    "Work utilizing SMOG should cite:\n\n"+
    "Noel JK, Levi M, Raghunathan M, Lammert H, Hayes R, Onuchic JN, and Whitford PC. (2016)\n"+
	"SMOG V2, A Versatile Software Package for Generating Structure-Based Models.\n"+
	"PLoS Comput Biol 12(3): e1004794. doi:10.1371/journal.pcbi.1004794\n"+
    "*****************************************************************";

	//Configuration parameters
	static IntHolder numDim,numFile,maxIterations,ntemps,ntempsF,ntempsC,setToZero,numThreads,numUmb,numBinProb;
	static DoubleHolder kb,tolerance,startT,startTF,startTC,deltaT,deltaTF,deltaTC,lowQ2,highQ2,startProb,stepProb;
	static Vector<LongHolder> numBin;
	static Vector<StringHolder> fileNames, umbrellaType;
	static Vector<DoubleHolder> fileTemps, start, step, umbrella_k, umbrella_sigma, umbrella_0, umbrella_0_2;
	static StringHolder config, dosFileName, run_cv_out, run_free_out, run_coord_out, free_energy_out, run_prob_out, convergenceType, convergenceCriteria;
	static BooleanHolder run_free,run_cv,run_coord,run_wham,readFreeEnergy,reweighting,run_FEP,run_prob,overwriting;
	
	//Useful global declarations
	static Vector<SmartBin> shist;
	static double[] scaledTemps;
	static int actualDim = 0;
	static int[] numBins;
	static Real[] nconf;
	
	public static void main(String[] args) {
        setUpConfiguration(args);
        int numU = numUmb.value;
        int numD = numDim.value;
        if(reweighting.value) { numD--; } //last column is the weights, forget about it
  	    if(run_wham.value) {
	        if(numU > 3) {
	            System.out.println("Only supports at most 3 umbrellas at a time at the moment, sorry!  numUmb = "+numU);
                System.exit(1);
	        }
	        if(numD-numU > 4) {
	            System.out.println("You really want to run wham with 4 order parameters (numDim >4 w/o umbrellas)? Not gonna happen.");
                System.exit(1);
	        }
	        System.out.println("Running high precision WHAM calculation.");
	        run_wham(numD-1-numU,numU);
	    }
        if(run_cv.value) { 
            SparseRealArray dos = readDOS(numD-numU);
            System.out.println("Printing specific heat.\n");
            for(int i = 0; i < (numD-1-numU); i++) {
                dos = SparseRealArray.sumOverLast(dos);
            }
            run_cv(dos);
        }
		if(run_FEP.value) {
            SparseRealArray dos = readDOS(numD-numU);
            System.out.println("Printing free energy perturbation. Will print Delta F.\n");
            if(numD - numU == 3) {
                run_FEP(dos);
            } else {
                System.out.println("Not sure what to do with >3 or <3 dimensions (not including umbrella) and FEP.  Need E,Q1,Q2.  Quitting.");
                System.exit(1);
            }
		}
        if(run_free.value) { 
            SparseRealArray dos = readDOS(numD-numU);
            System.out.println("Printing free energy vs order parameter(s).\n");
            if(numD - numU == 3) {
                run_free_2D(dos);
            } else if(numD - numU == 2) { 
                run_free(dos);
            } else {
                System.out.println("numD - numU = "+(numD - numU)+". Not sure what to do with >3 dimensions (not including umbrella) and free energy calculation.  Only need E,Q1,Q2.  Quitting.");
                System.exit(1);
            }
        }
        if(run_coord.value) {
            SparseRealArray dos = readDOS(numD-numU);
            System.out.println("Printing coordinate values.\n");
            if(numD - numU == 2) {
                run_coord(dos);
            } else if(numD - numU == 3) {
                run_coord_2D(dos);
            } else if(numD - numU == 4) {
				if(lowQ2.value == -314) { //default null value
					run_full_coord_3D(dos);
				} else { //user set lowQ2 and highQ2
					run_coord_3D(dos,lowQ2.value,highQ2.value);                    
				}
            } else {
                System.out.println("Not sure what to do with (numDim - numUmb != 2,3,4) and coordinate calculation. Quitting.");
                System.exit(1);
            }
        }
		if(run_prob.value) {
            SparseRealArray dos = readDOS(numD-numU);
			System.out.println("\nRegurgitating reweighted histograms.");
            if(numD - numU == 2) {
                run_prob(dos);
            } else if(numD - numU == 3) {
				if(lowQ2.value == -314) { //default null value
					if(startProb.value != -1 && stepProb.value != -1 && numBinProb.value != -1) {
						run_prob_heiko_2D(dos);
					} else {
		                System.out.println("Must set a Q2 range, lowQ2/highQ2, for run_prob and 3 dimensions, or startProb,stepProb,numBinProb.");
		                System.exit(1);
					}
				} else {
	            	run_prob_2D(dos,lowQ2.value,highQ2.value);
				}
			}
				
		}
        if(distribution) System.out.println(distString+"\nThank you, wham again.");
        else System.out.println("\nThank you, wham again.");
	}
	
	//Wrapper method for creating a FileIO in order to do some error checking
	private static FileIO getOutputFileIO(String name) {
		if(overwriting.value) return new FileIO(name,FileIO.WRITING);
		else { //will die if file already exists
			FileIO file = new FileIO(name,FileIO.NO_OVER_WRITING);
			if(file.alreadyExisted) {
				System.out.println("Error: A file named '"+name+ "' already exists. "+
					"Delete file, or add option 'overwriting' to config file or to command line.");
				System.exit(1);
			}
			return file;
		}
	}	

    /*********************************************
    * WHAM algorithm (outputs density of states)
    **********************************************/
		
	private static void run_wham2(int numOrderParams, int numUmbrella) {
        actualDim = numOrderParams + 1 + numUmbrella; //+1 because of energy
        shist = readSkinnyHistograms(actualDim);
        //Run WHAM
        int nfile = numFile.value;
        numBins = new int[actualDim];
        for (int i = 0; i < actualDim; i++) numBins[i] = (int)numBin.elementAt(i).value;
        //in histogram: energy 0, then orders, then umbrella last
        nconf = Utilities.initializeWithZeros(new Real[nfile]);
        Vector intIndices = new Vector<int[][]>();	    
        Vector doubleValues = new Vector<double[]>();
        for(int i = 0; i < nfile; i++) {
            intIndices.add(shist.elementAt(i).getIntegerIndex());
            doubleValues.add(shist.elementAt(i).getIndexValues());
        }
        System.out.println("\nSum total configurations.");
        for(int i = 0; i < nfile; i++) {
            for (int k = 0; k < ((int[][])intIndices.elementAt(i)).length; k++) {
                nconf[i].add( new Real ( new Double( ((double[])doubleValues.elementAt(i))[k]).toString() )); 
            }
        }
        for(int i = 0; i < nfile; i++) {
            System.out.println(fileNames.elementAt(i).value+" "+nconf[i]);
        }
        boolean done = false;
        System.out.println("\nBuilding histograms...");
        int iteration = 1;
		//N(i,j,k) is something, come back to this
		//lets make E(i,j,k)
		//E(i,j,k) = E_i + bias(i,j,k)
		
	}	    
    private static void run_wham(int numOrderParams, int numUmbrella) {
        actualDim = numOrderParams + 1 + numUmbrella; //+1 because of energy
        shist = readSkinnyHistograms(actualDim);
        //Run WHAM
        int nfile = numFile.value;
        numBins = new int[actualDim];
        for (int i = 0; i < actualDim; i++) numBins[i] = (int)numBin.elementAt(i).value;
        //in histogram: energy 0, then orders, then umbrella last
        nconf = Utilities.initializeWithZeros(new Real[nfile]);
        Vector intIndices = new Vector<int[][]>();	    
        Vector doubleValues = new Vector<double[]>();
        for(int i = 0; i < nfile; i++) {
            intIndices.add(shist.elementAt(i).getIntegerIndex());
            doubleValues.add(shist.elementAt(i).getIndexValues());
        }
        System.out.println("\nSum total configurations.");
        for(int i = 0; i < nfile; i++) {
            for (int k = 0; k < ((int[][])intIndices.elementAt(i)).length; k++) {
                nconf[i].add( new Real ( new Double( ((double[])doubleValues.elementAt(i))[k]).toString() )); 
            }
        }
        for(int i = 0; i < nfile; i++) {
            System.out.println(fileNames.elementAt(i).value+" "+nconf[i]);
        }
        boolean done = false;
        System.out.println("\nBuilding histograms...");
        int iteration = 1;
        double[] zold = Utilities.initializeWithZeros(new double[nfile]);
        if(readFreeEnergy.value) readFree(zold); //read in past energies
        double[] znew = Utilities.initializeWithZeros(new double[nfile]);
        double[] dz0 = Utilities.initializeWithZeros(new double[nfile]);
        double[] dz1 = Utilities.initializeWithZeros(new double[nfile]);
        double[] z00 = Utilities.initializeWithZeros(new double[nfile]);

        //fill energies
        double[] energy = Utilities.initializeWithZeros(new double[numBins[0]]);
        for(int j = 0; j < numBins[0]; j++) { //first column is the energy
            energy[j] = start.elementAt(0).value + step.elementAt(0).value * (j+0.5);
        }
		//This will contain precomputed bias energies along each umbrella coordinate
		//To determine total bias for a particular histogram bin, one must sum the bias energies for each umbrella
		//indexed by the umbrella coordinates for the bin
		double [][][] biases = new double[nfile][numUmbrella][];
		//FOR NOW WE ASSUME EVERY FILE HAS SAME NUMBER OF UMBRELLAS
		try {
			for(int i = 0; i < nfile; i++) {
				double strength;
				double center;
				int umbrellaColumn;
				for(int j = 0; j < numUmbrella; j++) {
					strength = umbrella_k.elementAt(i*numUmbrella+j).value;
					center = umbrella_0.elementAt(i*numUmbrella+j).value;	
					umbrellaColumn = actualDim-(numUmbrella-j);		
					biases[i][j] = new double[numBins[umbrellaColumn]];
					if(umbrellaType.elementAt(j).value.compareTo("harmonic")==0) {
						for( int k = 0; k < biases[i][j].length; k++) {
							double umbrellaCoordValue = start.elementAt(umbrellaColumn).value + step.elementAt(umbrellaColumn).value * (k+0.5);
							biases[i][j][k] = 0.5 * strength * Math.pow(umbrellaCoordValue-center,2);
						}
					} else if(umbrellaType.elementAt(j).value.compareTo("linear")==0) {
						for( int k = 0; k < biases[i][j].length; k++) {
							double umbrellaCoordValue = start.elementAt(umbrellaColumn).value + step.elementAt(umbrellaColumn).value * (k+0.5);
							biases[i][j][k] = strength * (umbrellaCoordValue-center);
						}
					} else {
		                System.out.println("Need to have an 'umbrellaType TYPE' for each umbrella column (i.e. numUmbrella). TYPE=[harmonic,linear] is currently allowed.\n");
		                System.exit(1);
					}
				}
			}
		} catch (Exception e) {
            e.printStackTrace();
            System.out.println("Something wrong with your umbrella input.");
            for(int i = 0; i < numUmbrella; i++) {
                System.out.println("You chose 'umbrellaType "+umbrellaType.elementAt(i).value+"'. ");
            }
            System.out.println("Did you define umbrella_0 and umbrella_k for each umbrella in each input histogram?");
            System.out.println("Quitting.");
            System.exit(1);
        }
        //Sum super histogram. This is just summing up equivalent bins from all histograms
        //straight up (no reweighting by energy or umbrellas).  It accomplishes this by going sequentially
        //through all the histogram indices, grabbing the next bin, adds them all together.
        SparseDoubleArray summed = new SparseDoubleArray(actualDim);
        for(int i = 0; i < nfile; i++) {
            //System.out.println(i);
            int[][] index = (int[][])intIndices.elementAt(i); //list of histogram i
            double[] values = (double[])doubleValues.elementAt(i);
            for (int k = 0; k < index.length; k++) {
                summed.increment(index[k],values[k]);
            }
        }
        int[][] summedIndex = summed.getIntegerIndex();
        double[] summedValues = summed.getIndexValues();
        Real[] summedValuesReal = Utilities.initializeWithZeros(new Real[summedValues.length]);
        for (int k = 0; k < summedValues.length; k++) {
            summedValuesReal[k] = new Real ( (new Double ( summedValues[k] )).toString() );
        }
        shist = null;
        intIndices = null;
        System.out.println("\nStarting iterations...");
        while (!done) {
            Real[] z0 = Utilities.initializeWithZeros(new Real[nfile]); //partition sum holder
            int threadLength = summedIndex.length;
            Real[][] z0part = Utilities.initializeWithZeros(new Real[nfile][numThreads.value]); //partition sum holder
			
            //*************** THREAD ************************
            for(int th = 0; th < numThreads.value; th++) {
                (new Thread(new InnerThread(summedIndex,summedValuesReal,energy,z0part,biases,zold,nfile,numUmbrella,(threadLength*(th))/numThreads.value,
                                                (threadLength*(th+1))/numThreads.value-1,th))).start();   
            }
			boolean notDone = true;
            while(notDone) {
				notDone = false;
                for(int th = 0; th < numThreads.value; th++) {
					notDone = notDone || !InnerThread.threads[th];
				} 
                try { Thread.sleep(3); }
                catch (InterruptedException e) { e.printStackTrace(); }
            }
            //*************** END THREAD ********************
            
            for(int i = 0; i < nfile; i++) {
                for(int th = 0; th < numThreads.value; th++) {
                    z0[i].add(z0part[i][th]); //add the parts back together (solves concurrency issues)
                }
            }
            //calculate new free energies
            for(int i = 0; i < nfile; i++) {
                if(z0[i].greaterThan(0)) {
          	        z0[i].ln();
          	        znew[i] = -scaledTemps[i]*(new Double(z0[i].toString())).doubleValue();
      	        } else { 
          	        System.out.println("z0 = "+z0[i]+" in file "+i); znew[i] = 0; 
          	    }
            }
            if(setToZero.value > 0 && setToZero.value < fileNames.size()) znew[setToZero.value] = 0; //pin this free energy to zero (in order the files list)
            System.out.println("******* ITERATION NUMBER "+iteration+" ************");
            for(int i = 0; i < nfile; i++) { 
                System.out.println(znew[i]+" "+zold[i]+" "+(znew[i]-zold[i]));
            }
            //test convergence
            double delta = 0;
			double maxdelta = 0;
            for(int i = 0; i < nfile; i++) { 
                dz1[i] = znew[i] - zold[i];
                delta += Math.abs(dz1[i]);
				if(Math.abs(dz1[i]) > maxdelta) {
					maxdelta = Math.abs(dz1[i]);
				}
            }	        
            if(convergenceType.value.compareTo("simple")==0) {
                for(int i = 0; i < nfile; i++) { 
                    zold[i] = znew[i];
                }
            } else if(convergenceType.value.compareTo("quick")==0){
                double cnum = 0, cden = 0, c = 0;
                 for(int i = 0; i < nfile; i++) { 
                     cnum += dz1[i]*(dz1[i] - dz0[i]);
                     cden += Math.pow(dz1[i] - dz0[i],2);
                 }
                 c = cnum / cden;
                 for(int i = 0; i < nfile; i++) { 
                     if(iteration <= 2) zold[i] = znew[i];                
                     else zold[i] = (1-c)*znew[i] + c*z00[i];
                     dz0[i] = dz1[i];
                     z00[i] = znew[i];
                 }
            } else {
                for(int i = 0; i < nfile; i++) { 
                    zold[i] = znew[i];
                }
            }
            iteration++;
			if(iteration >= maxIterations.value) {done = true;}
			else if(iteration > 2) {
				if(convergenceCriteria.value.compareTo("average")==0) {
		            System.out.println(delta/nfile+" is attempting to converge below "+tolerance.value);
					if(delta/nfile < tolerance.value) {done = true;}
				} else if(convergenceCriteria.value.compareTo("max")==0) {
		            System.out.println(maxdelta+" is attempting to converge below "+tolerance.value);
					if(maxdelta < tolerance.value) {done = true;}
				} else { //use average
		            System.out.println(delta/nfile+" is attempting to converge below "+tolerance.value);
					if(delta/nfile < tolerance.value) {done = true;}
				}
			}	
        }
        //output results!
        System.out.println("FINISHED ITERATIONS.  Preparing DOS."); 
        if(readFreeEnergy.value) { writeFree(zold); }
        SparseRealArray dos = new SparseRealArray(numOrderParams+1,true);

		int dim = numUmbrella + 1;
		int[] indices = new int[dim];
		int[] umbrellaIndices = new int[1]; //intialization for compilation
		if(numUmbrella>0) umbrellaIndices = new int[numUmbrella];
		int[] lengths = new int[dim];
		lengths[0] = numBins[0];
		for(int i = 1; i < dim; i++ ) { lengths[i] = numBins[actualDim-(dim-i)]; }
		double[] biasEnergy = Utilities.initializeWithZeros(new double[nfile]);
		HashMap zzRecipH = new HashMap<Long,Real>();
		
        for(int k = 0; k < summedIndex.length; k++) {
            int energyIndex = summedIndex[k][0];
            int umbrellaIndex1 = summedIndex[k][actualDim-2];
            int umbrellaIndex2 = summedIndex[k][actualDim-1];
			indices[0] = summedIndex[k][0];
			for(int i = 0; i < numUmbrella; i++ ) { 
				int umbrellaColumn = actualDim-(dim-1-i);
				indices[i+1] = summedIndex[k][umbrellaColumn]; 
				umbrellaIndices[i] = summedIndex[k][umbrellaColumn]; 
			}
			long hashindex = get1Dindex(indices,lengths,dim);
			if(numUmbrella > 0) { //there are umbrellas
				for(int i = 0; i < nfile; i++) {
					biasEnergy[i] = getBias(biases,i,umbrellaIndices);
				}
			}
            if(!zzRecipH.containsKey(hashindex)) { //dont calculate more than once!         
				zzRecipH.put(hashindex,new Real(0));      
                for(int i = 0; i < nfile; i++) {
                    Real exponent = new Real( (new Double( ( -energy[energyIndex] - biasEnergy[i] + zold[i] )/scaledTemps[i] )).toString() );
                                 exponent.exp();
                                 exponent.mul(nconf[i]);
								 ((Real)zzRecipH.get(hashindex)).add(exponent);
                }
				((Real)zzRecipH.get(hashindex)).recip();
            }
            int[] point = new int[numOrderParams+1]; 
            for(int p = 0; p < point.length; p++) point[p] = summedIndex[k][p]; //grab a subarray
            Real temp = new Real(summedValuesReal[k]);
            temp.mul((Real)zzRecipH.get(hashindex));
            dos.increment(point,temp);
        }
        writeDOS(dos,numOrderParams+1);
    } //end run_wham

    /*********************************************
    * Thermodynamic calculations (reads a dos)
    **********************************************/
    private static void run_free(SparseRealArray dos) {
            //problem: we have a dos(E,Q) but we need dos(Q,E) for this integral.  So switch it!
            dos = SparseRealArray.switchIndices(dos,0,1);
            int[] olist = dos.getTopList();
            //Utilities.print(olist,false); //abc
            for(double t = startTF.value ; t < startTF.value + ntempsF.value*deltaTF.value; t+=deltaTF.value) {
            // System.out.println("\nPrinting free energy vs order parameter."+startTF.value+" "+ntempsF.value+" "+deltaTF.value+" "+(ntempsF.value*deltaTF.value));
                FileIO free = getOutputFileIO(run_free_out.value+((int)(t*10+SMALL_SHIFT)));
                double tt = t*kb.value;
                Real ttR = new Real( (new Double(-tt)).toString() );
                Real min = new Real();
                //Do it once to get the minimum value of free energy so you can write it out
                for(int k = 0; k < olist.length; k++) {
                    Real pp = new Real(); 
                    Real ep = new Real();
                    Real pmf = new Real();
                    int[] p = {olist[k]};
                    int[] elist = dos.getSparseList(p);                    
                    for(int j = 0; j < elist.length; j++) {
                        int[] p1 = {olist[k],elist[j]};
                        double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[j]));
                        Real exponent = new Real( (new Double(-energy/tt)).toString() );
                        exponent.exp();
                        exponent.mul(dos.probeValue(p1));                    
                        pp.add(exponent);
                    }
                    pmf.assign(pp);
                    pmf.ln();
                    pmf.mul(ttR); 
                    if(k==0 ) min.assign(pmf);
                    if(min.greaterThan(pmf)) min.assign(pmf); //grabs the minimum value
                } 
                for(int k = 0; k < olist.length; k++) {
                    double order = start.elementAt(1).value + step.elementAt(1).value * (olist[k]);
                    Real pp = new Real(); 
                    Real ep = new Real();
                    Real pmf = new Real();
                    Real enth = new Real();
                    Real ts = new Real();
                    int[] p = {olist[k]};
                    int[] elist = dos.getSparseList(p);                    
                    for(int j = 0; j < elist.length; j++) {
                        double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[j]));
                        int[] p1 = {olist[k],elist[j]};
                        Real energyR = new Real( (new Double(energy)).toString() );
                        Real exponent = new Real( (new Double(-energy/tt)).toString() );
                        exponent.exp();
                        exponent.mul(dos.probeValue(p1));                    
                        pp.add(exponent);
                        exponent.mul(energyR);
                        ep.add(exponent);
                    }
                    if(pp.greaterThan(0)) {
                        pmf.assign(pp);
                        pmf.ln();
                        pmf.mul(ttR);
                        enth.assign(ep);
                        enth.div(pp);
                        pmf.sub(min);
                        ts.assign(enth);
                        ts.sub(pmf);
                        free.write(order+" "+pmf+" "+enth+" "+ts+"\n");
                    }
                }               
                free.close();
            }      
    }//end run_free
    private static void run_free_2D(SparseRealArray dos) {
	    int nfile = numFile.value;
        int nenergy = (int)numBin.elementAt(0).value;
        int norder1 = (int)numBin.elementAt(1).value;            
        int norder2 = (int)numBin.elementAt(2).value;
        //problem: we have a dos(E,Q1,Q2) but we need dos(Q1,Q2,E) for this integral.  So switch it!
        dos = SparseRealArray.switchIndices(dos,0,1); // dos(Q1,E,Q2)
        dos = SparseRealArray.switchIndices(dos,1,2);
        for(double t = startTF.value ; t < startTF.value + ntempsF.value*deltaTF.value; t+=deltaTF.value) {
    	   // System.out.println("\nPrinting free energy vs order parameter."+startTF.value+" "+ntempsF.value+" "+deltaTF.value+" "+(ntempsF.value*deltaTF.value));
            FileIO free = getOutputFileIO(run_free_out.value+((int)(t*10+SMALL_SHIFT)));
            double tt = t*kb.value;
            Real ttR = new Real( (new Double(-tt)).toString() );
            Real min = new Real();
            //Do it once to get the minimum value of free energy so you can write it out
            int[] o1list = dos.getTopList();
            for(int k = 0; k < o1list.length; k++) {
                int[] p = {o1list[k]};
                int[] o2list = dos.getSparseList(p);
                for(int l = 0; l < o2list.length; l++) {                
                    Real pp = new Real(); 
                    Real ep = new Real();
                    Real pmf = new Real();
                    int[] p1 = {o1list[k],o2list[l]};
                    int[] elist = dos.getSparseList(p1);                    
                    for(int j = 0; j < elist.length; j++) {
                        int[] p2 = {o1list[k],o2list[l],elist[j]};
                        double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[j]));
                        Real exponent = new Real( (new Double(-energy/tt)).toString() );
                        exponent.exp();
                        exponent.mul(dos.probeValue(p2));                    
                        pp.add(exponent);
                    }
                    pmf.assign(pp);
                    pmf.ln();
                    pmf.mul(ttR); 
                    if(k==0 ) min.assign(pmf);
                    if(min.greaterThan(pmf)) min.assign(pmf); //grabs the minimum value
                }
            } 
            for(int k = 0; k < o1list.length; k++) {
                double order1 = start.elementAt(1).value + step.elementAt(1).value * (o1list[k]+0.5);
                int[] p = {o1list[k]};
                int[] o2list = dos.getSparseList(p);
                for(int l = 0; l < o2list.length; l++) {     
                    double order2 = start.elementAt(2).value + step.elementAt(2).value * (o2list[l]+0.5);           
                    Real pp = new Real(); 
                    Real ep = new Real();
                    Real pmf = new Real();
                    Real enth = new Real();
                    Real ts = new Real();
                    int[] p1 = {o1list[k],o2list[l]};
                    int[] elist = dos.getSparseList(p1);                    
                    for(int j = 0; j < elist.length; j++) {
                        double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[j]+0.5));
                        int[] p2 = {o1list[k],o2list[l],elist[j]};
                        Real energyR = new Real( (new Double(energy)).toString() );
                        Real exponent = new Real( (new Double(-energy/tt)).toString() );
                        exponent.exp();
                        exponent.mul(dos.probeValue(p2));                    
                        pp.add(exponent);
                        exponent.mul(energyR);
                        ep.add(exponent);
                    }
                    pmf.assign(pp);
                    pmf.ln();
                    pmf.mul(ttR);
                    enth.assign(ep);
                    enth.div(pp);
                    pmf.sub(min);
                    ts.assign(enth);
                    ts.sub(pmf);
                    free.write(order1+" "+order2+" "+pmf+" "+enth+" "+ts+"\n");
                }
            }               
            free.close();
        }
	}//end run_free_2D
    private static void run_FEP(SparseRealArray dos) { 
		//Q2 is the perturbing coord, assuming perturbation is actually Q2 in energy units
        //problem: we have a dos(E,Q1,Q2) but we need dos(Q1,Q2,E) for this integral.  So switch it!
        dos = SparseRealArray.switchIndices(dos,0,1); // dos(Q1,E,Q2)
        dos = SparseRealArray.switchIndices(dos,1,2);
        for(double t = startTF.value ; t < startTF.value + ntempsF.value*deltaTF.value; t+=deltaTF.value) {
    	   // System.out.println("\nPrinting free energy vs order parameter."+startTF.value+" "+ntempsF.value+" "+deltaTF.value+" "+(ntempsF.value*deltaTF.value));
            FileIO free = getOutputFileIO(run_free_out.value+((int)(t*10+SMALL_SHIFT)));
            double tt = t*kb.value;
            Real ttR = new Real( (new Double(-tt)).toString() );
            int[] o1list = dos.getTopList();
            for(int k = 0; k < o1list.length; k++) {
                double order1 = start.elementAt(1).value + step.elementAt(1).value * (o1list[k]+0.5);
                int[] p = {o1list[k]};
                int[] o2list = dos.getSparseList(p);
                Real pert = new Real(); 
				Real pertNorm = new Real();
                for(int l = 0; l < o2list.length; l++) { //our perturbing coordinate
                    double order2 = start.elementAt(2).value + step.elementAt(2).value * (o2list[l]+0.5);           
                    //Real ep = new Real();
                    //Real pmf = new Real();
                    //Real enth = new Real();
                    //Real ts = new Real();
                    int[] p1 = {o1list[k],o2list[l]};
                    int[] elist = dos.getSparseList(p1);                    
                    for(int j = 0; j < elist.length; j++) {
                        double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[j]+0.5));
                        int[] p2 = {o1list[k],o2list[l],elist[j]};
                        Real energyR = new Real( (new Double(energy)).toString() );
                        Real exponent = new Real( (new Double( (- energy - order2)/tt )).toString() );
						Real exponentNorm = new Real( (new Double( (- energy )/tt )).toString() );
                        exponent.exp();
						exponentNorm.exp();
                        exponent.mul(dos.probeValue(p2));
						exponentNorm.mul(dos.probeValue(p2));
                        pert.add(exponent);
						pertNorm.add(exponentNorm);
                    }
                }
				pert.ln();
				pert.mul(ttR);
				pertNorm.ln();
				pertNorm.mul(ttR);
				pert.sub(pertNorm);
                free.write(order1+" "+pert+"\n");
            }               
            free.close();
        }
	}//end run_free_2D
    private static void run_cv(SparseRealArray dos) {
        //calculate CV
        FileIO cv = getOutputFileIO(run_cv_out.value);
        for(double t = startT.value ; t < startT.value + ntemps.value*deltaT.value; t+=deltaT.value) {
            Real enth = new Real(); 
            Real psum = new Real();
            Real sigma_e = new Real();
            double tt = t*kb.value;
            int[] list = dos.getTopList();
            for(int j = 0; j < list.length; j++) {
                double energy = start.elementAt(0).value + step.elementAt(0).value * (list[j]+0.5);
                Real energyR = new Real( (new Double(energy)).toString() );
                Real exponent = new Real( (new Double(-energy/tt)).toString() );
                exponent.exp();
                Real pp = new Real();
                int[] p = {list[j]};
                pp.add(dos.probeValue(p));
                pp.mul(exponent);
                Real energyHolder = new Real(energyR);
                energyR.mul(pp);
                enth.add(energyR);
                energyR.assign(energyHolder);
                energyR.sqr();
                energyR.mul(pp);
                sigma_e.add(energyR);
                psum.add(pp);
            }
            enth.div(psum);
            sigma_e.div(psum);
            Real enth2 = new Real(enth);
            enth2.sqr();
            sigma_e.sub(enth2);
            psum.ln();
            Real ttR = new Real( (new Double(-tt)).toString() );
            psum.mul(ttR);
            sigma_e.div(ttR);
            sigma_e.div(ttR);
            Real ent = new Real(enth);
            ent.sub(psum);
            ent.div(ttR);
            cv.write(t+" "+sigma_e+" "+enth+" "+psum+" "+ent+"\n");
        }
	}
    /****
    * computes Q(T) from a density of states assuming a dos(E,Q)
    */
    private static void run_coord(SparseRealArray dos) {
        FileIO coord = getOutputFileIO(run_coord_out.value);
        //problem: we have a dos(E,Q) but we need dos(Q,E) for this integral.  So switch it!
        dos = SparseRealArray.switchIndices(dos,0,1); // dos(Q,E)
        for(double t = startTC.value ; t < startTC.value + ntempsC.value*deltaTC.value; t+=deltaTC.value) {
            double tt = t*kb.value;
            int[] olist = dos.getTopList();
            Real pp = new Real();
            Real ep = new Real();
            for(int j = 0; j < olist.length; j++) {
                double order = start.elementAt(1).value + step.elementAt(1).value * (olist[j]+0.5);
                Real orderR = new Real( (new Double(order)).toString() );
                int[] p = {olist[j]};
                int[] elist = dos.getSparseList(p);
                for(int k = 0; k < elist.length; k++) {
                    int[] p1 = {olist[j],elist[k]};
                    double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[k]+0.5));
                    Real energyR = new Real( (new Double(energy)).toString() );
                    Real exponent = new Real( (new Double(-energy/tt)).toString() );
                    exponent.exp();
                    exponent.mul(dos.probeValue(p1));
                    pp.add(exponent);
                    exponent.mul(orderR);
                    ep.add(exponent);
                }
            }
            ep.div(pp);
            coord.write(t+" "+ep+"\n");
        }
    } //end run_coord
    /****
    * computes prob[Q,T](Q) from a density of states assuming a dos(E,Q). Basically regurgitates the underlying histogram.
    */
    private static void run_prob(SparseRealArray dos) {
        //problem: we have a dos(E,Q) but we need dos(Q,E) for this integral.  So switch it!
        dos = SparseRealArray.switchIndices(dos,0,1); // dos(Q,E)

        FileIO coord = getOutputFileIO(run_prob_out.value+startTC.value);
        double tt = startTC.value*kb.value;
        int[] olist = dos.getTopList();
        //Real ep = new Real();
        for(int j = 0; j < olist.length; j++) {
	        Real pp = new Real();
            double order = start.elementAt(1).value + step.elementAt(1).value * (olist[j]+0.5);
            Real orderR = new Real( (new Double(order)).toString() );
            int[] p = {olist[j]};
            int[] elist = dos.getSparseList(p);
            for(int k = 0; k < elist.length; k++) {
                int[] p1 = {olist[j],elist[k]};
                double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[k]+0.5));
                Real energyR = new Real( (new Double(energy)).toString() );
                Real exponent = new Real( (new Double(-energy/tt)).toString() );
                exponent.exp();
                exponent.mul(dos.probeValue(p1));
                pp.add(exponent);
                //exponent.mul(orderR);
                //ep.add(exponent);
            }
			coord.write(order+" "+pp+"\n");
        }
        //ep.div(pp);
            
    } //end run_coord
    /****
    * computes Q2(Q1,T) from a density of states assuming a dos(E,Q1,Q2)
    */
    private static void run_coord_2D(SparseRealArray dos) {
        Real pp,ep;
        //problem: we have a dos(E,Q1,Q2) but we need dos(Q1,Q2,E) for this integral.  So switch it!
        dos = SparseRealArray.switchIndices(dos,0,1); // dos(Q1,E,Q2)
        dos = SparseRealArray.switchIndices(dos,1,2); //dos(Q1,Q2,E)
        double t = startTC.value;
        double tt = t * kb.value;
        FileIO coord = getOutputFileIO(run_coord_out.value+((int)(t*10+SMALL_SHIFT)));
        int[] o1list = dos.getTopList();
        for(int k = 0; k < o1list.length; k++) { //Q1
            double order1 = start.elementAt(1).value + step.elementAt(1).value * (o1list[k]);
            int[] p = {o1list[k]};
            int[] o2list = dos.getSparseList(p);
            pp = new Real(0);
            ep = new Real(0);
            for(int l = 0; l < o2list.length; l++) { //Q2    
                double order2 = start.elementAt(2).value + step.elementAt(2).value * (o2list[l]);
                Real order2R = new Real( (new Double(order2)).toString() );           
                int[] p1 = {o1list[k],o2list[l]};
                int[] elist = dos.getSparseList(p1);                    
                for(int j = 0; j < elist.length; j++) {
                    double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[j]+0.5));
                    //Real energyR = new Real( (new Double(energy)).toString() );
                    Real exponent = new Real( (new Double(-energy/tt)).toString() );
                    exponent.exp();
                    int[] p2 = {o1list[k],o2list[l],elist[j]};
                    exponent.mul(dos.probeValue(p2));
                    pp.add(exponent);
                    exponent.mul(order2R);
                    ep.add(exponent);
                }
            }            
            ep.div(pp);
            coord.write(order1+" "+ep+"\n");
        }               
        coord.close();
    }//end run_coord_2D
    /****
    * computes prob[Q1,Q2*,T](Q1) from a density of states assuming a dos(E,Q1,Q2)
    */
    private static void run_prob_2D(SparseRealArray dos, double q2low, double q2high) {
        Real pp,ep;
        //problem: we have a dos(E,Q1,Q2) but we need dos(Q1,Q2,E) for this integral.  So switch it!
        dos = SparseRealArray.switchIndices(dos,0,1); // dos(Q1,E,Q2)
        dos = SparseRealArray.switchIndices(dos,1,2); //dos(Q1,Q2,E)
        double t = startTC.value;
        double tt = t * kb.value;
        FileIO coord = getOutputFileIO(run_prob_out.value+((int)(t*10+SMALL_SHIFT)));
        int[] o1list = dos.getTopList();
        for(int k = 0; k < o1list.length; k++) { //Q1
            double order1 = start.elementAt(1).value + step.elementAt(1).value * (o1list[k]+0.5);
            int[] p = {o1list[k]};
            int[] o2list = dos.getSparseList(p);
            pp = new Real(0);
            //ep = new Real(0);
            for(int l = 0; l < o2list.length; l++) { //Q2
                double order2 = start.elementAt(2).value + step.elementAt(2).value * (o2list[l]+0.5);
				if(order2 >= q2low && order2 <= q2high) {				    
	                Real order2R = new Real( (new Double(order2)).toString() );           
	                int[] p1 = {o1list[k],o2list[l]};
	                int[] elist = dos.getSparseList(p1);                    
	                for(int j = 0; j < elist.length; j++) {
	                    double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[j]+0.5));
	                    //Real energyR = new Real( (new Double(energy)).toString() );
	                    Real exponent = new Real( (new Double(-energy/tt)).toString() );
	                    exponent.exp();
	                    int[] p2 = {o1list[k],o2list[l],elist[j]};
	                    exponent.mul(dos.probeValue(p2));
	                    pp.add(exponent);
	                    //exponent.mul(order2R);
	                    //ep.add(exponent);
	                }
				}
            }            
            coord.write(order1+" "+pp+"\n");
        }               
        coord.close();
    }//end run_prob_2D
    /****
    * computes prob[Q1,Q2*,T](Q1) from a density of states assuming a dos(E,Q1,Q2)
    */
    private static void run_prob_heiko_2D(SparseRealArray dos) {
        Real pp,ep;
        //problem: we have a dos(E,Q1,Q2) but we need dos(Q1,Q2,E) for this integral.  So switch it!
        dos = SparseRealArray.switchIndices(dos,0,1); // dos(Q1,E,Q2)
        dos = SparseRealArray.switchIndices(dos,1,2); //dos(Q1,Q2,E)
        double t = startTC.value;
        double tt = t * kb.value;
		double q2prev = startProb.value;
		for(double q2 = startProb.value+stepProb.value; q2 < stepProb.value*numBinProb.value; q2 += stepProb.value) {
	        FileIO prob = getOutputFileIO(run_prob_out.value+((int)(t*10+SMALL_SHIFT))+"."+((int)(q2/stepProb.value+SMALL_SHIFT)));
	        int[] o1list = dos.getTopList();
	        for(int k = 0; k < o1list.length; k++) { //Q1
	            double order1 = start.elementAt(1).value + step.elementAt(1).value * (o1list[k]+0.5);
	            int[] p = {o1list[k]};
	            int[] o2list = dos.getSparseList(p);
	            pp = new Real(0);
	            //ep = new Real(0);
	            for(int l = 0; l < o2list.length; l++) { //Q2
	                double order2 = start.elementAt(2).value + step.elementAt(2).value * (o2list[l]+0.5);
					if(order2 >= q2prev && order2 <= q2) {				    
		                Real order2R = new Real( (new Double(order2)).toString() );           
		                int[] p1 = {o1list[k],o2list[l]};
		                int[] elist = dos.getSparseList(p1);                    
		                for(int j = 0; j < elist.length; j++) {
		                    double energy = (start.elementAt(0).value + step.elementAt(0).value * (elist[j]+0.5));
		                    //Real energyR = new Real( (new Double(energy)).toString() );
		                    Real exponent = new Real( (new Double(-energy/tt)).toString() );
		                    exponent.exp();
		                    int[] p2 = {o1list[k],o2list[l],elist[j]};
		                    exponent.mul(dos.probeValue(p2));
		                    pp.add(exponent);
		                    //exponent.mul(order2R);
		                    //ep.add(exponent);
		                }
					}
	            }            
	            prob.write(order1+" "+pp+"\n");
	        }               
	        prob.close();
			q2prev=q2;
		}
    }//end run_prob_2D
    /****
    * computes Q3(Q1,Q2*,T) from a density of states assuming a dos(E,Q1,Q2,Q3).  Q2* is a range of Q2 values.
    * 
    * @param q2low low limit of Q2 range
    * @param q2high high limit of Q2 range
    */
    private static void run_coord_3D(SparseRealArray dos, double q2low, double q2high) {
        Real pp,ep;
        //problem: we have a dos(E,Q1,Q2,Q3) but we need dos(Q1,Q2,Q3,E) for this integral.  So switch it!
        dos = SparseRealArray.switchIndices(dos,0,1);
        dos = SparseRealArray.switchIndices(dos,1,2);
        dos = SparseRealArray.switchIndices(dos,2,3);
        int[] olist = dos.getTopList();
        double t = startTC.value;
        double tt = t * kb.value;
        FileIO coord = getOutputFileIO(run_coord_out.value+((int)(t*10+SMALL_SHIFT)));
        int[] o1list = dos.getTopList();
        for(int k = 0; k < o1list.length; k++) { //Q1
            double order1 = start.elementAt(1).value + step.elementAt(1).value * (o1list[k]+0.5);
            int[] p = {o1list[k]};
            int[] o2list = dos.getSparseList(p);
            pp = new Real(0);
            ep = new Real(0);
            for(int l = 0; l < o2list.length; l++) { //Q2 only a range are valid, and valids are combined    
                double order2 = start.elementAt(2).value + step.elementAt(2).value * (o2list[l]+0.5);
                if(order2 >= q2low && order2 <= q2high) {   
                    int[] p1 = {o1list[k],o2list[l]};
                    int[] o3list = dos.getSparseList(p1);       
                    for(int m = 0; m < o3list.length; m++) { //Q3
                        double order3 = start.elementAt(3).value + step.elementAt(3).value * (o3list[m]+0.5);
                        Real order3R = new Real( (new Double(order3)).toString() );           
                        int[] p2 = {o1list[k],o2list[l],o3list[m]};
                        int[] elist = dos.getSparseList(p2);
                        for(int j = 0; j < elist.length; j++) { //E
                            double energy = start.elementAt(0).value + step.elementAt(0).value * (elist[j]+0.5);
                            Real exponent = new Real( (new Double(-energy/tt)).toString() );
                            int[] p3 = {o1list[k],o2list[l],o3list[m],elist[j]};
                            exponent.mul(dos.probeValue(p3));
                            pp.add(exponent);
                            exponent.mul(order3R);
                            ep.add(exponent);
                        }
                    }
                }
            }
            ep.div(pp);
            coord.write(order1+" "+ep+"\n");
        }               
        coord.close();
    }//end run_coord_3D
    /****
    * computes Q3(Q1,Q2,T) from a density of states assuming a dos(E,Q1,Q2,Q3).
    * 
    */
    private static void run_full_coord_3D(SparseRealArray dos) {
        Real pp,ep;
        //problem: we have a dos(E,Q1,Q2,Q3) but we need dos(Q1,Q2,Q3,E) for this integral.  So switch it!
        dos = SparseRealArray.switchIndices(dos,0,1);
        dos = SparseRealArray.switchIndices(dos,1,2);
        dos = SparseRealArray.switchIndices(dos,2,3);
        int[] olist = dos.getTopList();
        double t = startTC.value;
        double tt = t * kb.value;
        FileIO coord = getOutputFileIO(run_coord_out.value+((int)(t*10+SMALL_SHIFT)));
        int[] o1list = dos.getTopList();
        for(int k = 0; k < o1list.length; k++) { //Q1
            double order1 = start.elementAt(1).value + step.elementAt(1).value * (o1list[k]+0.5);
            int[] p = {o1list[k]};
            int[] o2list = dos.getSparseList(p);
            for(int l = 0; l < o2list.length; l++) { //Q2 only a range are valid, and valids are combined    
                double order2 = start.elementAt(2).value + step.elementAt(2).value * (o2list[l]+0.5);
//                if(order2 >= q2low && order2 <= q2high) {   
                int[] p1 = {o1list[k],o2list[l]};
                int[] o3list = dos.getSparseList(p1);       
	            pp = new Real(0);
	            ep = new Real(0);
                for(int m = 0; m < o3list.length; m++) { //Q3
                    double order3 = start.elementAt(3).value + step.elementAt(3).value * (o3list[m]+0.5);
                    Real order3R = new Real( (new Double(order3)).toString() );           
                    int[] p2 = {o1list[k],o2list[l],o3list[m]};
                    int[] elist = dos.getSparseList(p2);
                    for(int j = 0; j < elist.length; j++) { //E
                        double energy = start.elementAt(0).value + step.elementAt(0).value * (elist[j]+0.5);
                        Real exponent = new Real( (new Double(-energy/tt)).toString() );
                        int[] p3 = {o1list[k],o2list[l],o3list[m],elist[j]};
                        exponent.mul(dos.probeValue(p3));
                        pp.add(exponent);
                        exponent.mul(order3R);
                        ep.add(exponent);
                    }
                }
//                }
            	ep.div(pp);
				coord.write(order1+" "+order2+" "+ep+"\n");
			}
        }               
        coord.close();
    }//end run_full_coord_3D
    /***********************************
    * Data structures and input/output
    ***********************************/	
	/**
	* Bins the files.  Reads dim dimensions unless doing energy perturbation theory reads numDim+1.
	* The new weight for each snapshot should be put as the last column - this will not be read into the bin.
	*/
	private static Vector<SmartBin> readSkinnyHistograms(int dim) {
	    System.out.println("Reading in histograms.");
	    Vector<SmartBin> sbin = new Vector<SmartBin>();
	    Vector<Integer> lineCounts = new Vector<Integer>();
	    int fileCount = 0;
	    double weight = 1;
	    for(StringHolder sh : fileNames) {
	        sbin.add(new SmartBin(dim,numBin,start,step));
	        FileIO file = new FileIO(sh.value,FileIO.BUFFERED_READING);
	        String line = file.readLine();
	        double[] vals = new double[dim];
	        int count = 0;
	        while (line != null) {
	            String[] tokens = FileIO.getTokens(line," \t");
	            //first token is assumed to be energy
	            int i = 0;
	            try {
    	            for(i = 0; i < dim; i++) {
    	                vals[i] = Double.parseDouble(tokens[i]);
    	            }
    	            if(reweighting.value) {
    	                weight = Double.parseDouble(tokens[numDim.value-1]); //last column is the weight
    	            }
	            } catch (ArrayIndexOutOfBoundsException e) {
	                System.out.println("Died on column "+i+" of file "+sh.value+".  Probably missing data!");
	                System.exit(1);
	            }
	            sbin.elementAt(fileCount).add(vals,weight);
				
				//if(sbin.elementAt(fileCount).add(vals,weight) > 0) {System.out.println(sh.value);for(int j=0;j<vals.length;j++){System.out.println(vals[j]);}}
				//for(int j=0;j<vals.length;j++){System.out.println(vals[j]);}
				//System.exit(0);
	            line = file.readLine();
	            count++;
	        }
	        lineCounts.add(new Integer(count));
	        System.out.println(sh.value+" "+count+" lines.");
	        fileCount++;
	    }
	    return sbin;
	}

	/**
	* Represents an any dimensional histogram as a sparse array. 
	* If reweighting, takes the last column as the weight, otherwise assumes 1.
	*/
	private static class SmartBin {
	    SparseDoubleArray bin;
	    Vector<LongHolder> numBin;
	    Vector<DoubleHolder> start;
	    Vector<DoubleHolder> step;
	    int dim;
	    int[] indices;
	    public SmartBin(int dimension, Vector<LongHolder> numBin, Vector<DoubleHolder> start, Vector<DoubleHolder> step) {
	        bin = new SparseDoubleArray(dimension);  
	        this.dim = dimension;
	        this.numBin = numBin;
	        this.start = start;
	        this.step = step;  
	        indices = new int[dim];
	    }
	    /**
	    * Adds a point to the histogram of arbitrary weight (can be different than 1 if using energy perturbation theory).
	    * @param double[]
	    * @return int 0 if successful, 1 if not within stated bounds
	    */
	    public int add(double[] values, double weight) {
	        if(values.length == dim) {
	            for(int i = 0; i < dim; i++) {
	                indices[i] = (int) ((values[i] - start.elementAt(i).value)/step.elementAt(i).value);
					//System.out.println(values[i]+" "+start.elementAt(i).value+" "+step.elementAt(i).value+" "+indices[i]);
	            }
        	    
	            //check that this is a valid point (within bounds)
	            boolean valid = true;
	            for(int i = 0; i < dim; i++) {
	                if(indices[i] < 0 || indices[i] >= ((int)numBin.elementAt(i).value)) valid = false;
					//if(!valid){System.out.println(values[i]+" "+start.elementAt(i).value+" "+step.elementAt(i).value+" "+indices[i]);}
	            }
	            if(valid) { bin.increment(indices,weight); return 0; }
	            else return 1;
	        }
	        return 1;
	    }

	    public int[][] getIntegerIndex() {
	        return bin.getIntegerIndex();
	    }
	    public double[] getIndexValues() {
	        return bin.getIndexValues();
	    }
	    public SparseDoubleArray getBin() {
	        return bin;
	    }
	}
	

	private static void setUpConfiguration(String[] args) {
	    //set up configuration
	    numDim = new IntHolder();
	    numThreads = new IntHolder();
	    setToZero = new IntHolder();
	    numFile = new IntHolder();
	    maxIterations = new IntHolder();
	    numUmb = new IntHolder();
	    kb = new DoubleHolder();
	    tolerance = new DoubleHolder();
	    numBin = new Vector<LongHolder>(3);
	    start = new Vector<DoubleHolder>(3);
	    step = new Vector<DoubleHolder>(3);
	    umbrella_k = new Vector<DoubleHolder>(10);
	    umbrella_sigma = new Vector<DoubleHolder>(10);
	    umbrella_0 = new Vector<DoubleHolder>(10);
	    umbrella_0_2 = new Vector<DoubleHolder>(10);
	    fileNames = new Vector<StringHolder>(10);
	    fileTemps = new Vector<DoubleHolder>(10);
	    config = new StringHolder();
	    umbrellaType = new Vector<StringHolder>(10);
	    dosFileName = new StringHolder();
	    run_cv_out = new StringHolder();
	    run_free_out = new StringHolder();
	    convergenceType = new StringHolder();
	    convergenceCriteria = new StringHolder();
	    run_coord_out = new StringHolder();
	    run_prob_out = new StringHolder();
	    free_energy_out = new StringHolder();
	    run_free = new BooleanHolder();
	    run_cv = new BooleanHolder();
	    run_wham = new BooleanHolder();
	    run_FEP = new BooleanHolder();
	    run_coord = new BooleanHolder();
	    run_prob = new BooleanHolder();
	    readFreeEnergy = new BooleanHolder();
	    reweighting = new BooleanHolder();
		overwriting = new BooleanHolder();
	    startT = new DoubleHolder();
	    startTF = new DoubleHolder();
	    startTC = new DoubleHolder();
	    numBinProb = new IntHolder();
	    stepProb = new DoubleHolder();
	    startProb = new DoubleHolder();
	    lowQ2 = new DoubleHolder();
	    highQ2 = new DoubleHolder();
	    deltaT = new DoubleHolder();
	    deltaTF = new DoubleHolder();
	    deltaTC = new DoubleHolder();
	    ntemps = new IntHolder();
	    ntempsF = new IntHolder();
	    ntempsC = new IntHolder();
	    
	    ArgParser parser;
	    String usageString = "java -jar [[-Xmx1000m]] WHAM.jar --config [filename] [[other options]]\n\nVersion: "+versionNum+"\n";
	    if(distribution) parser = new ArgParser(usageString+distString);
	    else parser = new ArgParser(usageString);
	    
	    parser.addOption("--config %s #Configuration file",config);
	    parser.addOption("threads %d {[1,48]} #<int>#choose number of threads (1-48)",numThreads);
	    parser.addOption("numDimensions %d {[1,6]} #<int>#Number of histogram dimensions (columns)",numDim);
	    parser.addOption("numUmbrella %d {[0,4]} #<int>#Number of umbrellas included",numUmb);
	    parser.addOption("umbrellaType %s {harmonic,linear} #Type of umbrellas, name columns from left to right",umbrellaType);
	    parser.addOption("tolerance %f #Sets the required convergence level",tolerance);
	    parser.addOption("maxIterations %d #Sets the maximum iterations to attempt for convergence",maxIterations);
	    parser.addOption("convergenceType %s {simple,quick} #If convergence is jumpy switch to simple",convergenceType);	    
	    parser.addOption("convergenceCriteria %s {average,max} #either average (default) of deltas or maximum delta must go below tolerance",convergenceCriteria);	    
	    parser.addOption("dosFile %s #Sets the filename for the density of states",dosFileName);
	    parser.addOption("setToZero %d #<int>#Fix this file's free energy to zero",setToZero);
	    parser.addOption("kB %f #Boltzmann constant to use",kb);
	    parser.addOption("numBins %d #Number of bins in this histogram",numBin);
	    parser.addOption("numBinProb %d #Number of bins in this histogram",numBinProb);
	    parser.addOption("startProb %f #Starting value for this histogram",startProb);
	    parser.addOption("stepProb %f #Step between bins for this histogram",stepProb);
	    parser.addOption("start %f #Starting value for this histogram",start);
	    parser.addOption("step %f #Step between bins for this histogram",step);
	    parser.addOption("numFiles %d #Number of files to add to histogram",numFile);
	    parser.addOption("name %s #Names of the histogram files",fileNames);
	    parser.addOption("temp %f #Temperature of the histrogram files",fileTemps);
	    parser.addOption("umbrella_k %f #stiffness of harmonic/quartic umbrella",umbrella_k);
	    parser.addOption("umbrella_sigma %f #width of the gaussian umbrella",umbrella_sigma);
	    parser.addOption("umbrella_0 %f #center of restraint in this run",umbrella_0);
	    parser.addOption("umbrella_0_2 %f #center of Q_high harmonic restraint in this run",umbrella_0_2);
	    parser.addOption("reweighting %v #use the last column to reweight (for EPT)",reweighting);
	    parser.addOption("startT %f #temperature to start cv at",startT);
	    parser.addOption("startTF %f #temperature to start free plots at",startTF);
	    parser.addOption("startTC %f #temperature to start coord plots at",startTC);
	    parser.addOption("lowQ2 %f #lower limit for coordinate 2 avg",lowQ2);
	    parser.addOption("highQ2 %f #upper limit for coordinate 2 avg",highQ2);
	    parser.addOption("deltaT %f #how fine to make the cv",deltaT);
	    parser.addOption("deltaTF %f #how fine to make the free plots",deltaTF);
	    parser.addOption("deltaTC %f #how fine to make the coordinate plots",deltaTC);
	    parser.addOption("ntemps %d #number of cv points",ntemps);
	    parser.addOption("ntempsF %d #number of free energy plots",ntempsF);
	    parser.addOption("ntempsC %d #number of coordinate points",ntempsC);
		parser.addOption("overwriting %v #allows overwriting of existing files",overwriting);
	    parser.addOption("run_cv_out %s #name of cv file",run_cv_out);
	    parser.addOption("run_free_out %s #appended name of free file",run_free_out);
	    parser.addOption("run_coord_out %s #appended name of free file",run_coord_out);
	    parser.addOption("run_prob_out %s #appended name of free file",run_coord_out);
	    parser.addOption("run_cv %v #tells program to calculate cv",run_cv);
	    parser.addOption("run_free %v #tells program to calculate free vs orders",run_free);
	    parser.addOption("run_FEP %v #tells program to calculate FEP",run_FEP);
	    parser.addOption("run_prob %v #tells program to calculate reweighted histograms",run_prob);
	    parser.addOption("run_wham %v #tells program to calculate new dos",run_wham);	
	    parser.addOption("run_coord %v #tells program to calculate order vs order",run_coord);
	    parser.addOption("readFreeEnergy %v #tells program to read in past energies if file exists",readFreeEnergy);
	    parser.addOption("free_energy_out %s #name of output/input free energy file",free_energy_out);    
	    
	    /**** Set up default values for non-required values ****/
        maxIterations.value = 500;
        kb.value = 0.008314;
        tolerance.value = 0.001;
        dosFileName.value = "dos";
        run_cv_out.value = "cv";
        run_free_out.value = "free";
        run_coord_out.value = "coord";
        run_prob_out.value = "prob";
        free_energy_out.value = "last_free";        
        numUmb.value = 0;
        setToZero.value = -1;
        numThreads.value = 1;
        convergenceType.value = "quick";
        convergenceCriteria.value = "average";
		lowQ2.value = -314;
		startProb.value = -1;
		stepProb.value = -1;
		numBinProb.value = -1;
	    
        parser.matchAllArgs(args);
        if(config.value == null) {parser.printErrorAndExit("Configuration file is required!\nUsage: --config [file]\n");}
        
	    //read in configuration
        String[] configArgs = null;
    	try {
    	    System.out.println("Parsing configuration file: "+config.value);
			File configFile = new File(config.value);
			if(!configFile.exists()) { parser.printErrorAndExit("Some problem with reading configuration file; maybe it doesn't exist?"); }
			configArgs = parser.prependArgs(new File(config.value),configArgs);
	    } catch (IOException e) {
	    	parser.printErrorAndExit("Some problem with reading configuration file; maybe it doesn't exist?");
	    }
	    
        parser.matchAllArgs(configArgs);
		
		//make sure some options are ok
        if(numDim.value == 5 && run_coord.value) {
            if(lowQ2 == null || highQ2 == null) { parser.printErrorAndExit("4/5 dimensions and run_coord requested.  You are required to set lowQ2/highQ2.\n");}
		}
		if(numUmb.value > 0) {
			System.out.println(numUmb.value);
		}	
		//System.exit(0);        
        
        //System.out.println("Energy bins: "+((LongHolder)numBin.elementAt(0)).value);
        //print out args
        if(fileNames.size()!=fileTemps.size()) {
            System.out.println("Must be the same number of files and temperatures.  Exiting");
            System.exit(1);
        }
        System.out.println("\nListing files and temperatures");
        for(int i = 0; i < fileNames.size();i++) {
            System.out.println(fileNames.elementAt(i).value+" "+fileTemps.elementAt(i).value);
        }
        
        scaledTemps = new double[fileTemps.size()];
	    System.out.println("\nScaling temperatures * kB");
	    for(int i = 0; i < fileTemps.size(); i++) {
	        scaledTemps[i] = fileTemps.elementAt(i).value * kb.value;
	        //System.out.println(fileTemps.elementAt(i).value+" ---> "+scaledTemps[i]);
	    }
	}
    private static void writeDOS(SparseRealArray dos, int dim) {
	    System.out.println("Writing density of states to file: "+dosFileName.value);
		FileIO dosFile = getOutputFileIO(dosFileName.value);
        int nenergy = (int)numBin.elementAt(0).value;
	    int norder = (int)numBin.elementAt(1).value;
        int[][] index = dos.getIntegerIndex();
        Real[] values = dos.getIndexValues();
        Real.NumberFormat format = new Real.NumberFormat();
        format.fse = Real.NumberFormat.FSE_SCI;
        for(int i = 0; i < index.length; i++) {
            for(int j = 0; j < dim; j++) {
                float value = (float)(start.elementAt(j).value + step.elementAt(j).value * index[i][j]);
                dosFile.write(value+" ");
            }
            dosFile.write(values[i].toString(format)+"\n");
        }
	}
    private static SparseRealArray readDOS(int dim) {
        //assume its the same as the current configuration
        SparseRealArray dos = new SparseRealArray(dim,true);
        int nenergy = (int)numBin.elementAt(0).value;
	    int norder = (int)numBin.elementAt(1).value;
	    System.out.println("Reading density of states from file: "+dosFileName.value);
        FileIO dosFile = new FileIO(dosFileName.value,FileIO.BUFFERED_READING);
        String line = dosFile.readLine();
        while (line != null) {
            String[] tokens = FileIO.getTokens(line," ");
            int[] point = new int[dim];
            for(int i = 0; i < dim; i++) {
                point[i] = (int)((Double.parseDouble(tokens[i]) - start.elementAt(i).value + SMALL_SHIFT)/step.elementAt(i).value);
            }
            dos.increment(point,new Real(tokens[dim]));
            line = dosFile.readLine();
        }
        return dos;
    }
    
    private static void readFree(double[] z) {
        //look for file of name the value of free_energy_out 
        try {
            File f = new File(free_energy_out.value);
            System.out.println(free_energy_out.value);
            if(f.exists()) {
                //read in file (just a list of doubles)
                FileIO file = new FileIO(free_energy_out.value,FileIO.BUFFERED_READING);
    	        String line = file.readLine();
    	        int count = 0;
    	        while (line != null) {
    	            String[] tokens = FileIO.getTokens(line," ");
    	            z[count] = Double.parseDouble(tokens[0]);
    	            System.out.println(z[count]);
    	            count++;
    	            line = file.readLine();
                }
            }
        } catch (Exception e) {
            System.out.println("Error in reading initial free energies from file: "+free_energy_out.value+" \n Maybe something changed, like file count. Using zeros to start.");
            z = Utilities.initializeWithZeros(new double[numFile.value]);
        }
    }
    private static void writeFree(double[] z) {
        FileIO free = getOutputFileIO(free_energy_out.value);
        for (int i = 0; i < z.length; i++) {
            free.write(z[i]+"\n");
        }
        free.close();
    }
	// sums the right entries in double[][][] biases in order to get the bias for a 
	//particular histogram bin identified by file index and umbrella indicies
	private static double getBias(double[][][] biases, int file, int[] umbrellaIndices) {
		double bias = 0;
		for(int i = 0; i < umbrellaIndices.length; i++) {
			bias += biases[file][i][umbrellaIndices[i]];
		}
		return bias;
	}
    
	//Unique linear index into multidimensional array of dimension dim
	//with given indices a_i and given row lengths N_i
	//returns \sum_{i=0}^{dim} a_i * \product_{j=0}^{i-1} N_j
	private static long get1Dindex(int[] indices, int[] lengths, int dim) {
		long arrayIndex = 0;
		for(int i = 0; i < dim; i++) {
			long product = 1;
			for( int j = 0; j < i; j++) {
				product *= lengths[j];
			}
			arrayIndex += product * indices[i];
		}
		return arrayIndex;
	}
	
    /*********************************************
    * Thread for the run_wham algorithm
    **********************************************/

    private static class InnerThread implements Runnable {
        int START, END, index, nfile;
        double[] energy;
        double[][][] bias;
        Real[][] z0part;
        //Real[][][] zzRecip;
		HashMap zzRecipH;
        int[][] summedIndex;
        Real[] summedValuesReal;
        double[] zold;
		int dim;
		double[][][] biases;

  	    static boolean[] threads = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
								true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
								true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true };
								
        public InnerThread(int[][] summedIndex, Real[] summedValuesReal, double[] energy, Real[][] z0part, double[][][] biases, double zold[], int nfile, int numUmbrella, int START, int END, int index) {
            this.energy = energy;
            this.bias = bias;
            this.summedIndex = summedIndex;
            this.summedValuesReal = summedValuesReal;
			//(a shared zzRecip creates concurrency issues so each thread has its own)
			//this.zzRecip = Utilities.initializeWithZeros(new Real[numBins[0]][numBins[actualDim-2]][numBins[actualDim-1]]);   	
			this.zzRecipH = new HashMap<Long,Real>();
            this.z0part = z0part;
            this.nfile = numFile.value;
            threads[index] = false;
            this.index = index;
            this.START = START;
            this.END = END;
            this.zold = zold;
			this.dim = numUmbrella+1; //+1 because of energy
			this.biases = biases;

        }
        public void run() {
            //System.out.println("Hello, I am thread "+index+".");
			int[] indices = new int[dim];
			int[] umbrellaIndices = new int[dim]; //intialization for compilation
			if(dim>1) umbrellaIndices = new int[dim-1];
			int[] lengths = new int[dim];
			lengths[0] = numBins[0];
			for(int i = 1; i < dim; i++ ) { lengths[i] = numBins[actualDim-(dim-i)]; }
			double[] biasEnergy = Utilities.initializeWithZeros(new double[nfile]);
			
			for(int k = START; k <= END; k++) {
                int energyIndex = summedIndex[k][0];
                int umbrellaIndex1 = summedIndex[k][actualDim-2];
                int umbrellaIndex2 = summedIndex[k][actualDim-1];
				indices[0] = summedIndex[k][0];
				for(int i = 1; i < dim; i++ ) { 
					indices[i] = summedIndex[k][actualDim-(dim-i)]; 
					umbrellaIndices[i-1] = summedIndex[k][actualDim-(dim-i)]; 
				}
				long hashindex = get1Dindex(indices,lengths,dim);
				if(dim > 1) { //there are umbrellas
					for(int i = 0; i < nfile; i++) {
						biasEnergy[i] = getBias(biases,i,umbrellaIndices);
					}
				}
				//System.out.println(energyIndex+" "+umbrellaIndex1+" "+umbrellaIndex2+" "+hashindex+"\n");
				//Utilities.print(lengths,false);
				//System.exit(0);
				if(!zzRecipH.containsKey(hashindex)) { //don't do more than once
					zzRecipH.put(hashindex,new Real(0));
                    for(int i = 0; i < nfile; i++) {
						
                        //System.out.println(energyIndex+" "+umbrellaIndex1+" "+umbrellaIndex2+" "+biasEnergy+" "+bias[i][umbrellaIndex1][umbrellaIndex2]);
						Real exponent = new Real( (new Double( ( -energy[energyIndex] - biasEnergy[i] + zold[i] )/scaledTemps[i] )).toString() );
                        exponent.exp();
                        exponent.mul(nconf[i]);
						((Real)zzRecipH.get(hashindex)).add(exponent);
						//System.out.println("putting "+exponent+" into index "+hashindex)
                    }
					((Real)zzRecipH.get(hashindex)).recip();
                }
                for(int i = 0; i < nfile; i++) {
                    Real exponent = new Real( (new Double( ( -energy[energyIndex] - biasEnergy[i] ) / scaledTemps[i])).toString());
                    exponent.exp();
                    exponent.mul(summedValuesReal[k]);
					exponent.mul((Real)zzRecipH.get(hashindex));
                    z0part[i][index].add(exponent);
                }
            }
            threads[index] = true; //I am finished! 
        }
    }


}
