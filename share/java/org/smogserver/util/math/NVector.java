/*
 * NVector.java
 *
 * Created on November 9, 2003, 10:38 PM
 */

package org.smogserver.util.math;

import java.util.*;
/**
 *
 * @author  jknoel
 * @version
 */
public class NVector implements Cloneable {
    public double[] values;
    public int length;
    /** Creates new NVector */
    public NVector(int length) {
        this.length = length;
        this.values = new double[length];
    }
    public NVector(int length, boolean initialize) {
    	this.length = length;
        this.values = new double[length];
        for (int i = 0; i < length; i++) {
        	this.values[i] = 0;
        }
    }
    public NVector(double[] values) {
        this.length = values.length;
        this.values = values;
    }
	public NVector(Coordinates a, Coordinates b){
		this.length = 3;
		this.values = new double[length];
		values[0] = a.x - b.x;
		values[1] = a.y - b.y;
		values[2] = a.z - b.z;
	}
    public double[] returnArray() {
        return values;
    }
    public void setArray(double[] values) {
    	this.values = values;
    }
    public int length() {
        return length;
    }
    public double valueAt(int i) {
        return values[i];
    }
    public void setValueAt(int i, double v) {
        values[i] = v;
    }
    public void setNVector(double[] v) {
        values = v;
    }
    public NVector add(NVector p) {
        NVector n = new NVector(length);
        for (int i = 0; i < length; i++) {
            n.setValueAt(i, values[i]+p.valueAt(i));
        }
        return n;
    }
    public void addVoid(NVector p) {
        //NVector n = new NVector(length);
        for (int i = 0; i < length; i++) {
            setValueAt(i, values[i]+p.valueAt(i));
        }
    }
    public NVector sub(NVector p) {
        NVector n = new NVector(length);
        for (int i = 0; i < length; i++) {
            n.setValueAt(i, values[i]-p.valueAt(i));
        }
        return n;
    }
    public NVector mul(NVector p) {
        NVector n = new NVector(length);
        for (int i = 0; i < length; i++) {
            n.setValueAt(i, values[i]*p.valueAt(i));
        }
        return n;
    }
    public NVector mul(double p) {
        NVector n = new NVector(length);
        for (int i = 0; i < length; i++) {
            n.setValueAt(i, values[i]*p);
        }
        return n;
    }
    public void mulVoid(double p) {
        for (int i = 0; i < length; i++) {
            setValueAt(i, values[i]*p);
        }
    }
    public NVector add(double p) {
        NVector n = new NVector(length);
        for (int i = 0; i < length; i++) {
            n.setValueAt(i, values[i]+p);
        }
        return n;
    }
    public NVector sub(double p) {
        NVector n = new NVector(length);
        for (int i = 0; i < length; i++) {
            n.setValueAt(i, values[i]-p);
        }
        return n;
    }
    public Object clone() {
        NVector w = new NVector(this.length);
        for (int i = 0; i < length; i++) {
            w.setValueAt(i,values[i]);
        }
        return w;
    }
    public String toString() {
        StringBuffer s = new StringBuffer();
        for (int i = 0; i < length; i++) {
            s.append(values[i]);
            s.append(" ");
        }
        return s.toString();
    }
    public double distance(NVector a) {
        double sum = 0;
        double[] avalues = a.returnArray();
        for (int i = 0; i < length; i++) {
            double dist = avalues[i] - values[i];
            sum += dist*dist;
            //sum += Math.pow(avalues[i] - values[i],2); //slow!
        }
        return Math.sqrt(sum);
    }
	public NVector cross(NVector a){
        NVector n = new NVector(length);
		if(length != 3) return n;
		else {
			n.setValueAt(0,values[1]*a.valueAt(2)-values[2]*a.valueAt(1));
			n.setValueAt(1,values[2]*a.valueAt(0)-values[0]*a.valueAt(2));
			n.setValueAt(2,values[0]*a.valueAt(1)-values[1]*a.valueAt(0));
		}
		return n;
	}
	public double dot(NVector a){		
		double sum = 0;
		for(int i = 0; i < length; i++) {
			sum+=values[i]*a.valueAt(i);
		}
		return sum;
	}
	public double mag() {
		double mag = 0;
		for(int i = 0; i < length; i++) {
			mag+=values[i]*values[i];
		}
		return Math.sqrt(mag);
	}
	public void normalizeVoid() {
	    mulVoid(1/mag());
	}
	public NVector normalize() {
	   return mul(1/mag());
	}
	//returns an nvector with components uniformly between 0 and max
	public static NVector getRandomNVector(int dimension, long seed, double max) {
		Random rand = new Random(seed);
		double[] values = new double[dimension];
		for (int i = 0; i < dimension; i++) {
			values[i] = rand.nextDouble() * max;
		}
		return new NVector(values);
	}
}
