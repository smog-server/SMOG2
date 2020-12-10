package org.smogserver.util.math;


public class Coordinates {
    public double x;
    public double y;
    public double z;
    
    public Coordinates(double x, double y, double z) {
    	this.x = x;
    	this.y = y;
    	this.z = z;
    }
    
    /*
     * Gives the absolute distance between two coordinates
     */
    public static double distance(Coordinates c1, Coordinates c2) {
	    double diffx = c1.x - c2.x;
    	double diffy = c1.y - c2.y;
    	double diffz = c1.z - c2.z;

    	return Math.sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    }
    public static double diffx(Coordinates c1, Coordinates c2) {
	    return c1.x - c2.x;
    }
    public static double diffy(Coordinates c1, Coordinates c2) {
	    return c1.y - c2.y;
    }
    public static double diffz(Coordinates c1, Coordinates c2) {
	    return c1.z - c2.z;
    }
    public static double distanceSq(Coordinates c1, Coordinates c2) {
	    double diffx = c1.x - c2.x;
    	double diffy = c1.y - c2.y;
    	double diffz = c1.z - c2.z;
	
    	return diffx*diffx + diffy*diffy + diffz*diffz;
    }
	
    public String toString() {
	    return x+" "+y+" "+z;
    }
    public String toFormattedString(String format) {
	    return String.format(format,x,y,z);
    }
    public void update(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    public Coordinates shift(NVector v) {
        double x = this.x + v.valueAt(0);
        double y = this.y + v.valueAt(1);
        double z = this.z + v.valueAt(2);
        return new Coordinates(x,y,z);
    }
    public void shiftVoid(NVector v) {
        this.x += v.valueAt(0);
        this.y += v.valueAt(1);
        this.z += v.valueAt(2);
    }
}
