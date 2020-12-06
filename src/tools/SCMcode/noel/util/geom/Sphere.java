package noel.util.geom;
import noel.util.math.*;

public class Sphere {
	double radius;
	Coordinates coord; //center

	//Create a sphere with radius r and center c
	public Sphere(double r, Coordinates c) {
		radius = r;
		coord = c;
	}
	public Coordinates getCenter() {
		return coord;
	}
	public double getRadius() {
		return radius;
	}
	/**
	 * Gives the intersection volume between two spheres with radii r1 and r2,
	 * separated by a distance d
	 * @param d separation
	 * @param r1 radius sphere 1
	 * @param r2 radius sphere 2
	 */
	public static double intersectionVolume(double r1, double r2, double d) {
	    return Math.PI * (r1 + r2 - d)*(r1 + r2 -d)*(d*d + 2*d*(r1+r2) - 3*(r1*r1 + r2*r2) + 6*r1*r2) / d / 12;
	}
}
