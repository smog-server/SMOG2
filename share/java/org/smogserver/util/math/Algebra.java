package org.smogserver.util.math;

public class Algebra {
	//returns plus solution to quadratic equation with coeffs a,b,c
	public static double plusQuadratic(double a, double b, double c) {
		return ( (-b) + Math.sqrt(b*b-4*a*c) ) / (2 * a);
	}
	//returns minus solution to quadratic equation with coeffs a,b,c
	public static double minusQuadratic(double a, double b, double c) {
		return ( (-b) - Math.sqrt(b*b-4*a*c) ) / (2 * a);
	}
}

