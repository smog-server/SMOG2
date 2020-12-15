package org.smogserver.util.math;
import org.smogserver.util.*;

public class LinearAlgebra {
	
	/* tridag solves a tridiagonal system of linear equations
	 *
	 *@param double[] a - diagonal furthest left
	 *@param double[] b - diagonal on center
	 *@param double[] c - diagonal to the right
	 *@param double[] d - right hand sides
	 *@return double[] - the solution
	 */
	public static double[] tridag(double[] a, double[] b, double[] c, double[] d) {
		int num = a.length; //number of equations
		//double[] sol = new double[num];
		double m;
		
		for (int k = 1; k < num; k++) {
			if (b[k-1] == 0.0) {System.out.println("element is zero!"); }
			m = a[k] / b[k-1];
			b[k] = b[k] - m * c[k-1];
			d[k] = d[k] - m * d[k-1];
			System.out.println("d[k] = "+d[k]+" c[k] = "+c[k]);
			System.out.println("b[k] = "+b[k]+" m = "+m);
		}
		d[num-1] = d[num-1] / b[num-1];
		for (int k = num - 2; k >= 0; k--) {
			System.out.println("d[k] = "+d[k]+" d[k+1] = "+d[k+1]);
			d[k] = (d[k] - c[k]*d[k+1]) / b[k];
			System.out.println("d[k] = "+d[k]+" c[k] = "+c[k]+" b[k] = "+b[k]);
		}
		return d;
	}
	public static void main(String args[]) {
		LinearAlgebra.test();
	}
	
	public static void test() {
		int num = 10;
		double[] a = new double[num];
		double[] b = new double[num];
		double[] c = new double[num];
		double[] d = new double[num];
		for (int i = 0; i < num; i++) {
			a[i] = -1;
			b[i] = 2;
			c[i] = -1;
			d[i] = 0;
		}
		d[0] = 1;
		//Utilities.print(d,true);
		d = LinearAlgebra.tridag(a,b,c,d);
		Utilities.print(d,true);
	}
}
