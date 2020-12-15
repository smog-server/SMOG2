/*
 * Complex.java
 *
 * Represents a complex number--> a + ib
 *
 * Created on August 7, 2003, 10:57 PM
 */

/**
 * @author  jknoel
 * @version 
 */
package org.smogserver.util.math;

//actually the Comparable implementation sucks
//All arithmetic operations are static and return Complex objects
public class Complex implements Comparable {

    private double real;
    private double imag;
    
    private static double a, b, c, d, t1, t2;
    
    /** Creates new Complex */
    public Complex(double real, double imag) {
        this.real = real;
        this.imag = imag;
    }

    public double getRealPart() {
        return real;
    }
    
    public double getImaginaryPart() {
        return imag;
    }
    
    public double modulus()
    {
        return Math.sqrt(real*real + imag*imag);
    }
    
    public int compareTo(java.lang.Object obj) {
        return 0;
    }
    
    public String toString() {
        return "(" + String.valueOf(real) + ", " + String.valueOf(imag) + ")";
    }
    
    //STATIC METHODS
    
    public static Complex add(Complex n1, Complex n2) {
        a = n1.getRealPart();
        b = n1.getImaginaryPart();
        c = n2.getRealPart();
        d = n2.getImaginaryPart();
        return new Complex(a+c,b+d);
    }
    
    public static Complex sub(Complex n1, Complex n2) {
        a = n1.getRealPart();
        b = n1.getImaginaryPart();
        c = n2.getRealPart();
        d = n2.getImaginaryPart();
        return new Complex(a-c,b-d);
    }
    
    public static Complex mul(Complex n1, Complex n2) {
        a = n1.getRealPart();
        b = n1.getImaginaryPart();
        c = n2.getRealPart();
        d = n2.getImaginaryPart();
        t1 = a*c - b*d;
        t2 = b*c + d*a;
        return new Complex(t1,t2);
    }
    
    public static Complex div(Complex n1, Complex n2) {
        a = n1.getRealPart();
        b = n1.getImaginaryPart();
        c = n2.getRealPart();
        d = n2.getImaginaryPart();
        double square = Math.pow(c,2) + Math.pow(d,2);
        t1 = (a*c + b*d) / square;
        t2 = (b*c - d*a) / square;
        return new Complex(t1,t2);
    }
    //returns complex conjugate
    public static Complex conj(Complex n1) {
        a = n1.getRealPart();
        b = n1.getImaginaryPart();
        return new Complex(a, -1.0 * b);
    }
    // returns real part
     public static double real(Complex n1) {
        return n1.getRealPart();
    }
    //returns imaginary part
    public static double imag(Complex n1) {
        return n1.getImaginaryPart();
    }
    //return exponential of complex number
    public static Complex exp(Complex n1) {
        a = n1.getRealPart();
        b = n1.getImaginaryPart();
        t1 = Math.exp(a);
        t2 = t1 * Math.sin(b);
        t1 = t1 * Math.cos(b);
        return new Complex(t1,t2);
    }
        
}
