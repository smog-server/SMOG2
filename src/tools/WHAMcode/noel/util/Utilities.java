/*
 * Utilities.java
 *
 * Created on October 22, 2003, 10:46 PM
 */

package noel.util;

import java.awt.Color;
import java.util.*;
import ral.*;

/**
 *
 * @author  jknoel
 * @version
 */
public class Utilities {

    public static Random random = new Random();
    /** Creates new Utilities */
    public Utilities() {
    }

    public static Color getColor(int value) {
        return new Color(value);
    }
    public static Color getRandomColor() {
        Random rand = new Random((long)(random.nextDouble()*Long.MAX_VALUE));
        Color c;
        c = new Color((int)(255*rand.nextDouble()),(int)(255*rand.nextDouble()),
            (int)(255*rand.nextDouble()));
        return c;
    }
    public static int getRandomIntOf(int max) {
        return (int)(random.nextDouble() * max);
    }
    /**
     * Returns a string with up to 3 leading zeros - max number is 9999
     * @param int number
     * @return String number
     */
    public static String getLeadingZeros(int n) {
        if (n < 10) return "000"+n;
        if (n < 100) return "00"+n;
        if (n < 1000) return "0"+n;
        return (new Integer(n)).toString();
    }
    public static int findLastEntry(String[] o) {
        int i = 0;
        while (o[i] != null && i < o.length) i++;
        return i-1;
    }
    public static String[] enlarge(String[] f)
    {
        int numLinesInArray = f.length;
        String[] tempLines = new String[numLinesInArray * 2];
        for (int i=0; i < numLinesInArray; i++)
        {
            tempLines[i] = f[i];
        }
        f = tempLines;
        return f;
    }
    //give algorithm standard deviation and average
    public static double getRandomGaussian(double average, double sd) {
        Random rand = new Random((long)(random.nextDouble()*Long.MAX_VALUE));
        return random.nextGaussian() * sd + average;
    }
    /*use clone stupid - not this!*/
    public static double[][] copy2DDoubleArray(double[][] a) {
    	double[][] b = new double[a.length][a[0].length];
    	for (int i = 0; i < a.length; i++) {
    		for (int j = 0; j < a[i].length; j++) {
    			b[i][j] = a[i][j];
    		}
    	}
    	return b;
    }
    /*public static String[] trimArray(String[] f) {
        int numLinesInArray = f.length;
        for (int i = 0; i < numLinesInArray; i++)
        {
            tempLines[i] = f[i];
        }
        String[] tempLines = new String[numLinesInArray];
        f = tempLines;
        return f;
     } */

    /* Tons of print methods! */
    public static void print(double[] d, boolean columns) {
        if (columns) {
            for (int i = 0; i < d.length; i++) {
                System.out.println(d[i]);
            }
        } else {
            for (int i = 0; i < d.length; i++) {
                System.out.print(d[i] + " ");
            }
            System.out.println();
        }
    }
    public static String getPrint(double[] d, boolean columns) {
        StringBuffer buf = new StringBuffer();
        if (columns) {
            for (int i = 0; i < d.length; i++) {
                buf.append(d[i]+"\n");
            }
        } else {
            for (int i = 0; i < d.length; i++) {
                buf.append(d[i] + " ");
            }
            buf.append("\n");
        }
        return buf.toString();
    }
    public static void print(int[] d, boolean columns) {
        if (columns) {
            for (int i = 0; i < d.length; i++) {
                System.out.println(d[i]);
            }
        } else {
            for (int i = 0; i < d.length; i++) {
                System.out.print(d[i] + " ");
            }
            System.out.println();
        }
    }
    public static String getPrint(int[] d, boolean columns) {
        StringBuffer buf = new StringBuffer();
        if (columns) {
            for (int i = 0; i < d.length; i++) {
                buf.append(d[i]+"\n");
            }
        } else {
            for (int i = 0; i < d.length; i++) {
                buf.append(d[i] + " ");
            }
            buf.append("\n");
        }
        return buf.toString();
    }
    public static void print(double[][] d) {
        for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                System.out.print(d[i][j] + " ");
            }
            System.out.println();
        }
    }
    public static void print(int[][] d) {
        for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                System.out.print(d[i][j] + " ");
            }
            System.out.println();
        }
    }
    public static void print(String[] s) {
    	for (int i = 0; i < s.length; i++) {
    		System.out.print(s[i]);
    	}
    }
    public static void println(String[] s) {
    	for (int i = 0; i < s.length; i++) {
    		System.out.println(s[i]);
    	}
    }
    public static double[][] copy(double[][] d) {
        double[][] n = new double[d.length][d[0].length];
        for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                n[i][j] = d[i][j];
            }
        }
        return n;
    }
    //matrix initialization routines!
    public static double[][] initializeWithZeros(double[][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                d[i][j] = 0;
            }
        }
        return d;
    }
    public static Real[] initializeWithZeros(Real[] d) {
    	for (int i = 0; i < d.length; i++) {
                d[i] = new Real();                    
        }
        return d;
    }
    public static Real[][] initializeWithZeros(Real[][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                d[i][j] = new Real();                    
            }
        }
        return d;
    }
    public static Real[][][] initializeWithZeros(Real[][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    d[i][j][k] = new Real();                    
                }
            }
        }
        return d;
    }
    public static Real[][][][] initializeWithZeros(Real[][][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    for(int l = 0; l < d[i][j][k].length; l++) {
                        d[i][j][k][l] = new Real();                    
                    }
                }
            }
        }
        return d;
    }
    public static Real[][][][][] initializeWithZeros(Real[][][][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    for(int l = 0; l < d[i][j][k].length; l++) {
                        for(int m = 0; m < d[i][j][k][l].length; m++) {
                            d[i][j][k][l][m] = new Real();                    
                        }
                    }
                }
            }
        }
        return d;
    }
    public static double[][][] initializeWithZeros(double[][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    d[i][j][k] = 0;                    
                }
            }
        }
        return d;
    }
    public static double[][][][] initializeWithZeros(double[][][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    for(int l = 0; l < d[i][j][k].length; l++) {
                        d[i][j][k][l] = 0;  
                    }                  
                }
            }
        }
        return d;
    }
    public static double[][][][][] initializeWithZeros(double[][][][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    for(int l = 0; l < d[i][j][k].length; l++) {
                        for(int m = 0; m < d[i][j][k][l].length; m++) {
                            d[i][j][k][l][m] = 0;  
                        }
                    }                  
                }
            }
        }
        return d;
    }
    public static int[][][][] initializeWithZeros(int[][][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    for(int l = 0; l < d[i][j][k].length; l++) {
                        d[i][j][k][l] = 0;  
                    }                  
                }
            }
        }
        return d;
    }
    public static int[][][] initializeWithZeros(int[][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    d[i][j][k] = 0;                    
                }
            }
        }
        return d;
    }
    public static boolean[][] initializeWithZeros(boolean[][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                d[i][j] = false;
            }
        }
        return d;
    }
    public static boolean[][][] initializeWithZeros(boolean[][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    d[i][j][k] = false;                    
                }
            }
        }
        return d;
    }
    public static boolean[][][][] initializeWithZeros(boolean[][][][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for(int k = 0; k < d[i][j].length; k++) {
                    for(int l = 0; l < d[i][j][k].length; l++) {
                        d[i][j][k][l] = false;  
                    }                  
                }
            }
        }
        return d;
    }
    public static int[][] initializeWithZeros(int[][] d) {
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                d[i][j] = 0;
            }
        }
        return d;
    }
    public static double sum(int[][] d) {
        double sum = 0;
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                sum += d[i][j];
            }
        }
        return sum;
    }
    public static double sum(int[][][] d) {
        double sum = 0;
    	for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[i].length; j++) {
                for (int k = 0; k < d[i][j].length; k++) {
                    sum += d[i][j][k];
                }
            }
        }
        return sum;
    }
    public static int[] initializeWithZeros(int[] d) {
    	for (int i = 0; i < d.length; i++) {
            d[i] = 0;
        }
        return d;
    }
    public static double[] initializeWithZeros(double[] d) {
    	for (int i = 0; i < d.length; i++) {
            d[i] = 0;
        }
        return d;
    }
    public static Double[] wrapDoubleArray(double[] d) {
    	Double[] db = new Double[d.length];
    	for (int i = 0; i < d.length; i++) {
    		db[i] = new Double(d[i]);
    	}
    	return db;
    }
}
