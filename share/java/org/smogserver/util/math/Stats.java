package org.smogserver.util.math;
import org.smogserver.util.*;

public class Stats {
	public static double[] testNums = {1,2,3,4,5,6,7,8,9,10};
	public static double[] testNums2 = {4,9,11,12,17,5,8,12,14};
	//test this class!


	public static void main(String[] args) {
		//System.out.println("The standard devation is: "+Stats.stdDev(Stats.testNums2));
		/*double[] p = testNums;
		  p[2] = 8;
		  Utilities.print(testNums,false);
		  Utilities.print(p,false);*/
		Stats p = new Stats();
		// String[] s = new String[2];
		// s[1] = "b";
		// System.out.print(s[1]);
		// changeS(s);
		// System.out.print(s[1]);
		int total = 100;
		double testCorr[] = new double[total];
		for(int i = 0; i < total; i++) {
			testCorr[i] = Math.sin(i*(Math.PI*2/20));
		}
		Utilities.print(testCorr,true);
		System.out.println("yep");
		Utilities.print(p.getAutoCorrelationTime(testCorr,20),true);
	}

	public Stats() {	}
	public static void changeS(String[] p) {
		p[1] = "a";
	}
	public static double stdDev(double[] d) {
		double sum = 0, sumsq = 0;
		for (int i = 0; i < d.length; i++) {
			sum += d[i];
			sumsq += Math.pow(d[i],2);
		}
		double top = sumsq*d.length - sum*sum;
		double bottom = d.length*(d.length-1);
		return Math.sqrt(top / bottom);
	}
	public static double var(double[] d) {
		double sum = 0, sumsq = 0;
		for (int i = 0; i < d.length; i++) {
			sum += d[i];
			sumsq += Math.pow(d[i],2);
		}
		double top = sumsq*d.length - sum*sum;
		double bottom = d.length*(d.length-1);
		return top / bottom;
	}
	public static double stdDev(int[] d) {
		int sum = 0, sumsq = 0;
		for (int i = 0; i < d.length; i++) {
			sum += d[i];
			sumsq += (int)Math.pow(d[i],2);
		}
		int top = sumsq*d.length - sum*sum;
		int bottom = d.length*(d.length-1);
		return Math.sqrt(top / bottom);
	}
	public static double getPearsonCorrelation(double[] scores1,double[] scores2){
		double result = 0;
		double sum_sq_x = 0;
		double sum_sq_y = 0;
		double sum_coproduct = 0;
		double mean_x = scores1[0];
		double mean_y = scores2[0];
		for(int i=2;i<scores1.length+1;i+=1){
			double sweep =Double.valueOf(i-1)/i;
			double delta_x = scores1[i-1]-mean_x;
			double delta_y = scores2[i-1]-mean_y;
			sum_sq_x += delta_x * delta_x * sweep;
			sum_sq_y += delta_y * delta_y * sweep;
			sum_coproduct += delta_x * delta_y * sweep;
			mean_x += delta_x / i;
			mean_y += delta_y / i;
		}
		double pop_sd_x = (double) Math.sqrt(sum_sq_x/scores1.length);
		double pop_sd_y = (double) Math.sqrt(sum_sq_y/scores1.length);
		double cov_x_y = sum_coproduct / scores1.length;
		result = cov_x_y / (pop_sd_x*pop_sd_y);
		return result;
	}

	public static double getPearsonCorrelation(Double[] scores1,Double[] scores2){
		double result = 0;
		double sum_sq_x = 0;
		double sum_sq_y = 0;
		double sum_coproduct = 0;
		double mean_x = scores1[0].doubleValue();
		double mean_y = scores2[0].doubleValue();
		for(int i=2;i<scores1.length+1;i+=1){
			double sweep =Double.valueOf(i-1)/i;
			double delta_x = scores1[i-1].doubleValue()-mean_x;
			double delta_y = scores2[i-1].doubleValue()-mean_y;
			sum_sq_x += delta_x * delta_x * sweep;
			sum_sq_y += delta_y * delta_y * sweep;
			sum_coproduct += delta_x * delta_y * sweep;
			mean_x += delta_x / i;
			mean_y += delta_y / i;
		}
		double pop_sd_x = (double) Math.sqrt(sum_sq_x/scores1.length);
		double pop_sd_y = (double) Math.sqrt(sum_sq_y/scores1.length);
		double cov_x_y = sum_coproduct / scores1.length;
		result = cov_x_y / (pop_sd_x*pop_sd_y);
		return result;
	}
	public static double getOverlap(Double[] X,Double[] Y){
		double sumProducts=0;
		double sumxx=0;
		double sumyy=0;
		double x,y;
		for(int i=0;i<X.length;i++){
			x=X[i].doubleValue();
			y=Y[i].doubleValue();
			sumProducts+= x*y;
			sumxx+= x*x;
			sumyy+= y*y;
		}
		return sumProducts/Math.sqrt(sumxx*sumyy);
	}
	//assumes x(t)
	public static double[] getAutoCorrelationTime(double[] x, int maxTime) {
		//int counts = new int[maxTime];
		double[] corr = new double[maxTime];
		double avg, avgSq, sum = 0, sumsqrs = 0;
		for(int i = 0; i < x.length; i++) { //average x,std dev x
			sum+=x[i];
			sumsqrs+=x[i]*x[i];
		}
		avg = sum/x.length;
		avgSq = sumsqrs/x.length;
		double variance = avgSq - avg*avg;
		corr[0] = 1;
		for(int i = 1; i < maxTime; i++) { //tau
			for(int j = 0; j < x.length - i; j++) {
				corr[i] += (x[j]-avg)*(x[j+i]-avg);
			}
			corr[i] = (corr[i] / (x.length - i)) / (variance);
		}
		//System.out.println(variance+" "+avg+" "+avgSq);
		return corr;
	}
}
