package noel.util.math;


public class SimpleDoubleBin {
	
	double[] bins;
	double min;
	double max;
	double spread;
	double increment;
	long totalPoints;
	
    public SimpleDoubleBin (int num, double min, double max) {
        this.min = min;
        this.max = max;
        this.bins = new double[num];
        spread = (max-min);
        increment = spread/num;
        totalPoints = 0;
    }
    
	public void bin(double num) {
	    int index = (int)((num - min) / increment);
	    bins[index]++;
	    totalPoints++;
	    //System.out.println("here");
	}
	
	public String toString() {
	    String s = "";
	    for(int i = 0; i < bins.length; i++) {
	        s = s+""+i+" "+i*increment+" "+bins[i]+"\n";
	    }
	    return s;
	}
	
	public double probabilityOfBin(int num) {
	    return bins[num] / totalPoints;
	}
	
	public double valueOfBin(int num) {
	    return num*increment+increment/2+min;
	}
	
	//assumes the min in each bin is 0
	public double maximum() {
	    double max = 0;
	    for(int i=0;i<bins.length;i++){
	        if(bins[i]>max) max=bins[i];
        }
        return max;
	}
}
