package noel.folding;
import noel.util.math.*;
import noel.util.geom.*;

/**
 * Contact between atom left and atom right.  left.index < right.index
 */
public class Contact {
	public Atom left;
	public Atom right;
	public double distance;
	public double strength;
	public int leftIndex;
	public int rightIndex;

	public Contact(Atom left, Atom right, double distance, double strength) {
		this.left = left;
		this.right = right;
		this.distance = distance;
		this.strength = strength;
		this.leftIndex = left.position;
		this.rightIndex = right.position;
	}

	public Contact(int leftIndex, int rightIndex, double distance, double strength) {
		this.leftIndex = leftIndex;
		this.rightIndex = rightIndex;
		this.distance = distance;
		this.strength = strength;
	}

	public String skinnyToString() {
	                        return leftIndex+" "+rightIndex;
}
	public String toString() {
			return left.getChain()+" "+leftIndex+" "+right.getChain()+" "+rightIndex+" "+distance+" "+strength;
	}
}
