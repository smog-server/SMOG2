package org.smogserver.util.geom;   
import org.smogserver.util.math.*;

public class Triangle {

	Coordinates[] coords;

	public Triangle(Coordinates[] coords) {
		if (coords.length != 3) { 
			System.out.println("Triangle has "+coords.length+" coordinates!");
		}
		this.coords = new Coordinates[3];
		for (int i = 0; i < 3; i++) {
			this.coords[i] = new Coordinates(coords[i].x,coords[i].y,coords[i].z);
		}
	}
	public Triangle(Coordinates a, Coordinates b, Coordinates c) {
		this.coords = new Coordinates[3];
		coords[0] = new Coordinates(a.x,a.y,a.z);
		coords[1] = new Coordinates(b.x,b.y,b.z);
		coords[2] = new Coordinates(c.x,c.y,c.z); 
	}
	public String toString() {
		String string = "";
		for (int i = 0; i < 3; i++) {
			string+=coords[i].toString()+"\n";
		}
		return string;
	}
	public Coordinates[] getCoords() {
		return coords;
	}
	public NVector getNormal() {
	    NVector side1 = new NVector(coords[0],coords[1]);
	    NVector side2 = new NVector(coords[0],coords[2]);
        NVector perp = side1.cross(side2);
        return perp.mul(1/perp.mag());
	}
}
