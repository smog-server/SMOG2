package noel.util.geom;
import noel.util.math.*;

// A collection of geometry routines
public class Geometry {
	static double PI = Math.PI;
	/** 
* Takes three coordinates and returns the angle opposite of the last two coordinates.
*/
	public static double getOppositeAngle(Coordinates a, Coordinates b, Coordinates c) {
		double ab = Coordinates.distance(a,b);
		double bc = Coordinates.distance(b,c);
		double ca = Coordinates.distance(c,a);
		double theta = Math.acos(-(bc*bc-ab*ab-ca*ca)/(2*ab*ca));
		if (b.y > c.y) {
			return (2*PI-theta)*180./PI;
		} else return theta*180./PI;	
	}
	/**
	* Returns a congruent triangle in 2D, first point is the origin.
*/
	public static Triangle get2DTriangle(Triangle t) {
		double a,b,c;
		Coordinates[] coords = new Coordinates[3];
		Coordinates[] coord = t.getCoords();
		a = Coordinates.distance(coord[0],coord[1]);
		b = Coordinates.distance(coord[1],coord[2]);
		c = Coordinates.distance(coord[2],coord[0]);
		coords[0] = new Coordinates(0,0,0);
		coords[1] = new Coordinates(a,0,0);
		double c1 = (c*c-b*b+a*a)/(2*a); 
		coords[2] = new Coordinates(c1,Math.sqrt(c*c-c1*c1),0);
		return new Triangle(coords);
	}
/**
* Returns whether Sphere b is shadowed by Sphere a from the perspective of point p.
*/
	public static boolean isShadowed(Coordinates p, Sphere a, Sphere b) {
		Coordinates ap = a.getCenter();
		Coordinates bp = b.getCenter();
		Triangle d3 = new Triangle(p,ap,bp);
		Coordinates[] d2 = get2DTriangle(d3).getCoords();
		double thetaA = Math.atan(d2[1].y / d2[1].x);
		double thetaB = Math.atan(d2[2].y / d2[2].x);
		double DthetaA = Math.asin(a.getRadius()/Coordinates.distance(d2[0],d2[1]));
		double DthetaB = Math.asin(b.getRadius()/Coordinates.distance(d2[0],d2[2]));
		if (thetaA < 0) thetaA = PI+thetaA;
		if (thetaB < 0) thetaB = PI+thetaB;
		double AL = thetaA - DthetaA;
		double AU = thetaA + DthetaA;
		double BL = thetaB - DthetaB;
		double BU = thetaB + DthetaB;
		if (BL > AL && BL < AU) return true; //is my lower bound sandwiched?
		if (BU > AL && BU < AU) return true; //is my upper bound snadwiched?
		if (BL < AL && BU > AU) return true; //am i enveloping the shadower?
		return false;
	}
	
/**
* LEGACY METHOD - there is an error where indicated, but allowing this code to
* remain in case reproducibility is important in the future.
*
* Returns whether Sphere b is shadowed by Sphere a from the perspective of point p.
*/
	public static boolean isShadowedLegacy(Coordinates p, Sphere a, Sphere b) {
		Coordinates ap = a.getCenter();
		Coordinates bp = b.getCenter();
		Triangle d3 = new Triangle(p,ap,bp);
		Coordinates[] d2 = get2DTriangle(d3).getCoords();
		double thetaA = Math.atan(d2[1].y / d2[1].x);
		double thetaB = Math.atan(d2[2].y / d2[2].x);
		//Next two lines should be asin!!
		double DthetaA = Math.atan(a.getRadius()/Coordinates.distance(d2[0],d2[1]));
		double DthetaB = Math.atan(b.getRadius()/Coordinates.distance(d2[0],d2[2]));
		if (thetaA < 0) thetaA = PI+thetaA;
		if (thetaB < 0) thetaB = PI+thetaB;
		double AL = thetaA - DthetaA;
		double AU = thetaA + DthetaA;
		double BL = thetaB - DthetaB;
		double BU = thetaB + DthetaB;
		if (BL > AL && BL < AU) return true; //is my lower bound sandwiched?
		if (BU > AL && BU < AU) return true; //is my upper bound snadwiched?
		if (BL < AL && BU > AU) return true; //am i enveloping the shadower?
		return false;
	}
	
	/**
	*   This routine returns whether the two line segments appear to intersect when 
	*	viewed (from far away) along a given vector (c). Actual intersection is
	*	not tested for because we are in 3D. Routine created for testing whether two dynamin dimer links
	*	are intersecting. Consider segments a and b, with points (a1,a2) and (b1,b2). 
	*
	*	Algorithm:
	*	1. (vec(a1,a2) x vec(a2,b1)) * c has a different sign than (vec(a1,a2) x vec(a2,b2)) * c
	*	2. (vec(b1,b2) x vec(b2,a1)) * c has a different sign than (vec(b1,b2) x vec(b2,a2)) * c
	*	3. if(1 and 2) then true else false.
	*	4. special cases: zero length, parallel should be handled without issue
	*/
	public static boolean isApparentIntersection(Coordinates a1, Coordinates a2, Coordinates b1, Coordinates b2, 
		NVector c) {
			//(vec(a1,a2) x vec(a2,b1)) * c has a different sign than (vec(a1,a2) x vec(a2,b2)) * c
			NVector e = new NVector(a1,a2);
			NVector f = new NVector(a2,b1);
			NVector g = new NVector(a2,b2);
			double sgn1 = Math.signum(e.cross(f).dot(c));
			double sgn2 = Math.signum(e.cross(g).dot(c));			
			if(sgn1 == sgn2) return false;
			//(vec(b1,b2) x vec(b2,a1)) * c has a different sign than (vec(b1,b2) x vec(b2,a2)) * c
			e = new NVector(b1,b2);
			f = new NVector(b2,a1);
			g = new NVector(b2,a2);
			sgn1 = Math.signum(e.cross(f).dot(c));
			sgn2 = Math.signum(e.cross(g).dot(c));			
			if(sgn1 == sgn2) return false;
			return true;
	}
		
		

}
