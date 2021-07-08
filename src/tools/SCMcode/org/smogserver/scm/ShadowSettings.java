package org.smogserver.scm;
import java.util.*;
import org.smogserver.io.*;

class ShadowSettings {
    public static boolean SMOG_ERRORS_ON = false;    
    public static boolean SHOW_PROGRESS = false;
    public static boolean IGNORE_HYDROGEN = false;
	public static boolean BIF_PARSING = false;
	public static boolean SMOG2_OUTPUT_ON = false;
	public static boolean USE_SHADOW_MAP = false;
	public static boolean USE_CUTOFF_MAP = false;
	public static boolean LEGACY_SHADOW = true;
	public static boolean FREE_FORM_COORDINATES = false;
	
	public static double SHADOW_RADIUS = 1.0;
	public static double BONDED_RADIUS = 0.5;
	public static double DISTANCE_CUTOFF = 6.0;
	public static int NUCLEIC_DELTA = 0; 
	public static int PROTEIN_DELTA = 3;
	
	public static int GRO_PRECISION = 3;
	
	
	public static void throwError(int type) {
		System.out.println(" !! SCM ERROR !! ");
		switch (type) {
			case 1:	System.out.println("Unsure of what type of contact map to produce. Both ShadowSettings.USE_SHADOW_MAP or ShadowSettings.USE_CUTOFF_MAP are false. Quitting.");
					System.exit(1);
					break;
			case 2: System.out.println("Both ShadowSettings.USE_SHADOW_MAP and ShadowSettings.USE_CUTOFF_MAP are true. Quitting.");
					System.exit(1);
					break;
			default: System.out.println("SCM threw an error but type was unknown.");
		}
	}
}