package org.smogserver.scm;

public interface Topology {
public BondedList getBondedList();
public Atom[] getAtoms();
public int[][] getBonds();
public int[][] getAngles();
public int[][] getDihedrals();
}
