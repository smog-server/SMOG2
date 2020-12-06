package noel.folding;

public interface Topology {
public BondedList getBondedList();
public Atom[] getAtoms();
public int[][] getBonds();
public int[][] getAngles();
public int[][] getDihedrals();
}
