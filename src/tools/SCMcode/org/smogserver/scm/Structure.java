package org.smogserver.scm;
import java.util.*;

public interface Structure {
    public Atom[] getAtoms();
    //public Atom[] getCA();
    public void setChains(String chainfile);
    public Vector getChains();
    public int getNumChains();
}