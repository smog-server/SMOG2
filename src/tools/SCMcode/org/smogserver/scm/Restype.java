package org.smogserver.scm;

/**
 * Defines all the valid residue names and whether they are Nucleic Acids, Proteins or Ligands.
 * Use this class like this: this.type = Restype.restype(residueName) 
 */
public enum Restype {

ADE (Residue.NUCLEIC_ACID),
CYT (Residue.NUCLEIC_ACID),
URA (Residue.NUCLEIC_ACID),
GUA (Residue.NUCLEIC_ACID),
THY (Residue.NUCLEIC_ACID),
A   (Residue.NUCLEIC_ACID),
C (Residue.NUCLEIC_ACID),
G (Residue.NUCLEIC_ACID),
U (Residue.NUCLEIC_ACID),
T  (Residue.NUCLEIC_ACID),
DA (Residue.NUCLEIC_ACID),
DC (Residue.NUCLEIC_ACID),
DG  (Residue.NUCLEIC_ACID),
DT  (Residue.NUCLEIC_ACID),
A5   (Residue.NUCLEIC_ACID),
C5 (Residue.NUCLEIC_ACID),
G5 (Residue.NUCLEIC_ACID),
U5 (Residue.NUCLEIC_ACID),
T5  (Residue.NUCLEIC_ACID),
DA5 (Residue.NUCLEIC_ACID),
DC5 (Residue.NUCLEIC_ACID),
DG5  (Residue.NUCLEIC_ACID),
DT5  (Residue.NUCLEIC_ACID),
A3   (Residue.NUCLEIC_ACID),
C3 (Residue.NUCLEIC_ACID),
G3 (Residue.NUCLEIC_ACID),
U3 (Residue.NUCLEIC_ACID),
T3  (Residue.NUCLEIC_ACID),
DA3 (Residue.NUCLEIC_ACID),
DC3 (Residue.NUCLEIC_ACID),
DG3  (Residue.NUCLEIC_ACID),
DT3  (Residue.NUCLEIC_ACID),
MIA  (Residue.NUCLEIC_ACID),
ALA (Residue.PROTEIN),
ARG (Residue.PROTEIN),
ASN (Residue.PROTEIN),
ASP (Residue.PROTEIN),
CYS (Residue.PROTEIN),
GLN (Residue.PROTEIN),
GLU (Residue.PROTEIN),
GLY (Residue.PROTEIN),
HIS (Residue.PROTEIN),
ILE (Residue.PROTEIN),
LEU (Residue.PROTEIN),
LYS (Residue.PROTEIN),
MET (Residue.PROTEIN),
PHE (Residue.PROTEIN),
PRO (Residue.PROTEIN),
SER (Residue.PROTEIN),
THR (Residue.PROTEIN),
TRP (Residue.PROTEIN),
TYR (Residue.PROTEIN),
VAL (Residue.PROTEIN),
DBZ (Residue.PROTEIN),
CBS (Residue.PROTEIN), //Sugar 1GYA_chitobiose.pdb
ACE (Residue.PROTEIN), //N-terminal acetylation
SAM (Residue.LIGAND),
GNP (Residue.LIGAND),
ATP (Residue.LIGAND),
ADP (Residue.LIGAND),
AMP (Residue.LIGAND),
FUA (Residue.LIGAND),
GTP (Residue.LIGAND),
GDP (Residue.LIGAND),
BMG (Residue.LIGAND),
B12 (Residue.LIGAND),
ZN  (Residue.LIGAND),
UNDEFINED ( Residue.LIGAND );
    
    private final int restype; 
    Restype(int restype) {
        this.restype = restype;
    }
    public static int restype(String str)   { 
	try {
		return valueOf(str).restype; 
	} catch (Exception ex) { return UNDEFINED.restype; }
   }
}


