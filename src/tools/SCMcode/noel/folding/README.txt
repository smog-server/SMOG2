Welcome to the shadow contact map (SCM) Java application.  This code is 
tailored for use with Paul Whitford's READ application.

Three main uses:
1. Convert a .gro to a PDB.
2. Get a cutoff contact map.
3. Get a shadow map.

Quick example for CI2:
1. Compile READ and run the AA version on CI2.pdb.  
	a. Choose contact map option 1 in AAsettings.dat to let READ 
	   choose contacts, we will overwrite them later.

2. Use the output .gro and .top file as input to SCM.  SCM will use the
   .gro to get coordinates and the .top to create a bonded list.  Let us
   assume that the output is called CI2.gro and CI2.top.

3. Run "java -jar $SCM_DIR/SCM.jar -help" and look at the options

4. Run SCM with the following command:

java -jar $SCM_DIR/SCM.jar -g CI2.gro -t CI2.top -m shadow -c 6.0 -s 1.0 \
                           -p CI2.pdb -o CI2.shadow

Output to standard out should be:
Reading grofile: CI2.gro.
Found 504 atoms.
Reading Topology.
Initializing connected list.
Found 481 pairs.
Found 511 bonds.
Found 693 angles.
Found 1522 dihedrals.
Calculating contacts.
Printing PDB to file: CI2.pdb.
Finished.  Printing contacts to file: CI2.shadow.

5. If your PDB has multiple chains, READ should output an index file.  You
   can include this information in SCM with the option "--chain indexfile"

6. Now change contact map setting in AAsettings.dat to 2.  Use the PDB file
   outputted by SCM as input to READ (to keep the indexing constant).  Put
   the new contact map in Readsettings.dat.  Run READ again and the output 
   is ready for GROMACS!

HOW SCM WORKS:
   * Determines which atoms are in contact using the following algorithm:
   * 1. Checks if two atoms are within cutoff (-c) angstroms.
   * 2. Checks if two atoms are greater than the specified number of
   *    residues apart or in separate chains.
   * 3. If so then checks to see if any atoms are
   *    between the two.  This is determined by
   *    putting a light bulb at the center of each of
   *    the atoms and seeing if any of the
   *    other atoms cast a shadow on the potential contacting atom.
   *    If so then the contact is thrown out.
   * 3b. The size of the atoms is
   *     determined by size (-s).
   * 3c. A potentially shadowing
   *     atom is given a size of 0.5
   *     angstrom when it is bonded
   *     to the shadowee.  This is
   *     because the C-C bond
   *     length is 1.5 angstrom
   *     and bonded neighbors should not shadow.
