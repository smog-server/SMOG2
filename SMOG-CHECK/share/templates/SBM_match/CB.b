<?xml version='1.0'?>
<b>
<!-- BONDS -->
<bonds>
	<bond func="bond_harmonic(0.12,20000)">
	<bType>X</bType>
	<bType>*</bType>
	</bond>
	<bond func="bond_harmonic(0.15,50000)">
	<bType>X1</bType>
	<bType>X1</bType>
	</bond>
	<bond func="bond_harmonic(0.1,2000)">
	<bType>*</bType>
	<bType>*</bType>
	</bond>
        
</bonds>

<!-- ANGLES -->
<angles>
	<angle func="angle_harmonic(0.5,40)">
	<bType>*</bType>
	<bType>*</bType>
	<bType>X</bType>
	</angle>
	<angle func="angle_harmonic(0.6,10)">
	<bType>*</bType>
	<bType>X</bType>
	<bType>X</bType>
	</angle>
	<angle func="angle_harmonic(0.7,4)">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</angle>
</angles>

<!-- DIHEDRALS -->
<dihedrals>																															
<!-- AMINO DIHEDRALS -->
	<dihedral func="dihedral_cosine(1,0.1,1)+dihedral_cosine(1,0.05,3)" energyGroup="bb">
	<bType>X1</bType>
	<bType>X</bType>
	<bType>X</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_cosine(10,7,1)+dihedral_cosine(10,3.5,3)" energyGroup="bb">
	<bType>X</bType>
	<bType>X</bType>
	<bType>X</bType>
	<bType>X</bType>
	</dihedral>
	<dihedral func="dihedral_cosine(1,100,1)+dihedral_cosine(1,50,3)" energyGroup="bb">
	<bType>*</bType>
	<bType>X1</bType>
	<bType>X1</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_cosine(1,1,1)+dihedral_cosine(1,0.5,3)" energyGroup="bb">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>


	<dihedral func="dihedral_cosine(1,1,1)+dihedral_cosine(1,0.5,3)" energyGroup="sc">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
</dihedrals>
</b>
