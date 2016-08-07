<?xml version='1.0'?>
<b>
<!-- BONDS -->
<bonds>
	<bond func="bond_harmonic(?,20000)">
	<bType>*</bType>
	<bType>*</bType>
	</bond>
        
</bonds>

<!-- ANGLES -->
<angles>
	<angle func="angle_harmonic(?,40)">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</angle>
</angles>

<!-- DIHEDRALS -->

<dihedrals>
	<!-- AMINO DIHEDRALS -->
	<dihedral func="dihedral_cosine(?,1,1)+dihedral_cosine(?,0.5,3)" energyGroup="bb_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_cosine(?,1,1)+dihedral_cosine(?,0.5,3)" energyGroup="sc_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_harmonic(?,40)" energyGroup="pr_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_harmonic(?,10)" energyGroup="r_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<!-- dummy dihedrals for replacement -->
	<dihedral func="dihedral_free()" energyGroup="sc_ah">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_free()" energyGroup="bb_ah">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>

</dihedrals>

<!-- IMPROPERS -->
<impropers>
	<improper func="dihedral_harmonic(?,10)">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</improper>
</impropers>



</b>
