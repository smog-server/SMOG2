<?xml version='1.0'?>
<b>
<!-- BONDS -->
<bonds>
	<bond func="bond_harmonic(?,20000)">
	<bType>*</bType>
	<bType>*</bType>
        </bond>
	<bond func="bond_type6(?,200)">
	<bType>*</bType>
    	<bType>MG</bType>
        </bond>
</bonds>

<!-- ANGLES -->
<angles>
	<angle func="angle_free()">
		<bType>*</bType>
        <bType>*</bType>
        <bType>*</bType>
	</angle>

</angles>

<!-- DIHEDRALS -->
<dihedrals>	
	<dihedral func="dihedral_free()" energyGroup="free">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<!-- NUCLEIC DIHEDRALS -->
	<dihedral func="dihedral_harmonic(?,40)" energyGroup="pr_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_harmonic(?,10)" energyGroup="r_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>

	<!-- AMINO DIHEDRALS -->
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
	<!-- LIGAND DIHEDRALS -->
    <dihedral func="dihedral_harmonic(?,10)" energyGroup="lig">
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
