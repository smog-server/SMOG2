<?xml version='1.0'?>
<b>
<!-- BONDS -->
<bonds>
	<bond func="bond_harmonic(?*0.6,10000)">
	<bType>*</bType>
	<bType>*</bType>
        </bond>
</bonds>

<!-- ANGLES -->
<angles>
	<angle func="angle_harmonic(1.35*?,80)">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</angle>
	<angle func="angle_free()">
	<bType>B_3</bType>
	<bType>B_3</bType>
	<bType>B_3</bType>
	</angle>
</angles>

<!-- DIHEDRALS -->
<dihedrals>	
	<!-- NUCLEIC DIHEDRALS -->
	<dihedral func="dihedral_cosine(?,1,1)+dihedral_cosine(?,0.5,6)" energyGroup="bb_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_cosine(?,1,1)+dihedral_cosine(?,0.5,6)" energyGroup="sc_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_harmonic(?,40)" energyGroup="pr_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<!-- AMINO DIHEDRALS -->
	<dihedral func="dihedral_cosine(?,1,1)+dihedral_cosine(?,0.5,6)" energyGroup="bb_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_cosine(?,1,1)+dihedral_cosine(?,0.5,6)" energyGroup="sc_a">
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

	<dihedral func="dihedral_free()" energyGroup="free">
	<bType>*</bType>
	<bType>B_4</bType>
	<bType>B_4</bType>
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
