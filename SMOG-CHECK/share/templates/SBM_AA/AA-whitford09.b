<?xml version='1.0'?>
<b>
<!-- BONDS -->
<bonds>
	<bond func="bond_harmonic(2.0*(?/2),10000)">
	<bType>*</bType>
	<bType>*</bType>
	</bond>
        <bond func="bond_harmonic(0.21,10000)">
	<bType>B_FES</bType>
	<bType>B_FES</bType>
        </bond>
</bonds>

<!-- ANGLES -->
<angles>
	<angle func="angle_harmonic((?)/(1.0*1.0),80)">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</angle>
</angles>

<!-- DIHEDRALS -->
<dihedrals>	
	<!-- NUCLEIC DIHEDRALS -->
	<dihedral func="dihedral_cosine(?*(3.0**0),1,1)+dihedral_cosine(?/1.0,0.5,3)" energyGroup="bb_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_cosine(1.0*?,1,1)+dihedral_cosine(1*?,0.5,3)" energyGroup="sc_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="dihedral_harmonic(1.0*?/1.0,40)" energyGroup="pr_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<!-- AMINO DIHEDRALS -->
	<dihedral func="dihedral_cosine(?*1.0**2,1,1)+dihedral_cosine(?,0.5,3)" energyGroup="bb_a">
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

	<!-- LIGAND DIHEDRALS -->
        <dihedral func="dihedral_harmonic(?,40)" energyGroup="lig">
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
