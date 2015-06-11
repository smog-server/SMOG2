<?xml version='1.0'?>
<b>
<!-- BONDS -->
<bonds>
	<bond func="sbm_bonds(?,20000)">
	<bType>*</bType>
	<bType>*</bType>
        </bond>
	<bond func="sbm_bonds_6(?,200)">
	<bType>*</bType>
    	<bType>MG</bType>
        </bond>
</bonds>

<!-- ANGLES -->
<angles>
	<angle func="sbm_angles(?,40)">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</angle>
</angles>

<!-- DIHEDRALS -->
<dihedrals>	
	<!-- NUCLEIC DIHEDRALS -->
	<dihedral func="sbm_dihedrals(?,?,1)+sbm_dihedrals(?,?*0.5,3)" energyGroup="bb_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="sbm_dihedrals(?,?,1)+sbm_dihedrals(?,?*0.5,3)" energyGroup="sc_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="sbm_planarRigid(?,40)" energyGroup="pr_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="sbm_rigid(?,10)" energyGroup="r_n">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<!-- AMINO DIHEDRALS -->
	<dihedral func="sbm_dihedrals(?,?,1)+sbm_dihedrals(?,?*0.5,3)" energyGroup="bb_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="sbm_dihedrals(?,?,1)+sbm_dihedrals(?,?*0.5,3)" energyGroup="sc_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="sbm_planarRigid(?,40)" energyGroup="pr_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="sbm_rigid(?,10)" energyGroup="r_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<!-- LIGAND DIHEDRALS -->
        <dihedral func="sbm_rigid(?,10)" energyGroup="lig">
        <bType>*</bType>
        <bType>*</bType>
        <bType>*</bType>
        <bType>*</bType>
        </dihedral>
</dihedrals>

<!-- IMPROPERS -->
<impropers>
	<improper func="sbm_improper(?,10)">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</improper>
</impropers>

</b>
