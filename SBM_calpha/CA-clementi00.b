<?xml version='1.0'?>
<b>
<!-- BONDS -->
<bonds>
	<bond func="sbm_bonds(?,20000)">
	<bType>*</bType>
	<bType>*</bType>
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
<!-- AMINO DIHEDRALS -->
	<dihedral func="sbm_dihedrals(?,1,1)+sbm_dihedrals(?*3,0.5,3)" energyGroup="bb_a">
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	<bType>*</bType>
	</dihedral>
	<dihedral func="sbm_dihedrals(?,1,1)+sbm_dihedrals(?*3,0.5,3)" energyGroup="sc_a">
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
