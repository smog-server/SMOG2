<?xml version="1.0"?>
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
    <dihedral func="dihedral_cosine(?,1,1)+dihedral_cosine(?,0.5,3)" energyGroup="bb">
      <bType>*</bType>
      <bType>*</bType>
      <bType>*</bType>
      <bType>*</bType>
    </dihedral>
  </dihedrals>
</b>
