# EPIC
Electron-phonon coupling in narrow bands intermetallics

EPIC contains a suite of scripts that can be used to compute the electron-phonon coupling (EPC) and parametrize it into deformation potential values, for narrow band systems. EPIC extends the methodology contained in the EPHOS repository for organic semiconductor. The EPC is computed as a bandwidth modulation over the phonon eigenmode coordinates, thus it can be applied when the bands are narrow. Mode-resolved EPC values are then DOS-weighted averaged over the BZ points. 

The current version of EPIC is under development and testing.

To be used:
- After having run the phonon calculations, for each q-point of interest (typically Gamma and zone boundary) and phonon branch or mode, modulate the structure and run, for each q-point/mode pair run a scf and nscf calc. and then save the bandstructure in the bxsf format. You find example files in the VASP_phonopy_examples directory.
- Copy in the modulations directory the scripts:  analyze_extract_bands_DefPot.sh bands_mod_extract_DefPot_narrowbands_shifted_script_L.m  bxsf_to_ELECTRA_shifting_option.m  shifting_bands.m
- Now, you fill the analyze_extract_bands_DefPot.sh with the proper data, current lines 14-22. You must edit
      n_q → number of q-points considered
      n_modes → nu,ber of modes you are interested in
      n_atoms → number of atoms in the POSCAR or other structure file
      modes_to_skip → if you want to avoid some modes
      fileName_bxsf → the name of the bxsf file, use the form as in the provided file "'name'"
      material_name → name of the material you will use for subsequent ELECTRA run
     E_band_window → you can change the energy range, in ev
     labels_q → you can list the labels of the q-point to further save in table format
     force_metallic → to force a semimetal or very narrow gap bandstructure to be treated as metallic, to allow scattering transition through the narrow gap if the phonon energy is larger than the gap - the system will must be set as metallic in ELECTRA calc.. In case, it must be "'yes'", otherwise "'no'"
  Finally launch it. The saved file contains the Deformation Potential in structures which are indexed for the number of q-points and modes.
- The script overall_EPC_narrow_bands_wADP extracts these values for each mode and perform a mode-resolved DOS-weighted average of the values and save it in a format suitable for the ELECTRA scattering parameter file. The lines 5-17 contains the relevant informations, files path, eventual cut at a certain frequency, data for ADP. Line 21 allows the user to enter the indexes of the bands that are actually used by ELECTRA, the variable bands_transp. This is important for metals with many bands because it is likely that the bands considered to parametrize the EPC, and saved in the DeformationPotentials.mat file, are more than the ones considered in the transport and since the inter-band scattering parameters are sensitive to that, it is important to restrict them to the range used. Moreover, you can freely change it without re-computing the EPC. It is not relevant in semiconductors.




