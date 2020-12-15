# dipolar.tcl
VMD script for computing dipolar couplings from MD simulations and PDB files.

## Features

* Calculate dipolar coupling from PDB file or MD trajectories. 
* The supports common spin-half nuclei (1H, 13C, 15N, 19F, 31P).
* Factors in angular term to assess dipolar coupling strength in oriented systems.
* Can ignore angular term for predicting dipolar couplings for MAS samples.

## Example usage

To run the script, load a structure into VMD (i.e., PDB or MD trajectory) and open the Tk Console (under Extensions menu). With the script copied into the working directory, load it by: 

	source dipolar.tcl

### Example 1: Simple dipolar coupling calculation

The very simplest use of this script is compute a dipolar couping from an internuclear distance and angle term using the "coupling" function. For example, let's compute the dipolar coupling of an amide N-H bond (1.04 angstroms), parallel (0 deg), perpendicular (90 deg) and at the magic angle (54.74 deg) with respect to the magnetic field:

	coupling 1.04 0.0 1H 15N
	> 21647.696742495056
	
	coupling 1.04 90.0 1H 15N
	> 10823.848348741496
	
	coupling 1.04 54.74 1H 15N
	> 1.8552543205173964

Note that the dipolar coupling is always reported in Hz. Now for some other nuclei:

	coupling 3.52 0.0 1H 1H
	> 5508.304047530693
	
	coupling 10.0 0.0 19F 19F
	> 212.59917189747392
	
	coupling 1.09 0.0 13C 1H
	> 46656.33823222464

Note that setting the angle to 0.0 essentially negates the angular term and reports the full dipolar coupling.

### Example 2: Dipolar coupling between two atoms in a structure

Dipolar couplings for two atoms can be computed using the "pair" function. To use this, simply make two VMD selections for the atoms of interest. Distances, angles and dipolar couplings will be reported as the average and standard deviation over all frames in the simulation. Standard deviations will be "0.0" if only one frame is used. For example, to compute the dipolar coupling of the backbone nitrogens of residue 25 and 26:

	# Make VMD selections
	set s1 [atomselect top "resid 25 and name N"]
	set s2 [atomselect top "resid 26 and name N"]
	
	# Compute dipolar coupling
	pair $s1 15N $s2 15N
	> 2.7953279636930186 0.007089910965287658 94.17377211294986 0.8613440759735665 55.57010945149467 0.513435325952776 55.60290751176622
	
	# Can also compute the coupling without the angle
	pair $s1 15N $s2 15N -noangle
	> 2.7953279636930186 0.007089910965287658 0.0 0.0 113.00414695517762 0.8646612831915333 113.00020317610571

The function will return a list of 7 numbers, which for the first calculations is:

* 2.7953279636930186 0.007089910965287658: Mean and stdev of distance (in angstroms)
* 94.17377211294986 0.8613440759735665: Mean and stdev of angle (in degrees)
* 55.57010945149467 0.513435325952776: Mean and stdev of dipolar coupling (in Hz)
* 55.60290751176622: Dipolar couplings computed from final averaged distance and angle (in Hz)

The dipolar coupling is calculated in two ways. For the first method, the dipolar coupling is calculated at each frame using the instantaneous distances and angles and reported as the mean and standard deviation. In the second method, the dipolar coupling is calculated at the end once the average distance and angle has been determined. The first method is likely more accurate for heterogeous MD ensembles since strong dipolar couplings (non-linear with respect to distance) associated with sporadic atomic contact can be factored into the final value. The second method may be more applicable for NMR ensembles that are restrained by harmonic terms.

Other considerations:

* Note that the atom name still needs to be explicitly defined in the function as it is not automatically recognized. 
* Be sure that the selections each only contain one atom, else the center of mass will be used (which would not really mean much). If unsure, check the number of atoms in the selection by typing "$s1 num" into the console.
* The angular term is only meaningful if structures are correctly oriented with respect to the Z-axis (i.e., a membrane protein measured by oriented solid-state NMR or residual dipolar couplings).
* Order parameters are not factored into the equations.

### Example 3: Dipolar coupling from two sets of atoms

The function "pairlist" builds on the others by automatically computing dipolar couplings of all possible pairwise combinations of atoms included in two selection (homo or heteronuclear). This can take a while depending on how many atoms are included in each selection. Note that the same selection can be used for both lists. For example, to compute the homonuclear 15N couplings for assessing an oriented 15N-PDSD spectrum:

	# Select all backbone nitrogens
	set s1 [atomselect top "name N"]
	
	# Compute dipolar coupling with angular terms
	pairlist $s1 15N $s1 15N -filter 1.0 -out dipolar.dat

	# Without anglular dependance
	pairlist $s1 15N $s1 15N -filter 1.0 -noangle -out dipolar_noangle.dat

 
