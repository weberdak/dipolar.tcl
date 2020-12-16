# dipolar.tcl
VMD script for computing dipolar couplings from MD simulations and PDB files.

## Features

* Calculate dipolar coupling from PDB file or MD trajectories. 
* Supports common spin-half nuclei (1H, 13C, 15N, 19F and 31P).
* Factors in angular term to assess dipolar coupling strength in oriented samples.
* Can ignore the angular term. Useful for predicting dipolar couplings for MAS samples.

## Example usage

To run the script, load a structure into VMD (i.e., PDB or MD trajectory) and open the Tk Console (under Extensions menu). With the script copied into the working directory, load it by: 

	source dipolar.tcl

### Example 1: Simple dipolar coupling calculation

The very simplest use of this script is to compute a dipolar couping from an input internuclear distance and angle using the "coupling" function. For example, let's compute the dipolar coupling of an amide N-H bond (1.04 angstroms), and parallel (0 deg), perpendicular (90 deg) and magic angle (54.74 deg) orientations with respect to the magnetic field (Z-axis of the frame):

	coupling 1.04 0.0 1H 15N
	> 21647.696742495056
	
	coupling 1.04 90.0 1H 15N
	> 10823.848348741496
	
	coupling 1.04 54.74 1H 15N
	> 1.8552543205173964

Note that the absolute dipolar coupling is reported in Hz. Now for some other nuclei:

	coupling 3.52 0.0 1H 1H
	> 5508.304047530693
	
	coupling 10.0 0.0 19F 19F
	> 212.59917189747392
	
	coupling 1.09 0.0 13C 1H
	> 46656.33823222464

Note that setting the angle to 0.0 essentially negates the angular term and reports the full dipolar coupling.

### Example 2: Dipolar coupling between two atoms in a structure

Dipolar couplings of two atoms can be computed using the "pair" function. To use this, simply make two VMD selections of the atoms of interest. Distances, angles and dipolar couplings will be reported as thier average and standard deviation over all frames of the simulation. Standard deviations will be "0.0" if only one frame is used. For example, to compute the dipolar coupling of backbone nitrogens of residue 25 and 26 in an oriented alpha-helix:

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

The dipolar coupling is calculated in two ways. For the first method, the dipolar coupling is calculated at each frame using the instantaneous distances and angles and reported as the mean and standard deviation at the end. In the second method, the dipolar coupling is calculated only at the end once the average distance and angle has been determined. The first method is likely to be more accurate for heterogeous MD ensembles since strong dipolar couplings (non-linear with respect to distance) associated with sporadic atomic contact will be better factored into the final value. The second method may be more applicable for NMR ensembles that are restrained by harmonic terms in which the average distance and angles have more meaning.

Other considerations:

* Note that the atom name still needs to be explicitly defined in the function as it is not automatically recognized.  You may also use another atom type. I.e., select a hydrogen and specify it as 19F.
* Be sure that the selections each only contain one atom, else the center of mass of the selection will be used (which would not really mean much). If unsure, check the number of atoms in the selection by typing "$s1 num" into the console.
* The angular term is only meaningful if structures are correctly oriented with respect to the Z-axis (i.e., a membrane protein measured by oriented solid-state NMR or residual dipolar couplings).
* Order parameters are not factored into the equations.

### Example 3: Dipolar coupling from two sets of atoms

The function "pairlist" builds on the other two function above by automatically computing dipolar couplings of all possible pairwise combinations of atoms included in two selections (homo or heteronuclear). This can take a while depending on how many atoms are included in each selection. Note that the same selection can be used for both lists (combinations involving the same atom will be skipped). For example, to compute homonuclear 15N couplings expected produce peaks in an oriented 15N-PDSD spectrum:

	# Select all backbone nitrogens
	set s1 [atomselect top "name N"]
	
	# Compute dipolar coupling with angular terms
	pairlist $s1 15N $s1 15N -filter 1.0 -out dipolar.dat
	
	# Without anglular dependance
	pairlist $s1 15N $s1 15N -filter 1.0 -noangle -out dipolar_noangle.dat

In the first example, "-filter 1.0" specifies that only dipolar couplings greater the 1.0 Hz will be recorded in the output file. Results or each pair are written to "dipolar.dat" specified by the "-out" flag. Like Example 2, the second run negates the angular term using the "-noangle" flag. For a 15N-PDSD spectrum of TM helical segment, this data could be useful for assessing the effect of orientation  on the crosspeaks used to make sequential assignments. Note that all flags "-out" (default: dipolar.dat), "-filter" (default: 0.0), and "-noangle" (default: off) are optional. 

Dipolar couplings between amide N and H (>500 Hz) are determined in the next example:

	# Select amide protons and nitrogen
	set s1 [atomselect top "name HN"]
	set s2 [atomselect top "name N"]
	pairlist $s1 1H $s2 15N -filter 500

This could be useful for calculation dipolar waves observed by either SLF of DIPSHIFT. However, large errors can be associated with the bond lengths assigned by the forcefield. The [MDSLF](https://github.com/weberdak/md2slf/blob/master/slf.tcl) script would be a better option in this case.

