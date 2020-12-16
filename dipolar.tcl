# dipolar.tcl
# ---------------
# Compute pairwise dipolar couplings for oriented and unoriented proteins.
# Written by D. K. Weber (Veglia Lab)
# Last revised: Dec 15 2020
#
# Features
# --------
# - proc coupling: Compute homo/heteronuclear dipolar coupling from distance, angle and gyros.
# - proc pair: Compute dipolar couplings for atom pair over trajectory.
# - proc pairlist: Compute pair-wise dipolar couplings for two atom lists over trajectory.

proc coupling { dist angl type1 type2 } {
    # Compute dipolar coupling from distance and angle
    #
    # Parameters
    # ----------
    # dist: float
    #    Internuclear distance in angstroms.
    # angl: float
    #    Internuclear angle from Z axis in degrees.
    # type1: str
    #    Nucleus of atom 1 (1H, 13C, 15N, 19F or 31P)
    # type2: str
    #    Nucleus of atom 2 (1H, 13C, 15N, 19F or 31P)
    #
    # Returns
    # -------
    # dc: float
    #    Dipolar coupling in Hz
    #
    # Example - Amide N-H coupling perpendicular to B0
    # ------------------------------------------------
    # > dipolar_coupling 1.04 90.0 1H 15N
    # 10823.848348741496

    # Constants
    set pi 3.14159265358979
    set mu [expr {4*$pi/10000000}]; # Vacuum permeability (N/A^2) = 4*pi/10000000
    set h 6.62607004E-034; # Planck constant (J S)
    
    # Gyromagnetic ratios (rad s-1 T-1) https://en.wikipedia.org/wiki/Gyromagnetic_ratio
    set gyro(1H)  267522187.44
    set gyro(13C)  67282800.00
    set gyro(15N) -27116000.00
    set gyro(19F) 251662000.00
    set gyro(31P) 108291000.00

    # Convert inputs
    set angl_rad [expr {$angl*0.017453}]; # degrees to radians
    set dist_met [expr {$dist*(10E-011)}]; # distance in meters

    # Compute dipolar coupling from distance and anglular components
    set c [expr {($mu*$h*-1)/(8*($pi**3))}]
    set d [expr {($gyro($type1)*$gyro($type2))/$dist_met**3}]
    set dc_dist [expr {$c*$d}]
    set dc_angl [expr {(3*(cos($angl_rad)**2)-1)/2}]
    set dc [expr {abs($dc_dist*$dc_angl)}]; # Output only positive values
    #puts $dist_met
    #puts $c
    #puts $d
    #puts $dc_dist
    #puts $dc_angl
    #puts $dc
    return $dc 
}


proc pair { atom1 type1 atom2 type2 { args } } {
    # Compute dipolar coupling from average distance and angle over trajectory
    #
    # Parameters
    # ----------
    # atom1: obj
    #    VMD selection of atom 1 (one atom only, otherwise uses center-of-mass)
    # type1: str
    #    Nucleus of atom 1 (1H, 13C, 15N, 19F or 31P)
    # atom2: obj
    #    VMD selection of atom 2 (one atom only, otherwise uses center-of-mass)
    # type2: str
    #    Nucleus of atom 2 (1H, 13C, 15N, 19F or 31P)
    #
    # Optional
    # --------
    # -noangle
    #    Ignore angular term (sets to 0.0 degrees) in dipolar coupling
    #
    # Returns
    # -------
    # result: list of floats
    #    0: avg distance (angstroms)
    #    1: std distances (angstroms)
    #    2: avg angle (degrees)
    #    3: std angles (degrees)
    #    4: avg dipolar coupling (Hz)
    #    5: std dipolar couplings (Hz)
    #    6: dipolar coupling (Hz) - from avg distance and angle

    # Handle optional arguments
    set noangle [ switch_on $args "-noangle" off ]
    
    # Intialize results lists
    set dists []
    set angls []
    set dipcs []
       
    # Loop through each frame and record distance and angle
    set num_frames [molinfo top get numframes]
    for {set i 0} {$i < $num_frames} {incr i} {
	$atom1 frame $i
	$atom2 frame $i

	# Get coordinates of atoms
	set a1_c [measure center $atom1]
	set a2_c [measure center $atom2]

	# Compute internuclear distance
	set vs [vecsub $a2_c $a1_c]
	set vd [veclength $vs]
	
	# Compute internuclear angle (set to 0.0 if -noangle set)
	if { $noangle == on } {
	    set va 0.0
	} else {
	    set va [expr acos([lindex [vecnorm $vs] 2]) * 57.2957795130823]
	}
	
	# Compute dipolar coupling for frame
	set dc [coupling $vd $va $type1 $type2]

	# Record results
	lappend dists $vd
	lappend angls $va
	lappend dipcs $dc
    }
    
    # Compute averages and standard deviations
    set stats_dists [stats $dists]
    set stats_angls [stats $angls]
    set stats_dipcs [stats $dipcs]

    # Compute dipolar coupling from average dist and angle
    set dc_f [coupling [lindex $stats_dists 0] [lindex $stats_angls 0] $type1 $type2]

    # Return results list
    set da [lindex $stats_dists 0]
    set ds [lindex $stats_dists 1]
    set aa [lindex $stats_angls 0]
    set as [lindex $stats_angls 1]
    set dca [lindex $stats_dipcs 0]
    set dcs [lindex $stats_dipcs 1]
    return [list $da $ds $aa $as $dca $dcs $dc_f]
}


proc pairlist { selection1 type1 selection2 type2 { args } } {
    # Centers helix selection to origin. Aligns along Z. Measures tilt and azimuthal angle.
    #
    # Parameters
    # ----------
    # selection1: obj
    #    VMD selections of first set of atoms - should be of the same nucleus type.
    # type1: str
    #    Nucleus type of selection1 (1H, 13C, 15N, 19F or 31P).
    # selection2: obj
    #    VMD selections of second set of atoms - should be of the same nucleus type.
    # type2: str
    #    Nucleus type of selection2 (1H, 13C, 15N, 19F or 31P).
    #
    # Optional
    # --------
    # -noangle
    #    Ignore angular term (sets to 0.0 degrees) in dipolar coupling.
    # -filter <coupling>
    #    <coupling>: Minimum dipolar coupling (in Hz) to be recorded in output file. Default: null.
    # -out: <file name>
    #    <file name>: Name of file to output tilt and azimuthal angles. Default: null.
    
    
    # Handle optional arguments
    set noangle [ switch_on $args "-noangle" off ]
    set out [argparse $args "-out" 1 "dipolar.dat"]
    set filter [argparse $args "-filter" 1 "null"]
    
    # Initialize output file
    set outf [open $out w]
    puts $outf "# Atom1\tAtom2\t\tDistance\tAngle\t\tCoupling\tCouplingF"
    
    # Get list of atom indices for each selection
    set ind1 [$selection1 get index]
    set ind2 [$selection2 get index]

    # Intialize Progress meter
    set n1 [llength $ind1]
    puts "$n1 atoms in selection 1 assigned type $type1"
    set n2 [llength $ind2]
    puts "$n2 atoms in selection 2 assigned type $type2"
    set tot 0
    foreach i $ind1 { foreach j $ind2 { if { $i != $j } { incr tot }}}
    set num_frames [molinfo top get numframes]
    puts "Computing dipolar couplings for $tot internuclear pairs over $num_frames frames"

    # Iterate through each pairwise combination
    set c 0
    set ignored 0
    set recorded 0
    foreach i $ind1 {
	foreach j $ind2 {
	    if { $i != $j } {
		set a1 [atomselect top "index $i"]
		set a2 [atomselect top "index $j"]

		if { $noangle == on } {
		    set r [pair $a1 $type1 $a2 $type2 -noangle]
		} else {
		    set r [pair $a1 $type1 $a2 $type2]
		}
		
		    
		# Clean up results for output
		set rd1 [lindex [$a1 get resid] 0]
		set rd2 [lindex [$a2 get resid] 0]
		set rn1 [lindex [$a1 get resname] 0]
		set rn2 [lindex [$a2 get resname] 0]
		set n1 [lindex [$a1 get name] 0]
		set n2 [lindex [$a2 get name] 0]
		set da [format "%.2f" [lindex $r 0]]
		set ds [format "%.2f" [lindex $r 1]]
		set aa [format "%.1f" [lindex $r 2]]
		set as [format "%.1f" [lindex $r 3]]
		set dca [format "%.1f" [lindex $r 4]]
		set dcs [format "%.1f" [lindex $r 5]]
		set dcf [format "%.1f" [lindex $r 6]]

		if { $filter != "null" } {
		    if { [lindex $r 4] >= $filter } {
			puts $outf "$rn1-$rd1-$n1\t$rn2-$rd2-$n2\t$da +/- $ds\t$aa +/- $as\t$dca +/- $dcs\t$dcf"
			incr recorded
		    } else { incr ignored }
		} else {
		    puts $outf "$rn1-$rd1-$n1\t$rn2-$rd2-$n2\t$da +/- $ds\t$aa +/- $as\t$dca +/- $dcs\t$dcf"
		    incr recorded
		}
		    
		# Clear memory from temporary selections
		$a1 delete
		$a2 delete
		unset a1 a2

		# Update progress meter
		incr c
		if { $tot >= 10 } { progress $c $tot }
	    }
	}
    }
    puts ""
    puts "Output $recorded couplings to $out"
    puts "Ignored $ignored couplings less than $filter Hz"
    close $outf
}


proc argparse { args flag input_index default } {
    # Crude argument parser that does the trick.
    #
    # Parameters
    # ----------
    # args: list of str and values
    #    List of arguments and values
    # flag: str
    #    <args> are searched for the presence of this flag
    # input_index: int
    #    Position that value occurs after <flag>. I.e., if "1", the
    #    value immediately following <flag> will be returned.
    # default: anything
    #    If <flag> isn't found in <args>, then this value is returned.
    #
    # Returns
    # -------
    # value: anything
    #    The value parsed in from <args> of default.
   
    set value $default
    if { [ lsearch $args $flag ] != -1 } { 
	set value [lindex $args [expr ([lsearch $args $flag ] + $input_index)]]
    } 
    return $value
}


# Proc to handle switched
proc switch_on { args flag default } {
    set value $default
    if { [ lsearch $args $flag ] != -1 } { 
	set value on
    }
    return $value
}


proc stats { data } {
    # Compute mean and standard deviation of list
    #
    # Parameters
    # ----------
    # data: list of floats
    #    List of numbers.
    #
    # Returns
    # -------
    # [mean, std]: list of two floats
    #    Mean and standard deviation.
    
    # Comute mean
    set csum 0
    set csum2 0
    set num [llength $data]
    foreach i $data {
	set csum [expr $csum + $i]
	set csum2 [expr $csum2 + pow($i,2)] 
    }
    set mean [expr $csum / $num]
    if { $num < 2 } {
	set std 0
    } else {
	set std [expr sqrt((($num*$csum2)-pow($csum,2))/($num*($num-1)))]
    }
    return [list $mean $std]
}


proc progress {cur tot} {
    # Progress meter (Modified from http://wiki.tcl.tk/16939)
    #
    # Parameters
    # ----------
    # cur: int
    #    Current interation
    # tot: tot
    #    Total interations
    
    if {$cur % ($tot/10)} { return }
    # set to total width of progress bar
    set total 100
    set percent [expr {100.*$cur/$tot}]
    set val (\ [format "%6.2f%%" $percent]\ )
    set str "[expr {round($percent*$total/100)}]% "
    puts -nonewline $str
}
