#!/usr/bin/tclsh
# Comprehensive MD Analysis Script for protein complexes
# Designed to work with any segment/chain naming convention

# ===== Procedure Definitions =====
# Center of mass calculation
proc center_of_mass {selection} {
    if {[$selection num] <= 0} {
        error "center_of_mass: needs a selection with atoms"
    }
    set com [veczero]
    set mass 0
    foreach coord [$selection get {x y z}] m [$selection get mass] {
        set mass [expr {$mass + $m}]
        set com [vecadd $com [vecscale $m $coord]]
    }
    if {$mass == 0} {
        error "center_of_mass: total mass is zero"
    }
    return [vecscale [expr {1.0/$mass}] $com]
}

# Radius of gyration calculation
proc gyr_radius {sel} {
    if {[$sel num] <= 0} {
        error "gyr_radius: must have at least one atom in selection"
    }
    set com [center_of_mass $sel]
    set sum 0
    foreach coord [$sel get {x y z}] {
        set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
    }
    return [expr {sqrt($sum / ([$sel num] + 0.0))}]
}

# Periodic boundary modulo function
proc pmodulo {n m} {
    return [expr {$n-$m * floor(1.0 * $n/$m)}]
}

# Clean up selections to prevent memory leaks
proc cleanup_selections {selections} {
    foreach sel $selections {
        catch {$sel delete}
    }
}

# ===== Main Script =====
proc main {} {
    global errorInfo
    
    # Load required packages
    catch {package require pbctools}
    catch {package require hbonds}

    # Set default file names
    set name "step3_input"
    set namepdb "step5_production_noPBC"
    
    # Ask for input files if needed
    puts "\nDefault input files: $name.psf and $namepdb.xtc"
    puts "Use default files? (y/n)"
    flush stdout
    set use_default [gets stdin]
    
    if {![string match -nocase "y*" $use_default]} {
        puts "Enter PSF file name (without extension):"
        flush stdout
        set name [gets stdin]
        
        puts "Enter trajectory file name (without extension):"
        flush stdout
        set namepdb [gets stdin]
    }
    
    # Load molecule
    puts "\nLoading molecule: $name.psf"
    if {[catch {mol new $name.psf waitfor all} err]} {
        puts "Error loading PSF: $err"
        return
    }
    
    # Frame limit option to prevent memory issues
    puts "\nHow many frames would you like to analyze? (0 for all frames)"
    flush stdout
    set max_frames_input [gets stdin]
    set max_frames 0
    if {$max_frames_input != ""} {
        set max_frames [expr int($max_frames_input)]
    }
    
    # Load trajectory
    puts "\nLoading trajectory: $namepdb.xtc"
    if {$max_frames > 0} {
        puts "Will load only the first $max_frames frames to avoid memory issues"
        if {[catch {mol addfile $namepdb.xtc first 0 last [expr $max_frames-1] waitfor all} err]} {
            puts "Error loading trajectory: $err"
            return
        }
    } else {
        if {[catch {mol addfile $namepdb.xtc waitfor all} err]} {
            puts "Error loading trajectory: $err"
            return
        }
    }
    
    # Get number of frames
    set nf [molinfo top get numframes]
    puts "Loaded $nf frames"
    
    # ===== Display Structure Information =====
    set all [atomselect top "all" frame 0]
    
    puts "\n=== STRUCTURE INFORMATION ==="
    
    # Check chains
    set chains [lsort -unique [$all get chain]]
    puts "Chains: $chains"
    
    # Check segments
    set segments [lsort -unique [$all get segname]]
    puts "Segments: $segments"
    
    # Check protein segments
    set prot [atomselect top "protein" frame 0]
    set prot_segments [lsort -unique [$prot get segname]]
    puts "Protein segments: $prot_segments"
    
    # Information about each segment
    foreach seg $prot_segments {
        set sel [atomselect top "segname $seg" frame 0]
        set atomcount [$sel num]
        set rescount [llength [lsort -unique [$sel get resid]]]
        puts "Segment $seg: $atomcount atoms, $rescount residues"
        $sel delete
    }
    
    $all delete
    $prot delete
    
    # ===== User Input for Selection Method =====
    puts "\nHow would you like to select your protein groups?"
    puts "1) By segment ID (recommended, e.g., PROA PROB)"
    puts "2) By chain ID (e.g., A B)"
    flush stdout
    set select_method [gets stdin]
    
    if {$select_method == "2"} {
        # Selection by Chain ID
        puts "\nEnter chains for first protein group (space separated):"
        flush stdout
        set sel1_parts [gets stdin]
        
        puts "Enter name for this group:"
        flush stdout
        set group1_name [gets stdin]
        
        puts "Enter chains for second protein group (space separated):"
        flush stdout
        set sel2_parts [gets stdin]
        
        puts "Enter name for this group:"
        flush stdout
        set group2_name [gets stdin]
        
        # Create selection strings
        set group1_selection ""
        foreach part [split $sel1_parts] {
            if {$group1_selection != ""} {
                append group1_selection " or "
            }
            append group1_selection "chain $part"
        }
        
        set group2_selection ""
        foreach part [split $sel2_parts] {
            if {$group2_selection != ""} {
                append group2_selection " or "
            }
            append group2_selection "chain $part"
        }
    } else {
        # Default: Selection by Segment ID
        puts "\nEnter segments for first protein group (space separated):"
        flush stdout
        set sel1_parts [gets stdin]
        
        puts "Enter name for this group:"
        flush stdout
        set group1_name [gets stdin]
        
        puts "Enter segments for second protein group (space separated):"
        flush stdout
        set sel2_parts [gets stdin]
        
        puts "Enter name for this group:"
        flush stdout
        set group2_name [gets stdin]
        
        # Create selection strings
        set group1_selection ""
        foreach part [split $sel1_parts] {
            if {$group1_selection != ""} {
                append group1_selection " or "
            }
            append group1_selection "segname $part"
        }
        
        set group2_selection ""
        foreach part [split $sel2_parts] {
            if {$group2_selection != ""} {
                append group2_selection " or "
            }
            append group2_selection "segname $part"
        }
    }
    
    # Test selections
    puts "\nTesting selections..."
    set test1 [atomselect top "protein and ($group1_selection)" frame 0]
    set test2 [atomselect top "protein and ($group2_selection)" frame 0]
    
    puts "Group 1 ($group1_name): [$test1 num] atoms selected with: protein and ($group1_selection)"
    puts "Group 2 ($group2_name): [$test2 num] atoms selected with: protein and ($group2_selection)"
    
    if {[$test1 num] == 0} {
        puts "Error: No atoms found for group 1. Please check your selection."
        cleanup_selections [list $test1 $test2]
        return
    }
    
    if {[$test2 num] == 0} {
        puts "Error: No atoms found for group 2. Please check your selection."
        cleanup_selections [list $test1 $test2]
        return
    }
    
    cleanup_selections [list $test1 $test2]
    
    # Create output directory with timestamp
    set timestamp [clock format [clock seconds] -format "%Y%m%d_%H%M%S"]
    set outdir "${group1_name}_${group2_name}_${timestamp}"
    catch {file mkdir $outdir}
    puts "Output will be saved to directory: $outdir"
    
    # ===== Select Analyses to Run =====
    puts "\nSelect analyses to run (1=yes, 0=no):"
    puts "1) Alignment (recommended before other analyses)"
    puts "2) RMSD"
    puts "3) RMSF"
    puts "4) Radius of Gyration"
    puts "5) SASA"
    puts "6) Hydrogen Bonds"
    puts "7) COM-COM Distance"
    puts "8) Interface SASA"
    puts "9) Contact Analysis"
    puts "A) Run ALL analyses"
    puts "0) Run NONE (exit)"
    flush stdout
    
    # Default to run all
    set run_alignment 1
    set run_rmsd 1
    set run_rmsf 1
    set run_rog 1
    set run_sasa 1
    set run_hbonds 1
    set run_comcom 1
    set run_interface_sasa 1
    set run_contacts 1
    
    set choice [gets stdin]
    switch -nocase $choice {
        "1" {
            # Just alignment
            set run_rmsd 0
            set run_rmsf 0
            set run_rog 0
            set run_sasa 0
            set run_hbonds 0
            set run_comcom 0
            set run_interface_sasa 0
            set run_contacts 0
        }
        "2" { 
            # Just RMSD
            set run_alignment 1  # Need alignment for RMSD
            set run_rmsf 0
            set run_rog 0
            set run_sasa 0
            set run_hbonds 0
            set run_comcom 0
            set run_interface_sasa 0
            set run_contacts 0
        }
        "3" { 
            # Just RMSF
            set run_alignment 1  # Need alignment for RMSF
            set run_rmsd 0
            set run_rog 0
            set run_sasa 0
            set run_hbonds 0
            set run_comcom 0
            set run_interface_sasa 0
            set run_contacts 0
        }
        "4" { 
            # Just ROG
            set run_alignment 0
            set run_rmsd 0
            set run_rmsf 0
            set run_sasa 0
            set run_hbonds 0
            set run_comcom 0
            set run_interface_sasa 0
            set run_contacts 0
        }
        "5" { 
            # Just SASA
            set run_alignment 0
            set run_rmsd 0
            set run_rmsf 0
            set run_rog 0
            set run_hbonds 0
            set run_comcom 0
            set run_interface_sasa 0
            set run_contacts 0
        }
        "6" { 
            # Just H-bonds
            set run_alignment 0
            set run_rmsd 0
            set run_rmsf 0
            set run_rog 0
            set run_sasa 0
            set run_comcom 0
            set run_interface_sasa 0
            set run_contacts 0
        }
        "7" { 
            # Just COM-COM
            set run_alignment 0
            set run_rmsd 0
            set run_rmsf 0
            set run_rog 0
            set run_sasa 0
            set run_hbonds 0
            set run_interface_sasa 0
            set run_contacts 0
        }
        "8" { 
            # Just Interface SASA
            set run_alignment 0
            set run_rmsd 0
            set run_rmsf 0
            set run_rog 0
            set run_sasa 0
            set run_hbonds 0
            set run_comcom 0
            set run_contacts 0
        }
        "9" { 
            # Just Contacts
            set run_alignment 0
            set run_rmsd 0
            set run_rmsf 0
            set run_rog 0
            set run_sasa 0
            set run_hbonds 0
            set run_comcom 0
            set run_interface_sasa 0
        }
        "a" -
        "A" { 
            # All analyses (default)
        }
        "0" { 
            puts "Exiting without analyses"
            return
        }
        default { 
            puts "Invalid choice, using defaults (all analyses)"
        }
    }
    
    # ===== EXECUTE ANALYSES =====
    
    # ===== Alignment =====
    if {$run_alignment} {
        puts "\nAligning trajectory based on $group1_name backbone..."
        if {[catch {
            set reference [atomselect top "protein and ($group1_selection) and backbone" frame 0]
            set compare [atomselect top "protein and ($group1_selection) and backbone"]
            set all [atomselect top "all"]
            
            for {set frame 0} {$frame < $nf} {incr frame} {
                $compare frame $frame
                $all frame $frame
                
                set trans_mat [measure fit $compare $reference]
                $all move $trans_mat
                
                if {$frame % 100 == 0} {
                    puts "  Aligned frame $frame / $nf"
                }
            }
            
            cleanup_selections [list $reference $compare $all]
            puts "Alignment complete"
        } err]} {
            puts "Error during alignment: $err"
        }
    }
    
    # ===== RMSD =====
    if {$run_rmsd} {
        # Group 1 RMSD
        puts "\nCalculating RMSD of $group1_name..."
        if {[catch {
            set reference [atomselect top "protein and ($group1_selection) and backbone" frame 0]
            set compare [atomselect top "protein and ($group1_selection) and backbone"]
            set outfile [open "$outdir/RMSD_${group1_name}.dat" w]
            
            for {set frame 0} {$frame < $nf} {incr frame} {
                $compare frame $frame
                set rmsd [measure rmsd $compare $reference]
                puts $outfile "$frame $rmsd"
                
                if {$frame % 100 == 0} {
                    puts "  $group1_name RMSD frame $frame: $rmsd"
                }
            }
            
            close $outfile
            cleanup_selections [list $reference $compare]
            puts "$group1_name RMSD calculation complete"
        } err]} {
            puts "Error calculating $group1_name RMSD: $err"
        }
        
        # Group 2 RMSD
        puts "\nCalculating RMSD of $group2_name..."
        if {[catch {
            set reference [atomselect top "protein and ($group2_selection) and backbone" frame 0]
            set compare [atomselect top "protein and ($group2_selection) and backbone"]
            set outfile [open "$outdir/RMSD_${group2_name}.dat" w]
            
            for {set frame 0} {$frame < $nf} {incr frame} {
                $compare frame $frame
                set rmsd [measure rmsd $compare $reference]
                puts $outfile "$frame $rmsd"
                
                if {$frame % 100 == 0} {
                    puts "  $group2_name RMSD frame $frame: $rmsd"
                }
            }
            
            close $outfile
            cleanup_selections [list $reference $compare]
            puts "$group2_name RMSD calculation complete"
        } err]} {
            puts "Error calculating $group2_name RMSD: $err"
        }
        
        # Complex RMSD
        puts "\nCalculating complex RMSD..."
        if {[catch {
            set reference [atomselect top "protein and (($group1_selection) or ($group2_selection)) and backbone" frame 0]
            set compare [atomselect top "protein and (($group1_selection) or ($group2_selection)) and backbone"]
            set outfile [open "$outdir/RMSD_complex.dat" w]
            
            for {set frame 0} {$frame < $nf} {incr frame} {
                $compare frame $frame
                set rmsd [measure rmsd $compare $reference]
                puts $outfile "$frame $rmsd"
                
                if {$frame % 100 == 0} {
                    puts "  Complex RMSD frame $frame: $rmsd"
                }
            }
            
            close $outfile
            cleanup_selections [list $reference $compare]
            puts "Complex RMSD calculation complete"
        } err]} {
            puts "Error calculating complex RMSD: $err"
        }
    }
    
    # ===== RMSF =====
    if {$run_rmsf} {
        puts "\nCalculating RMSF..."
        
        # Group 1 RMSF
        if {[catch {
            set sel [atomselect top "protein and ($group1_selection) and name CA"]
            set outfile [open "$outdir/RMSF_${group1_name}.dat" w]
            
            set rmsf [measure rmsf $sel]
            set resids [$sel get resid]
            set chains [$sel get chain]
            set resnames [$sel get resname]
            
            puts $outfile "# Chain ResID ResName RMSF"
            for {set i 0} {$i < [$sel num]} {incr i} {
                puts $outfile "[lindex $chains $i] [lindex $resids $i] [lindex $resnames $i] [lindex $rmsf $i]"
            }
            
            close $outfile
            $sel delete
            puts "RMSF calculation for $group1_name complete"
        } err]} {
            puts "Error calculating $group1_name RMSF: $err"
        }
        
        # Group 2 RMSF
        if {[catch {
            set sel [atomselect top "protein and ($group2_selection) and name CA"]
            set outfile [open "$outdir/RMSF_${group2_name}.dat" w]
            
            set rmsf [measure rmsf $sel]
            set resids [$sel get resid]
            set chains [$sel get chain]
            set resnames [$sel get resname]
            
            puts $outfile "# Chain ResID ResName RMSF"
            for {set i 0} {$i < [$sel num]} {incr i} {
                puts $outfile "[lindex $chains $i] [lindex $resids $i] [lindex $resnames $i] [lindex $rmsf $i]"
            }
            
            close $outfile
            $sel delete
            puts "RMSF calculation for $group2_name complete"
        } err]} {
            puts "Error calculating $group2_name RMSF: $err"
        }
    }
    
    # ===== Radius of Gyration =====
    if {$run_rog} {
        puts "\nCalculating Radius of Gyration..."
        
        # Group 1 RoG
        if {[catch {
            set outfile [open "$outdir/rog_${group1_name}.dat" w]
            puts $outfile "frame rad_of_gyr"
            set prot [atomselect top "protein and ($group1_selection)"]
            
            for {set i 0} {$i < $nf} {incr i} {
                $prot frame $i
                $prot update
                set rog [gyr_radius $prot]
                puts $outfile "$i $rog"
                
                if {$i % 100 == 0} {
                    puts "  RoG $group1_name frame $i: $rog"
                }
            }
            
            close $outfile
            $prot delete
            puts "Radius of Gyration calculation for $group1_name complete"
        } err]} {
            puts "Error calculating $group1_name RoG: $err"
        }
        
        # Group 2 RoG
        if {[catch {
            set outfile [open "$outdir/rog_${group2_name}.dat" w]
            puts $outfile "frame rad_of_gyr"
            set prot [atomselect top "protein and ($group2_selection)"]
            
            for {set i 0} {$i < $nf} {incr i} {
                $prot frame $i
                $prot update
                set rog [gyr_radius $prot]
                puts $outfile "$i $rog"
                
                if {$i % 100 == 0} {
                    puts "  RoG $group2_name frame $i: $rog"
                }
            }
            
            close $outfile
            $prot delete
            puts "Radius of Gyration calculation for $group2_name complete"
        } err]} {
            puts "Error calculating $group2_name RoG: $err"
        }
    }
    
    # ===== SASA =====
    if {$run_sasa} {
        puts "\nCalculating Solvent Accessible Surface Area..."
        
        # Group 1 SASA
        if {[catch {
            set sel [atomselect top "protein and ($group1_selection)"]
            set output [open "$outdir/SASA_${group1_name}.dat" w]
            
            for {set i 0} {$i < $nf} {incr i} {
                molinfo top set frame $i
                set sasa [measure sasa 1.4 $sel]
                
                if {$i % 100 == 0} {
                    puts "  SASA $group1_name progress: $i/$nf"
                }
                puts $output "$i $sasa"
            }
            
            close $output
            $sel delete
            puts "SASA calculation for $group1_name complete"
        } err]} {
            puts "Error calculating $group1_name SASA: $err"
        }
        
        # Group 2 SASA
        if {[catch {
            set sel [atomselect top "protein and ($group2_selection)"]
            set output [open "$outdir/SASA_${group2_name}.dat" w]
            
            for {set i 0} {$i < $nf} {incr i} {
                molinfo top set frame $i
                set sasa [measure sasa 1.4 $sel]
                
                if {$i % 100 == 0} {
                    puts "  SASA $group2_name progress: $i/$nf"
                }
                puts $output "$i $sasa"
            }
            
            close $output
            $sel delete
            puts "SASA calculation for $group2_name complete"
        } err]} {
            puts "Error calculating $group2_name SASA: $err"
        }
    }
    
    # ===== Hydrogen Bonds =====
    if {$run_hbonds} {
        puts "\nAnalyzing Hydrogen Bonds between $group1_name and $group2_name..."
        if {[catch {
            package require hbonds
            
            set sel1 [atomselect top "protein and ($group1_selection)"]
            set sel2 [atomselect top "protein and ($group2_selection)"]
            
            set hfile "$outdir/hbonds_${group1_name}_${group2_name}.dat"
            hbonds -sel1 $sel1 -sel2 $sel2 -writefile yes -outfile $hfile -type all
            
            cleanup_selections [list $sel1 $sel2]
            puts "Hydrogen bonds analysis complete"
        } err]} {
            puts "Error analyzing hydrogen bonds: $err"
        }
    }
    
    # ===== COM-COM Distance =====
    if {$run_comcom} {
        puts "\nCalculating center of mass distance..."
        if {[catch {
            set sel1 [atomselect top "protein and ($group1_selection)"]
            set sel2 [atomselect top "protein and ($group2_selection)"]
            
            set out [open "$outdir/comcom_${group1_name}_${group2_name}.dat" w]
            
            for {set i 0} {$i < $nf} {incr i} {
                if {$i % 100 == 0} {
                    puts "  Processing COM-COM Frame $i/$nf"
                }
                
                $sel1 frame $i
                $sel1 update
                $sel2 frame $i
                $sel2 update
                
                set a [measure center $sel1 weight mass]
                set b [measure center $sel2 weight mass]
                
                # Calculate distance with PBC consideration
                if {[catch {set c [lindex [pbc get -first $i -last $i] 0]} err]} {
                    # If PBC calculation fails, use direct distance
                    set dx [expr {[lindex $a 0] - [lindex $b 0]}]
                    set dy [expr {[lindex $a 1] - [lindex $b 1]}]
                    set dz [expr {[lindex $a 2] - [lindex $b 2]}]
                    set d [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
                } else {
                    # With PBC
                    set x1 [lindex $a 0]; set y1 [lindex $a 1]; set z1 [lindex $a 2]
                    set x2 [lindex $b 0]; set y2 [lindex $b 1]; set z2 [lindex $b 2]
                    set boxx [lindex $c 0]; set boxy [lindex $c 1]; set boxz [lindex $c 2]
                    
                    set dx [expr {[pmodulo [expr {$x1-$x2+($boxx/2.0)}] $boxx] - ($boxx/2.0)}]
                    set dy [expr {[pmodulo [expr {$y1-$y2+($boxy/2.0)}] $boxy] - ($boxy/2.0)}]
                    set dz [expr {[pmodulo [expr {$z1-$z2+($boxz/2.0)}] $boxz] - ($boxz/2.0)}]
                    set d [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
                }
                
                puts $out "$i $d"
            }
            
            cleanup_selections [list $sel1 $sel2]
            close $out
            puts "COM-COM distance calculation complete"
        } err]} {
            puts "Error calculating COM-COM distance: $err"
        }
    }
    
    # ===== Interface SASA =====
    if {$run_interface_sasa} {
        puts "\nCalculating Interface SASA..."
        if {[catch {
            set output [open "$outdir/interface_SASA_${group1_name}_${group2_name}.dat" w]
            
            set frame_step 10
            if {$nf < 100} {set frame_step 1}
            
            for {set i 0} {$i < $nf} {incr i $frame_step} {
                molinfo top set frame $i
                
                # Select the groups for this frame
                set sel1 [atomselect top "protein and ($group1_selection)" frame $i]
                set sel2 [atomselect top "protein and ($group2_selection)" frame $i]
                set both [atomselect top "protein and (($group1_selection) or ($group2_selection))" frame $i]
                
                # Calculate SASA for each group separately and together
                set sasa1 [measure sasa 1.4 $sel1]
                set sasa2 [measure sasa 1.4 $sel2]
                set sasaBoth [measure sasa 1.4 $both]
                
                # Calculate interface SASA (buried surface area)
                set interfaceSASA [expr {$sasa1 + $sasa2 - $sasaBoth}]
                
                puts $output "$i $interfaceSASA"
                
                if {$i % 100 == 0 || $i == 0} {
                    puts "  Interface SASA progress: $i/$nf"
                }
                
                cleanup_selections [list $sel1 $sel2 $both]
            }
            
            close $output
            puts "Interface SASA calculation complete"
        } err]} {
            puts "Error calculating interface SASA: $err"
        }
    }
    
    # ===== Contact Analysis =====
    if {$run_contacts} {
        puts "\nPerforming contact analysis..."
        if {[catch {
            # Define contact cutoffs to analyze
            set cutoffs {3.5 5.0 8.0}
            
            foreach cutoff $cutoffs {
                puts "  Analyzing contacts within ${cutoff}Å..."
                set output [open "$outdir/contacts_${group1_name}_${group2_name}_${cutoff}A.dat" w]
                puts $output "# Frame ContactCount"
                
                # Frame step for efficiency
                set frame_step 10
                if {$nf < 100} {set frame_step 1}
                
                for {set i 0} {$i < $nf} {incr i $frame_step} {
                    if {$i % 100 == 0 || $i == 0} {
                        puts "    Frame $i/$nf"
                    }
                    
                    set sel1 [atomselect top "protein and ($group1_selection) and noh" frame $i]
                    set sel2 [atomselect top "protein and ($group2_selection) and noh" frame $i]
                    
                    # Get contact atom pairs
                    set contacts [measure contacts $cutoff $sel1 $sel2]
                    set indices1 [lindex $contacts 0]
                    set indices2 [lindex $contacts 1]
                    set contact_count [llength $indices1]
                    
                    puts $output "$i $contact_count"
                    
                    cleanup_selections [list $sel1 $sel2]
                }
                
                close $output
            }
            
            # Detailed contact analysis for middle frame
            set detail_frame [expr {int($nf / 2)}]
            puts "\nGenerating detailed contact map for frame $detail_frame..."
            
           set cutoff 5.0  # Use 5Å for detailed analysis
            set detail_output [open "$outdir/contact_details_${group1_name}_${group2_name}_${cutoff}A_frame${detail_frame}.dat" w]
            puts $detail_output "# Resid1 Chain1 Resname1 Atom1 Resid2 Chain2 Resname2 Atom2 Distance"
            
            set sel1 [atomselect top "protein and ($group1_selection) and noh" frame $detail_frame]
            set sel2 [atomselect top "protein and ($group2_selection) and noh" frame $detail_frame]
            
            # Get contact atom pairs for detailed analysis
            set contacts [measure contacts $cutoff $sel1 $sel2]
            set indices1 [lindex $contacts 0]
            set indices2 [lindex $contacts 1]
            set contact_count [llength $indices1]
            
            for {set j 0} {$j < $contact_count} {incr j} {
                set idx1 [lindex $indices1 $j]
                set idx2 [lindex $indices2 $j]
                
                set atom1 [atomselect top "index $idx1" frame $detail_frame]
                set atom2 [atomselect top "index $idx2" frame $detail_frame]
                
                set resid1 [$atom1 get resid]
                set chain1 [$atom1 get chain]
                set resname1 [$atom1 get resname]
                set name1 [$atom1 get name]
                
                set resid2 [$atom2 get resid]
                set chain2 [$atom2 get chain]
                set resname2 [$atom2 get resname]
                set name2 [$atom2 get name]
                
                # Calculate actual distance
                set coord1 [lindex [$atom1 get {x y z}] 0]
                set coord2 [lindex [$atom2 get {x y z}] 0]
                set distance [veclength [vecsub $coord1 $coord2]]
                
                puts $detail_output "$resid1 $chain1 $resname1 $name1 $resid2 $chain2 $resname2 $name2 $distance"
                
                $atom1 delete
                $atom2 delete
            }
            
            close $detail_output
            cleanup_selections [list $sel1 $sel2]
            puts "Contact analysis complete"
        } err]} {
            puts "Error performing contact analysis: $err"
        }
    }
    
    # ===== Cleanup and Summary =====
    puts "\n============================================="
    puts "All analyses complete!"
    puts "Output files saved to: $outdir/"
    puts "============================================="
    puts "\nAnalyzed:"
    puts "  Group 1 ($group1_name): $group1_selection"
    puts "  Group 2 ($group2_name): $group2_selection"
    puts "\nFrames analyzed: $nf"
    puts "\nOutput files:"
    if {$run_rmsd} {
        puts "  RMSD: RMSD_${group1_name}.dat, RMSD_${group2_name}.dat, RMSD_complex.dat"
    }
    if {$run_rmsf} {
        puts "  RMSF: RMSF_${group1_name}.dat, RMSF_${group2_name}.dat"
    }
    if {$run_rog} {
        puts "  Radius of Gyration: rog_${group1_name}.dat, rog_${group2_name}.dat"
    }
    if {$run_sasa} {
        puts "  SASA: SASA_${group1_name}.dat, SASA_${group2_name}.dat"
    }
    if {$run_hbonds} {
        puts "  Hydrogen Bonds: hbonds_${group1_name}_${group2_name}.dat"
    }
    if {$run_comcom} {
        puts "  Center of Mass Distance: comcom_${group1_name}_${group2_name}.dat"
    }
    if {$run_interface_sasa} {
        puts "  Interface SASA: interface_SASA_${group1_name}_${group2_name}.dat"
    }
    if {$run_contacts} {
        puts "  Contacts: contacts_${group1_name}_${group2_name}_*A.dat, contact_details_*.dat"
    }
    puts "============================================="
}

# Run the main procedure with error handling
if {[catch {main} err]} {
    puts "ERROR: $err"
    puts "Error Info: $errorInfo"
}