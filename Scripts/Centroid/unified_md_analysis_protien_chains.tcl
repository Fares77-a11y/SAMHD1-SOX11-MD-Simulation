#!/usr/bin/tclsh
# Fixed Complete MD Analysis Script
# Addresses selection issues and provides full functionality

# ===== Global Variables =====
set SCRIPT_VERSION "3.1"
set PROGRESS_INTERVAL 100

# ===== Utility Procedures =====
proc log_message {level message} {
    set timestamp [clock format [clock seconds] -format "%H:%M:%S"]
    puts "\[$timestamp\] \[$level\] $message"
    flush stdout
}

proc progress_bar {current total {width 40}} {
    if {$total <= 0} {return}
    set percent [expr {double($current) / $total}]
    set filled [expr {int($percent * $width)}]
    set bar [string repeat "=" $filled]
    set empty [string repeat " " [expr {$width - $filled}]]
    puts -nonewline "\r\[${bar}${empty}\] [format "%.1f" [expr {$percent * 100}]]% ($current/$total)"
    flush stdout
    if {$current >= $total} {puts ""}
}

# ===== Core Calculation Procedures =====
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

proc gyr_radius {sel} {
    if {[$sel num] <= 0} {
        error "gyr_radius: must have at least one atom in selection"
    }
    set com [center_of_mass $sel]
    set sum 0
    foreach coord [$sel get {x y z}] {
        set sum [expr {$sum + [veclength2 [vecsub $coord $com]]}]
    }
    return [expr {sqrt($sum / ([$sel num] + 0.0))}]
}

proc pmodulo {n m} {
    return [expr {$n - $m * floor(1.0 * $n/$m)}]
}

proc cleanup_selections {selections} {
    foreach sel $selections {
        catch {$sel delete}
    }
}

# ===== Fixed Safe Selection Function =====
proc safe_selection {selection_text {frame 0}} {
    # Try multiple approaches to create atom selection
    set methods [list]
    
    # Method 1: Standard top molecule
    lappend methods "atomselect top \"$selection_text\" frame $frame"
    
    # Method 2: Use explicit molecule ID if available
    if {[info exists ::top_mol]} {
        lappend methods "atomselect $::top_mol \"$selection_text\" frame $frame"
    }
    
    # Method 3: Try molecule 0
    lappend methods "atomselect 0 \"$selection_text\" frame $frame"
    
    foreach method $methods {
        if {[catch {
            set sel [eval $method]
            set count [$sel num]
            return [list $count $sel]
        } err]} {
            continue
        }
    }
    
    return [list -1 "All selection methods failed for: '$selection_text'"]
}

# ===== File Loading =====
proc load_system {} {
    log_message "INFO" "Starting Fixed MD Analysis Script v$::SCRIPT_VERSION"
    
    # Load packages
    foreach pkg {pbctools hbonds} {
        if {[catch {package require $pkg} err]} {
            log_message "WARN" "Could not load $pkg: $err"
        }
    }
    
    # Get file names
    set default_psf "step3_input"
    set default_traj "step5_production_noPBC"
    
    puts "\n=== FILE INPUT ==="
    puts "Default: ${default_psf}.psf & ${default_traj}.xtc"
    puts "Use defaults? (y/n)"
    flush stdout
    set choice [gets stdin]
    
    if {[string match -nocase "y*" $choice]} {
        set psf_file $default_psf
        set traj_file $default_traj
    } else {
        puts "PSF file (no extension):"
        flush stdout
        set psf_file [gets stdin]
        puts "Trajectory file (no extension):"
        flush stdout
        set traj_file [gets stdin]
    }
    
    puts "Frames to analyze? (0=all)"
    flush stdout
    set mf [gets stdin]
    if {[string is integer $mf]} {
        set max_frames $mf
    } else {
        set max_frames 0
    }
    
    # Load PSF
    log_message "INFO" "Loading ${psf_file}.psf"
    if {[catch {set mol_id [mol new ${psf_file}.psf waitfor all]} err]} {
        error "PSF load failed: $err"
    }
    
    # Set global molecule ID
    set ::top_mol $mol_id
    
    # Load trajectory
    log_message "INFO" "Loading ${traj_file}.xtc"
    if {$max_frames > 0} {
        if {[catch {mol addfile ${traj_file}.xtc molid $mol_id first 0 last [expr $max_frames-1] waitfor all} err]} {
            error "Trajectory load failed: $err"
        }
    } else {
        if {[catch {mol addfile ${traj_file}.xtc molid $mol_id waitfor all} err]} {
            error "Trajectory load failed: $err"
        }
    }
    
    set nf [molinfo $mol_id get numframes]
    log_message "INFO" "Loaded $nf frames"
    return $nf
}

# ===== Structure Analysis =====
proc analyze_structure {} {
    log_message "INFO" "Analyzing system structure..."
    
    # Test basic selection
    set res_all [safe_selection "all"]
    if {[lindex $res_all 0] <= 0} {
        error "No atoms found in system: [lindex $res_all 1]"
    }
    set sel_all [lindex $res_all 1]
    puts "Total atoms: [lindex $res_all 0]"
    
    # Test protein selection
    set res_prot [safe_selection "protein"]
    if {[lindex $res_prot 0] <= 0} {
        log_message "WARN" "No 'protein' atoms found, trying alternatives..."
        
        # Try alternative protein selections
        set alt_methods {
            "name CA" "Alpha carbons"
            "backbone" "Backbone atoms"
            "resname ALA VAL LEU ILE PHE TRP TYR MET GLY PRO SER THR CYS ASN GLN ASP GLU HIS LYS ARG" "Standard amino acids"
        }
        
        set found_protein 0
        foreach {method desc} $alt_methods {
            set alt_res [safe_selection $method]
            if {[lindex $alt_res 0] > 0} {
                puts "$desc: [lindex $alt_res 0] atoms"
                if {!$found_protein} {
                    set sel_prot [lindex $alt_res 1]
                    set found_protein 1
                } else {
                    cleanup_selections [list [lindex $alt_res 1]]
                }
            }
        }
        
        if {!$found_protein} {
            cleanup_selections [list $sel_all]
            error "No protein atoms could be identified"
        }
    } else {
        set sel_prot [lindex $res_prot 1]
    }
    
    puts "Protein atoms: [$sel_prot num]"
    
    # Get chains and segments
    set chains [lsort -unique [$sel_prot get chain]]
    set segs [lsort -unique [$sel_prot get segname]]
    
    puts "Chains: $chains"
    puts "Segments: $segs"
    
    # Show detailed information
    puts "\n--- CHAIN DETAILS ---"
    foreach chain $chains {
        if {$chain == ""} {set chain "EMPTY"}
        set chain_res [safe_selection "protein and chain \"$chain\""]
        if {[lindex $chain_res 0] > 0} {
            set chain_sel [lindex $chain_res 1]
            set residues [llength [lsort -unique [$chain_sel get resid]]]
            puts "Chain '$chain': [lindex $chain_res 0] atoms, $residues residues"
            $chain_sel delete
        }
    }
    
    puts "\n--- SEGMENT DETAILS ---"
    foreach seg $segs {
        if {$seg == ""} {set seg "EMPTY"}
        set seg_res [safe_selection "protein and segname \"$seg\""]
        if {[lindex $seg_res 0] > 0} {
            set seg_sel [lindex $seg_res 1]
            set residues [llength [lsort -unique [$seg_sel get resid]]]
            puts "Segment '$seg': [lindex $seg_res 0] atoms, $residues residues"
            $seg_sel delete
        }
    }
    
    cleanup_selections [list $sel_all $sel_prot]
    return [list $chains $segs]
}

# ===== User Selection Input =====
proc get_user_selections {structure_info} {
    set chains [lindex $structure_info 0]
    set segs [lindex $structure_info 1]
    
    puts "\n=== PROTEIN GROUP SELECTION ==="
    
    # Determine selection method
    set has_chains [expr {[llength $chains] > 1 && [lindex $chains 0] != ""}]
    set has_segs [expr {[llength $segs] > 1 && [lindex $segs 0] != ""}]
    
    if {$has_chains && $has_segs} {
        puts "Available options:"
        puts "1) Use segment names: $segs"
        puts "2) Use chain IDs: $chains"
        puts "Choose (1 or 2):"
        flush stdout
        set choice [gets stdin]
    } elseif {$has_chains} {
        puts "Using chain IDs: $chains"
        set choice 2
    } elseif {$has_segs} {
        puts "Using segment names: $segs"
        set choice 1
    } else {
        error "Cannot determine valid group identifiers"
    }
    
    if {$choice == 1} {
        # Segment-based selection
        puts "\nEnter segments for Group 1 (space separated):"
        flush stdout
        set parts1 [gets stdin]
        puts "Name for Group 1:"
        flush stdout
        set name1 [gets stdin]
        puts "Enter segments for Group 2 (space separated):"
        flush stdout
        set parts2 [gets stdin]
        puts "Name for Group 2:"
        flush stdout
        set name2 [gets stdin]
        
        # Build selections
        set sel1_parts {}
        foreach p [split $parts1] {
            lappend sel1_parts "segname $p"
        }
        set sel1 [join $sel1_parts " or "]
        
        set sel2_parts {}
        foreach p [split $parts2] {
            lappend sel2_parts "segname $p"
        }
        set sel2 [join $sel2_parts " or "]
        
    } else {
        # Chain-based selection
        puts "\nEnter chains for Group 1 (space separated):"
        flush stdout
        set parts1 [gets stdin]
        puts "Name for Group 1:"
        flush stdout
        set name1 [gets stdin]
        puts "Enter chains for Group 2 (space separated):"
        flush stdout
        set parts2 [gets stdin]
        puts "Name for Group 2:"
        flush stdout
        set name2 [gets stdin]
        
        # Build selections
        set sel1_parts {}
        foreach p [split $parts1] {
            lappend sel1_parts "chain $p"
        }
        set sel1 [join $sel1_parts " or "]
        
        set sel2_parts {}
        foreach p [split $parts2] {
            lappend sel2_parts "chain $p"
        }
        set sel2 [join $sel2_parts " or "]
    }
    
    # Validate selections
    puts "\n=== VALIDATING SELECTIONS ==="
    set r1 [safe_selection "protein and ($sel1)"]
    set r2 [safe_selection "protein and ($sel2)"]
    
    puts "Group 1 ($name1): [lindex $r1 0] atoms - protein and ($sel1)"
    puts "Group 2 ($name2): [lindex $r2 0] atoms - protein and ($sel2)"
    
    if {[lindex $r1 0] <= 0} {
        error "Group 1 selection failed: [lindex $r1 1]"
    }
    if {[lindex $r2 0] <= 0} {
        error "Group 2 selection failed: [lindex $r2 1]"
    }
    
    # Clean up test selections
    if {[lindex $r1 0] > 0} {cleanup_selections [list [lindex $r1 1]]}
    if {[lindex $r2 0] > 0} {cleanup_selections [list [lindex $r2 1]]}
    
    puts "✓ Selections validated successfully!"
    
    return [list $sel1 $name1 $sel2 $name2]
}

# ===== Analysis Selection =====
proc select_analyses {} {
    puts "\n=== ANALYSIS SELECTION ==="
    puts "1) Alignment"
    puts "2) RMSD"
    puts "3) RMSF"
    puts "4) Radius of Gyration"
    puts "5) SASA"
    puts "6) Hydrogen Bonds"
    puts "7) COM-COM Distance"
    puts "8) Interface SASA"
    puts "9) Contact Analysis"
    puts ""
    puts "A) All analyses"
    puts "Q) Quick (Alignment, RMSD, COM-COM)"
    puts "C) Custom selection"
    flush stdout
    set choice [gets stdin]
    
    array set a {alignment 0 rmsd 0 rmsf 0 rog 0 sasa 0 hbonds 0 comcom 0 interface_sasa 0 contacts 0}
    
    switch -nocase $choice {
        "A" {
            foreach k [array names a] {set a($k) 1}
            puts "Selected: All analyses"
        }
        "Q" {
            set a(alignment) 1; set a(rmsd) 1; set a(comcom) 1
            puts "Selected: Quick analyses"
        }
        "C" {
            puts "Enter analysis numbers (space separated, e.g., '1 2 7'):"
            flush stdout
            set nums [gets stdin]
            foreach n [split $nums] {
                switch $n {
                    1 {set a(alignment) 1}
                    2 {set a(rmsd) 1; set a(alignment) 1}
                    3 {set a(rmsf) 1; set a(alignment) 1}
                    4 {set a(rog) 1}
                    5 {set a(sasa) 1}
                    6 {set a(hbonds) 1}
                    7 {set a(comcom) 1}
                    8 {set a(interface_sasa) 1}
                    9 {set a(contacts) 1}
                }
            }
            puts "Selected: Custom analyses"
        }
        default {
            foreach k [array names a] {set a($k) 1}
            puts "Invalid choice, defaulting to all analyses"
        }
    }
    
    return [array get a]
}

# ===== Core Analysis Functions =====
proc perform_alignment {nf reference_selection reference_name} {
    log_message "INFO" "Aligning trajectory based on $reference_name backbone"
    
    if {[catch {
        set ref_result [safe_selection "protein and ($reference_selection) and backbone"]
        if {[lindex $ref_result 0] <= 0} {
            error "No backbone atoms found for alignment"
        }
        
        set reference [lindex $ref_result 1]
        set compare_result [safe_selection "protein and ($reference_selection) and backbone"]
        set compare [lindex $compare_result 1]
        set all_result [safe_selection "all"]
        set all [lindex $all_result 1]
        
        for {set frame 0} {$frame < $nf} {incr frame} {
            $compare frame $frame
            $all frame $frame
            
            set trans_mat [measure fit $compare $reference]
            $all move $trans_mat
            
            if {$frame % $::PROGRESS_INTERVAL == 0} {
                progress_bar $frame $nf
            }
        }
        progress_bar $nf $nf
        
        cleanup_selections [list $reference $compare $all]
        log_message "INFO" "Alignment complete"
    } err]} {
        log_message "ERROR" "Alignment failed: $err"
    }
}

proc calculate_rmsd {nf outdir selection_command group_name} {
    log_message "INFO" "Calculating RMSD for $group_name"
    
    if {[catch {
        set ref_result [safe_selection "protein and ($selection_command) and backbone"]
        if {[lindex $ref_result 0] <= 0} {
            error "No backbone atoms found for RMSD"
        }
        
        set reference [lindex $ref_result 1]
        set compare_result [safe_selection "protein and ($selection_command) and backbone"]
        set compare [lindex $compare_result 1]
        set outfile [open "$outdir/RMSD_${group_name}.dat" w]
        
        puts $outfile "# RMSD data for $group_name"
        puts $outfile "# Generated by Fixed MD Analysis Script v$::SCRIPT_VERSION"
        puts $outfile "# Frame RMSD(Angstrom)"
        
        for {set frame 0} {$frame < $nf} {incr frame} {
            $compare frame $frame
            set rmsd [measure rmsd $compare $reference]
            puts $outfile "$frame $rmsd"
            
            if {$frame % $::PROGRESS_INTERVAL == 0} {
                progress_bar $frame $nf
            }
        }
        progress_bar $nf $nf
        
        close $outfile
        cleanup_selections [list $reference $compare]
        log_message "INFO" "RMSD calculation complete for $group_name"
    } err]} {
        log_message "ERROR" "RMSD calculation failed for $group_name: $err"
    }
}

proc calculate_complex_rmsd {nf outdir sel1 sel2 name1 name2} {
    log_message "INFO" "Calculating complex RMSD for ${name1}_${name2}"
    
    if {[catch {
        set ref_result [safe_selection "protein and (($sel1) or ($sel2)) and backbone"]
        if {[lindex $ref_result 0] <= 0} {
            error "No backbone atoms found for complex RMSD"
        }
        
        set reference [lindex $ref_result 1]
        set compare_result [safe_selection "protein and (($sel1) or ($sel2)) and backbone"]
        set compare [lindex $compare_result 1]
        set outfile [open "$outdir/RMSD_complex_${name1}_${name2}.dat" w]
        
        puts $outfile "# Complex RMSD data for ${name1} + ${name2}"
        puts $outfile "# Generated by Fixed MD Analysis Script v$::SCRIPT_VERSION"
        puts $outfile "# Frame RMSD(Angstrom)"
        
        for {set frame 0} {$frame < $nf} {incr frame} {
            $compare frame $frame
            set rmsd [measure rmsd $compare $reference]
            puts $outfile "$frame $rmsd"
            
            if {$frame % $::PROGRESS_INTERVAL == 0} {
                progress_bar $frame $nf
            }
        }
        progress_bar $nf $nf
        
        close $outfile
        cleanup_selections [list $reference $compare]
        log_message "INFO" "Complex RMSD calculation complete"
    } err]} {
        log_message "ERROR" "Complex RMSD calculation failed: $err"
    }
}

proc calculate_rmsf {outdir selection_command group_name} {
    log_message "INFO" "Calculating RMSF for $group_name"
    
    if {[catch {
        set ca_result [safe_selection "protein and ($selection_command) and name CA"]
        if {[lindex $ca_result 0] <= 0} {
            error "No CA atoms found for RMSF"
        }
        
        set sel [lindex $ca_result 1]
        set outfile [open "$outdir/RMSF_${group_name}.dat" w]
        
        puts $outfile "# RMSF data for $group_name (CA atoms)"
        puts $outfile "# Generated by Fixed MD Analysis Script v$::SCRIPT_VERSION"
        puts $outfile "# Chain ResID ResName RMSF(Angstrom)"
        
        set rmsf [measure rmsf $sel]
        set resids [$sel get resid]
        set chains [$sel get chain]
        set resnames [$sel get resname]
        
        for {set i 0} {$i < [$sel num]} {incr i} {
            puts $outfile "[lindex $chains $i] [lindex $resids $i] [lindex $resnames $i] [lindex $rmsf $i]"
        }
        
        close $outfile
        $sel delete
        log_message "INFO" "RMSF calculation complete for $group_name"
    } err]} {
        log_message "ERROR" "RMSF calculation failed for $group_name: $err"
    }
}

proc calculate_rog {nf outdir selection_command group_name} {
    log_message "INFO" "Calculating Radius of Gyration for $group_name"
    
    if {[catch {
        set prot_result [safe_selection "protein and ($selection_command)"]
        if {[lindex $prot_result 0] <= 0} {
            error "No protein atoms found for RoG"
        }
        
        set prot [lindex $prot_result 1]
        set outfile [open "$outdir/rog_${group_name}.dat" w]
        
        puts $outfile "# Radius of Gyration data for $group_name"
        puts $outfile "# Generated by Fixed MD Analysis Script v$::SCRIPT_VERSION"
        puts $outfile "# Frame RadiusOfGyration(Angstrom)"
        
        for {set i 0} {$i < $nf} {incr i} {
            $prot frame $i
            $prot update
            set rog [gyr_radius $prot]
            puts $outfile "$i $rog"
            
            if {$i % $::PROGRESS_INTERVAL == 0} {
                progress_bar $i $nf
            }
        }
        progress_bar $nf $nf
        
        close $outfile
        $prot delete
        log_message "INFO" "Radius of Gyration calculation complete for $group_name"
    } err]} {
        log_message "ERROR" "Radius of Gyration calculation failed for $group_name: $err"
    }
}

proc calculate_sasa {nf outdir selection_command group_name} {
    log_message "INFO" "Calculating SASA for $group_name"
    
    if {[catch {
        set sel_result [safe_selection "protein and ($selection_command)"]
        if {[lindex $sel_result 0] <= 0} {
            error "No protein atoms found for SASA"
        }
        
        set sel [lindex $sel_result 1]
        set output [open "$outdir/SASA_${group_name}.dat" w]
        
        puts $output "# SASA data for $group_name"
        puts $output "# Generated by Fixed MD Analysis Script v$::SCRIPT_VERSION"
        puts $output "# Frame SASA(Angstrom^2)"
        
        for {set i 0} {$i < $nf} {incr i} {
            molinfo top set frame $i
            set sasa [measure sasa 1.4 $sel]
            puts $output "$i $sasa"
            
            if {$i % $::PROGRESS_INTERVAL == 0} {
                progress_bar $i $nf
            }
        }
        progress_bar $nf $nf
        
        close $output
        $sel delete
        log_message "INFO" "SASA calculation complete for $group_name"
    } err]} {
        log_message "ERROR" "SASA calculation failed for $group_name: $err"
    }
}

proc calculate_comcom_distance {nf outdir sel1 sel2 name1 name2} {
    log_message "INFO" "Calculating COM-COM distance between $name1 and $name2"
    
    if {[catch {
        set selection1_result [safe_selection "protein and ($sel1)"]
        set selection2_result [safe_selection "protein and ($sel2)"]
        
        if {[lindex $selection1_result 0] <= 0 || [lindex $selection2_result 0] <= 0} {
            error "Invalid selections for COM-COM distance"
        }
        
        set selection1 [lindex $selection1_result 1]
        set selection2 [lindex $selection2_result 1]
        set out [open "$outdir/comcom_${name1}_${name2}.dat" w]
        
        puts $out "# Center of Mass distance between $name1 and $name2"
        puts $out "# Generated by Fixed MD Analysis Script v$::SCRIPT_VERSION"
        puts $out "# Frame Distance(Angstrom)"
        
        for {set i 0} {$i < $nf} {incr i} {
            $selection1 frame $i
            $selection1 update
            $selection2 frame $i
            $selection2 update
            
            set com1 [measure center $selection1 weight mass]
            set com2 [measure center $selection2 weight mass]
            
            set dx [expr {[lindex $com1 0] - [lindex $com2 0]}]
            set dy [expr {[lindex $com1 1] - [lindex $com2 1]}]
            set dz [expr {[lindex $com1 2] - [lindex $com2 2]}]
            set distance [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
            
            puts $out "$i $distance"
            
            if {$i % $::PROGRESS_INTERVAL == 0} {
                progress_bar $i $nf
            }
        }
        progress_bar $nf $nf
        
        cleanup_selections [list $selection1 $selection2]
        close $out
        log_message "INFO" "COM-COM distance calculation complete"
    } err]} {
        log_message "ERROR" "COM-COM distance calculation failed: $err"
    }
}

proc analyze_hbonds {outdir sel1 sel2 name1 name2} {
    log_message "INFO" "Analyzing hydrogen bonds between $name1 and $name2"
    
    if {[catch {
        package require hbonds
        
        set selection1_result [safe_selection "protein and ($sel1)"]
        set selection2_result [safe_selection "protein and ($sel2)"]
        
        if {[lindex $selection1_result 0] <= 0 || [lindex $selection2_result 0] <= 0} {
            error "Invalid selections for hydrogen bond analysis"
        }
        
        set selection1 [lindex $selection1_result 1]
        set selection2 [lindex $selection2_result 1]
        
        set hfile "$outdir/hbonds_${name1}_${name2}.dat"
        hbonds -sel1 $selection1 -sel2 $selection2 -writefile yes -outfile $hfile -type all
        
        cleanup_selections [list $selection1 $selection2]
        log_message "INFO" "Hydrogen bonds analysis complete"
    } err]} {
        log_message "ERROR" "Hydrogen bonds analysis failed: $err"
    }
}

# ===== Main Execution =====
proc execute_analysis {nf sel1 name1 sel2 name2 analyses outdir} {
    array set run_analyses $analyses
    
    log_message "INFO" "Starting analysis execution..."
    
    # Alignment
    if {$run_analyses(alignment)} {
        perform_alignment $nf $sel1 $name1
    }
    
    # RMSD Analysis
    if {$run_analyses(rmsd)} {
        calculate_rmsd $nf $outdir $sel1 $name1
        calculate_rmsd $nf $outdir $sel2 $name2
        calculate_complex_rmsd $nf $outdir $sel1 $sel2 $name1 $name2
    }
    
    # RMSF Analysis
    if {$run_analyses(rmsf)} {
        calculate_rmsf $outdir $sel1 $name1
        calculate_rmsf $outdir $sel2 $name2
    }
    
    # Radius of Gyration
    if {$run_analyses(rog)} {
        calculate_rog $nf $outdir $sel1 $name1
        calculate_rog $nf $outdir $sel2 $name2
    }
    
    # SASA Analysis
    if {$run_analyses(sasa)} {
        calculate_sasa $nf $outdir $sel1 $name1
        calculate_sasa $nf $outdir $sel2 $name2
    }
    
    # Hydrogen Bonds
    if {$run_analyses(hbonds)} {
        analyze_hbonds $outdir $sel1 $sel2 $name1 $name2
    }
    
    # COM-COM Distance
    if {$run_analyses(comcom)} {
        calculate_comcom_distance $nf $outdir $sel1 $sel2 $name1 $name2
    }
}

proc generate_summary {outdir name1 name2 analyses nf} {
    array set run_analyses $analyses
    
    set summary_file [open "$outdir/analysis_summary.txt" w]
    
    puts $summary_file "MD Analysis Summary"
    puts $summary_file "=================="
    puts $summary_file "Generated: [clock format [clock seconds]]"
    puts $summary_file "Script Version: $::SCRIPT_VERSION"
    puts $summary_file ""
    puts $summary_file "System Information:"
    puts $summary_file "  Group 1: $name1"
    puts $summary_file "  Group 2: $name2"
    puts $summary_file "  Frames analyzed: $nf"
    puts $summary_file ""
    puts $summary_file "Analyses Performed:"
    
    if {$run_analyses(alignment)} {puts $summary_file "  ✓ Trajectory Alignment"}
    if {$run_analyses(rmsd)} {puts $summary_file "  ✓ RMSD Analysis"}
    if {$run_analyses(rmsf)} {puts $summary_file "  ✓ RMSF Analysis"}
    if {$run_analyses(rog)} {puts $summary_file "  ✓ Radius of Gyration"}
    if {$run_analyses(sasa)} {puts $summary_file "  ✓ SASA Analysis"}
    if {$run_analyses(hbonds)} {puts $summary_file "  ✓ Hydrogen Bond Analysis"}
    if {$run_analyses(comcom)} {puts $summary_file "  ✓ COM-COM Distance"}
    if {$run_analyses(interface_sasa)} {puts $summary_file "  ✓ Interface SASA"}
    if {$run_analyses(contacts)} {puts $summary_file "  ✓ Contact Analysis"}
    
    puts $summary_file ""
    puts $summary_file "Output Files:"
    puts $summary_file "-------------"
    
    set files [glob -nocomplain "$outdir/*.dat"]
    foreach file $files {
        set filename [file tail $file]
        set filesize [file size $file]
        puts $summary_file "  $filename ([expr {$filesize / 1024}] KB)"
    }
    
    close $summary_file
    log_message "INFO" "Analysis summary generated: $outdir/analysis_summary.txt"
}

proc main {} {
    global errorInfo
    
    # Load system and trajectory
    set nf [load_system]
    
    # Analyze structure
    set structure_info [analyze_structure]
    
    # Get user selections
    set selection_info [get_user_selections $structure_info]
    set sel1 [lindex $selection_info 0]
    set name1 [lindex $selection_info 1]
    set sel2 [lindex $selection_info 2]
    set name2 [lindex $selection_info 3]
    
    # Select analyses
    set analyses [select_analyses]
    
    # Create output directory
    set timestamp [clock format [clock seconds] -format "%Y%m%d_%H%M%S"]
    set outdir "${name1}_${name2}_${timestamp}"
    if {[catch {file mkdir $outdir} err]} {
        error "Cannot create output directory: $err"
    }
    log_message "INFO" "Output directory created: $outdir"
    
    # Execute analyses
    execute_analysis $nf $sel1 $name1 $sel2 $name2 $analyses $outdir
    
    # Generate summary
    generate_summary $outdir $name1 $name2 $analyses $nf
    
    # Final summary
    puts "\n🎉 ANALYSIS COMPLETE! 🎉"
    puts "==============================="
    puts "Groups analyzed:"
    puts "  • $name1: $sel1"
    puts "  • $name2: $sel2"
    puts ""
    puts "Frames processed: $nf"
    puts "Output directory: $outdir/"
    puts ""
    puts "Check analysis_summary.txt for details"
    puts "==============================="
    
    log_message "INFO" "All analyses completed successfully!"
}

# ===== Error Handling =====
proc cleanup_and_exit {{exit_code 0}} {
    log_message "INFO" "Cleaning up resources"
    
    # Close any open files
    foreach channel [file channels] {
        if {$channel != "stdin" && $channel != "stdout" && $channel != "stderr"} {
            catch {close $channel}
        }
    }
    
    if {$exit_code == 0} {
        log_message "INFO" "Analysis completed successfully"
    } else {
        log_message "ERROR" "Analysis terminated with errors"
    }
    
    # DO NOT EXIT - let user see any errors
}

# ===== Script Execution =====
if {[catch {main} err]} {
    log_message "ERROR" "Fatal error in main execution: $err"
    puts "\n=== DETAILED ERROR INFORMATION ==="
    puts "Error message: $err"
    puts "\nFull error trace:"
    puts $errorInfo
    puts "\n=== TROUBLESHOOTING TIPS ==="
    puts "1. Check that your PSF and XTC files exist and are readable"
    puts "2. Verify VMD can load your files manually:"
    puts "   mol new step3_input.psf"
    puts "   mol addfile step5_production_noPBC.xtc"
    puts "3. Test basic selections:"
    puts "   set sel \[atomselect top \"all\"\]"
    puts "   puts \[\$sel num\]"
    puts "   \$sel delete"
    puts "4. Check if protein atoms exist:"
    puts "   set prot \[atomselect top \"protein\"\]"
    puts "   puts \[\$prot num\]"
    puts "   \$prot delete"
    puts "\n=== CURRENT SYSTEM STATE ==="
    catch {
        puts "Number of molecules loaded: [molinfo num]"
        if {[molinfo num] > 0} {
            puts "Top molecule ID: [molinfo top]"
            puts "Frames in top molecule: [molinfo top get numframes]"
            puts "Atoms in top molecule: [molinfo top get numatoms]"
        }
    }
    puts "\n=== SELECTION TEST ==="
    catch {
        puts "Testing 'all' selection..."
        set test_all [atomselect top "all"]
        puts "All atoms: [$test_all num]"
        $test_all delete
        
        puts "Testing 'protein' selection..."
        set test_prot [atomselect top "protein"]
        puts "Protein atoms: [$test_prot num]"
        $test_prot delete
    }
    cleanup_and_exit 1
} else {
    cleanup_and_exit 0
}

# ===== Additional Functions for Comprehensive Analysis =====

# Add individual chain analysis capability
proc analyze_individual_chains {nf outdir samhd1_chains sox11_chains analyses} {
    array set run_analyses $analyses
    
    log_message "INFO" "Starting individual chain analysis..."
    
    # Analyze each SAMHD1 chain individually
    foreach chain $samhd1_chains {
        set chain_sel "chain $chain"
        set chain_name "SAMHD1_${chain}"
        
        log_message "INFO" "Analyzing SAMHD1 chain $chain"
        
        if {$run_analyses(rmsd)} {
            calculate_rmsd $nf $outdir $chain_sel $chain_name
        }
        
        if {$run_analyses(rmsf)} {
            calculate_rmsf $outdir $chain_sel $chain_name
        }
        
        if {$run_analyses(rog)} {
            calculate_rog $nf $outdir $chain_sel $chain_name
        }
        
        if {$run_analyses(sasa)} {
            calculate_sasa $nf $outdir $chain_sel $chain_name
        }
    }
    
    # Analyze each SOX11 chain individually
    foreach chain $sox11_chains {
        set chain_sel "chain $chain"
        set chain_name "SOX11_${chain}"
        
        log_message "INFO" "Analyzing SOX11 chain $chain"
        
        if {$run_analyses(rmsd)} {
            calculate_rmsd $nf $outdir $chain_sel $chain_name
        }
        
        if {$run_analyses(rmsf)} {
            calculate_rmsf $outdir $chain_sel $chain_name
        }
        
        if {$run_analyses(rog)} {
            calculate_rog $nf $outdir $chain_sel $chain_name
        }
        
        if {$run_analyses(sasa)} {
            calculate_sasa $nf $outdir $chain_sel $chain_name
        }
    }
    
    # Inter-chain analyses
    if {[llength $samhd1_chains] > 1} {
        log_message "INFO" "Analyzing inter-SAMHD1 chain interactions"
        
        for {set i 0} {$i < [llength $samhd1_chains]} {incr i} {
            for {set j [expr $i+1]} {$j < [llength $samhd1_chains]} {incr j} {
                set chain1 [lindex $samhd1_chains $i]
                set chain2 [lindex $samhd1_chains $j]
                
                if {$run_analyses(comcom)} {
                    calculate_comcom_distance $nf $outdir "chain $chain1" "chain $chain2" "SAMHD1_${chain1}" "SAMHD1_${chain2}"
                }
                
                if {$run_analyses(hbonds)} {
                    analyze_hbonds $outdir "chain $chain1" "chain $chain2" "SAMHD1_${chain1}" "SAMHD1_${chain2}"
                }
            }
        }
    }
    
    # Individual SAMHD1 chains vs SOX11
    foreach samhd1_chain $samhd1_chains {
        foreach sox11_chain $sox11_chains {
            if {$run_analyses(comcom)} {
                calculate_comcom_distance $nf $outdir "chain $samhd1_chain" "chain $sox11_chain" "SAMHD1_${samhd1_chain}" "SOX11_${sox11_chain}"
            }
            
            if {$run_analyses(hbonds)} {
                analyze_hbonds $outdir "chain $samhd1_chain" "chain $sox11_chain" "SAMHD1_${samhd1_chain}" "SOX11_${sox11_chain}"
            }
        }
    }
}

# Function to run comprehensive analysis (individual + combined + complex)
proc run_comprehensive_analysis {} {
    puts "\n=== COMPREHENSIVE ANALYSIS MODE ==="
    puts "This will analyze:"
    puts "1. Individual SAMHD1 chains"
    puts "2. Individual SOX11 chains"  
    puts "3. Combined SAMHD1"
    puts "4. Combined SOX11"
    puts "5. SAMHD1-SOX11 complex"
    puts "6. All pairwise interactions"
    puts ""
    puts "Proceed? (y/n)"
    flush stdout
    set proceed [gets stdin]
    
    if {![string match -nocase "y*" $proceed]} {
        puts "Returning to standard analysis..."
        return 0
    }
    
    return 1
}

# Usage instructions
proc print_usage {} {
    puts "\n=== USAGE INSTRUCTIONS ==="
    puts "This script provides comprehensive MD analysis for protein complexes."
    puts ""
    puts "To run the analysis:"
    puts "1. Load this script in VMD: source fixed_complete_md_analysis.tcl"
    puts "2. The script will automatically start and guide you through:"
    puts "   - File selection (PSF and trajectory)"
    puts "   - System structure analysis"
    puts "   - Protein group identification"
    puts "   - Analysis selection"
    puts "   - Execution and results"
    puts ""
    puts "Supported analyses:"
    puts "- Trajectory alignment"
    puts "- RMSD (Root Mean Square Deviation)"
    puts "- RMSF (Root Mean Square Fluctuation)"
    puts "- Radius of Gyration"
    puts "- SASA (Solvent Accessible Surface Area)"
    puts "- Hydrogen bond analysis"
    puts "- Center of mass distances"
    puts "- Interface SASA"
    puts "- Contact analysis"
    puts ""
    puts "For comprehensive analysis including individual chains:"
    puts "- The script will automatically detect and offer individual chain analysis"
    puts "- This includes SAMHD1 chains individually, SOX11 individually, and all interactions"
    puts ""
    puts "Output:"
    puts "- All results saved to timestamped directory"
    puts "- Summary file with complete analysis report"
    puts "- Individual data files for each analysis"
    puts "==============================================="
}

# Print usage when script is loaded
print_usage