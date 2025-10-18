#!/usr/bin/tclsh
# Comprehensive MD Analysis Script for SAMHD1-SOX11 Complex
# Modified for automated, chain-specific analysis with flexible selection method

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

    # --- Configuration ---
    # Set input file names (modify if needed)
    set psf_file "step3_input.psf"
    set traj_file "step5_production_noPBC.xtc"

    # Set maximum frames to analyze (0 for all frames)
    set max_frames 0

    # Choose selection method: "segname" or "chain"
    # Based on diagnostics, "segname" is required for the current PSF
    set selection_method "segname"

    # Define base names and corresponding identifiers
    set samhd1_base_names {A B C D}
    set sox11_base_name "E"
    set samhd1_identifiers [list "PROA" "PROB" "PROC" "PROD"]
    set sox11_identifier "PROE"

    # Build selection strings dynamically
    set selections [dict create]
    set samhd1_combined_sel_parts {}
    set complex_sel_parts {}

    # SAMHD1 individual chains
    for {set i 0} {$i < [llength $samhd1_base_names]} {incr i} {
        set base_name [lindex $samhd1_base_names $i]
        set identifier [lindex $samhd1_identifiers $i]
        if {$selection_method eq "segname"} {
            set sel_string "protein and segname $identifier"
        } else {
            set sel_string "protein and chain $base_name"
        }
        dict set selections "SAMHD1_$base_name" $sel_string
        lappend samhd1_combined_sel_parts "($sel_string)"
        lappend complex_sel_parts "($sel_string)"
    }

    # SOX11
    if {$selection_method eq "segname"} {
        set sox11_sel_string "protein and segname $sox11_identifier"
    } else {
        set sox11_sel_string "protein and chain $sox11_base_name"
    }
    dict set selections "SOX11_$sox11_base_name" $sox11_sel_string
    lappend complex_sel_parts "($sox11_sel_string)"

    # Combined selections
    dict set selections "SAMHD1_combined" [join $samhd1_combined_sel_parts " or "]
    dict set selections "Complex" [join $complex_sel_parts " or "]

    # Define analyses to run (1=yes, 0=no)
    set run_alignment 1
    set run_rmsd 1
    set run_rmsf 1
    set run_rog 1
    set run_sasa 1
    set run_hbonds 1
    set run_comcom 1
    set run_interface_sasa 1
    set run_contacts 1
    set contact_dist 5.0 ;# Contact distance in Angstroms
    # --- End Configuration ---

    # Load molecule
    puts "\nLoading molecule: $psf_file"
    if {[catch {mol new $psf_file waitfor all} err]} {
        puts "Error loading PSF: $err"
        return
    }

    # Load trajectory
    puts "\nLoading trajectory: $traj_file"
    if {$max_frames > 0} {
        puts "Will load only the first $max_frames frames"
        if {[catch {mol addfile $traj_file first 0 last [expr $max_frames-1] waitfor all} err]} {
            puts "Error loading trajectory: $err"
            return
        }
    } else {
        if {[catch {mol addfile $traj_file waitfor all} err]} {
            puts "Error loading trajectory: $err"
            return
        }
    }

    # Get number of frames
    set nf [molinfo top get numframes]
    puts "Loaded $nf frames"

    # Verify selections
    puts "\nVerifying selections using method: $selection_method..."
    set all_selections_valid 1
    dict for {name sel_string} $selections {
        set test_sel [atomselect top "$sel_string" frame 0]
        puts "  Selection 	'$name':	[$test_sel num] atoms selected with: '$sel_string'"
        if {[$test_sel num] == 0} {
            puts "Error: Selection '$name' resulted in 0 atoms. Check identifiers or selection method."
            set all_selections_valid 0
        }
        $test_sel delete
    }
    if {!$all_selections_valid} {
        puts "Exiting due to invalid selections."
        return
    }

    # Create output directory with timestamp
    set timestamp [clock format [clock seconds] -format "%Y%m%d_%H%M%S"]
    set outdir "SAMHD1_SOX11_comprehensive_${timestamp}"
    catch {file mkdir $outdir}
    puts "Output will be saved to directory: $outdir"

    # Open summary file
    set summary_file [open "$outdir/comprehensive_analysis_summary.txt" w]
    puts $summary_file "Comprehensive MD Analysis Summary"
    puts $summary_file "Timestamp: $timestamp"
    puts $summary_file "PSF: $psf_file"
    puts $summary_file "Trajectory: $traj_file"
    puts $summary_file "Frames Analyzed: $nf"
    puts $summary_file "Selection Method: $selection_method"
    puts $summary_file "\nSelections Defined:"
    dict for {name sel_string} $selections {
        puts $summary_file "  $name: $sel_string"
    }
    puts $summary_file "\nAnalyses Performed:"
    if {$run_alignment} {puts $summary_file "  - Alignment (based on SAMHD1_combined backbone)"}
    if {$run_rmsd} {puts $summary_file "  - RMSD (backbone)"}
    if {$run_rmsf} {puts $summary_file "  - RMSF (CA atoms)"}
    if {$run_rog} {puts $summary_file "  - Radius of Gyration (protein)"}
    if {$run_sasa} {puts $summary_file "  - SASA (protein)"}
    if {$run_hbonds} {puts $summary_file "  - Hydrogen Bonds (protein within 3.5A, 30 deg)"}
    if {$run_comcom} {puts $summary_file "  - COM-COM Distances"}
    if {$run_interface_sasa} {puts $summary_file "  - Interface SASA"}
    if {$run_contacts} {puts $summary_file "  - Contact Analysis (protein within ${contact_dist}A)"}
    puts $summary_file "\n--- Analysis Log ---"

    # ===== EXECUTE ANALYSES =====

    # ===== Alignment =====
    if {$run_alignment} {
        puts "\nAligning trajectory based on SAMHD1_combined backbone..."
        puts $summary_file "\nAlignment:"
        if {[catch {
            set align_sel_string [dict get $selections SAMHD1_combined]
            set reference [atomselect top "$align_sel_string and backbone" frame 0]
            set compare [atomselect top "$align_sel_string and backbone"]
            set all [atomselect top "all"]

            for {set frame 0} {$frame < $nf} {incr frame} {
                $compare frame $frame
                $all frame $frame
                set trans_mat [measure fit $compare $reference]
                $all move $trans_mat
                if {$frame % 100 == 0} {
                    puts -nonewline "\r  Aligned frame $frame / $nf"
                    flush stdout
                }
            }
            puts "\nAlignment complete"
            puts $summary_file "  Trajectory aligned successfully using SAMHD1_combined backbone."
            cleanup_selections [list $reference $compare $all]
        } err]} {
            puts "\nError during alignment: $err"
            puts $summary_file "  Error during alignment: $err"
        }
    }

    # ===== RMSD =====
    if {$run_rmsd} {
        puts "\nCalculating RMSD..."
        puts $summary_file "\nRMSD Analysis:"
        dict for {name sel_string} $selections {
            puts "  Calculating RMSD for $name..."
            if {[catch {
                set reference [atomselect top "$sel_string and backbone" frame 0]
                set compare [atomselect top "$sel_string and backbone"]
                set outfile [open "$outdir/RMSD_${name}.dat" w]
                puts $outfile "# Frame RMSD(A)"

                for {set frame 0} {$frame < $nf} {incr frame} {
                    $compare frame $frame
                    set rmsd [measure rmsd $compare $reference]
                    puts $outfile "$frame $rmsd"
                }

                close $outfile
                cleanup_selections [list $reference $compare]
                puts "    RMSD calculation for $name complete. Output: $outdir/RMSD_${name}.dat"
                puts $summary_file "  - $name RMSD calculated."
            } err]} {
                puts "    Error calculating $name RMSD: $err"
                puts $summary_file "  - Error calculating $name RMSD: $err"
            }
        }
    }

    # ===== RMSF =====
    if {$run_rmsf} {
        puts "\nCalculating RMSF..."
        puts $summary_file "\nRMSF Analysis:"
        dict for {name sel_string} $selections {
             # Skip RMSF for combined selections if desired, as it might be less meaningful
             # if {[string match "*combined*" $name] || [string match "*Complex*" $name]} { continue }

            puts "  Calculating RMSF for $name..."
            if {[catch {
                set sel [atomselect top "$sel_string and name CA"]
                if {[$sel num] == 0} {
                    puts "    Skipping RMSF for $name: No CA atoms found."
                    puts $summary_file "  - Skipped RMSF for $name: No CA atoms found."
                    $sel delete
                    continue
                }
                set outfile [open "$outdir/RMSF_${name}.dat" w]

                set rmsf_values [measure rmsf $sel]
                set resids [$sel get resid]
                set resnames [$sel get resname]
                # Get chain/segment info directly from the selection
                set identifiers [$sel get $selection_method] ; # Use 'chain' or 'segname'

                puts $outfile "# Identifier ResID ResName RMSF(A)"
                for {set i 0} {$i < [$sel num]} {incr i} {
                    puts $outfile "[lindex $identifiers $i] [lindex $resids $i] [lindex $resnames $i] [lindex $rmsf_values $i]"
                }

                close $outfile
                $sel delete
                puts "    RMSF calculation for $name complete. Output: $outdir/RMSF_${name}.dat"
                puts $summary_file "  - $name RMSF calculated."
            } err]} {
                puts "    Error calculating $name RMSF: $err"
                puts $summary_file "  - Error calculating $name RMSF: $err"
            }
        }
    }

    # ===== Radius of Gyration =====
    if {$run_rog} {
        puts "\nCalculating Radius of Gyration..."
        puts $summary_file "\nRadius of Gyration Analysis:"
        dict for {name sel_string} $selections {
            puts "  Calculating RoG for $name..."
            if {[catch {
                set sel [atomselect top "$sel_string"]
                set outfile [open "$outdir/RoG_${name}.dat" w]
                puts $outfile "# Frame RoG(A)"

                for {set frame 0} {$frame < $nf} {incr frame} {
                    $sel frame $frame
                    set rog [gyr_radius $sel]
                    puts $outfile "$frame $rog"
                }

                close $outfile
                $sel delete
                puts "    RoG calculation for $name complete. Output: $outdir/RoG_${name}.dat"
                puts $summary_file "  - $name RoG calculated."
            } err]} {
                puts "    Error calculating $name RoG: $err"
                puts $summary_file "  - Error calculating $name RoG: $err"
            }
        }
    }

    # ===== SASA =====
    if {$run_sasa} {
        puts "\nCalculating SASA..."
        puts $summary_file "\nSASA Analysis:"
        dict for {name sel_string} $selections {
            puts "  Calculating SASA for $name..."
            if {[catch {
                set sel [atomselect top "$sel_string"]
                set outfile [open "$outdir/SASA_${name}.dat" w]
                puts $outfile "# Frame SASA(A^2)"

                for {set frame 0} {$frame < $nf} {incr frame} {
                    $sel frame $frame
                    # Use the selection object itself for restriction
                    set sasa [measure sasa 1.4 $sel -restrict $sel]
                    puts $outfile "$frame $sasa"
                }

                close $outfile
                $sel delete
                puts "    SASA calculation for $name complete. Output: $outdir/SASA_${name}.dat"
                puts $summary_file "  - $name SASA calculated."
            } err]} {
                puts "    Error calculating $name SASA: $err"
                puts $summary_file "  - Error calculating $name SASA: $err"
            }
        }
    }

    # ===== Hydrogen Bonds =====
    if {$run_hbonds} {
        puts "\\nCalculating Hydrogen Bonds..."
        puts $summary_file "\\nHydrogen Bond Analysis:"

        # Define pairs for H-bond analysis (using base names)
        set hbond_pairs [list \\
            [list SAMHD1_combined SAMHD1_combined] \\
            [list SOX11_$sox11_base_name SOX11_$sox11_base_name] \\
            [list SAMHD1_combined SOX11_$sox11_base_name]
        ]

        # Also calculate H-bonds within each SAMHD1 chain
        foreach chain_base_name $samhd1_base_names {
            lappend hbond_pairs [list SAMHD1_${chain_base_name} SAMHD1_${chain_base_name}]
        }

        # Use foreach to directly unpack the pair
        foreach {name1 name2} $hbond_pairs {
            # name1 and name2 are now directly assigned by the loop
            set sel_string1 [dict get $selections $name1]
            set sel_string2 [dict get $selections $name2]
            set pair_name "${name1}_${name2}"
            if {$name1 eq $name2} {
                set pair_name "within_${name1}"
            }
            # Define the output filename clearly
            set output_filename "$outdir/HBonds_${pair_name}.dat"

            puts "  Calculating H-bonds for $pair_name..."
            if {[catch {
                # Define selections ONCE before calling hbonds, covering all frames
                # This is more efficient than redefining selections per frame
                set sel1 [atomselect top "($sel_string1)" frame all]
                set sel2 [atomselect top "($sel_string2)" frame all]

                # Use hbonds command to process all frames and write file directly
                # Set -writefile yes and specify the correct output file
                # Set -loghb no to prevent extra log file creation unless needed
                # The plugin will write the detailed format to the outfile
                # Format: Frame DonorIdx AcceptorIdx DonorResid AcceptorResid DonorAtom AcceptorAtom Distance Angle
                hbonds -sel1 $sel1 -sel2 $sel2 -dist 3.5 -ang 30 -writefile yes -plot no -outfile "$output_filename" -loghb no

                # Clean up the selections after the command finishes
                cleanup_selections [list $sel1 $sel2]

                puts "    H-bond calculation for $pair_name complete. Output: $output_filename"
                # Update summary log - note we don't have the total count easily anymore
                puts $summary_file "  - $pair_name H-bonds calculated and saved to $output_filename."
            } err]} {
                puts stderr "    Error calculating H-bonds for pair '$pair_name':"
                puts stderr "    >> $err"
                # Attempt cleanup even if error occurs
                catch {cleanup_selections [list $sel1 $sel2]}
                puts $summary_file "  - Error calculating $pair_name H-bonds: $err"
            }
        }
    }

    # ===== COM-COM Distance =====
    if {$run_comcom} {
        puts "\nCalculating COM-COM Distances..."
        puts $summary_file "\nCOM-COM Distance Analysis:"

        # Define pairs for COM distance analysis (using base names)
        set com_pairs {}
        # SAMHD1 inter-chain distances
        for {set i 0} {$i < [llength $samhd1_base_names]} {incr i} {
            for {set j [expr {$i + 1}]} {$j < [llength $samhd1_base_names]} {incr j} {
                lappend com_pairs [list SAMHD1_[lindex $samhd1_base_names $i] SAMHD1_[lindex $samhd1_base_names $j]]
            }
        }
        # SAMHD1 chains vs SOX11
        foreach chain_base_name $samhd1_base_names {
            lappend com_pairs [list SAMHD1_${chain_base_name} SOX11_$sox11_base_name]
        }
        # Combined SAMHD1 vs SOX11
        lappend com_pairs [list SAMHD1_combined SOX11_$sox11_base_name]

        foreach pair $com_pairs {
            set name1 [lindex $pair 0]
            set name2 [lindex $pair 1]
            set sel_string1 [dict get $selections $name1]
            set sel_string2 [dict get $selections $name2]
            set pair_name "${name1}_${name2}"

            puts "  Calculating COM distance between $name1 and $name2..."
            if {[catch {
                set outfile [open "$outdir/COM_${pair_name}.dat" w]
                puts $outfile "# Frame Distance(A)"

                for {set frame 0} {$frame < $nf} {incr frame} {
                    set sel1 [atomselect top "$sel_string1"]
                    set sel2 [atomselect top "$sel_string2"]
                    $sel1 frame $frame
                    $sel2 frame $frame

                    set com1 [center_of_mass $sel1]
                    set com2 [center_of_mass $sel2]
                    set dist [veclength [vecsub $com1 $com2]]

                    puts $outfile "$frame $dist"
                    cleanup_selections [list $sel1 $sel2]
                }
                close $outfile
                puts "    COM distance calculation for $pair_name complete. Output: $outdir/COM_${pair_name}.dat"
                puts $summary_file "  - $pair_name COM distance calculated."
            } err]} {
                puts "    Error calculating $pair_name COM distance: $err"
                puts $summary_file "  - Error calculating $pair_name COM distance: $err"
            }
        }
    }

    # ===== Interface SASA =====
    if {$run_interface_sasa} {
        puts "\nCalculating Interface SASA..."
        puts $summary_file "\nInterface SASA Analysis:"

        # Define pairs for interface SASA analysis (using base names)
        set interface_pairs {}
        # SAMHD1 inter-chain interfaces
        for {set i 0} {$i < [llength $samhd1_base_names]} {incr i} {
            for {set j [expr {$i + 1}]} {$j < [llength $samhd1_base_names]} {incr j} {
                lappend interface_pairs [list SAMHD1_[lindex $samhd1_base_names $i] SAMHD1_[lindex $samhd1_base_names $j]]
            }
        }
        # SAMHD1 chains vs SOX11
        foreach chain_base_name $samhd1_base_names {
            lappend interface_pairs [list SAMHD1_${chain_base_name} SOX11_$sox11_base_name]
        }
        # Combined SAMHD1 vs SOX11
        lappend interface_pairs [list SAMHD1_combined SOX11_$sox11_base_name]

        foreach pair $interface_pairs {
            set name1 [lindex $pair 0]
            set name2 [lindex $pair 1]
            set sel_string1 [dict get $selections $name1]
            set sel_string2 [dict get $selections $name2]
            set pair_name "${name1}_${name2}"

            puts "  Calculating Interface SASA between $name1 and $name2..."
            if {[catch {
                set outfile [open "$outdir/InterfaceSASA_${pair_name}.dat" w]
                puts $outfile "# Frame SASA1 SASA2 ComplexSASA InterfaceSASA(A^2)"

                for {set frame 0} {$frame < $nf} {incr frame} {
                    set sel1 [atomselect top "$sel_string1"]
                    set sel2 [atomselect top "$sel_string2"]
                    set sel_complex [atomselect top "($sel_string1) or ($sel_string2)"]
                    $sel1 frame $frame
                    $sel2 frame $frame
                    $sel_complex frame $frame

                    set sasa1 [measure sasa 1.4 $sel1 -restrict "$sel_string1"]
                    set sasa2 [measure sasa 1.4 $sel2 -restrict "$sel_string2"]
                    set sasa_complex [measure sasa 1.4 $sel_complex -restrict "($sel_string1) or ($sel_string2)"]
                    set interface_sasa [expr {($sasa1 + $sasa2 - $sasa_complex) / 2.0}]

                    puts $outfile "$frame $sasa1 $sasa2 $sasa_complex $interface_sasa"
                    cleanup_selections [list $sel1 $sel2 $sel_complex]
                }
                close $outfile
                puts "    Interface SASA calculation for $pair_name complete. Output: $outdir/InterfaceSASA_${pair_name}.dat"
                puts $summary_file "  - $pair_name Interface SASA calculated."
            } err]} {
                puts "    Error calculating $pair_name Interface SASA: $err"
                puts $summary_file "  - Error calculating $pair_name Interface SASA: $err"
            }
        }
    }

    # ===== Contact Analysis =====
    if {$run_contacts} {
        puts "\nCalculating Contacts..."
        puts $summary_file "\nContact Analysis (Distance < ${contact_dist}A):"

        # Define pairs for contact analysis (using base names)
        set contact_pairs {}
        # SAMHD1 inter-chain contacts
        for {set i 0} {$i < [llength $samhd1_base_names]} {incr i} {
            for {set j [expr {$i + 1}]} {$j < [llength $samhd1_base_names]} {incr j} {
                lappend contact_pairs [list SAMHD1_[lindex $samhd1_base_names $i] SAMHD1_[lindex $samhd1_base_names $j]]
            }
        }
        # SAMHD1 chains vs SOX11
        foreach chain_base_name $samhd1_base_names {
            lappend contact_pairs [list SAMHD1_${chain_base_name} SOX11_$sox11_base_name]
        }
        # Combined SAMHD1 vs SOX11
        lappend contact_pairs [list SAMHD1_combined SOX11_$sox11_base_name]

        foreach pair $contact_pairs {
            set name1 [lindex $pair 0]
            set name2 [lindex $pair 1]
            set sel_string1 [dict get $selections $name1]
            set sel_string2 [dict get $selections $name2]
            set pair_name "${name1}_${name2}"

            puts "  Calculating Contacts between $name1 and $name2..."
            if {[catch {
                set outfile [open "$outdir/Contacts_${pair_name}_${contact_dist}A.dat" w]
                puts $outfile "# Frame NumContacts"
                set contact_map_file [open "$outdir/ContactMap_${pair_name}_${contact_dist}A.dat" w]
                puts $contact_map_file "# Frame Res1_Identifier Res1_Resid Res2_Identifier Res2_Resid"

                for {set frame 0} {$frame < $nf} {incr frame} {
                    set sel1 [atomselect top "$sel_string1"]
                    set sel2 [atomselect top "$sel_string2"]
                    $sel1 frame $frame
                    $sel2 frame $frame

                    set contacts [measure contacts $contact_dist $sel1 $sel2]
                    set num_contacts [llength $contacts]
                    puts $outfile "$frame $num_contacts"

                    # Write contact map details
                    set resids1 [$sel1 get resid]
                    set resids2 [$sel2 get resid]
                    set identifiers1 [$sel1 get $selection_method] ; # Use 'chain' or 'segname'
                    set identifiers2 [$sel2 get $selection_method] ; # Use 'chain' or 'segname'

                    foreach contact_pair $contacts {
                        set idx1 [lindex $contact_pair 0]
                        set idx2 [lindex $contact_pair 1]
                        puts $contact_map_file "$frame [lindex $identifiers1 $idx1] [lindex $resids1 $idx1] [lindex $identifiers2 $idx2] [lindex $resids2 $idx2]"
                    }

                    cleanup_selections [list $sel1 $sel2]
                }
                close $outfile
                close $contact_map_file
                puts "    Contact calculation for $pair_name complete. Output: $outdir/Contacts_${pair_name}_${contact_dist}A.dat, $outdir/ContactMap_${pair_name}_${contact_dist}A.dat"
                puts $summary_file "  - $pair_name Contacts calculated."
            } err]} {
                puts "    Error calculating $pair_name Contacts: $err"
                puts $summary_file "  - Error calculating $pair_name Contacts: $err"
            }
        }
    }

    # ===== Final Summary =====
    puts $summary_file "\n--- End of Analysis Log ---"
    close $summary_file

    puts "\n============================================="
    puts "All analyses complete!"
    puts "Output files saved to: $outdir/"
    puts "Summary log: $outdir/comprehensive_analysis_summary.txt"
    puts "============================================="

    # Quit VMD
    quit
}

# ===== Execute Main Script =====
main

