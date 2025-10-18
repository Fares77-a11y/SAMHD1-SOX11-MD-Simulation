{\rtf1\ansi\ansicpg1252\cocoartf2821
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/usr/bin/tclsh\
# Comprehensive MD Analysis Script for protein complexes\
# Designed for SAMHD1-SOX11 analysis with flexible chain/segment identification\
\
# ===== Procedure Definitions =====\
# Center of mass calculation\
proc center_of_mass \{selection\} \{\
    if \{[$selection num] <= 0\} \{\
        error "center_of_mass: needs a selection with atoms"\
    \}\
    set com [veczero]\
    set mass 0\
    foreach coord [$selection get \{x y z\}] m [$selection get mass] \{\
        set mass [expr \{$mass + $m\}]\
        set com [vecadd $com [vecscale $m $coord]]\
    \}\
    if \{$mass == 0\} \{\
        error "center_of_mass: total mass is zero"\
    \}\
    return [vecscale [expr \{1.0/$mass\}] $com]\
\}\
\
# Radius of gyration calculation\
proc gyr_radius \{sel\} \{\
    if \{[$sel num] <= 0\} \{\
        error "gyr_radius: must have at least one atom in selection"\
    \}\
    set com [center_of_mass $sel]\
    set sum 0\
    foreach coord [$sel get \{x y z\}] \{\
        set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]\
    \}\
    return [expr \{sqrt($sum / ([$sel num] + 0.0))\}]\
\}\
\
# Periodic boundary modulo function\
proc pmodulo \{n m\} \{\
    return [expr \{$n-$m * floor(1.0 * $n/$m)\}]\
\}\
\
# Clean up selections to prevent memory leaks\
proc cleanup_selections \{selections\} \{\
    foreach sel $selections \{\
        catch \{$sel delete\}\
    \}\
\}\
\
# ===== Main Script =====\
proc main \{\} \{\
    global errorInfo\
    \
    # Load required packages\
    catch \{package require pbctools\}\
    catch \{package require hbonds\}\
\
    # Set default file names\
    set name "step3_input"\
    set namepdb "step5_production_noPBC"\
    \
    # Ask for input files if needed\
    puts "\\nDefault input files: $name.psf and $namepdb.xtc"\
    puts "Use default files? (y/n)"\
    flush stdout\
    set use_default [gets stdin]\
    \
    if \{![string match -nocase "y*" $use_default]\} \{\
        puts "Enter PSF file name (without extension):"\
        flush stdout\
        set name [gets stdin]\
        \
        puts "Enter trajectory file name (without extension):"\
        flush stdout\
        set namepdb [gets stdin]\
    \}\
    \
    # Load molecule\
    puts "\\nLoading molecule: $name.psf"\
    if \{[catch \{mol new $name.psf waitfor all\} err]\} \{\
        puts "Error loading PSF: $err"\
        return\
    \}\
    \
    # Frame limit option to prevent memory issues\
    puts "\\nHow many frames would you like to analyze? (0 for all frames)"\
    flush stdout\
    set max_frames_input [gets stdin]\
    set max_frames 0\
    if \{$max_frames_input != ""\} \{\
        set max_frames [expr int($max_frames_input)]\
    \}\
    \
    # Load trajectory\
    puts "\\nLoading trajectory: $namepdb.xtc"\
    if \{$max_frames > 0\} \{\
        puts "Will load only the first $max_frames frames to avoid memory issues"\
        if \{[catch \{mol addfile $namepdb.xtc first 0 last [expr $max_frames-1] waitfor all\} err]\} \{\
            puts "Error loading trajectory: $err"\
            return\
        \}\
    \} else \{\
        if \{[catch \{mol addfile $namepdb.xtc waitfor all\} err]\} \{\
            puts "Error loading trajectory: $err"\
            return\
        \}\
    \}\
    \
    # Get number of frames\
    set nf [molinfo top get numframes]\
    puts "Loaded $nf frames"\
    \
    # ===== Display Structure Information =====\
    set all [atomselect top "all" frame 0]\
    \
    puts "\\n=== STRUCTURE INFORMATION ==="\
    \
    # Check chains\
    set chains [lsort -unique [$all get chain]]\
    puts "Chains: $chains"\
    \
    # Check segments\
    set segments [lsort -unique [$all get segname]]\
    puts "Segments: $segments"\
    \
    # Check protein segments\
    set prot [atomselect top "protein" frame 0]\
    set prot_segments [lsort -unique [$prot get segname]]\
    puts "Protein segments: $prot_segments"\
    \
    # Information about each segment\
    foreach seg $prot_segments \{\
        set sel [atomselect top "segname $seg" frame 0]\
        set atomcount [$sel num]\
        set rescount [llength [lsort -unique [$sel get resid]]]\
        puts "Segment $seg: $atomcount atoms, $rescount residues"\
        $sel delete\
    \}\
    \
    $all delete\
    $prot delete\
    \
    # ===== Determine Selection Method =====\
    puts "\\nHow would you like to select protein parts?"\
    puts "1) By segment ID (e.g., PROA, PROB)"\
    puts "2) By chain ID (e.g., A, B)"\
    flush stdout\
    set select_method [gets stdin]\
    \
    # Setup for SAMHD1 and SOX11 selections\
    set samhd1_parts ""\
    set sox11_parts ""\
    \
    if \{$select_method == "2"\} \{\
        # Chain selection mode\
        puts "\\nEnter chain ID(s) for SAMHD1 (space separated):"\
        flush stdout\
        set samhd1_parts [gets stdin]\
        \
        puts "Enter chain ID(s) for SOX11:"\
        flush stdout\
        set sox11_parts [gets stdin]\
        \
        # Build selection commands\
        set samhd1_sel_command ""\
        foreach part [split $samhd1_parts] \{\
            if \{$samhd1_sel_command != ""\} \{\
                append samhd1_sel_command " or "\
            \}\
            append samhd1_sel_command "chain $part"\
        \}\
        \
        set sox11_sel_command ""\
        foreach part [split $sox11_parts] \{\
            if \{$sox11_sel_command != ""\} \{\
                append sox11_sel_command " or "\
            \}\
            append sox11_sel_command "chain $part"\
        \}\
    \} else \{\
        # Segment selection mode (default)\
        puts "\\nEnter segment ID(s) for SAMHD1 (space separated):"\
        flush stdout\
        set samhd1_parts [gets stdin]\
        \
        puts "Enter segment ID(s) for SOX11:"\
        flush stdout\
        set sox11_parts [gets stdin]\
        \
        # Build selection commands\
        set samhd1_sel_command ""\
        foreach part [split $samhd1_parts] \{\
            if \{$samhd1_sel_command != ""\} \{\
                append samhd1_sel_command " or "\
            \}\
            append samhd1_sel_command "segname $part"\
        \}\
        \
        set sox11_sel_command ""\
        foreach part [split $sox11_parts] \{\
            if \{$sox11_sel_command != ""\} \{\
                append sox11_sel_command " or "\
            \}\
            append sox11_sel_command "segname $part"\
        \}\
    \}\
    \
    # Test SAMHD1 and SOX11 selections\
    puts "\\nTesting selections..."\
    set samhd1_sel [atomselect top "protein and ($samhd1_sel_command)" frame 0]\
    set sox11_sel [atomselect top "protein and ($sox11_sel_command)" frame 0]\
    \
    puts "SAMHD1: [$samhd1_sel num] atoms selected with: protein and ($samhd1_sel_command)"\
    puts "SOX11: [$sox11_sel num] atoms selected with: protein and ($sox11_sel_command)"\
    \
    if \{[$samhd1_sel num] == 0\} \{\
        puts "Error: No atoms found for SAMHD1. Please check your selection."\
        cleanup_selections [list $samhd1_sel $sox11_sel]\
        return\
    \}\
    \
    if \{[$sox11_sel num] == 0\} \{\
        puts "Error: No atoms found for SOX11. Please check your selection."\
        cleanup_selections [list $samhd1_sel $sox11_sel]\
        return\
    \}\
    \
    # Get SAMHD1 chains (for tetramer analysis)\
    if \{$select_method == "2"\} \{\
        # Using chain IDs\
        set samhd1_chains [split $samhd1_parts]\
    \} else \{\
        # Using segments - need to determine corresponding chains\
        set samhd1_chains \{\}\
        foreach part [split $samhd1_parts] \{\
            set seg_sel [atomselect top "segname $part" frame 0]\
            set seg_chains [lsort -unique [$seg_sel get chain]]\
            set samhd1_chains [concat $samhd1_chains $seg_chains]\
            $seg_sel delete\
        \}\
        set samhd1_chains [lsort -unique $samhd1_chains]\
    \}\
    \
    cleanup_selections [list $samhd1_sel $sox11_sel]\
    \
    # ===== User Input for Analysis Type =====\
    puts "\\nThis script analyzes the SAMHD1-SOX11 complex."\
    puts "Select analysis type:"\
    puts "1) Analyze SAMHD1 as one unit vs SOX11"\
    puts "2) Analyze each SAMHD1 part individually vs SOX11"\
    puts "3) Analyze specific part pairs"\
    puts "4) Analyze SAMHD1 assembly (interactions between SAMHD1 parts)"\
    flush stdout\
    set analysis_type [gets stdin]\
    \
    switch $analysis_type \{\
        "1" \{\
            # SAMHD1 as one unit vs SOX11\
            set group1_selection $samhd1_sel_command\
            set group1_name "SAMHD1"\
            set group2_selection $sox11_sel_command\
            set group2_name "SOX11"\
        \}\
        "2" \{\
            # Each SAMHD1 part vs SOX11\
            puts "\\nSAMHD1 has multiple parts: $samhd1_parts"\
            puts "Which SAMHD1 part would you like to analyze against SOX11?"\
            flush stdout\
            set selected_part [gets stdin]\
            \
            # Build selection for the specific part\
            if \{$select_method == "2"\} \{\
                # Chain selection\
                set group1_selection "chain $selected_part"\
            \} else \{\
                # Segment selection\
                set group1_selection "segname $selected_part"\
            \}\
            \
            set group1_name "SAMHD1_$\{selected_part\}"\
            set group2_selection $sox11_sel_command\
            set group2_name "SOX11"\
        \}\
        "3" \{\
            # Custom part pairs\
            puts "\\nEnter custom selection for first group:"\
            puts "Available parts: [concat $samhd1_parts $sox11_parts]"\
            if \{$select_method == "2"\} \{\
                puts "Enter chain(s) for first group (space separated):"\
            \} else \{\
                puts "Enter segment(s) for first group (space separated):"\
            \}\
            flush stdout\
            set custom_parts1 [gets stdin]\
            \
            puts "Enter a name for this group:"\
            flush stdout\
            set group1_name [gets stdin]\
            \
            puts "\\nEnter custom selection for second group:"\
            if \{$select_method == "2"\} \{\
                puts "Enter chain(s) for second group (space separated):"\
            \} else \{\
                puts "Enter segment(s) for second group (space separated):"\
            \}\
            flush stdout\
            set custom_parts2 [gets stdin]\
            \
            puts "Enter a name for this group:"\
            flush stdout\
            set group2_name [gets stdin]\
            \
            # Build selection strings\
            set group1_selection ""\
            foreach part [split $custom_parts1] \{\
                if \{$group1_selection != ""\} \{\
                    append group1_selection " or "\
                \}\
                if \{$select_method == "2"\} \{\
                    append group1_selection "chain $part"\
                \} else \{\
                    append group1_selection "segname $part"\
                \}\
            \}\
            \
            set group2_selection ""\
            foreach part [split $custom_parts2] \{\
                if \{$group2_selection != ""\} \{\
                    append group2_selection " or "\
                \}\
                if \{$select_method == "2"\} \{\
                    append group2_selection "chain $part"\
                \} else \{\
                    append group2_selection "segname $part"\
                \}\
            \}\
        \}\
        "4" \{\
            # SAMHD1 tetramer/assembly analysis\
            puts "Analyzing SAMHD1 assembly interactions"\
            set group1_selection $samhd1_sel_command\
            set group1_name "SAMHD1_assembly"\
            set group2_selection ""\
            set group2_name ""\
        \}\
        default \{\
            puts "Invalid choice. Using SAMHD1 vs SOX11 analysis."\
            set group1_selection $samhd1_sel_command\
            set group1_name "SAMHD1"\
            set group2_selection $sox11_sel_command\
            set group2_name "SOX11"\
        \}\
    \}\
    \
    # Test final selections\
    if \{$analysis_type != "4"\} \{\
        puts "\\nVerifying final selections..."\
        set test1 [atomselect top "protein and ($group1_selection)" frame 0]\
        set test2 [atomselect top "protein and ($group2_selection)" frame 0]\
        \
        puts "Group 1 ($group1_name): [$test1 num] atoms selected with: protein and ($group1_selection)"\
        puts "Group 2 ($group2_name): [$test2 num] atoms selected with: protein and ($group2_selection)"\
        \
        if \{[$test1 num] == 0 || [$test2 num] == 0\} \{\
            puts "Error: One or both selections have no atoms. Please check your selections."\
            cleanup_selections [list $test1 $test2]\
            return\
        \}\
        \
        cleanup_selections [list $test1 $test2]\
    \} else \{\
        puts "\\nVerifying SAMHD1 assembly selection..."\
        set test1 [atomselect top "protein and ($group1_selection)" frame 0]\
        puts "SAMHD1 assembly: [$test1 num] atoms selected with: protein and ($group1_selection)"\
        \
        if \{[$test1 num] == 0\} \{\
            puts "Error: No atoms found for SAMHD1. Please check your selection."\
            $test1 delete\
            return\
        \}\
        \
        $test1 delete\
    \}\
    \
    # Create output directory with timestamp\
    set timestamp [clock format [clock seconds] -format "%Y%m%d_%H%M%S"]\
    if \{$analysis_type == "4"\} \{\
        set outdir "$\{group1_name\}_$\{timestamp\}"\
    \} else \{\
        set outdir "$\{group1_name\}_$\{group2_name\}_$\{timestamp\}"\
    \}\
    catch \{file mkdir $outdir\}\
    puts "Output will be saved to directory: $outdir"\
    \
    # ===== Select Analyses to Run =====\
    puts "\\nSelect analyses to run (1=yes, 0=no):"\
    puts "1) Alignment (recommended before other analyses)"\
    puts "2) RMSD"\
    puts "3) RMSF"\
    puts "4) Radius of Gyration"\
    puts "5) SASA"\
    puts "6) Hydrogen Bonds"\
    puts "7) COM-COM Distance"\
    puts "8) Interface SASA"\
    puts "9) Contact Analysis"\
    puts "A) Run ALL analyses"\
    puts "0) Run NONE (exit)"\
    flush stdout\
    \
    # Default to run all\
    set run_alignment 1\
    set run_rmsd 1\
    set run_rmsf 1\
    set run_rog 1\
    set run_sasa 1\
    set run_hbonds 1\
    set run_comcom 1\
    set run_interface_sasa 1\
    set run_contacts 1\
    \
    set choice [gets stdin]\
    switch -nocase $choice \{\
        "1" \{\
            # Just alignment\
            set run_rmsd 0\
            set run_rmsf 0\
            set run_rog 0\
            set run_sasa 0\
            set run_hbonds 0\
            set run_comcom 0\
            set run_interface_sasa 0\
            set run_contacts 0\
        \}\
        "2" \{ \
            # Just RMSD\
            set run_alignment 1  # Need alignment for RMSD\
            set run_rmsf 0\
            set run_rog 0\
            set run_sasa 0\
            set run_hbonds 0\
            set run_comcom 0\
            set run_interface_sasa 0\
            set run_contacts 0\
        \}\
        "3" \{ \
            # Just RMSF\
            set run_alignment 1  # Need alignment for RMSF\
            set run_rmsd 0\
            set run_rog 0\
            set run_sasa 0\
            set run_hbonds 0\
            set run_comcom 0\
            set run_interface_sasa 0\
            set run_contacts 0\
        \}\
        "4" \{ \
            # Just ROG\
            set run_alignment 0\
            set run_rmsd 0\
            set run_rmsf 0\
            set run_sasa 0\
            set run_hbonds 0\
            set run_comcom 0\
            set run_interface_sasa 0\
            set run_contacts 0\
        \}\
        "5" \{ \
            # Just SASA\
            set run_alignment 0\
            set run_rmsd 0\
            set run_rmsf 0\
            set run_rog 0\
            set run_hbonds 0\
            set run_comcom 0\
            set run_interface_sasa 0\
            set run_contacts 0\
        \}\
        "6" \{ \
            # Just H-bonds\
            set run_alignment 0\
            set run_rmsd 0\
            set run_rmsf 0\
            set run_rog 0\
            set run_sasa 0\
            set run_comcom 0\
            set run_interface_sasa 0\
            set run_contacts 0\
        \}\
        "7" \{ \
            # Just COM-COM\
            set run_alignment 0\
            set run_rmsd 0\
            set run_rmsf 0\
            set run_rog 0\
            set run_sasa 0\
            set run_hbonds 0\
            set run_interface_sasa 0\
            set run_contacts 0\
        \}\
        "8" \{ \
            # Just Interface SASA\
            set run_alignment 0\
            set run_rmsd 0\
            set run_rmsf 0\
            set run_rog 0\
            set run_sasa 0\
            set run_hbonds 0\
            set run_comcom 0\
            set run_contacts 0\
        \}\
        "9" \{ \
            # Just Contacts\
            set run_alignment 0\
            set run_rmsd 0\
            set run_rmsf 0\
            set run_rog 0\
            set run_sasa 0\
            set run_hbonds 0\
            set run_comcom 0\
            set run_interface_sasa 0\
        \}\
        "a" -\
        "A" \{ \
            # All analyses (default)\
        \}\
        "0" \{ \
            puts "Exiting without analyses"\
            return\
        \}\
        default \{ \
            puts "Invalid choice, using defaults (all analyses)"\
        \}\
    \}\
    \
    # ===== EXECUTE ANALYSES =====\
    \
    # ===== Analysis modules =====\
    # Choose between standard analysis or SAMHD1 assembly analysis\
    if \{$analysis_type == "4"\} \{\
        # Execute SAMHD1 assembly analysis\
        puts "\\nExecuting SAMHD1 assembly analysis..."\
        if \{[catch \{\
            execute_assembly_analysis $nf $outdir $samhd1_sel_command $samhd1_parts $samhd1_chains $select_method $run_alignment $run_rmsd $run_rmsf $run_rog $run_sasa $run_hbonds $run_comcom $run_interface_sasa $run_contacts\
        \} err]\} \{\
            puts "Error in assembly analysis: $err"\
        \}\
    \} else \{\
        # Execute standard analysis\
        puts "\\nExecuting standard analysis..."\
        if \{[catch \{\
            execute_standard_analysis $nf $outdir $group1_selection $group1_name $group2_selection $group2_name $run_alignment $run_rmsd $run_rmsf $run_rog $run_sasa $run_hbonds $run_comcom $run_interface_sasa $run_contacts\
        \} err]\} \{\
            puts "Error in standard analysis: $err"\
        \}\
    \}\
    \
    # ===== Summary =====\
    puts "\\n============================================="\
    puts "All analyses complete!"\
    puts "Output files saved to: $outdir/"\
    puts "============================================="\
\}\
\
# Procedure for standard analysis between two groups\
proc execute_standard_analysis \{nf outdir group1_selection group1_name group2_selection group2_name run_alignment run_rmsd run_rmsf run_rog run_sasa run_hbonds run_comcom run_interface_sasa run_contacts\} \{\
    # ===== Alignment =====\
    if \{$run_alignment\} \{\
        puts "\\nAligning trajectory based on $group1_name backbone..."\
        if \{[catch \{\
            set reference [atomselect top "protein and ($group1_selection) and backbone" frame 0]\
            set compare [atomselect top "protein and ($group1_selection) and backbone"]\
            set all [atomselect top "all"]\
            \
            for \{set frame 0\} \{$frame < $nf\} \{incr frame\} \{\
                $compare frame $frame\
                $all frame $frame\
                \
                set trans_mat [measure fit $compare $reference]\
                $all move $trans_mat\
                \
                if \{$frame % 100 == 0\} \{\
                    puts "  Aligned frame $frame / $nf"\
                \}\
            \}\
            \
            cleanup_selections [list $reference $compare $all]\
            puts "Alignment complete"\
        \} err]\} \{\
            puts "Error during alignment: $err"\
        \}\
    \}\
    \
    # ===== RMSD =====\
    if \{$run_rmsd\} \{\
        # Group 1 RMSD\
        puts "\\nCalculating RMSD of $group1_name..."\
        if \{[catch \{\
            set reference [atomselect top "protein and ($group1_selection) and backbone" frame 0]\
            set compare [atomselect top "protein and ($group1_selection) and backbone"]\
            set outfile [open "$outdir/RMSD_$\{group1_name\}.dat" w]\
            \
            for \{set frame 0\} \{$frame < $nf\} \{incr frame\} \{\
                $compare frame $frame\
                set rmsd [measure rmsd $compare $reference]\
                puts $outfile "$frame $rmsd"\
                \
                if \{$frame % 100 == 0\} \{\
                    puts "  $group1_name RMSD frame $frame: $rmsd"\
                \}\
            \}\
            \
            close $outfile\
            cleanup_selections [list $reference $compare]\
            puts "$group1_name RMSD calculation complete"\
        \} err]\} \{\
            puts "Error calculating $group1_name RMSD: $err"\
        \}\
        \
        # Group 2 RMSD\
        puts "\\nCalculating RMSD of $group2_name..."\
        if \{[catch \{\
            set reference [atomselect top "protein and ($group2_selection) and backbone" frame 0]\
            set compare [atomselect top "protein and ($group2_selection) and backbone"]\
            set outfile [open "$outdir/RMSD_$\{group2_name\}.dat" w]\
            \
            for \{set frame 0\} \{$frame < $nf\} \{incr frame\} \{\
                $compare frame $frame\
                set rmsd [measure rmsd $compare $reference]\
                puts $outfile "$frame $rmsd"\
                \
                if \{$frame % 100 == 0\} \{\
                    puts "  $group2_name RMSD frame $frame: $rmsd"\
                \}\
            \}\
            \
            close $outfile\
            cleanup_selections [list $reference $compare]\
            puts "$group2_name RMSD calculation complete"\
        \} err]\} \{\
            puts "Error calculating $group2_name RMSD: $err"\
        \}\
        \
        # Complex RMSD\
        puts "\\nCalculating complex RMSD..."\
        if \{[catch \{\
            set reference [atomselect top "protein and (($group1_selection) or ($group2_selection)) and backbone" frame 0]\
            set compare [atomselect top "protein and (($group1_selection) or ($group2_selection)) and backbone"]\
            set outfile [open "$outdir/RMSD_complex.dat" w]\
            \
            for \{set frame 0\} \{$frame < $nf\} \{incr frame\} \{\
                $compare frame $frame\
                set rmsd [measure rmsd $compare $reference]\
                puts $outfile "$frame $rmsd"\
                \
                if \{$frame % 100 == 0\} \{\
                    puts "  Complex RMSD frame $frame: $rmsd"\
                \}\
            \}\
            \
            close $outfile\
            cleanup_selections [list $reference $compare]\
            puts "Complex RMSD calculation complete"\
        \} err]\} \{\
            puts "Error calculating complex RMSD: $err"\
        \}\
    \}\
    \
    # ===== RMSF =====\
    if \{$run_rmsf\} \{\
        puts "\\nCalculating RMSF..."\
        \
        # Group 1 RMSF\
        if \{[catch \{\
            set sel [atomselect top "protein and ($group1_selection) and name CA"]\
            set outfile [open "$outdir/RMSF_$\{group1_name\}.dat" w]\
            \
            set rmsf [measure rmsf $sel]\
            set resids [$sel get resid]\
            set chains [$sel get chain]\
            set resnames [$sel get resname]\
            \
            puts $outfile "# Chain ResID ResName RMSF"\
            for \{set i 0\} \{$i < [$sel num]\} \{incr i\} \{\
                puts $outfile "[lindex $chains $i] [lindex $resids $i] [lindex $resnames $i] [lindex $rmsf $i]"\
            \}\
            \
            close $outfile\
            $sel delete\
            puts "RMSF calculation for $group1_name complete"\
        \} err]\} \{\
            puts "Error calculating $group1_name RMSF: $err"\
        \}\
        \
        # Group 2 RMSF\
        if \{[catch \{\
            set sel [atomselect top "protein and ($group2_selection) and name CA"]\
            set outfile [open "$outdir/RMSF_$\{group2_name\}.dat" w]\
            \
            set rmsf [measure rmsf $sel]\
            set resids [$sel get resid]\
            set chains [$sel get chain]\
            set resnames [$sel get resname]\
            \
            puts $outfile "# Chain ResID ResName RMSF"\
            for \{set i 0\} \{$i < [$sel num]\} \{incr i\} \{\
                puts $outfile "[lindex $chains $i] [lindex $resids $i] [lindex $resnames $i] [lindex $rmsf $i]"\
            \}\
            \
            close $outfile\
            $sel delete\
            puts "RMSF calculation for $group2_name complete"\
        \} err]\} \{\
            puts "Error calculating $group2_name RMSF: $err"\
        \}\
    \}\
    \
    # ===== Radius of Gyration =====\
    if \{$run_rog\} \{\
        puts "\\nCalculating Radius of Gyration..."\
        \
        # Group 1 RoG\
        if \{[catch \{\
            set outfile [open "$outdir/rog_$\{group1_name\}.dat" w]\
            puts $outfile "frame rad_of_gyr"\
            set prot [atomselect top "protein and ($group1_selection)"]\
            \
            for \{set i 0\} \{$i < $nf\} \{incr i\} \{\
                $prot frame $i\
                $prot update\
                set rog [gyr_radius $prot]\
                puts $outfile "$i $rog"\
                \
                if \{$i % 100 == 0\} \{\
                    puts "  RoG $group1_name frame $i: $rog"\
                \}\
            \}\
            \
            close $outfile\
            $prot delete\
            puts "Radius of Gyration calculation for $group1_name complete"\
        \} err]\} \{\
            puts "Error calculating $group1_name RoG: $err"\
        \}\
        \
        # Group 2 RoG\
        if \{[catch \{\
            set outfile [open "$outdir/rog_$\{group2_name\}.dat" w]\
            puts $outfile "frame rad_of_gyr"\
            set prot [atomselect top "protein and ($group2_selection)"]\
            \
            for \{set i 0\} \{$i < $nf\} \{incr i\} \{\
                $prot frame $i\
                $prot update\
                set rog [gyr_radius $prot]\
                puts $outfile "$i $rog"\
                \
                if \{$i % 100 == 0\} \{\
                    puts "  RoG $group2_name frame $i: $rog"\
                \}\
            \}\
            \
            close $outfile\
            $prot delete\
            puts "Radius of Gyration calculation for $group2_name complete"\
        \} err]\} \{\
            puts "Error calculating $group2_name RoG: $err"\
        \}\
    \}\
    \
    # ===== SASA =====\
    if \{$run_sasa\} \{\
        puts "\\nCalculating Solvent Accessible Surface Area..."\
        \
        # Group 1 SASA\
        if \{[catch \{\
            set sel [atomselect top "protein and ($group1_selection)"]\
            set output [open "$outdir/SASA_$\{group1_name\}.dat" w]\
            \
            for \{set i 0\} \{$i < $nf\} \{incr i\} \{\
                molinfo top set frame $i\
                set sasa [measure sasa 1.4 $sel]\
                \
                if \{$i % 100 == 0\} \{\
                    puts "  SASA $group1_name progress: $i/$nf"\
                \}\
                puts $output "$i $sasa"\
            \}\
            \
            close $output\
            $sel delete\
            puts "SASA calculation for $group1_name complete"\
        \} err]\} \{\
            puts "Error calculating $group1_name SASA: $err"\
        \}\
        \
        # Group 2 SASA\
        if \{[catch \{\
            set sel [atomselect top "protein and ($group2_selection)"]\
            set output [open "$outdir/SASA_$\{group}