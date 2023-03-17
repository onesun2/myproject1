#!/usr/bin/tclsh
package require Tcl 8
package require math::statistics
#starting/ending monomer number
set sm 2
set em 9 
########################
set filename "angle${sm}to${em}.txt"
set fileId [open $filename "w"]
set num_steps [molinfo top get numframes]
set data {}
set thetasum 0.0
set ttt {}
for {set i $sm} {$i <= $em} {incr i} {
                set monomer($i) [atomselect top "(name C01 C02 C03 C04 C05 C06) and (segname PLY) and (resid $i)"]
}

proc angle {a1 b1 a2 b2} {
        set x1 [expr $a1/(sqrt($a1*$a1+$b1*$b1))]
        set y1 [expr $b1/(sqrt($a1*$a1+$b1*$b1))]
        set x2 [expr $a2/(sqrt($a2*$a2+$b2*$b2))]
        set y2 [expr $b2/(sqrt($a2*$a2+$b2*$b2))]   
        set M_PI 3.14159265358979323846
        set theta [expr (180/$M_PI)*acos(($x2*$x1+$y2*$y1)/($x1*$x1+$y1*$y1))]
        set sign [expr $x1*$y2 - $x2*$y1]
        if {$sign >= 0} {
                set temp 1
        } elseif {$sign < 0} {
                set temp -1
        }
        return [expr $temp * $theta]
}

proc geom_center {selection} {
        set gc [veczero]
        foreach coord [$selection get {x y z}] {
                set gc [vecadd $gc $coord]
        }
        return [vecscale [expr 1.0/[$selection num]] $gc]
}

proc Bestline {l} {
	set sumx 0.0
	set sumy 0.0
	set sum2x 0.0
	set sum2y 0.0
	set sumxy 0.0
	set num 0

	for {set i 0} {$i < [llength $l]/2} {incr i} {
		set x [lindex $l [expr 2*$i]]
		set y [lindex $l [expr 2*$i+1]]
		set sumx [expr $x + $sumx]
		set sum2x [expr $x * $x + $sum2x]
		set sumy [expr $y + $sumy]
		set sum2y [expr $y * $y + $sum2y]
		set sumxy [expr $x * $y + $sumxy]
		incr num
	}

	set denx [expr $num * $sum2x - $sumx * $sumx]
	set deny [expr $num * $sum2y - $sumy * $sumy]
	set top  [expr $num * $sumxy - $sumx * $sumy]

	set corr [expr $top / sqrt($denx * $deny)]
	set slp [expr $top / $denx]
	set inter [expr ($sumy - $sumx * $slp) / $num]

	return "$slp $inter $corr"
}
###################################################################
###           MAIN                                             ####
###################################################################
for {set frame 0} {$frame < $num_steps} {incr frame} {
#	puts "frame = $frame"
	for {set i $sm} {$i <= $em} {incr i} {
                $monomer($i) frame $frame
        }
	for {set i $sm} {$i <= $em} {incr i} {
                set cm($i) [geom_center $monomer($i)]
        }
	for {set i $sm} {$i < $em} {incr i} {
		set x1 [lindex $cm($i) 0]
		set y1 [lindex $cm($i) 1]
		set x2 [lindex $cm([expr $i+1]) 0]
		set y2 [lindex $cm([expr $i+1]) 1]
		set z2 [lindex $cm([expr $i+1]) 2]

                set theta [angle $x1 $y1 $x2 $y2] 
#		puts $fileId "angle($i,[expr $i+1]) = $theta ; z = $z2" 
#		puts "angle($i,[expr $i+1]) = $theta ; z = $z2" 
		set thetasum [expr $thetasum + $theta]
		lappend data $thetasum 
		lappend data $z2
        }
		set fit [Bestline $data]
		set slp [lindex $fit 0]
		set pitch [expr 360*$slp]
		lappend ttt $pitch
		puts $fileId "spl = $slp, pitch = $pitch, inter = [lindex $fit 1], corr = [lindex $fit 2]"
		set data {}
		set thetasum 0.0
}
puts $fileId "pitch_av2 = [::math::statistics::mean $ttt], stdev = [::math::statistics::stdev $ttt]"
puts "pitch_av2 = [::math::statistics::mean $ttt], stdev = [::math::statistics::stdev $ttt]"
set ttt {} 
close $fileId
