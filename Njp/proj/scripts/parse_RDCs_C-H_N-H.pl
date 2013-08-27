#!/usr/bin/perl

%res_num_to_rdc = ();


#open(MR, "< CH_unparsed.mr");
open(MR, "< 1C05_rdcs.txt");
while ($line = <MR>) {
#    if ($line =~ m/^! RDC\((.+)\),/) {
#    if ($line =~ m/^! RDCs/) {
    if ($line =~ m/^!! Dipolar coupling restraints/) {
	if (keys(%res_num_to_rdc) > 0) {
	    if (-e "parsedResonances.txt") {
		open(OUT, "> $rdc_type.m");
		open(HSQC, "< parsedResonances.txt");
		while ($line = <HSQC>) {
		    if ($line =~ m/^(\d+)\s+/) {
			$res_num = $1;
			if (exists($res_num_to_rdc{$res_num})) {
			    $rdc = $res_num_to_rdc{$res_num};
#			    if ($rdc_type =~ m/C-H/) {
#				    $rdc = $rdc * 0.491;
#				}
			} else {
			    $rdc = -999;
			}
			
			print OUT "$res_num\t$rdc\n";
		    }
		}
		close(HSQC);
		close(OUT);
	    }
	}
	%res_num_to_rdc = ();
#	$rdc_type = $1;
	$rdc_type = "NH";
    }

    if ($line =~ /resid\s+(\d+)\s+and name\s+(\w+)\s*\)\s+(\-*\d+.\d+).+$/) {
	$res_num   = $1;
	$atom_type = $2;
	$rdc       = $3;
	print "individual RDC read. atomType = $atom_type.\n";

	
	if ($atom_type eq "HA") {
	    $res_num_to_rdc{$res_num} = $rdc;
	} elsif ($atom_type eq "HN") {
	    $res_num_to_rdc{$res_num} = $rdc;
	}
    }
}
close(MR);
