#!/usr/bin/perl

%res_num_to_rdc = ();

print "started computation.\n";
open(MR, "<CH_unparsed.m");
while ($line = <MR>) {
    print "line read 1.\n";
    if ($line =~ /^! RDCs/) {
	print "line read 2.\n";
	if (keys(%res_num_to_rdc) > 0) {
	    if (-e "parsedResonances.txt") {
		open(OUT, "> $rdc_type.m");
		open(HSQC, "< parsedResonances.txt");
		while ($line = <HSQC>) {
		    if ($line =~ m/^(\d+)\s+/) {
			$res_num = $1;
			if (exists($res_num_to_rdc{$res_num})) {
			    $rdc = $res_num_to_rdc{$res_num};
			    $rdc = $rdc * 0.491;
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
	$rdc_type = "CH";
    }

    if ($line =~ /^\s+(\d+)\s+\w\w\w\s+\w\w\s+\d+\s+\w\w\w\s+\w\w\W*\s+(\-*\d+.\d+)/) {
	$res_num   = $1;
	$rdc       = $2;
	print "individual residue read. resnum = $res_num.\n";
	print "individual RDC read. rdc = $rdc.\n";

	
#	if ($atom_type eq "HA") {
	$res_num_to_rdc{$res_num} = $rdc;
#	} elsif ($atom_type eq "HN") {
#	    $res_num_to_rdc{$res_num} = $rdc;
#	}
    }
}
close(MR);
