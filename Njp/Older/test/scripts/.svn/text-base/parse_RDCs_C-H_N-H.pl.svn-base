#!/usr/bin/perl

%res_num_to_rdc = ();

open(MR, "< 1D3Z.mr");
while ($line = <MR>) {
	if ($line =~ m/^!!! DipolarCouplings\.(.+)/) {
		if (keys(%res_num_to_rdc) > 0) {
			if (-e "hsqcdata.m") {
				open(OUT, "> $rdc_type.m");
				open(HSQC, "< hsqcdata.m");
				$peak = 1;
				while ($line = <HSQC>) {
					if ($line =~ m/^(\d+)\s+/) {
						$res_num = $1;
		
						if (exists($res_num_to_rdc{$res_num})) {
							$rdc = $res_num_to_rdc{$res_num};
							if ($rdc_type =~ m/C-H/) {
								$rdc = $rdc * 0.491;
							}
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
		$rdc_type = $1;
	}
	if ($line =~ m/resid\s+(\d+) and name\s+(\w+)\)\s+(\-*\d+.\d+).+$/) {
		$res_num = $1;
		$atom_type = $2;
		$rdc = $3;

		if ($atom_type eq "HA") {
			$res_num_to_rdc{$res_num} = $rdc;
		} elsif ($atom_type eq "HN") {
			$res_num_to_rdc{$res_num} = $rdc;
		}
	}
}
close(MR);