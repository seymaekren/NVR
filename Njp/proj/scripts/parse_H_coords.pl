#!/usr/bin/perl

################################################
# parse_H_coords.pl                            #
# parse amide proton coordinates from PDB file #
#                                              #
# Nick Patrick                                 #
################################################

%res_num_to_coords = ();
$max_res = -1;

open(PDB, "< 1UBI.pdb") or die "can't open file for reading";
while ($line = <PDB>) {
	if ($line =~ m/^ATOM/) {
		$atom_type = substr($line, 13, 2);
		$res_type = substr($line, 17, 3);
		$res_num = substr($line, 22, 5);
		$res_num =~ s/^\s+//;
		$res_num =~ s/\s+$//;
		$atom_x = substr($line, 30, 8);
		$atom_y = substr($line, 38, 8);
		$atom_z = substr($line, 46, 8);

		if ($res_type eq "HOH") {
			last;
		}

		if ($atom_type eq "H ") {
            if ($res_num > $max_res) {
				$max_res = $res_num;
			}
			$coords = [$atom_x, $atom_y, $atom_z];
			$res_num_to_coords{$res_num} = $coords;
		}
	}
}
close(PDB);

#################################
# write data to file H_coords.m #
#################################

open(OUT, "> H_coords.m") or die "can't open file for writing";
for ($i = 1; $i <= $max_res; $i++) {
	if (exists($res_num_to_coords{$i})) {
		$coord_ref = $res_num_to_coords{$i};
		@coords = @$coord_ref;
		print OUT "$i\t$coords[0]\t$coords[1]\t$coords[2]\n";
	} else {
		print OUT "$i\t-999\t-999\t-999\n";
	}
}
close(OUT);