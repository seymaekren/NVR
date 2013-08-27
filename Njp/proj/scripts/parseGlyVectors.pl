#!/usr/bin/perl

###############################################
# parse_vectors_C-H.pl                        #
# calculate CA-HA bond vectors from .pdb file #
#                                             #
# Nick Patrick                                #
###############################################

%res_num_to_vector_ref1 = ();
%res_num_to_vector_ref2 = ();
$max_res = -1;

open(PDB, "< 1UBI.pdb") or die "can't open file for reading";
while ($line = <PDB>) {
	if ($line =~ m/^ATOM/) {
		$atom_type = substr($line, 12, 3);
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

		if ($atom_type eq " CA") {
			$ca_x = $atom_x;
			$ca_y = $atom_y;
			$ca_z = $atom_z;
			
			if ($res_num > $max_res) {
				$max_res = $res_num;
			}
		} elsif ($atom_type eq "1HA") {
			#normalize, store the vector
			$x_comp = $atom_x - $ca_x;
			$y_comp = $atom_y - $ca_y;
			$z_comp = $atom_z - $ca_z;

			$length = sqrt($x_comp * $x_comp + $y_comp * $y_comp + $z_comp * $z_comp);
			$vector_ref = [$x_comp / $length, $y_comp / $length, $z_comp / $length];

			$res_num_to_vector_ref1{$res_num} = $vector_ref;
		} elsif ($atom_type eq "2HA") {
			#normalize, store the vector
			$x_comp = $atom_x - $ca_x;
			$y_comp = $atom_y - $ca_y;
			$z_comp = $atom_z - $ca_z;

			$length = sqrt($x_comp * $x_comp + $y_comp * $y_comp + $z_comp * $z_comp);
			$vector_ref = [$x_comp / $length, $y_comp / $length, $z_comp / $length];

			$res_num_to_vector_ref2{$res_num} = $vector_ref;
		} 
	}
}
close(PDB);

######################################
# write data to file CA-HA_vectors.m #
######################################

open(HSQC, "< myinput.m");
open(OUT, "> glycineC-H_vectors.m") or die "can't open file for writing";

while ($line = <HSQC>) {
    if ($line =~ m/^(\d+)\s+/) {
	$res_num = $1;
	
	$i = $res_num; #for some reason doesn't work otherwise

	if (exists($res_num_to_vector_ref1{$i})) {
		$vector_ref = $res_num_to_vector_ref1{$i};
		@vector = @$vector_ref;
		print OUT "$i\t$vector[0]\t$vector[1]\t$vector[2]\n";
	} 
	if (exists($res_num_to_vector_ref2{$i})) {
	    $vector_ref = $res_num_to_vector_ref2{$i};
	    @vector = @$vector_ref;
	    print OUT "$i\t$vector[0]\t$vector[1]\t$vector[2]\n";
	}

#else {
	#	print OUT "$i\t-999\t-999\t-999\n";
	#}
    }
}
#for ($i = 1; $i <= $max_res; $i++) {
#	if (exists($res_num_to_vector_ref{$i})) {
#		$vector_ref = $res_num_to_vector_ref{$i};
#		@vector = @$vector_ref;
#		print OUT "$i\t$vector[0]\t$vector[1]\t$vector[2]\n";
#	} else {
#		print OUT "$i\t-999\t-999\t-999\n";
#	}
#}
close(OUT);
