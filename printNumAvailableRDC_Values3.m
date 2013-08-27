function   printNumAvailableRDC_Values3(CCa_RDCS, NH_RDCS, NC_RDCS);
	   
nonAvailableCCa_RDCs = find(CCa_RDCS == -999);
nonAvailableNH_RDCs  = find(NH_RDCS == -999);
nonAvailableNC_RDCs  = find(NC_RDCS == -999);

fprintf(1, 'numAvailable CCa_RDCs = %d\n', length(CCa_RDCS) ...
	- length(nonAvailableCCa_RDCs));
fprintf(1, 'numAvailable NH_RDCs = %d\n', length(NH_RDCS) ...
	- length(nonAvailableNH_RDCs));
fprintf(1, 'numAvailable NC_RDCs = %d\n', length(NC_RDCS) ...
	- length(nonAvailableNC_RDCs));

fprintf(1, 'dont forget to remove the removed peaks from these.\n');
%WHAT I mean above is that in order to compute the number
%of available CH or NH RDCs we have to substract from the
%number of available CH or NH RDCs a number corresponding
%to the removed peaks.
keyboard