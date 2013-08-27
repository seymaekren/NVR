function   printNumAvailableRDC_Values(CH_RDCS, NH_RDCS);
	   
nonAvailableCH_RDCs = find(CH_RDCS == -999);
nonAvailableNH_RDCs = find(NH_RDCS == -999);

fprintf(1, 'numAvailable CH_RDCs = %d\n', length(CH_RDCS) ...
	- length(nonAvailableCH_RDCs));
fprintf(1, 'numAvailable NH_RDCs = %d\n', length(NH_RDCS) ...
	- length(nonAvailableNH_RDCs));

fprintf(1, 'dont forget to remove the removed peaks from these.\n');
%WHAT I mean above is that in order to compute the number
%of available CH or NH RDCs we have to substract from the
%number of available CH or NH RDCs a number corresponding
%to the removed peaks.
keyboard