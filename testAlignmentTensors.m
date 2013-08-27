load allPeaksAfterInitialAssignments--For1UBI.mat

fprintf(1, 'there are %d assignments in MASTER.\n', sum(sum(MASTER)));


S1_partial                               = updateTen(MASTER,RDC1,VECTORS);
S2_partial                               = updateTen(MASTER,RDC2,VECTORS);



fullMASTER = MASTER;
for peakIndex = 1:size(MASTER,1)
  fullMASTER(peakIndex,peakIndex) = 1;
  assert (sum(fullMASTER(peakIndex,:)) == 1);
end

S1_full                                       = updateTen(fullMASTER, ...
						  RDC1, VECTORS);

S2_full                                       = updateTen(fullMASTER, ...
						  RDC2, VECTORS);

percentile1                               = ...
    NVR_COMP_TEN(S1_partial,S1_full);

percentile2                               =     NVR_COMP_TEN(S2_partial,S2_full);

fprintf(1, 'the percentile difference between the partial and full');
fprintf(1, ' is %f and %f\n', percentile1, percentile2);

figure; plot(percentile1, percentile2, 'r*'), hold on;

maxNumIter = 1000; numPartialAssignments = 5;
percentiles1 = []; percentiles2= [];

for iterIndex = 1:maxNumIter
  randomMASTER = MASTER*0;
  randomOrder = randperm(size(MASTER,1));
  for relPeakIndex  = 1:numPartialAssignments
    peakIndex = randomOrder(relPeakIndex);
    randomMASTER(peakIndex,peakIndex) = 1;
  end

  S1_partial_random = updateTen(randomMASTER,RDC1,VECTORS);
  S2_partial_random = updateTen(randomMASTER,RDC2,VECTORS);
  
  percentile1  = NVR_COMP_TEN(S1_partial_random,S1_full);
  percentile2  = NVR_COMP_TEN(S2_partial_random, S2_full);
  percentiles1 = [percentiles1 percentile1];
  percentiles2 = [percentiles2 percentile2];
  
end

plot(percentiles1, percentiles2, '*');