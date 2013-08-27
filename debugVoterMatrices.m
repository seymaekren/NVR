load allPeaksAfterInitialAssignments--For1UBI.mat
SSCP_1UBI                = SSCP;
unassignedResidueIndices = COLIN;
unassignedPeakIndices     = ROWIN;

load allPeaksNewAssignmentEnvironment-forBayesianScoring.mat
SSCP_allPeaks            = SSCP;


assert (twoMatricesAreClose(SSCP_allPeaks(unassignedRowIndices, ...
					  unassignedResidueIndices), ...
			    SSCP_1UBI, EPSILON) == 1);

fprintf(1, 'assert passed.\n');

keyboard
load ('allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat');

SSCP_original = SSCP;



EPSILON       = 1E-7;

assert (twoMatricesAreClose(SSCP_original, SSCP_1UBI, EPSILON) == 1);





load ('allPeaksAssignmentEnvironment-for1UBI-forBayesianScoring-Take2.mat');

SSCP_newBayesianScoring_allPeaks = SSCP;

if (twoMatricesAreClose(SSCP_newBayesianScoring_allPeaks, SSCP_allPeaks, ...
			EPSILON) == 1)
  fprintf(1, 'matrices from allPeaksAssignmentEnvironment-for1UBI-forBayesianScoring-Take2.mat');
  fprintf(1, ' and allPeaksNewAssignmentEnvironment-forBayesianScoring.mat\n');
  fprintf(1, ' are similar.\n');
end

assert (twoMatricesAreClose(SSCP_1UBI, ...
			    SSCP_newBayesianScoring_allPeaks(unassignedPeakIndices,unassignedResidueIndices),EPSILON) == 1);