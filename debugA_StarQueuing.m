load insertQueueAssertFailure.mat
load bayesianMatrix.mat


% $$$ relPeakIndex = find(peakIndicess{nodeIndex} == 43);
% $$$ relResidueIndex = find(residueIndicess{nodeIndex} == 43);
% $$$ 
% $$$ assert (~isempty(relPeakIndex));
% $$$ assert (~isempty(relResidueIndex));
% $$$ 
% $$$ fprintf(1, 'Bayesian score for the assignment of 23 to 23 is %f',voters{nodeIndex}{1}(relPeakIndex,relResidueIndex));
% $$$ 
% $$$ fprintf(1, 'ASSIGNTABLE for 23 to 23 is %d\n', ...
% $$$ 	ASSIGNTABLES{nodeIndex}(relPeakIndex,relResidueIndex));
% $$$ 
% $$$ relResidueIndices = find(ASSIGNTABLES{nodeIndex}(relPeakIndex,:));
% $$$ fprintf(1, 'possible residues for peak 23 are: ');
% $$$ for relResidueIndicesCounter = 1:length(relResidueIndices)
% $$$   fprintf(1, '%d ', residueIndicess{nodeIndex}(relResidueIndices(relResidueIndicesCounter)));
% $$$   fprintf(1, 'the corresponding score is %f\n', voters{nodeIndex}{1}(relPeakIndex,relResidueIndices(relResidueIndicesCounter)));
% $$$ end
% $$$ 
% $$$ fprintf(1, 'print keyboard to continue.\n');
% $$$ keyboard

g = computeG            (newMaster, fullBayesianMatrix);
h = computeBetterH_Value(newAssigntable, newVoter, newNumVoters, ...
			 newUnassignedPeakIndices, newUnassignedResidueIndices);

fprintf(1, 'the new node has a combined g of %f ',g);
fprintf(1, 'and h of %f\n',h);

fprintf(1, 'for comparison, the function call had these parameters.\n');

fprintf(1, 'scoreSoFar = %f h_value = %f\n', newScoreSoFar, ...
	newH_Value);

g_in_function = newScoreSoFar; h_in_function = newH_Value;


fprintf(1, 'now the g and h of the parent node.\n');
g_parent = computeG(MASTERS{nodeIndex}, fullBayesianMatrix);
h_parent = computeBetterH_Value(ASSIGNTABLES{nodeIndex}, voters{nodeIndex}, ...
				numVoterss(nodeIndex), peakIndicess(nodeIndex),residueIndicess(nodeIndex));

fprintf(1, 'g_parent = %f, h_parent = %f\n',g_parent, h_parent);
fprintf(1, 'for comparison, these values were set to the following.\n');

g_parent_in_function = scoreSoFars(nodeIndex); 
h_parent_in_function = h_values(nodeIndex);
fprintf(1, 'g_parent_in_function = %f, h_parent_in_function = %f\n',scoreSoFars(nodeIndex),h_values(nodeIndex));



EPSILON = 1E-4;

fprintf(1, 'g_parent + assignment of new peak-residue pairs = ');
fprintf(1, '%f\n', -fullBayesianMatrix(67,67)-fullBayesianMatrix(68,68) - ...
	fullBayesianMatrix(69,69) + scoreSoFars(nodeIndex));

fprintf(1, 'difference between g_parent+score of 3 assignments - g_child=%f\n', g_parent - fullBayesianMatrix(67,67)- ...
	fullBayesianMatrix(68,68) - fullBayesianMatrix(69,69) - g);


if (g_parent_in_function+h_parent_in_function>g_in_function+h_in_function)
  fprintf(1, 'ERROR: Margin:%f\n', g_parent_in_function+h_parent_in_function-g_in_function-h_in_function);
end

fprintf(1, 'the score of assigning 67-67:%f\n', ...
	fullBayesianMatrix(67,67));
fprintf(1, 'the score of assigning 68-68:%f\n', ...
	fullBayesianMatrix(68,68));
fprintf(1, 'the score of assigning 69-69:%f\n',fullBayesianMatrix(69,69));



assert (newMaster(67,67) == 1);
assert (newMaster(68,68) == 1);
assert (newMaster(69,69) == 1);

assert (MASTERS{nodeIndex}(67,67) == 0);
assert (MASTERS{nodeIndex}(68,68) == 0);
assert (MASTERS{nodeIndex}(69,69) == 0);

numAssignmentsInParent = sum(sum(MASTERS{nodeIndex}));
numAssignmentsInChild  = sum(sum(newMaster));

assert (numAssignmentsInChild - numAssignmentsInParent == 3);

fprintf(1, 'again, g_child = %f\n', newScoreSoFar);

assert (abs(-fullBayesianMatrix(67,67)-fullBayesianMatrix(68,68) - ...
	    fullBayesianMatrix(69,69) + scoreSoFars(nodeIndex) - ...
	    newScoreSoFar) < EPSILON);




