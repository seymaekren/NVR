clear all;
%this one checks whether the constraints are satisfied by the
%soln. found.


%BinaryDistancesFilename = '/home2/apaydin/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/WithRDC_SecondRound/BinaryDistances.txt';
%NOE_List_Filename       = '/home2/apaydin/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/WithRDC_SecondRound/NOE_List.txt';


BinaryDistancesFilename = '/home2/apaydin/Workdir/OptimizationFiles/MBP/WithRDCs/BinaryDistances.txt';
NOE_List_Filename       = '/home2/apaydin/Workdir/OptimizationFiles/MBP/WithRDCs/NOE_List.txt';
scoreMatrixFilename     = '/home2/apaydin/Workdir/OptimizationFiles/MBP/WithRDCs/combinedScoringMatrix_3RDC_Matrices.txt';

%AssignmentMatrixFilename = '/home2/apaydin/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/TruncatingWithSpecialCoefficients/WithRDC_SecondRound/hSRI.txt';

%AssignmentMatrix = load (AssignmentMatrixFilename);
B                = load (BinaryDistancesFilename);

%AssignmentMatrix           = AssignmentMatrix(1:size(AssignmentMatrix,1),1:size(AssignmentMatrix,2)-3);
useOrigData                = 1;
NOE_List                   = textread(NOE_List_Filename);

%[HSQCDATA, NOES]           = readNMR_Data2(useOrigData);
%FullCorrectAssigmentMatrix = eye(size(AssignmentMatrix,1));
%AssignmentMatrix           = FullCorrectAssigmentMatrix;

%numPeaks         = size(AssignmentMatrix,1);
scoreMatrix      = load(scoreMatrixFilename);
numPeaks         = size(scoreMatrix,1);
for i = 1:numPeaks
  peak1    = i;
%  residue1 = find(AssignmentMatrix(peak1,:));
  residue1 = peak1;
if (isempty(residue1))
    continue;
  end
  entryIndexInNOE_List = find(NOE_List(:,2) == peak1);
%  peaksHavingAnNOE_WithThisOne = find(NOES(peak1, :));
  peaksHavingAnNOE_WithThisOne = NOE_List(entryIndexInNOE_List, 3:NOE_List(entryIndexInNOE_List,1)+3-1);
  for j = 1:length(peaksHavingAnNOE_WithThisOne)
    peak2    = peaksHavingAnNOE_WithThisOne(j);
%    residue2 = find(AssignmentMatrix(peak2, :));
    residue2 = peak2;
    if (isempty(residue2))
      continue;
    end
    if (B(residue1,residue2) == 0)
      fprintf(1, 'violation of constraint. Here the details:\n');
      fprintf(1, 'peak1 %d peak2 %d residue1 %d residue2 %d\n', ...
	      peak1, peak2, residue1, residue2);
      keyboard
    else
      fprintf(1, 'peak1 %d peak2 %d have an NOE and their assignment',peak1,peak2);
      fprintf(1, ' to %d and %d satisfies NOE constraints.\n',residue1, ...
	      residue2);
%      keyboard
    end
    
  end
end
