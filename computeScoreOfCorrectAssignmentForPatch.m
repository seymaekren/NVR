relPeakIndicesToBeAllowed = 61:65; %load 'patch6-10-allPeaks.mat'


%for i = 1:numVoters
%  voter{i} = voter{i}(relPeakIndicesToBeAllowed,:);
%end
%ASSIGNTABLE     = ASSIGNTABLE(relPeakIndicesToBeAllowed,:);
figure; plot(finalScores(1:totalNumAssignments),'*');
score = 0;
for i = 1:numVoters
  for j = 1:length(relPeakIndicesToBeAllowed)
    assert (voter{i}(relPeakIndicesToBeAllowed(j),relPeakIndicesToBeAllowed(j))~=0);
    score = score - log(voter{i}(relPeakIndicesToBeAllowed(j),relPeakIndicesToBeAllowed(j)));
  end
end


fprintf(1, 'score of opt assignment is %f\n',score);