function filterScores
load correctScores.txt
load incorrectScores.txt
load correctScoreIndices.txt
load incorrectScoreIndices.txt
correctnessLabel   =  1;
incorrectnessLabel = -1;
filterScoresAndPrintNewOnes(correctScores  , correctScoreIndices, ...
			    'filteredCorrectScores.txt', ...
			    'filteredCorrectScoreIndices.txt', correctnessLabel);
filterScoresAndPrintNewOnes(incorrectScores, incorrectScoreIndices, ...
			    'filteredIncorrectScores.txt', ...
			    'filteredIncorrectScoreIndices.txt', incorrectnessLabel);

function filterScoresAndPrintNewOnes(scores, scoreIndices, scoreFilename, ...
				     indicesFilename, correctnessLabel);
fid_scoreFile = fopen(scoreFilename, 'w');
fid_indexFile = fopen(indicesFilename, 'w');
fprintf(1, 'check out %s and %s\n', scoreFilename, indicesFilename);

[numScoreVectors,numScoringSources] = size(scores);
assert (size(scoreIndices,1) == numScoreVectors);
numDiscarded = 0;
for i = 1:numScoreVectors
  discardThisEntry = 0;
  for j = 1:numScoringSources
    entry = scores(i,j);
    if (entry == 1E+9)
      discardThisEntry = 1;
      numDiscarded = numDiscarded + 1;
      break;
    end
  end
  if (discardThisEntry)
    continue;
  end
  fprintf(fid_scoreFile, '%d ', correctnessLabel);
  for j = 1:numScoringSources   
    if (j == 5)
      %we are omitting the HDE
      %entry which is always 0 in the
      %filtered scores
      continue;
    end
    if (j < 5)
      fprintf(fid_scoreFile, '%d:%f ', j,scores(i,j));
    else
      fprintf(fid_scoreFile, '%d:%f ', j-1,scores(i,j));
    end
  end
  fprintf(fid_indexFile, '%d %d\n', scoreIndices(i,1), scoreIndices(i,2));
  fprintf(fid_scoreFile, '\n');
end
fclose(fid_scoreFile);
fclose(fid_indexFile);
fprintf(1,'remaining %d entries out of %d entries.\n',numScoreVectors-numDiscarded, ...
	numScoreVectors);