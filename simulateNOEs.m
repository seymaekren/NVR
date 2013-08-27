function simulateNOEs;

%parsedPdbFilename       = '1dmb.10.model6.parsedPDB';
parsedPdbFilename        = '1SY9_model1.parsedPDB';
%parsedPdbFilename       = '1EZA.parsedPDB';
[protonName aaName residueIndex HX HY HZ] ...
    = textread(parsedPdbFilename,'%s %s %d %f %f %f');



ambiguousResidueIndices = load('ambiguousHSQC_ResidueIndices.txt');
order                   = load('order.m');




ALLDISTS=zeros(length(HX));
for(i=1:size(ALLDISTS,1))
   for(j=1:size(ALLDISTS,1))
      ALLDISTS(i,j)  = sqrt((HX(i)-HX(j)).^2+(HY(i)-HY(j)).^2+(HZ(i)-HZ(j)).^2);
   end   
end

threshold = 5;

[residue1,residue2] = find(ALLDISTS<threshold); 

numDistancesLessThanThreshold = (length(residue1) - size(ALLDISTS,1))/2;

fprintf(1, 'there are %d distances less than threshold in this protein.\n',numDistancesLessThanThreshold);

%keyboard


outputFilename = 'simulatedNOEs.txt';
fid = fopen(outputFilename, 'w');
fprintf(1, 'check out %s\n', outputFilename);

numResidues = size(ALLDISTS,1);



for i = 1:numResidues
  for j = i+1:numResidues
    if (ALLDISTS(i,j) < threshold)
      [relativeIndex1 relativeIndex2] = find(ambiguousResidueIndices == order(i)); 
      if (isempty(relativeIndex1))
	[relativeIndex3 relativeIndex4] = find(ambiguousResidueIndices == order(j));
	if (isempty(relativeIndex3))
	  fprintf(fid, '%d %d\n', order(i),order(j));
	end
      end
    end
  end
end

fprintf(1, 'printed all NOEs that are closer than %f A.\n', ...
	threshold);
fprintf(1, 'note. these are residue indices.\n');
fclose(fid);
keyboard

numSimulatedNOEs = 250;
NTH = 6;




for noeIndex = 1:numSimulatedNOEs

  foundA_PeakResiduePair   = 0;

  while (~foundA_PeakResiduePair)
    randomRelResidue1Index   = randi(length(HX),1,1);
    relResidueIndices        = find(ALLDISTS(randomRelResidue1Index, :) < NTH);
    if (~isempty(relResidueIndices))
      relResidue2Index       = relResidueIndices(randi(length(relResidueIndices),1,1));
      if (relResidue2Index ~= randomRelResidue1Index)
	foundA_PeakResiduePair = 1;
	fprintf(1, 'found a pair at a distance of %f ', ...
		ALLDISTS(randomRelResidue1Index, relResidue2Index));
      end
    end
  end
  
  fprintf(1,   'for residues %d and %d\n', residueIndex(randomRelResidue1Index), residueIndex(relResidue2Index));
  fprintf(fid, '%d %d\n', residueIndex(randomRelResidue1Index), residueIndex(relResidue2Index));
end

fclose(fid);