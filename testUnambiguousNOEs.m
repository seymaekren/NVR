%[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
% ignoredHSQCDATA] = loaddata('myinput.m');

%noeList = load  ('ubqParsedUnambiguousNOEs.txt');
%noeList = load  ('parsedDnns.txt');

[RESNUMS, ALLDISTS] = loadHSRI_Data('combinedResonancesAndProtonCoordinates.txt');

%noeList = load  ('parsedBackboneNOEs_renumbered.txt'); for ff2 it seems
%noeList = load  ('parsedBackboneNOEs.txt'); 
noeList  = load  ('simulatedNOEs.txt'); 

numNOEs   = size(noeList,1);
distances = zeros(numNOEs,1);
numNOEs   = size(noeList,1);
assert (size(ALLDISTS,1) == length(RESNUMS));
noeMatrix = zeros(length(RESNUMS),length(RESNUMS));
%keyboard
numUniqueRestraints = 0;
for noeIndex = 1:size(noeList,1)
  aa1Index    = noeList(noeIndex,1);
  aa2Index    = noeList(noeIndex,2);
%  if ((aa1Index == 736) & (aa2Index == 725))
%    keyboard
%  end
  relAA1Index = find(RESNUMS == aa1Index);
  if (isempty(relAA1Index))
    continue;
  end
  relAA2Index = find(RESNUMS == aa2Index);
  if (isempty(relAA2Index))
    continue;
  end  
  distances(noeIndex) = ALLDISTS(relAA1Index,relAA2Index);

  if (noeMatrix(relAA1Index,relAA2Index) == 0)
    noeMatrix(relAA1Index,relAA2Index) = 1;
    noeMatrix(relAA2Index,relAA1Index) = 1;
    numUniqueRestraints = numUniqueRestraints+1;
  end
  
end
figure; plot(1:noeIndex,distances,'*');
fprintf(1, 'the max. distance is %f\n',max(distances));
fprintf(1, 'numUniqueRestraints = %f\n',numUniqueRestraints);