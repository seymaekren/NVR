function [noeScore,numHN_NOES] = computeMB_NoeScore(HSQCDATA,ALLDISTS,NTH,MASTER);
%[h_ppm n_ppm h_ppm2 ] = load hn-hn-noes.txt
 
%[h_ppm n_ppm h_ppm2] =
%textread('hnNoeManualPeak.txt.parsedNOE','%f %f %f');

%could also skip passing HSQCDATA,ALLDISTS,NTH SB

persistent h_ppm n_ppm h_ppm2 numNOES hn_closePeaksCell
  
persistent secondProtonClosePeaksCell intensities distances

persistent smallDistanceThreshold largeDistanceThreshold

if (isempty(h_ppm))
 
  fprintf(1, 'first call, reading NOEs...\n');
  [h_ppm n_ppm h_ppm2 intensities] = textread('dnns.txt','%f %f %f %f');

  numNOES                    = length(h_ppm);

  hn_closePeaksCell          = cell(numNOES,1);
  secondProtonClosePeaksCell = cell(numNOES,1);
  distances                  = zeros(numNOES,1);
  
  for noeIndex = 1:length(h_ppm)
    hn_closePeaksCell{noeIndex}          = findCloseToHN(h_ppm(noeIndex),n_ppm(noeIndex),  HSQCDATA);
    secondProtonClosePeaksCell{noeIndex} = findCloseToH (h_ppm2(noeIndex), HSQCDATA);
    distances(noeIndex)                  = ((intensities(noeIndex) + 1.11)/1442.04)^(1/-4.5);
  end
  
%  smallDistanceThreshold = 0.25; 
%  largeDistanceThreshold = 1.0;

  smallDistanceThreshold = 1; 
  largeDistanceThreshold = 4;
  
end
  
%keyboard

noeScore   = 0; numHN_NOES = numNOES; 
distanceDiffs = []; 
structureDistances = []; correspondingNoeDistances = [];

for noeIndex = 1:length(h_ppm)
%  hn_closePeaks                               = findCloseToHN(h_ppm(i),n_ppm(i),  HSQCDATA);
%  secondProtonClosePeaks                      = findCloseToH(h_ppm2(i), HSQCDATA);
  [rowIndices,hn_correspondingLikelyResidues] = find(MASTER(hn_closePeaksCell{noeIndex},:));
  [rowIndex, secondProtonLikelyResidues]      = find(MASTER(secondProtonClosePeaksCell{noeIndex},:));
  foundVeryCloseNoe                           = 0;
  bestScoreAddition                           = 0;
  
  for relResidue1Index = 1:length(hn_correspondingLikelyResidues)
    
    residueIndex    = hn_correspondingLikelyResidues(relResidue1Index);

    foundDebugPeakIndex = 0;

    for debugRelPeakIndex = 1:length(hn_closePeaksCell{noeIndex})
      debugPeakIndex = hn_closePeaksCell{noeIndex}(debugRelPeakIndex);
      if (MASTER(debugPeakIndex,residueIndex) == 1)
	foundDebugPeakIndex = 1;
	break;
      end
    end
      
    for relResidue2Index = 1:length(secondProtonLikelyResidues)

      
      secondResidueIndex = ...
	  secondProtonLikelyResidues(relResidue2Index);
      
      foundSecondDebugPeakIndex = 0;
      
      for debugRelPeakIndex = 1:length(secondProtonClosePeaksCell{noeIndex})
	debugSecondPeakIndex = ...
	    secondProtonClosePeaksCell{noeIndex}(debugRelPeakIndex);
	if (MASTER(debugSecondPeakIndex,secondResidueIndex) == 1)
	  foundSecondDebugPeakIndex = 1;
	  break;
	end
      end

      assert ((foundDebugPeakIndex == 1) & (foundSecondDebugPeakIndex ...
					    == 1));
      
      %it seems this is a min operation that could be accelerated SB.
      if (ALLDISTS(hn_correspondingLikelyResidues(relResidue1Index),secondProtonLikelyResidues(relResidue2Index)) ...
	  > NTH)
	continue;
%     else
%	foundVeryCloseNoe = 1;%for going back to the old NOE
                              %scoring function.
%	break;
      end
       distance = ...
	   ALLDISTS(hn_correspondingLikelyResidues(relResidue1Index),secondProtonLikelyResidues(relResidue2Index));
       distanceDiff = abs(distance - distances(noeIndex));

       
       if (distanceDiff < smallDistanceThreshold)

% $$$ 	 fprintf(1, 'noe # %d, the distance between residue # %d and # %d is very close_to_noe_distance_computed_from_intensity.\n',noeIndex,hn_correspondingLikelyResidues(relResidue1Index),secondProtonLikelyResidues(relResidue2Index));
% $$$ 	 fprintf(1, 'corresponding HSQC peaks are: %d and %d\n', ...
% $$$ 		 debugPeakIndex, debugSecondPeakIndex);
% $$$ 	
% $$$ 	 if (bestScoreAddition == 0)
% $$$ 	   structureDistances = [structureDistances distance];
% $$$ 	   correspondingNoeDistances = [correspondingNoeDistances ...
% $$$ 		    distances(noeIndex)];
% $$$ 	 else
% $$$ 	   structureDistances(length(structureDistances)) =  distance;
% $$$ 	 end
% $$$ 	
%	if (foundCloseNoe == 0)
         bestScoreAddition = 1;
	 foundVeryCloseNoe = 1;
	 
%	else
%	  hn_closePeaks
%	  secondProtonClosePeaks
%	  figure
%	  keyboard
%	end
        break;
       elseif (distanceDiff < largeDistanceThreshold)
% $$$ 	 
% $$$  	 fprintf(1, 'noe # %d, the distance computed between residue # %d and # %d is somewhat close_to_the_distance_of_the_noe_computed_from_intensity.\n',noeIndex,hn_correspondingLikelyResidues(relResidue1Index),secondProtonLikelyResidues(relResidue2Index));
% $$$       	 fprintf(1, 'corresponding HSQC peaks are: %d and %d\n', ...
% $$$       		 debugPeakIndex, debugSecondPeakIndex);
	 
	 scoreAddition     = (largeDistanceThreshold - distanceDiff)/ ...
	     (largeDistanceThreshold - smallDistanceThreshold);
	 if (scoreAddition > bestScoreAddition)
% $$$       	   if (bestScoreAddition ~= 0)
% $$$       	     structureDistances(length(structureDistances)) = distance;
% $$$ % $$$ %	     keyboard
% $$$       	   else
% $$$       	     structureDistances                             = [structureDistances distance];
% $$$       	     correspondingNoeDistances                      = [correspondingNoeDistances distances(noeIndex)];
% $$$       	   end
	   bestScoreAddition                                = scoreAddition;
% $$$ 	   fprintf(1, 'best score addition changed to %f\n', bestScoreAddition);
						      %keyboard
	 end
% $$$       	 distanceDiffs = [distanceDiffs distanceDiff];
       elseif (distance < NTH)
% $$$ 	 fprintf(1, 'distance diff is too large. It is %f\n',distanceDiff);
% $$$ 	 fprintf(1, 'the noe # is %d\n',noeIndex);
% $$$ 	 fprintf(1, 'residues are %d and %d\n', hn_correspondingLikelyResidues(relResidue1Index),secondProtonLikelyResidues(relResidue2Index));
% $$$ 	 fprintf(1, 'corresponding HSQC peaks are: %d and %d\n', ...
% $$$ 		 debugPeakIndex, debugSecondPeakIndex);
% $$$ 	 fprintf(1, 'the distance in the structure is %f\n', distance );
% $$$ 	 fprintf(1, 'the distance in NOE intensity %f\n', distances(noeIndex));
% $$$ 	 distanceDiffs = [distanceDiffs distanceDiff];
%	 keyboard
       end
      end
      if (foundVeryCloseNoe)
	break;
      end
  end
  assert ((bestScoreAddition <= 1) & (bestScoreAddition >= 0));
  noeScore = noeScore + bestScoreAddition;
end
% $$$ figure; plot(distanceDiffs,'*');
% $$$ figure; plot(structureDistances, correspondingNoeDistances,'*');
% $$$ xlabel('structure distance')
% $$$ ylabel('corresponding NOE distance');