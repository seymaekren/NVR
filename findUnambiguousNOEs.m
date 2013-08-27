function findUnambiguousNOEs

NOESY_Filename = 'Peak-CAM13C_Dieckmann-NNOESY-freq-Refined.dat';
HSQC_Filename  = 'Peak-CAM13C_Dieckmann-NHSQC-freq-Refined.dat';

%NOESY_Filename = 'ff2_3D_hhn.xpk';
%HSQC_Filename  = 'hsqcData_ff2.txt';
outFilename = 'unambiguousNOEs.txt';

%NOESY_Filename = 'ff2_n15_sansHeader.parsedXpk';
%HSQC_Filename  = 'hsqcData_ff2.txt';

H_EPS = 0.02; N_EPS = 0.2;

[dummyString NOE_H2_ppm NOE_N1_ppm NOE_H1_ppm intensity] = ...
    textread(NOESY_Filename, '%s %f %f %f %f');

[dummyString HSQC_N_ppm HSQC_HN_ppm intensity] = ...
    textread(HSQC_Filename, '%s %f %f %f');
%ORDER_Filename = 'order.m.ff2';
%[ORDER] = textread(ORDER_Filename, '%f');

%[NOE_H1_ppm NOE_H2_ppm NOE_N1_ppm] = textread(NOESY_Filename, '%f %f %f');

%[HSQC_HN_ppm HSQC_N_ppm ] = textread(HSQC_Filename, '%f %f');


numNOEs  = length(NOE_H2_ppm);

HSQCDATA = zeros(length(HSQC_N_ppm), 3);

for i = 1:length(HSQC_N_ppm)
  HSQCDATA(i,3) = HSQC_N_ppm(i);
  HSQCDATA(i,2) = HSQC_HN_ppm(i);
end

%the following is debugging code
%newHSQCDATA = zeros(2,3);
%%hsqcIndices = [119 126];
%hsqcIndices = [37 148];
%for i = 1:2
%  newHSQCDATA(i,3) = HSQC_N_ppm (hsqcIndices(i));
%  newHSQCDATA(i,2) = HSQC_HN_ppm(hsqcIndices(i));
%end

%HSQCDATA = newHSQCDATA;
%end debugging code

%why is HSQCDATA(i,1) empty?compatibility i guess.

NOEs = zeros(numNOEs, 3);
for i = 1:numNOEs
  NOEs(i, 1) = NOE_H2_ppm(i);
  NOEs(i, 2) = NOE_N1_ppm (i);
  NOEs(i, 3) = NOE_H1_ppm (i);
end


fid         = fopen(outFilename, 'w');
fprintf(1, 'check out unambiguousNOEs.txt\n');


%noeIndices = [25 218];

for noeIndex = 1:numNOEs
%for noeRelIndex = 1:length(noeIndices)
%  noeIndex = noeIndices(noeRelIndex);

  [noeIndex2 p1 p2] = findCorrespondingHSQC_Peak(noeIndex, NOE_H2_ppm, ...
						 NOE_N1_ppm, NOE_H1_ppm, ...
						 HSQCDATA, ...
						 H_EPS, N_EPS);
%  if (noeIndex == 37)
%    fprintf(1, 'noeIndex = %d\n', noeIndex);
%    keyboard
%  end


%    keyboard
    %in the above, what is p1? p2? noeIndex2?
    %why does this function pass noeIndex, and NOE_H2_ppm,
    %NOE_N1_ppm, NOEs etc. separately? 
    %i guess p1 and p2 are the HSQC indices. 
    %also noeIndex2 is the other NOE which is explained by these
    %HSQC peaks. i.e., p1(H),p1(N),p2(H) form an NOE, and 
    %p2(H),p2(N),p1(H) form another NOE and these are indexed
    %noeIndex and noeIndex2.
    
    
  if ((~isempty(p1)) & (noeIndex < noeIndex2))
    %why the check for noeIndex < noeIndex2?
    %probably to avoid double counting.
    
    
    [noeIndex3 p3 p4] = findCorrespondingHSQC_Peak(noeIndex2,NOE_H2_ppm, ...
						 NOE_N1_ppm, NOE_H1_ppm, ...
						   HSQCDATA, ...
						   H_EPS, N_EPS);
    if ((~isempty(p3)) & ((noeIndex3 == noeIndex) & (min(p1,p2) == min(p3,p4)) & ...
	  (max(p1,p2) == max(p3,p4))))

      %what does this check do, min(p1,p2) == min(p3,p4) and
      %max(p1,p2) == max(p3,p4)? Essentially whether (p1,p2) and
      %(p3,p4) are the same pair. whether the NOE # noeIndex is
      %explained by the HSQC pair p1 and p2.
      
%      fprintf(fid, '%d %d\n', ORDER(p1), ORDER(p2));
      fprintf(fid, '%d %d\n', p1, p2);
%      if ((min(p1,p2) == 17) & (max(p1,p2) == 62))
%	noeIndex
%	keyboard
%      end
      %      keyboard
    end
  end
end


%this function returns the HSQC indices corresponding to the
%NOE(noeIndex). It also returns the other NOE index noeIndex2 which
%is explained by this HSQC peak (returnedP1, returnedP2).
function [noeIndex2 returnedP1 returnedP2] = findCorrespondingHSQC_Peak(noeIndex, ...
						  NOE_H2_ppm, ...
						  NOE_N1_ppm, NOE_H1_ppm,...
						  HSQCDATA,H_EPS, ...
						  N_EPS)

returnedP1 = []; returnedP2 = []; noeIndex2 = -1;

hnClosePeaks = findCloseToHN(NOE_H1_ppm(noeIndex), NOE_N1_ppm(noeIndex), HSQCDATA,H_EPS, ...
			     N_EPS);
%  find the closest HSQC peak p1 TO HN, N.
hClosePeaks  = findCloseToH(NOE_H2_ppm(noeIndex), HSQCDATA,H_EPS, ...
			    N_EPS);
%  find the closest HSQC peaks p2 to H.

%   what if there are multiple peaks p1 and p2?

if ((length(hnClosePeaks) > 0) & (length(hClosePeaks) > 0))
  fprintf(1, '******************************************\n');
  fprintf(1, 'found close HSQC peaks. They will be tested ');
  fprintf(1, 'for symmetric NOEs and then printed.\n');
  fprintf(1, '******************************************\n');
  %    keyboard
end

[foundHSQCPeaks,noe2Indices] = testSymmetricNOEs(noeIndex,hnClosePeaks, hClosePeaks, HSQCDATA, NOE_H1_ppm, ...
				   NOE_H2_ppm, NOE_N1_ppm, ...
				   H_EPS, N_EPS);

if ((length(foundHSQCPeaks) == 2) & (length(noe2Indices) == 1))
  %fprintf(fid, '%d %d\n', foundHSQCPeaks(1), foundHSQCPeaks(2));
  returnedP1 = foundHSQCPeaks(1);
  returnedP2 = foundHSQCPeaks(2);
  noeIndex2 = noe2Indices(1);
elseif (length(foundHSQCPeaks) >= 2)
  fprintf(1, 'found more than a pair of hsqc peaks\n');
  %    keyboard
elseif (length(noe2Indices)> 1)
  fprintf(1, 'found more than 1 NOE.\n');
end


function [foundHSQCPeaks,noe2Indices] = testSymmetricNOEs(noeIndex,hnClosePeaks, hClosePeaks, ...
						  HSQCDATA, NOE_H1_ppm, ...
						  NOE_H2_ppm, NOE_N1_ppm, ...
						  H_EPS, N_EPS);
%for each pair of HSQC peaks, checking whether there is another NOE
%that is also explained by this pair of HSQC peaks.

foundHSQCPeaks = [];noe2Indices = [];
  


for hnClosePeaksIndex = 1:length(hnClosePeaks)
  for hClosePeaksIndex = 1:length(hClosePeaks)

    p1   = hnClosePeaks(hnClosePeaksIndex);
    p2   = hClosePeaks (hClosePeaksIndex);
    
    if (p1 == p2)
      %I assume that if there is no NOE that is between the same
      %proton. If such a case should happen, it should be
      %detectable by checking the proton1 and proton2 chemical
      %shift values of the NOE.
      assert ( NOE_H1_ppm(noeIndex)~= NOE_H2_ppm(noeIndex)); 
      %if they are equal this is a weird
      %NOE that should be removed manually.
      
      continue;
    end
    
    p2_N = HSQCDATA(hClosePeaks(hClosePeaksIndex),3);
    p2_H = HSQCDATA(hClosePeaks(hClosePeaksIndex),2);
    p1_H = HSQCDATA(hnClosePeaks(hnClosePeaksIndex),2);
    p1_N = HSQCDATA(hnClosePeaks(hnClosePeaksIndex),3);

    [tecn, closeNoeIndices] = thereExistsCloseNOEs(p2_N, p2_H, ...
						   p1_H, NOE_N1_ppm, ...
						   NOE_H1_ppm, ...
						   NOE_H2_ppm, ...
						   H_EPS, N_EPS);
    if (tecn == 1)
      fprintf(1, 'unambiguous NOE from HSQC peaks #%d to #%d\n', p1, p2);
      fprintf(1, 'HSQC peak #%d CS are %f %f\n', p1, p1_N, p1_H);
      fprintf(1, 'HSQC peak #%d CS are %f %f\n', p2, p2_N, p2_H);
      fprintf(1, 'NOE  peak #%d is %f %f %f\n', noeIndex, ...
	      NOE_H2_ppm(noeIndex), NOE_N1_ppm(noeIndex), NOE_H1_ppm(noeIndex));
      for foundSimilarNOEsIndex = 1:length(closeNoeIndices)
	noeIndex2 = closeNoeIndices(foundSimilarNOEsIndex);
	fprintf(1, 'NOE peak #%d is %f %f %f\n', noeIndex2, NOE_H2_ppm(noeIndex2), NOE_N1_ppm(noeIndex2), NOE_H1_ppm(noeIndex2));
	noe2Indices = [noe2Indices noeIndex2];
      end
      
      %	fprintf(fid, '%d %d\n', p1+aaOffset, p2+aaOffset);
      foundHSQCPeaks = [foundHSQCPeaks p1 p2];
    end
  end
end  


%this is checking whether there is a close NOE to the given HSQC peaks.
function [tecn, closeNoeIndices] = thereExistsCloseNOEs(p2_N, p2_H, ...
						  p1_H, NOE_N1_ppm, NOE_H1_ppm, NOE_H2_ppm, H_EPS, ...
						  N_EPS);

tecn            = 0; %There Exists Close NOEs
%H_EPS           = 0.03; N_EPS = 0.3;  %these thresholds are also
                                      %set in findCloseToHN and
                                      %findCloseToH and findUnambiguousNOEs.m
closeNoeIndices = [];

for i = 1:size(NOE_H1_ppm,1)
  HN_CS      = NOE_H1_ppm(i);
  N_CS       = NOE_N1_ppm(i);
  H_CS       = NOE_H2_ppm(i);
  HN_CS_DIFF = abs(HN_CS - p2_H);
  N_CS_DIFF  = abs(N_CS  - p2_N);
  H_CS_DIFF  = abs(H_CS  - p1_H);
  if ((HN_CS_DIFF < H_EPS) & (N_CS_DIFF < N_EPS) & (H_CS_DIFF < H_EPS))
    tecn            = 1;
    closeNoeIndices = [closeNoeIndices i];
  end
end