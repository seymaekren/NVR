[R T RDC1 RDC2 X Y Z H_CS NH_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
    =textread('myinput.m','%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

%[extrN pnt1N ppm1 pnt2 ppm2 pnt3 ppm3 magn] = ...
%    textread('hnNoe.txt','%d %f %f %f %f %f %f %f');

%[ppm1 ppm3 ppm2] = textread('hnNoeManualPeak.txt.parsedNOE','%f %f %f');
unadjustedNOEs = load ('NOES.txt');
ppm1          = unadjustedNOEs(:,3);
ppm2          = unadjustedNOEs(:,4);
ppm3          = unadjustedNOEs(:,5);
residue1Index = unadjustedNOEs(:,1);
residue2Index = unadjustedNOEs(:,2);
%[resID1 resID2 ppm1 ppm2 ppm3 largeNumber anotherNumber] = textread('NOES.txt','%d %d %f %f %f %d %f');

h1_CS_Difference = []; h2_CS_Difference = []; n_CS_Difference = []; 

for noeIndex = 1:length(ppm1)
  hsqc1Index  = find(R == residue1Index(noeIndex));
  hsqc2Index  = find(R == residue2Index(noeIndex));
  hsqc_h1_CS  = H_CS(hsqc1Index);
  hsqc_h2_CS  = H_CS(hsqc2Index);
  hsqc_n_CS   = NH_CS(hsqc1Index);
  h1_CS_Difference(noeIndex) = ppm1(noeIndex)-hsqc_h1_CS;
  h2_CS_Difference(noeIndex) = ppm2(noeIndex)-hsqc_h2_CS;
  n_CS_Difference(noeIndex)  = ppm3(noeIndex)-hsqc_n_CS;
end

figure; plot(1:length(ppm1), h1_CS_Difference,'*');
title('H1 CS difference');
figure; plot(1:length(ppm1), h2_CS_Difference,'*');
title('H2 CS difference');
figure; plot(1:length(ppm1), n_CS_Difference,'*');
title('N CS difference');
keyboard





[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ignoredHSQCDATA] = loaddata('myinput.m');


H_EPS = 0.04;
N_EPS = 0.2;


load allPeaksAfterInitialAssignment.mat %why is this loaded?

%useOrigData = 1;
%[HSQCDATA, NOES] = loadUbqData(useOrigData);
%filename = sprintf('myinput.m')
%[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, HSQCDATA] = loaddata(filename);
case1 = 0; case2 = 0;
warning1 = 0; warning2 = 0;

for i = 1:length(ppm1)
  numFoundCloseHSQC_Peaks = 0;
  for j = 1:length(H_CS)
%    if (numFoundCloseHSQC_Peaks > 1)
%      break;
%    end
    if (abs(ppm1(i) - H_CS(j)) < H_EPS)
      for k = 1:length(H_CS)
	if (k == j)
	  continue;
	end

	if (abs(ppm2(i) - H_CS(k)) < H_EPS)
	  if (abs(ppm3(i) - NH_CS(j)) < N_EPS)

	    numFoundCloseHSQC_Peaks = numFoundCloseHSQC_Peaks + 1;
	    if (numFoundCloseHSQC_Peaks > 1)
	      break;
	      
% $$$ 	      caseNumber           = 1;
% $$$ 	      noeIndex             = i;
% $$$ 	      secondH_CS_1_peakNum = j;
% $$$ 	      secondH_CS_2_peakNum = k;
% $$$ 	      secondN_CS_peakNum   = j;
% $$$ 
% $$$ 	      secondH1_ppm_diff    = abs(ppm1(i) - H_CS(j));
% $$$ 	      secondH2_ppm_diff    = abs(ppm2(i) - H_CS(k));
% $$$ 	      secondN_ppm_diff     = abs(ppm3(i) - NH_CS(j));
% $$$ 	      
% $$$ 	      if (secondH1_ppm_diff < H1_ppm_diff) & (secondH2_ppm_diff ...
% $$$ 						      < H2_ppm_diff) ...
% $$$ 		    & (secondN_ppm_diff < N_ppm_diff)
% $$$ 		
% $$$ 		H_CS_1_peakNum = secondH_CS_1_peakNum;
% $$$ 		H_CS_2_peakNum = secondH_CS_2_peakNum;
% $$$ 		N_CS_peakNum   = secondN_CS_peakNum;
% $$$ 		fprintf(1, 'replaced the closest HSQC peak\n');
% $$$ 		numFoundCloseHSQC_Peaks = 1;
% $$$ 	      elseif (secondH1_ppm_diff > H1_ppm_diff) & (secondH2_ppm_diff ...
% $$$ 						      > H2_ppm_diff) ...
% $$$ 		    & (secondN_ppm_diff > N_ppm_diff)
% $$$ 		numFoundCloseHSQC_Peaks = 1;
% $$$ 	      else
% $$$ 		
% $$$ 	      end
	    end

	    caseNumber     = 1;
	    noeIndex       = i;
	    H_CS_1_peakNum = j;
	    H_CS_2_peakNum = k;
	    N_CS_peakNum   = j;
	    H1_ppm_diff    = abs(ppm1(i) - H_CS(j));
	    H2_ppm_diff    = abs(ppm2(i) - H_CS(k));
	    N_ppm_diff     = abs(ppm3(i) - NH_CS(j));
	    
	  end
	  
% $$$ 	  if (abs(ppm3(i) - NH_CS(k)) < N_EPS)
% $$$ 	    numFoundCloseHSQC_Peaks = numFoundCloseHSQC_Peaks + 1;
% $$$ 	    if (numFoundCloseHSQC_Peaks > 1)
% $$$ 	      break;
% $$$ 	    end
% $$$ 	    caseNumber     = 2;
% $$$ 	    noeIndex       = i;
% $$$ 	    H_CS_1_peakNum = j;
% $$$ 	    H_CS_2_peakNum = k;
% $$$ 	    N_CS_peakNum   = k;
% $$$ 	  end
	end
      end
    end
  end
  
% $$$   if (numFoundCloseHSQC_Peaks > 1)
% $$$     fprintf(1, 'found more than 1 close HSQC peaks\n');
% $$$     fprintf(1,'%d th noe is close to peak %d H_CS %d H_CS %d N_CS\n',noeIndex,R(H_CS_1_peakNum),R(H_CS_2_peakNum),R(N_CS_peakNum));
% $$$     fprintf(1,'NOE %d %d\n',R(H_CS_1_peakNum),R(H_CS_2_peakNum));
% $$$     fprintf(1, '%f %f %f\n',ppm1(noeIndex),ppm2(noeIndex),ppm3(noeIndex));
% $$$ %    fprintf(1, 'noe between residues %d and %d in the NOES.txt\n',resID1(noeIndex),resID2(noeIndex));
% $$$     fprintf(1, '%f %f %f\n',H_CS(H_CS_1_peakNum),H_CS(H_CS_2_peakNum),NH_CS(N_CS_peakNum));
% $$$     fprintf(1, 'case number is %d\n',caseNumber);
% $$$     
% $$$     
% $$$     fprintf(1,'%d th noe is close to peak %d H_CS %d H_CS %d N_CS\n',noeIndex,R(secondH_CS_1_peakNum),R(secondH_CS_2_peakNum),R(secondN_CS_peakNum));
% $$$     fprintf(1,'NOE %d %d\n',R(secondH_CS_1_peakNum),R(secondH_CS_2_peakNum));
% $$$     fprintf(1, '%f %f %f\n',ppm1(noeIndex),ppm2(noeIndex),ppm3(noeIndex));
% $$$ %    fprintf(1, 'noe between residues %d and %d in the NOES.txt\n',resID1(noeIndex),resID2(noeIndex));
% $$$     fprintf(1, '%f %f %f\n',H_CS(secondH_CS_1_peakNum),H_CS(secondH_CS_2_peakNum),NH_CS(secondN_CS_peakNum));
% $$$     fprintf(1, 'case number is %d\n',caseNumber);
% $$$     
% $$$     if (ALLDISTS(H_CS_1_peakNum,H_CS_2_peakNum) < NTH)
% $$$ 	    
% $$$       %	    if (j < size(NOES,1)) & (k < size(NOES,1))
% $$$       %	      if (NOES(j,k) == 1)
% $$$       
% $$$       fprintf(1, 'there is an NOE in the input file between these residues.\n');
% $$$     
% $$$       %	      end
% $$$       %	    end
% $$$     else
% $$$     
% $$$       fprintf(1, 'the distance is %f\n', ALLDISTS(H_CS_1_peakNum,H_CS_2_peakNum));
% $$$     end
    
    
    
% $$$     if (ALLDISTS(secondH_CS_1_peakNum,secondH_CS_2_peakNum) < NTH)
% $$$ 	    
% $$$       %	    if (j < size(NOES,1)) & (k < size(NOES,1))
% $$$       %	      if (NOES(j,k) == 1)
% $$$       
% $$$       fprintf(1, 'there is an NOE in the input file between these residues.\n');
% $$$     
% $$$       %	      end
% $$$       %	    end
% $$$     else
% $$$     
% $$$       fprintf(1, 'the distance is %f\n', ALLDISTS(secondH_CS_1_peakNum,secondH_CS_2_peakNum));
% $$$     end
% $$$     
% $$$     if (secondH1_ppm_diff < H1_ppm_diff) & (secondH2_ppm_diff ...
% $$$ 					    < H2_ppm_diff) ...
% $$$ 	  & (secondN_ppm_diff < N_ppm_diff)
% $$$ 
% $$$ 		H_CS_1_peakNum = secondH_CS_1_peakNum;
% $$$ 		H_CS_2_peakNum = secondH_CS_2_peakNum;
% $$$ 		N_CS_peakNum   = secondN_CS_peakNum;
% $$$ 		fprintf(1, 'replaced the closest HSQC peak\n');
% $$$ 		
% $$$        end
% $$$     
% $$$     
% $$$     keyboard
% $$$   end
  
  
  if (numFoundCloseHSQC_Peaks == 1)
    fprintf(1,'%d th noe is close to peak %d H_CS %d H_CS %d N_CS\n',noeIndex,R(H_CS_1_peakNum),R(H_CS_2_peakNum),R(N_CS_peakNum));
    fprintf(1,'NOE %d %d\n',R(H_CS_1_peakNum),R(H_CS_2_peakNum));
    fprintf(1, 'noe between residues %d and %d in the NOES.txt\n',residue1Index(noeIndex),residue2Index(noeIndex));
    fprintf(1, '%f %f %f\n',ppm1(noeIndex),ppm2(noeIndex),ppm3(noeIndex));
%    fprintf(1, 'noe between residues %d and %d in the NOES.txt\n',resID1(noeIndex),resID2(noeIndex));
    fprintf(1, '%f %f %f\n',H_CS(H_CS_1_peakNum),H_CS(H_CS_2_peakNum),NH_CS(N_CS_peakNum));
    correctPeak2Index = find(R == residue2Index(noeIndex));
    fprintf(1, '2nd Hydrogen correct H ppm is %f, diff from NOE ppm %f\n', ...
		H_CS(correctPeak2Index), H_CS(correctPeak2Index) - ppm2(noeIndex));
    fprintf(1, 'case number is %d\n',caseNumber);
    if (ALLDISTS(H_CS_1_peakNum,H_CS_2_peakNum) < NTH)
	    
      %	    if (j < size(NOES,1)) & (k < size(NOES,1))
      %	      if (NOES(j,k) == 1)
      
      fprintf(1, 'there is an NOE in the input file between these residues.\n');
      if (caseNumber == 1)
	case1 = case1 + 1;
      else
	assert (caseNumber == 2);
	case2 = case2 + 1;
      end
      %	      end
      %	    end
    else
      if (caseNumber == 1)
	warning1 = warning1 + 1;
      else
	assert (caseNumber == 2);
	warning2 = warning2 + 1;
      end
      fprintf(1, 'the distance is %f\n', ALLDISTS(H_CS_1_peakNum,H_CS_2_peakNum));
    end
    keyboard
  end
end

fprintf(1, 'case1 = %d, case2 = %d\n',case1,case2);
fprintf(1, 'warning1 = %d warning2 = %d\n', warning1, warning2);