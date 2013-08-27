useOrigData = 1;
 [HSQCDATA, NOES] = loadUbqData(useOrigData);
 [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
  ignoredHSQCDATA] = loaddata('myinput.m');

 H_CS   = HSQCDATA(:,2);
 N_CS  = HSQCDATA(:,3);

  unadjustedNOEs = load ('NOES.txt');
  order = load ('order.m');
 
[ppm1s ppm3s ppm2s] = textread('hnNoeManualPeak.txt.parsedNOE','%f %f %f');
numCorrespondingNOEs = zeros(size(HSQCDATA,1));
 
assert (size(ALLDISTS,1) == 72);
correct = 0; incorrect = 0; missingNOE = 0;absentNOE = 0; foundNOE ...
 	  = 0;
 
HN_EPS = 0.015; N_EPS = 0.35; H_EPS = 0.05; NOE_TH = 4; LARGER_NOE_TH ...
	 = 8;

for residue1Index = 1:size(HSQCDATA,1)
  for residue2Index = residue1Index+1:size(HSQCDATA,1)
%      noe1 = [H_CS(residue1Index) N_CS(residue1Index) H_CS(residue2Index)];
%      noe2 = ;
      for noeIndex = 1:length(ppm1s)
	if (abs(ppm1s(noeIndex)-H_CS(residue1Index))    < HN_EPS) & ...
	      (abs(ppm3s(noeIndex)-N_CS(residue1Index)) < N_EPS) & ...
	      (abs(ppm2s(noeIndex)-H_CS(residue2Index)) < H_EPS)
	  %fprintf(1, 'the close pair of peaks %d and %d have a corresponding NOE between them in the NOE file.\n',residue1Index,residue2Index);
	  numCorrespondingNOEs(residue1Index,residue2Index) = numCorrespondingNOEs(residue1Index,residue2Index) + 1;
	end
      
      	if (abs(ppm1s(noeIndex)-H_CS(residue2Index)) < HN_EPS) & ...
	      (abs(ppm3s(noeIndex)-N_CS(residue2Index))<N_EPS) & ...
	      (abs(ppm2s(noeIndex)-H_CS(residue1Index)) < H_EPS)
	  %fprintf(1, 'the close pair of peaks %d and %d have a corresponding NOE between them in the NOE file.\n',residue1Index,residue2Index);
	  numCorrespondingNOEs(residue1Index,residue2Index) = numCorrespondingNOEs(residue1Index,residue2Index) + 1;
	end
      end
      
      if (ALLDISTS(residue1Index,residue2Index) < NOE_TH)
	if (numCorrespondingNOEs(residue1Index,residue2Index) > 0)
	  foundNOE = foundNOE + 1;
	  fprintf(1, 'the pair of atoms are close and they have a corresponding NOE in the NOE file.\n');
	elseif (numCorrespondingNOEs(residue1Index,residue2Index) == 0)
	  fprintf(1, 'the pair of atoms are close and they DONT have a corresponding NOE in the NOE file.\n');
	  missingNOE= missingNOE + 1;
%	else
%	  fprintf(1, 'the pair of atoms are close and they have %d corresponding NOEs in the NOE file.\n',numCorrespondingNOEs(residue1Index,residue2Index));
%	  incorrect = incorrect + 1;
%	  keyboard
	end
      elseif (ALLDISTS(residue1Index,residue2Index) > LARGER_NOE_TH)
	if (numCorrespondingNOEs(residue1Index,residue2Index)  == 0)
	  absentNOE = absentNOE + 1;
	  fprintf(1, 'the pair of atoms are distant and they dont have a corresponding NOE in the NOE file.\n');
	else
	  fprintf(1, 'the pair of atoms are distant and they have %d corresponding NOEs in the NOE file.\n',numCorrespondingNOEs(residue1Index,residue2Index));
	  incorrect = incorrect + 1;
%%	  keyboard
	end
      end
      %    end
%      keyboard
  end
end

fprintf(1, 'total %d foundNOE %d absentNOE %d spurious NOE cases.\n',foundNOE, absentNOE, incorrect);
fprintf(1, 'total %d missing NOEs\n',missingNOE);

% $$$ maxH_ppm = 0;max_N_ppm_diff_case1 = 0; max_N_ppm_diff_case2 = 0;
% $$$ v_H1_diff = []; v_N_diff_case1 = [];v_N_diff_case2 = [];
% $$$ v_H2_diff = []; 
% $$$   
% $$$   
% $$$   
% $$$   
% $$$   
% $$$   
% $$$ for i = 1:size(unadjustedNOEs,1)
% $$$ %for i = 1:length(ppm1s)
% $$$   residue1Index = unadjustedNOEs(i,1);
% $$$   residue2Index = unadjustedNOEs(i,2);
% $$$   ppm1          = unadjustedNOEs(i,3);
% $$$   ppm2          = unadjustedNOEs(i,4);
% $$$   ppm3          = unadjustedNOEs(i,5);
% $$$   peak1Index    = find(order == residue1Index);
% $$$   peak2Index    = find(order == residue2Index);
% $$$   hsqc_ppm1     = H_CS(peak1Index);
% $$$   hsqc_ppm2     = H_CS(peak2Index);
% $$$   hsqc_ppm3_case1     = N_CS(peak1Index);
% $$$   hsqc_ppm3_case2     = N_CS(peak2Index);
% $$$   fprintf(1, 'NOE between %d and %d, ppm1 = %f ppm2 = %f, ppm = %f\n', residue1Index,residue2Index,ppm1,ppm2,ppm3);
% $$$   fprintf(1, 'H ppm of residue %d = %f, H ppm of residue %d = %f\n',residue1Index,hsqc_ppm1,residue2Index,hsqc_ppm2);
% $$$   fprintf(1, 'N ppm of residue %d = %f, %d = %f\n', residue1Index,hsqc_ppm3_case1,residue2Index,hsqc_ppm3_case2);
% $$$   fprintf(1, '%f %f %f %f\n',ppm1-hsqc_ppm1, ppm2-hsqc_ppm2, ppm3 ...
% $$$ 	   - hsqc_ppm3_case1, ppm3 - hsqc_ppm3_case2);
% $$$   if (abs(ppm1 - hsqc_ppm1) > maxH_ppm)
% $$$       maxH_ppm = abs(ppm1 - hsqc_ppm1);
% $$$   end
% $$$   if (abs(ppm2 - hsqc_ppm2) > maxH_ppm)
% $$$      maxH_ppm = abs(ppm2 - hsqc_ppm2);
% $$$   end
% $$$   
% $$$   v_H1_diff= [v_H1_diff ppm1 - hsqc_ppm1];
% $$$   v_H2_diff= [v_H2_diff    ppm2 - hsqc_ppm2];
% $$$   v_N_diff_case1= [v_N_diff_case1 ppm3 - hsqc_ppm3_case1] ;
% $$$   v_N_diff_case2= [v_N_diff_case2 ppm3 - hsqc_ppm3_case2] ;
% $$$     
% $$$   if (abs(ppm3 - hsqc_ppm3_case1) > max_N_ppm_diff_case1)
% $$$     max_N_ppm_diff_case1 = abs(ppm3 - hsqc_ppm3_case1);
% $$$   end
% $$$ 
% $$$    if (abs(ppm3 - hsqc_ppm3_case2) > max_N_ppm_diff_case2)
% $$$      max_N_ppm_diff_case2 = abs(ppm3 - hsqc_ppm3_case2);
% $$$    end
% $$$     
% $$$     
% $$$     
% $$$ % $$$ %   keyboard
% $$$ end
% $$$  
% $$$ fprintf(1, 'max H difference ppm is %f\n',maxH_ppm);
% $$$ fprintf(1, 'max N ppm difference is %f (case1), %f (case2)\n',max_N_ppm_diff_case1,max_N_ppm_diff_case2);
% $$$ 
% $$$ plot(v_H1_diff, '*');
% $$$ figure;
% $$$ plot(v_H2_diff, '*');
% $$$ figure;
% $$$ plot(v_N_diff_case1, 'r*')
% $$$ figure;
% $$$ % $$$ % plot(v_H_diff, '*');
% $$$ % $$$ % hold on
% $$$ plot(v_N_diff_case2,'r*');
% $$$  
% $$$  
% $$$  % for i = 1:size(NOES,1)
% $$$ %   peakIndices = find(NOES(i,:));
% $$$ %   for j = 1:length(peakIndices)
% $$$ %     if (peakIndices(j) > 
% $$$ %     assert (ALLDISTS(i, residueIndices(j)) < 6.5);
% $$$ %     fprintf(1, 'assert passed.\n');
% $$$ %   end
% $$$ % end
% $$$  