function [M,differenceMatrixH, differenceMatrixN,alphaHelixMatrix,betaStrandMatrix,coilMatrix] = ...
    NVR_SHIFTS2PROB(TABLE,H,N,TYPES,SSTRUCT,NOES,ALLDISTS, ...
		    NTH,ROWIN,COLIN, SHIFTS_Filename, SHIFTX_Filename, ...
		    truncateProbabilities);

%NVR_SHIFTS2PROB: This computes assignment probabilities based on the program SHIFTS, it is not meant
%             to be called by the user


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_SHIFTS2PROB.m
%//
%//  Version:		0.1
%//
%//  Description:	 This computes assignment probabilities based on BMRB statistics
%//
%// authors:
%//    initials    name            organization 					email
%//   ---------   --------------  ------------------------    ------------------------------
%//     CJL         Chris Langmead  Dartmouth College         langmead@dartmouth.edu
%//
%//
%// history:
%//     when        who     what
%//     --------    ----    ----------------------------------------------------------
%//     12/02/03    CJL 	 initial version for publication [Langmead et al, J Biomol NMR 2004]
%//
%////////////////////////////////////////////////////////////////////////////////////////////

%    NVR_SHIFTS2PROB
%    This library is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 2.1 of the License, or (at your option) any later version.

%    This library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.

%    You should have received a copy of the GNU Lesser General Public
%    License along with this library; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 		Contact Info:
%							Bruce Randall Donald
%							HB 6211
%							Dartmouth College
%							Hanover, NH 03755
%							brd@cs.dartmouth.edu

% 		If you use publish any results derived from the use of this program please cite:
%		"An Expectation/Maximization Nuclear Vector Replacement Algorithm for Automated NMR Resonance Assignments," 
%		C. J. Langmead and B. R. Donald, 
%		Journal of Biomolecular NMR, 2004 (in press)


%  Copyright (C) 2003  Christopher James Langmead and Bruce R. Donald
%
%  <signature of Bruce Donald>, 2 December 2003
%  Bruce Donald, Professor of Computer Science

FILTERS=load('~/NVR/trunk/CHEMSHIFTSTATS/SHIFTSFILTERS');FILTERS=FILTERS.FILTERS;

%[rn TY SS ha hn nf cb ca co]= textread('SHIFTX.m.backup','%f %s %s %f %f %f %f %f %f \n');
%[rn hn nf]= textread('SHIFTS.m.backup','%f %f %f \n');

[rn TY SS ha hn nf cb ca co]= textread(SHIFTX_Filename,'%f %s %s %f %f %f %f %f %f ');
[rn hn nf]                  = textread(SHIFTS_Filename,'%f %f %f ');

PRED                        = [rn  hn nf];

%size(PRED)

%get the scores for the two predictions
M=TABLE*0+1/size(TABLE,1);

differenceMatrixH = zeros(size(M,1),size(M,2));
differenceMatrixN = zeros(size(M,1),size(M,2));
alphaHelixMatrix  = zeros(size(M,1),size(M,2));
betaStrandMatrix  = zeros(size(M,1),size(M,2));
coilMatrix        = zeros(size(M,1),size(M,2));


% $$$ fid = fopen('shiftsCS_Differences.txt', 'w');
% $$$ fprintf(1, 'check out shiftsCS_Differences.txt\n');
% $$$ for(i=1:size(TABLE,1))
% $$$   j = i;
% $$$   h = H(i)-PRED(j,2);
% $$$   n = N(i)-PRED(j,3);
% $$$   fprintf(fid, '%f %f\n', h, n);
% $$$ end
% $$$ 
% $$$ for(i=1:size(TABLE,1))
% $$$   for(j=1:length(rn))
% $$$     h = H(i)-PRED(j,2);
% $$$     n = N(i)-PRED(j,3);
% $$$     fprintf(fid, '%f %f\n', h, n);
% $$$   end
% $$$ end
% $$$ fclose(fid);
for(i=1:size(TABLE,1))
   for(j=1:length(rn))
      
      T = TY(j);
      S = SS(j);
      
      h = H(i)-PRED(j,2);
      n = N(i)-PRED(j,3);

      differenceMatrixH(i,j) = 1/(1+exp(abs(h)));
      differenceMatrixN(i,j) = 1/(1+exp(abs(n)));
      
      if(strcmp(S,'C')==1)
	coilMatrix(i,j) = 1;
      elseif(strcmp(S,'B')==1)
	betaStrandMatrix(i,j) = 1;
      else
	alphaHelixMatrix(i,j) = 1;
      end

       if (i == j)
       printDetails = 1;
     else
       printDetails = 0;
     end
     
      if(strcmp(T,'A')==1)
         M(i,j)=getProb(h,n, S,FILTERS,1, truncateProbabilities,printDetails);
      elseif(strcmp(T,'C')==1)
         M(i,j)=getProb(h,n,S,FILTERS,2, truncateProbabilities,printDetails);
      elseif(strcmp(T,'D')==1)
         M(i,j)=getProb(h,n,S,FILTERS,3, truncateProbabilities,printDetails);
      elseif(strcmp(T,'E')==1)
         M(i,j)=getProb(h,n,S,FILTERS,4, truncateProbabilities,printDetails);
      elseif(strcmp(T,'F')==1)
         M(i,j)=getProb(h,n,S,FILTERS,5, truncateProbabilities,printDetails);
      elseif(strcmp(T,'G')==1)
         M(i,j)=getProb(h,n,S,FILTERS,6, truncateProbabilities,printDetails);
      elseif(strcmp(T,'H')==1)
         M(i,j)=getProb(h,n,S,FILTERS,7, truncateProbabilities,printDetails);
      elseif(strcmp(T,'I')==1)
         M(i,j)=getProb(h,n,S,FILTERS,8, truncateProbabilities,printDetails);
      elseif(strcmp(T,'K')==1)
         M(i,j)=getProb(h,n,S,FILTERS,9, truncateProbabilities,printDetails);
      elseif(strcmp(T,'L')==1)
         M(i,j)=getProb(h,n,S,FILTERS,10, truncateProbabilities,printDetails);
      elseif(strcmp(T,'M')==1)
         M(i,j)=getProb(h,n,S,FILTERS,11, truncateProbabilities,printDetails);
      elseif(strcmp(T,'N')==1)
         M(i,j)=getProb(h,n,S,FILTERS,12, truncateProbabilities,printDetails);
      elseif(strcmp(T,'Q')==1)
         M(i,j)=getProb(h,n,S,FILTERS,14, truncateProbabilities,printDetails);
      elseif(strcmp(T,'R')==1)
         M(i,j)=getProb(h,n,S,FILTERS,15, truncateProbabilities,printDetails);
      elseif(strcmp(T,'S')==1)
         M(i,j)=getProb(h,n,S,FILTERS,16, truncateProbabilities,printDetails);
      elseif(strcmp(T,'T')==1)
         M(i,j)=getProb(h,n,S,FILTERS,17, truncateProbabilities,printDetails);
      elseif(strcmp(T,'V')==1)
         M(i,j)=getProb(h,n,S,FILTERS,18, truncateProbabilities,printDetails);
      elseif(strcmp(T,'W')==1)
         M(i,j)=getProb(h,n,S,FILTERS,19, truncateProbabilities,printDetails);
      elseif(strcmp(T,'Y')==1)
         M(i,j)=getProb(h,n,S,FILTERS,20, truncateProbabilities,printDetails);
      else
         PROBLEM = TYPES(j)   
      end
      
    
   end
   
   if (sum(M(i,:)) == 0)
   %  keyboard
     M(i,:) = 1/size(M,2);

   else
     M(i,:) = M(i,:)/sum(M(i,:));%re-normalize
   end
   
end

M = M .* TABLE;
%renornmalize
for(i=1:size(M,1))
  if(sum(M(i,:))==0)
      M(i,:)=1;
   end
   %fprintf(1, 'shifts bpg construction removed all entries for peak #%d\n',i);
   %in that case don't you have to remove that row from further
   %consideration?
   %how can it be that a complete row is eliminated?
%  else
   M(i,:)=M(i,:)/sum(M(i,:));
%end
end

% $$$ nlast = sum(sum(TABLE));
% $$$ for(i=1:100)
% $$$    NP = NOE_PRUNE(TABLE(1:size(M,1),:),NOES,ALLDISTS,NTH,ROWIN,COLIN);
% $$$    TABLE(1:size(M,1),:)=and(TABLE(1:size(M,1),:),NP);
% $$$    if(sum(sum(TABLE)) == nlast)
% $$$       break;
% $$$    end
% $$$    nlast = sum(sum(TABLE));
% $$$ end
% $$$ 
% $$$ M = M.*TABLE;
% $$$ 
% $$$ %renornmalize
% $$$ for(i=1:size(M,1))
% $$$    if(sum(M(i,:))==0)
% $$$       M(i,:)=1;
% $$$    end
% $$$    %fprintf(1, 'shifts bpg construction removed all entries for peak #%d\n',i);
% $$$    %   else
% $$$    M(i,:)=M(i,:)/sum(M(i,:));
% $$$    %   end
% $$$ end



function p = getProb(h,n,SSTYPE,FILTERS,TY, truncateProbabilities, printDetails)%


persistent maxCoeff1 maxCoeff2 maxCoeff3 maxCoeff4;
persistent firstCall coeff1_ExtraMultiplier coeff2_ExtraMultiplier;
persistent coeff3_ExtraMultiplier coeff4_ExtraMultiplier;

if (isempty(maxCoeff1))
  maxCoeff1 = 1;
end
if (isempty(maxCoeff2))
  maxCoeff2 = 1;
end
if (isempty(maxCoeff3))
  maxCoeff3 = 1;
end
if (isempty(maxCoeff4))
  maxCoeff4 = 1;
end


  
%if (isempty(firstCall))
%   fprintf(1, 'will increase maxCoefficients for MBP.\n');
%   firstCall = 1;
%   coeff1_ExtraMultiplier = 15.93;
%   coeff2_ExtraMultiplier = 1.80;
%   coeff3_ExtraMultiplier = 1.76;
%   coeff4_ExtraMultiplier = 1;
%end

%if (isempty(firstCall))
%  fprintf(1, 'will increase maxCoefficients for EIN.\n');
%  firstCall = 1;
%  coeff1_ExtraMultiplier = 1.62;
%  coeff2_ExtraMultiplier = 1.22;
%  coeff3_ExtraMultiplier = 4.18;
%  coeff4_ExtraMultiplier = 1.55;
%end

if(strcmp(SSTYPE,'C')==1)
   
   xH = FILTERS(3,TY,1,1);
   xHS = FILTERS(3,TY,1,2);
   xN = FILTERS(3,TY,2,1);
   xNS = FILTERS(3,TY,2,2);
   xHMX = FILTERS(3,TY,1,4);
   xHMN = FILTERS(3,TY,1,3);
   xNMX = FILTERS(3,TY,2,4);
   xNMN = FILTERS(3,TY,2,3);
   
elseif(strcmp(SSTYPE,'B')==1)
   
   xH = FILTERS(2,TY,1,1);
   xHS = FILTERS(2,TY,1,2);
   xN = FILTERS(2,TY,2,1);
   xNS = FILTERS(2,TY,2,2);
   xHMX = FILTERS(2,TY,1,4);
   xHMN = FILTERS(2,TY,1,3);
   xNMX = FILTERS(2,TY,2,4);
   xNMN = FILTERS(2,TY,2,3);
   
else
   xH = FILTERS(1,TY,1,1);
   xHS = FILTERS(1,TY,1,2);
   xN = FILTERS(1,TY,2,1);
   xNS = FILTERS(1,TY,2,2);
   xHMX = FILTERS(1,TY,1,4);
   xHMN = FILTERS(1,TY,1,3);
   xNMX = FILTERS(1,TY,2,4);
   xNMN = FILTERS(1,TY,2,3);
   
end

%correction factor

%xHMX =xHMX*1.6*2.6*coeff1_ExtraMultiplier;
%xHMN =xHMN *1.6*4.5*coeff2_ExtraMultiplier;
%xNMX =xNMX *1.6*7.8*coeff3_ExtraMultiplier;
%xNMN =xNMN *1.6*6.6*coeff4_ExtraMultiplier;

xHMX =xHMX*1.6*2.6;
xHMN =xHMN *1.6*4.5;
xNMX =xNMX *1.6*7.8;
xNMN =xNMN *1.6*6.6;


%xHMX=xHMX*1.6;
%xHMN =xHMN *1.6;
%xNMX =xNMX *1.6;
%xNMN =xNMN *1.6;

if (truncateProbabilities)

 if(h>xH)
    if((h-xH)/xHS > xHMX)
      
      if (printDetails == 1)
       %fprintf(1, 'h: numStandardDeviationsAway : %f Limit : %f\n',(h-xH)/xHS,xHMX);
       fprintf(1, 'Ratio : %f\n',(h-xH)/xHS/xHMX);
       coeff1 = (h-xH)/xHS/xHMX;
       if (isempty(maxCoeff1))
	 maxCoeff1 = coeff1
       elseif (coeff1 > maxCoeff1)
	 maxCoeff1 = coeff1
       end
     end
      
      
       p=0;
       p1 = 0;
       %return
    else
       p1=normpdf(h,0,xHS);
    end
 else
    if((xH-h)/xHS > xHMN)
      
      if (printDetails == 1)
	%	fprintf(1, 'h: numStandardDeviationsAway : %f Limit :
	%	%f\n',-(h-xH)/xHS,xHMN);
	fprintf(1, 'Ratio : %f\n',-(h-xH)/xHS/xHMN);
      
	coeff2 = -(h-xH)/xHS/xHMN;
	if (isempty(maxCoeff2))
	  maxCoeff2 = coeff2
	elseif (coeff2 > maxCoeff2)
	  maxCoeff2 = coeff2
	end
      end
      
      
       p=0;
       p1 = 0;
       %return
    else
       p1=normpdf(h,0,xHS);
    end
 end
 if(n>xN)
    if((n-xN)/xNS > xNMX)
      if (printDetails == 1)
	%	fprintf(1, 'h: numStandardDeviationsAway : %f Limit : %f\n',(n-xN)/xNS,xNMX);
        fprintf(1, 'Ratio : %f\n',(n-xN)/xNS/xNMX);
	
	coeff3 = (n-xN)/xNS/xNMX;
	if (isempty(maxCoeff3))
	  maxCoeff3 = coeff3
	elseif (coeff3 > maxCoeff3)
	  maxCoeff3 = coeff3
	end
      end
       p=0;
       p2 = 0;
       %return
    else
       p2=normpdf(n,0,xNS);
    end
 else
    if((xN-n)/xNS > xNMN)
      if (printDetails == 1)
	%	fprintf(1, 'h: numStandardDeviationsAway : %f Limit : %f\n',-(n-xN)/xNS,xNMN);
        fprintf(1, 'Ratio : %f\n',-(n-xN)/xNS/xNMN);
	
	coeff4 = -(n-xN)/xNS/xNMN;
	if (isempty(maxCoeff4))
	  maxCoeff4 = coeff4
	elseif (coeff4 > maxCoeff4)
	  maxCoeff4 = coeff4
	end
      end
      p=0;
      p2 = 0;
      %return
    else
       p2=normpdf(n,0,xNS);
    end
 end
else
    p1=normpdf(h,0,xHS);
    p2=normpdf(n,0,xNS);
end
%input('d');

p = p1*p2;

if ((p == 0) & (printDetails == 1))
  fid = fopen('maxCoefficients_SHIFTS.txt','w');
  fprintf(1, 'check out maxCoefficients_SHIFTS.txt\n');
  fprintf(fid, '%f %f %f %f\n', maxCoeff1,maxCoeff2,maxCoeff3,maxCoeff4);
  fclose(fid);
end
