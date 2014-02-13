function [M, differenceMatrixH, differenceMatrixN] = NVR_SHIFTX2PROB(TABLE,...
						  H,N,CA,TYPES,SSTRUCT,NOES,...
						  ALLDISTS, ...
						  NTH,ROWIN,COLIN, ...
						  SHIFTX_Filename,...
						  truncateProbabilities, ...
						  b_runningMBP, ...
						  b_runningEIN, ...
						  b_runningPoln, b_running1FQB)

%NVR_SHIFTX2PROB: This computes assignment probabilities based on the program SHIFTX, it is not meant
%             to be called by the user


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_SHIFTX2PROB.m
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

%    NVR_SHIFTX2PROB
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

FILTERS=load('/InputFiles/SHIFTX_FILTERS');
FILTERS=FILTERS.FILTERS;

[rn TY SS ha hn nf ca cb co]= textread(SHIFTX_Filename,'%f %s %s %f %f %f %f %f %f');
PRED = [rn  hn nf ca];

%get the scores for the two predictions
M=TABLE*0+1/size(TABLE,1);

differenceMatrixH = zeros(size(M,1),size(M,2));
differenceMatrixN = zeros(size(M,1),size(M,2));
differenceMatrixCA = zeros(size(M,1),size(M,2));

for i=1:size(TABLE,1)
   for j=1:length(COLIN)
        
     
     if (i == j)
       printDetails = 1;
     else
       printDetails = 0;
     end
      T = TY(COLIN(j));
      S = SS(COLIN(j));
      
      h = H(i)-PRED(COLIN(j),2);
      n = N(i)-PRED(COLIN(j),3);
      ca = CA(i)-PRED(COLIN(j),4);
      
      differenceMatrixH(i,j) = 1/(1+exp(abs(h)));
      differenceMatrixN(i,j) = 1/(1+exp(abs(n)));
      differenceMatrixCA(i,j) = 1/(1+exp(abs(ca)));
      
      
      if(strcmp(T,'A')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,1,truncateProbabilities, ...
			  printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'C')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,2,truncateProbabilities, ...
			  printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'D')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,3,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'E')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,4,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'F')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,5,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'G')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,6,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'H')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,7,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'I')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,8,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'K')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,9,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'L')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,10,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'M')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,11,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'N')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,12,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'Q')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,14,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'R')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,15,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'S')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,16,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'T')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,17,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'V')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,18,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'W')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,19,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      elseif(strcmp(T,'Y')==1)
         M(i,j)=getProb(h,n,ca, S,FILTERS,20,truncateProbabilities,printDetails,b_runningMBP, b_runningEIN, b_runningPoln, b_running1FQB);
      else
         PROBLEM = TYPES(COLIN(j)); 
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
for i=1:size(M,1)
   if(sum(M(i,:))==0)
      M(i,:)=1;
   end
   M(i,:)=M(i,:)/sum(M(i,:));
end



function p = getProb(h,n,ca, SSTYPE,FILTERS,TY, truncateProbabilities, ...
		       printDetails, b_runningMBP, b_runningEIN, ...
		       b_runningPoln, b_running1FQB)

persistent maxCoeff1 maxCoeff2 maxCoeff3 maxCoeff4 maxCoeff5 maxCoeff6;
persistent firstCall coeff1_ExtraMultiplier coeff2_ExtraMultiplier ...
    coeff3_ExtraMultiplier coeff4_ExtraMultiplier coeff5_ExtraMultiplier coeff6_ExtraMultiplier;

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
if (isempty(maxCoeff5))
  maxCoeff5 = 1;
end
if (isempty(maxCoeff6))
  maxCoeff6 = 1;
end



if (isempty(firstCall))

  if (b_runningMBP == 1)
    fprintf(1, 'will increase maxCoefficients for MBP.\n');
    firstCall              = 1;
    coeff1_ExtraMultiplier = 14.12;
    coeff2_ExtraMultiplier = 115.09;
    coeff3_ExtraMultiplier = 2.42;
    coeff4_ExtraMultiplier = 2.25;  
    coeff5_ExtraMultiplier = 1;
    coeff6_ExtraMultiplier = 1;
  
  elseif (b_runningEIN == 1)


    fprintf(1, 'will increase maxCoefficients for EIN.\n');
    firstCall              = 1;
    coeff1_ExtraMultiplier = 9.13;
    coeff2_ExtraMultiplier = 1.96;
    coeff3_ExtraMultiplier = 1;
    coeff4_ExtraMultiplier = 2.52;  
    coeff5_ExtraMultiplier = 1;
    coeff6_ExtraMultiplier = 1;
  
    
  elseif (b_running1FQB == 1)


    fprintf(1, 'will increase maxCoefficients for 1FQB.\n');
    firstCall              = 1;
    coeff1_ExtraMultiplier = 3.36;
    coeff2_ExtraMultiplier = 114.7;
    coeff3_ExtraMultiplier = 4.07;
    coeff4_ExtraMultiplier = 1.64;  
    coeff5_ExtraMultiplier = 1;
    coeff6_ExtraMultiplier = 1;
      
  elseif (b_runningPoln == 1)

    fprintf(1, 'will increase maxCoefficients for poln.\n');
    firstCall              = 1;
    coeff1_ExtraMultiplier = 1;
    coeff2_ExtraMultiplier = 1;
    coeff3_ExtraMultiplier = 3.90;
    coeff4_ExtraMultiplier = 1;  
    coeff5_ExtraMultiplier = 1;
    coeff6_ExtraMultiplier = 1;
    
  else
    
    fprintf(1, 'will set the maxCoefficients to 1 for proteins ');
    fprintf(1, 'other than EIN and MBP.\n');
    firstCall = 1;
    coeff1_ExtraMultiplier = 1;
    coeff2_ExtraMultiplier = 1;
    coeff3_ExtraMultiplier = 1;
    coeff4_ExtraMultiplier = 1;    
    coeff5_ExtraMultiplier = 1;
    coeff6_ExtraMultiplier = 1;
  end

end


if(strcmp(SSTYPE,'C')==1)
   
   xH = FILTERS(3,TY,1,1);
   xHS = FILTERS(3,TY,1,2);
   xN = FILTERS(3,TY,2,1);
   xNS = FILTERS(3,TY,2,2);
   xCA = FILTERS(3,TY,3,1);
   xCAS = FILTERS(3,TY,3,2);
   xHMX = FILTERS(3,TY,1,4);
   xHMN = FILTERS(3,TY,1,3);
   xNMX = FILTERS(3,TY,2,4);
   xNMN = FILTERS(3,TY,2,3);
   xCAMX = FILTERS(3,TY,3,4);
   xCAMN = FILTERS(3,TY,3,3);
   
elseif(strcmp(SSTYPE,'B')==1)
   
   xH = FILTERS(2,TY,1,1);
   xHS = FILTERS(2,TY,1,2);
   xN = FILTERS(2,TY,2,1);
   xNS = FILTERS(2,TY,2,2);
   xCA = FILTERS(3,TY,3,1);
   xCAS = FILTERS(3,TY,3,2);
   xHMX = FILTERS(2,TY,1,4);
   xHMN = FILTERS(2,TY,1,3);
   xNMX = FILTERS(2,TY,2,4);
   xNMN = FILTERS(2,TY,2,3);
   xCAMX = FILTERS(3,TY,3,4);
   xCAMN = FILTERS(3,TY,3,3);
   
else
   xH = FILTERS(1,TY,1,1);
   xHS = FILTERS(1,TY,1,2);
   xN = FILTERS(1,TY,2,1);
   xNS = FILTERS(1,TY,2,2);
   xCA = FILTERS(3,TY,3,1);
   xCAS = FILTERS(3,TY,3,2);
   xHMX = FILTERS(1,TY,1,4);
   xHMN = FILTERS(1,TY,1,3);
   xNMX = FILTERS(1,TY,2,4);
   xNMN = FILTERS(1,TY,2,3);
   xCAMX = FILTERS(3,TY,3,4);
   xCAMN = FILTERS(3,TY,3,3);
   
end

%fprintf(1, 'computing probability.\n');
%keyboard

%correction factor
% COMMENT OUT SEYMA
% xHMX = xHMX * coeff1_ExtraMultiplier;
% xHMN = xHMN * coeff2_ExtraMultiplier;
% xNMX = xNMX * coeff3_ExtraMultiplier;
% xNMN = xNMN * coeff4_ExtraMultiplier;
% xCAMX =xCAMX * coeff5_ExtraMultiplier; 
% xCAMN =xCAMN * coeff6_ExtraMultiplier; 

%xHMX=xHMX  *1.2*3.0;
%xHMN =xHMN *1.2*2.2;
%xNMX =xNMX *1.2*6.4;
%% $$$ %xNMX =xNMX *1.2*6.4*3.90;
%xNMN =xNMN *1.2*1.7;


%xHMX=xHMX*1.2;
%xHMN =xHMN *1.2;
%xNMX =xNMX *1.2;
%xNMN =xNMN *1.2;

%fprintf(1, 'changed the truncation back to original.\n');
%keyboard

if (truncateProbabilities)
    if(h>xH)
        if((h-xH)/xHS > xHMX)
      
            if (printDetails == 1)

                %fprintf(1, 'h: numStandardDeviationsAway : %f Limit : %f\n',(h-xH)/xHS,xHMX);
                fprintf(1, 'Ratio : %f\n',(h-xH)/xHS/xHMX);
                coeff1 = (h-xH)/xHS/xHMX;
                if (isempty(maxCoeff1))
                    maxCoeff1 = coeff1;
                elseif (coeff1 > maxCoeff1)
                    maxCoeff1 = coeff1;
                end
            end
            p=0;
            p1 = 0;%return
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
                    maxCoeff2 = coeff2;
                elseif (coeff2 > maxCoeff2)
                    maxCoeff2 = coeff2;
                end
            end
            p=0;
            p1=0;%return
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
                    maxCoeff3 = coeff3;
                elseif (coeff3 > maxCoeff3)
                    maxCoeff3 = coeff3;
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
                    maxCoeff4 = coeff4;
                elseif (coeff4 > maxCoeff4)
                    maxCoeff4 = coeff4;
                end
            end
            p=0;
            p2 = 0;
        else
            p2=normpdf(n,0,xNS);
        end
    end
    if(ca>xCA)
        if((ca-xCA)/xCAS > xCAMX)
            if (printDetails == 1)
                %fprintf(1, 'h: numStandardDeviationsAway : %f Limit : %f\n',(h-xH)/xHS,xHMX);
                fprintf(1, 'Ratio : %f\n',(ca-xCA)/xCAS/xCAMX);
                coeff5 = (-xCA)/xCAS/xCAMX;
                if (isempty(maxCoeff5))
                    maxCoeff5 = coeff5;
                elseif (coeff5 > maxCoeff5)
                    maxCoeff5 = coeff5;
                end
            end
            p3 = 0;
            %return
        else
            p3=normpdf(ca,0,xCAS);
        end
    else
        if((xCA-ca)/xCAS > xCAMN)
            if (printDetails == 1)
                %	fprintf(1, 'h: numStandardDeviationsAway : %f Limit :
                %	%f\n',-(h-xH)/xHS,xHMN);
                fprintf(1, 'Ratio : %f\n',-(ca-xCA)/xCAS/xCAMN);
                coeff6 = -(ca-xCA)/xCAS/xCAMN;
                if (isempty(maxCoeff6))
                    maxCoeff6 = coeff6;
                elseif (coeff6 > maxCoeff6)
                    maxCoeff6 = coeff6;
                end
            end
            p3 = 0;
            %return
        else
            p3=normpdf(ca,0,xCAS);
        end
    end
else
    p1=normpdf(h,xH,xHS);
    p2=normpdf(n,xN,xNS);
    p3=normpdf(ca,xCA,xCAS);
end

p = p1*p2*p3;

if ((p == 0) && (printDetails == 1))
    fid = fopen('maxCoefficients_SHIFTX.txt','w');
    fprintf(1, 'check out maxCoefficients_SHIFTX.txt\n');
    fprintf(fid, '%f %f %f %f %f %f\n', maxCoeff1,maxCoeff2,maxCoeff3,maxCoeff4,maxCoeff5,maxCoeff6);  fclose(fid);
    fclose(fid);
end

