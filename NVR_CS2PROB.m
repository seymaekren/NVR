function [M] = NVR_CS2PROB(TABLE,H,N,TYPES,SSTRUCT,NOES,ALLDISTS, ...
			   NTH,ROWIN,COLIN, truncateProbabilities);

FILTERS=load('NVR/CHEMSHIFTSTATS/FILTERS');FILTERS=FILTERS.FILTERS;

%the stats for the shifts were obtained from the wishart group webpage
%http://redpoll.pharmacy.ualberta.ca/RefDB/stat.html
M = TABLE*0;
[rows,cols]=find(TABLE==1);
if(length(rows)==0)
   rows=-1;
   cols=-1;
end
mH=mean(H);sH=std(H);
mN=mean(N);sN=std(N);

for(i=1:size(TABLE,1))
   h = H(i);
   n = N(i);
   p1 = normpdf(h,mH,sH);
   p2 = normpdf(n,mN,sN);
   prior=p1*p2;
   for(j=1:length(TYPES))
     
     if (i == j)
       printDetails = 1;
     else
       printDetails = 0;
     end
     
     
     
      if(strcmp(TYPES(j),'ALA')==1)
         M(i,j)=getProb(h,n, SSTRUCT(j),FILTERS,1,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'CYS')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,2,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'ASP')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,3,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'GLU')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,4,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'PHE')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,5,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'GLY')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,6,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'HIS')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,7,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'ILE')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,8,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'LYS')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,9,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'LEU')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,10,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'MET')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,11,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'ASN')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,12,truncateProbabilities,printDetails)/prior;%
      elseif(strcmp(TYPES(j),'GLN')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,14,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'ARG')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,15,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'SER')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,16,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'THR')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,17,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'VAL')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,18,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'TRP')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,19,truncateProbabilities,printDetails)/prior;
      elseif(strcmp(TYPES(j),'TYR')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,20,truncateProbabilities,printDetails)/prior;
      else
         PROBLEM = TYPES(j)   
      end
   end
   if (sum(M(i,:))==0)
      M(i,:)=1;
   end
   M(i,:) = M(i,:)/sum(M(i,:));%re-normalize
end
M = M .* TABLE;
%renornmalize
for(i=1:size(M,1))
   if(sum(M(i,:))==0)
      M(i,:)=1;
   end
   M(i,:)=M(i,:)/sum(M(i,:));
end


function p = getProb(h,n,SSTYPE,FILTERS,TY,truncateProbabilities, printDetails)%

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
%   coeff1_ExtraMultiplier = 3.49;
%   coeff2_ExtraMultiplier = 1.21;
%   coeff3_ExtraMultiplier = 1.20;
%   coeff4_ExtraMultiplier = 1;
%end

%if (isempty(firstCall))
%  fprintf(1, 'will increase maxCoefficients for EIN.\n');
%  firstCall = 1;
%  coeff1_ExtraMultiplier = 5.97;
%  coeff2_ExtraMultiplier = 1;
%  coeff3_ExtraMultiplier = 16.63;
%  coeff4_ExtraMultiplier = 1;
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

%xHMX=xHMX*1.1*7.5;

%xHMX=xHMX*1.1*7.5*coeff1_ExtraMultiplier;
%xHMN =xHMN *1.1*36.2*coeff2_ExtraMultiplier;
%xNMX =xNMX *1.1*8.2*coeff3_ExtraMultiplier;
%xNMN =xNMN *1.1*69.7*coeff4_ExtraMultiplier;

xHMX=xHMX*1.1*7.5;
xHMN =xHMN *1.1*36.2;
xNMX =xNMX *1.1*8.2;
xNMN =xNMN *1.1*69.7;

%xHMX=xHMX*1.1;
%xHMN =xHMN *1.1;
%xNMX =xNMX *1.1;
%xNMN =xNMN *1.1;




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
     p1=normpdf(h,xH,xHS);
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
      p1=normpdf(h,xH,xHS);
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
      p2=normpdf(n,xN,xNS);
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
      p2=0;
      %return
    else
      p2=normpdf(n,xN,xNS);
    end
  end
else
  p1=normpdf(h,xH,xHS);
  p2=normpdf(n,xN,xNS);
%  assert (p1 ~= 0);
%  assert (p2 ~= 0);
end
p = p1*p2;
if ((p == 0) & (printDetails == 1))
  fid = fopen('maxCoefficients_BMRB.txt','w');
  fprintf(1, 'check out maxCoefficients_BMRB.txt\n');
  keyboard
  fprintf(fid, '%f %f %f %f\n', maxCoeff1,maxCoeff2,maxCoeff3,maxCoeff4);
  fclose(fid);
end

function M = deleteUnlikely(V,M,S, th)

T = V*0;
for(i=1:size(S,1))
   for(j=1:size(V,2))
      if(M(i,j)==0)
         T(i,j)=max(max(V));
      elseif(V(i,j)==0)
         if(S(i,j)>0)
            T(i,j)=S(i,j);
         else
            T(i,j)=max(max(V));
         end
      else
         T(i,j)=V(i,j);
      end
   end
end


h = hungarian(T);
dum = S*0;
for(i=1:size(S,1))
   dum(i,h(i))=1;
end

dum = dum .*S; 
dum = uthresh(dum,min(nonzeros(dum))*1.26);

%%%%%%%%%%%%%%%%%%%DEBUGGING
bad = sum(diag(dum));
if(bad>0)
   for(i=1:size(S,1))
      if(dum(i,i)>0 )
         length(nonzeros(M(i,:)))
 %        pro2=[(mean(nonzeros(S))-S(i,h(i)))/std(nonzeros(S))]
      end
   end
 %  input('d');
end
%%%%%%%%%%%%%%%%%%%DEBUGGING

sin = sum(sum(M));
mk = xor(dum, ones(size(dum)));
M = M.*mk;

killing = [sin - sum(sum(M)) sum(sum(M))];

