function M = NVR_TOCSY2PROB(peakIDs,H,N,TYPES,SSTRUCT,NOES,ALLDISTS,ROWIN,COLIN);


[TH1,TH2,TN,TRN]=textread('InputFiles/TOCSY.m','%f %f %f %f');

M = zeros(length(H),length(TYPES));

for(i=1:size(M,1))
   
  %get the number of tocsy peaks
   rn=peakIDs(i);
   NUMTOCSYPEAKS=length(find(TRN==rn));
   TOCSYPEAKS=TH2(find(TRN==rn));
   
   for(j=1:length(TYPES))

     if(strcmp(TYPES(j),'ALA')==1)
         if(NUMTOCSYPEAKS>3 )
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'A',SSTRUCT(j),length(TYPES));
         end   
      elseif(strcmp(TYPES(j),'CYS')==1)
         if(NUMTOCSYPEAKS>5)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'C',SSTRUCT(j),length(TYPES));
         end   
      elseif(strcmp(TYPES(j),'ASP')==1)
         if(NUMTOCSYPEAKS>4)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'D',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'GLU')==1)
         if(NUMTOCSYPEAKS>6)
	   M(i,j)=0; 
         else
	   M(i,j)=getTOCSYProb(TOCSYPEAKS,'E',SSTRUCT(j),length(TYPES));
         end   
      elseif(strcmp(TYPES(j),'PHE')==1)
         
         if(NUMTOCSYPEAKS>9)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'F',SSTRUCT(j),length(TYPES));
         end   
      elseif(strcmp(TYPES(j),'GLY')==1)
         if(NUMTOCSYPEAKS>3)
	   M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'G',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'HIS')==1)
         if(NUMTOCSYPEAKS>8)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'H',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'ILE')==1)
         if(NUMTOCSYPEAKS>7)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'I',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'LYS')==1)
         if(NUMTOCSYPEAKS>11)
            M(i,j)=0; 
         else
	   M(i,j)=getTOCSYProb(TOCSYPEAKS,'K',SSTRUCT(j),length(TYPES));
         end      
      elseif(strcmp(TYPES(j),'LEU')==1)
         if(NUMTOCSYPEAKS>7)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'L',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'MET')==1)
         if(NUMTOCSYPEAKS>7)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'M',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'ASN')==1)
         if(NUMTOCSYPEAKS>6)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'N',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'GLN')==1)
         if(NUMTOCSYPEAKS>8)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'Q',SSTRUCT(j),length(TYPES));
         end
     elseif(strcmp(TYPES(j),'ARG')==1)
       if(NUMTOCSYPEAKS>13)%this was added by me (MSA) so that
                           %cases where numTOCSYPeaks > 13 does not
                           %cause an error.
	 M(i,j)=0; 
       else
	 M(i,j)=getTOCSYProb(TOCSYPEAKS,'R',SSTRUCT(j),length(TYPES));
       end
       elseif(strcmp(TYPES(j),'SER')==1)
         if(NUMTOCSYPEAKS>5)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'S',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'THR')==1)
         if(NUMTOCSYPEAKS>5)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'T',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'VAL')==1)
         if(NUMTOCSYPEAKS>5)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'V',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'TRP')==1)
         if(NUMTOCSYPEAKS>10)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'W',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'TYR')==1)
         if(NUMTOCSYPEAKS>9)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'Y',SSTRUCT(j),length(TYPES));
         end
      else
         PROBLEM = TYPES(j); 
      end
      
   end
   if (sum(M(i,:)) == 0)
     M(i,:) = 1;
   end
   M(i,:) = M(i,:)/sum(M(i,:));%re-normalize
end

%keyboard

persistent firstCall;

if (isempty(firstCall))
  for(i=1:min(size(M,1),size(M,2)))
    if (M(i,i) == 0)
      fprintf(1, 'in TOCSY prob. computation, before thresh, i = %d M(i,i) = 0\n',i);
      fprintf(1, 'will not print similar messages (here and below)');
      fprintf(1, ' for CAM protein.\n');
      keyboard
      firstCall = 0;
    end
  end
end



for(i=1:size(M,1))
   M(i,:)=M(i,:)/sum(M(i,:));
end

if (firstCall ~= 0)
  for(i=1:size(M,1))
    if (M(i,i) == 0)
      fprintf(1, 'in TOCSY prob. computation, after NOE2PROB, i = %d M(i,i) = 0\n',i);
      keyboard
    end
  end
end

function tp=getTOCSYProb(TOCSYPEAKS,AATYPE,SSTYPE,len)
nm=sprintf('/Users/student/Desktop/NVR_CA/trunk/HDB/%s%s.mat',char(AATYPE),char(SSTYPE));

SHIFTS=load(nm);
SHIFTS=SHIFTS.shifts;
MB=zeros(length(SHIFTS(1,:)))+999;
for i=1:length(TOCSYPEAKS)
   for j=1:size(SHIFTS,2)
      MB(i,j)=abs(TOCSYPEAKS(i)-SHIFTS(1,j));
   end
end

h=hungarian(MB');

tp=1;
for i=1:length(TOCSYPEAKS)
   pos=h(i);
   if(MB(i,pos)<100)
      if(SHIFTS(2,pos)==0)
         p=0;
	 fprintf(1, 'setting a variable to 0 in getTOCSYPROB.\n');
      else
         p=MB(i,pos)/SHIFTS(2,pos);
      end
      tp=tp*(2*(1-tcdf(p,len)));%i should probably have just used pdf here
   end
end
