function [ASSIGNTABLE]=getUACMB(A,B, ASSIGNTABLE,NOES,ALLDISTS,NTH, ...
				  THR,rep,ROWIN,COLIN, MASTER);

%A cycle of making unambiguous assignments, followed by NOE pruning.
%makes the unambiguous assignments returned by FINDUACMB. 
%rep is the number of times this cycle will be executed.

tot=0;
for(k=1:rep)
   %find unambiguous
   for(i=1:size(A,1))
      A(i,:)=A(i,:)/sum(A(i,:));
      B(i,:)=B(i,:)/sum(B(i,:));
   end
   UA = FINDUACMB(A,B,ASSIGNTABLE,THR);
   if(sum(sum(UA))==tot)
      break;
   else
      tot = sum(sum(UA));
   end
   %make unambiguous assignments
   for(i=1:size(UA,1))
      x = find(UA(i,:));
      if(length(x)>0)
         ASSIGNTABLE(i,:)=0;
         ASSIGNTABLE(:,x)=0;
         ASSIGNTABLE(i,x)=1;
      end
   end
   %apply noes
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN);
end
