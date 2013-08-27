function [MASTER,A]=updateMASTER(MASTER,A,ROWIN,COLIN);
%first, propagate contraints
recur=1;
rm=ones(1,size(A,1));
cm=ones(1,size(A,2));
ain=A;
while(recur==1)
   recur=0;
   for(i=1:size(A,1))
      x=find(A(i,:));
      if(length(x)==1)
         A(:,x)=0;
         A(i,x)=1;
         if(rm(i)==1)
            rm(i)=0;
            recur=1;
         end
      end
   end
   for(i=1:size(A,2))
      x=find(A(:,i));
      if(length(x)==1 & size(A,1)==size(A,2))
         A(x,:)=0;
         A(x,i)=1;
         if(cm(i)==1)
            cm(i)=0;
            recur=1;
         end
      end
   end
end
for(i=1:size(A,1))
   if(sum(A(i,:))==0)
      A(i,:)=ain(i,:);
   end
end

%next, make assignments
for(i=1:size(A,1))
   x = find(A(i,:));
   if(length(x)==1)
      MASTER(ROWIN(i),COLIN(x))=1;
   end
end
