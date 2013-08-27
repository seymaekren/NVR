function [ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN)
sct=1;
for(i=1:size(MASTER,1))
   x = find(MASTER(i,:));
   if(length(x)==0)
      ROWIN(sct)=i;
      sct=sct+1;
   end
end
sct=1;
for(i=1:size(MASTER,2))
   x = find(MASTER(:,i));
   if(length(x)==0)
      COLIN(sct)=i;
      sct=sct+1;
   end
end
