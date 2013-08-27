function h = simp(M)
for(i=1:size(M,1))
   x = find(M(i,:));
   if(range(x)>0)
      [m,pos]=max(M(i,:));
      h(i)=pos;
   else
      r =randperm(length(x));
      h(i)=r(1);
   end
end
%throw away any that have more than one peak voting for it
for(i=1:size(M,2))
   x = find(M(:,i));
   if(length(x)>1)
      M(:,i)=0;
   end
end



