function out = thresh(M,x)

[a b] = size(M);

out = zeros(a,b);

for(i=1:a)
   for(j=1:b)
      if(M(i,j)>=x)
         out(i,j) = M(i,j);
      else
         out(i,j) = 0;
      end
      
   end
   
end
