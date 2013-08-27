function NOES = adjustNOEs(HSQCDATA,noes);

NOES = zeros(size(HSQCDATA,1));

for (i = 1:size(HSQCDATA,1))
  
  rn = HSQCDATA(i,1); 
  
  x  = find(noes(:,1)==rn);%see if there are any NOEs for this spin
                           %system
			   
  for (j=1:length(x))
    
    rn2 = noes(x(j),2);    %get the 2nd spin systems
    y   = find(HSQCDATA(:,1)==rn2);
    
    if (length(y)>0)
      NOES(i,y)=1; 
      NOES(y,i)=1;
    end
    
  end
end
