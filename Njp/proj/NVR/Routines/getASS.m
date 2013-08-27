function M = getASS(X,Y,ASSIGNTABLE,THR,atype,top)
if(size(ASSIGNTABLE,1)==0)
   M=[];
   return
end
if(size(ASSIGNTABLE,1)-size(ASSIGNTABLE,2)~=0)
   F=ones(max(size(ASSIGNTABLE)));
   F(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=ASSIGNTABLE;
   ASSIGNTABLE=F;   
   for(i=1:size(ASSIGNTABLE,1))
      ASSIGNTABLE(i,:)=ASSIGNTABLE(i,:)/sum(ASSIGNTABLE(i,:));
   end
end
ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);
if(size(X,1)-size(X,2)~=0)
   F=ones(max(size(X)));
   F(1:size(X,1),1:size(X,2))=X;
   HM=F;   
   for(i=1:size(HM,1))
     sum_ithRow = sum(HM(i,:));
     if (sum_ithRow == 0)
       HM(i,:) = 1/size(HM,2);
     else
       HM(i,:)=HM(i,:)/sum_ithRow;
     end
   end
else
   HM=X;
end
M=HM*0;
h=1;
if(atype==1)
   h=hungarian(HM*-1);
elseif(atype==2)
   sz=top;
   for(i=1:size(HM,1))
      [x ind]= sort(HM(i,:));
      in2=find(x);
      ind=ind(in2);x=x(in2);
      if(range(x)>0)
         if(length(x)<sz+1)
            M(i,ind)=1;
         else
            M(i,ind(length(ind)-sz:length(ind)))=1;
         end
      else
         M(i,ind)=1;
      end
      h=[];   
   end
else
   h=simp(HM);
end
for(i=1:length(h))
   M(i,h(i)) = 1;
end

