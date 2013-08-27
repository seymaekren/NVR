function [SAUPE, flag]= updateTen_CH(TABLE,NHRDC,CHRDC,NHVECS,CHVECS)
%flag is 1 if there is a problem in tensor computation. 0 otherwise.

SAUPE = ones(3,3);

M = zeros(1,5);
vec = 0;
ct = 1;
for(i=1:size(TABLE,1))
   x = find(TABLE(i,:));
   if(length(x)==1)
      if(NHRDC(i)>-999 && NHVECS(x,1)>-999)
         M(ct,1) = NHVECS(x,1).^2-NHVECS(x,3).^2;
         M(ct,2) = 2*(NHVECS(x,1)*NHVECS(x,2));
         M(ct,3) = 2*(NHVECS(x,1)*NHVECS(x,3));
         M(ct,4) = NHVECS(x,2).^2-NHVECS(x,3).^2;
         M(ct,5) = 2*(NHVECS(x,2)*NHVECS(x,3));
         vec(ct) = NHRDC(i);
         ct = ct+1;
      end
      if(CHRDC(i)>-999 && CHVECS(x,1)>-999)
         M(ct,1) = CHVECS(x,1).^2-CHVECS(x,3).^2; 
         M(ct,2) = 2*(CHVECS(x,1)*CHVECS(x,2)); 
         M(ct,3) = 2*(CHVECS(x,1)*CHVECS(x,3)); 
         M(ct,4) = CHVECS(x,2).^2-CHVECS(x,3).^2; 
         M(ct,5) = 2*(CHVECS(x,2)*CHVECS(x,3));
         vec(ct) = CHRDC(i);
         ct = ct+1;
      end
   end
end

flag = 0;
if(size(M,1)<5)
   size(M,1);
   flag=1;
   return
end
%invert the matrix
if(size(M,1)==5 && size(M,2)==5)
   M = inv(M);   
else
   M = pinv(M);
end
%do the least squares fitting
s = M*vec';
SAUPE = zeros(3,3);
SAUPE(1,1) = s(1);
SAUPE(2,2) = s(4);
SAUPE(3,3) = -1*(s(1)+s(4));
SAUPE(1,2) = s(2);
SAUPE(2,1) = s(2);
SAUPE(1,3) = s(3);
SAUPE(3,1) = s(3);
SAUPE(2,3) = s(5);
SAUPE(3,2) = s(5);
