%this is to compute the weights from the SVM model file.
[alpha xi1s xi2s xi3s xi4s] = textread('half_data_libsvm.txt.scale.model.sansHeader','%f %s %s %s %s');

xis = [0 0 0 0];
for i = 1:length(alpha)

  [T,R] = strtok(xi1s(i),':');
  xi1   = strtok(R,':');
  xi1   = str2num(xi1{1});
  
  [T,R] = strtok(xi2s(i),':');
  xi2   = strtok(R,':');
  xi2   = str2num(xi2{1});
  
  
  [T,R] = strtok(xi3s(i),':');
  xi3   = strtok(R,':');
  xi3   = str2num(xi3{1});
  
  
  [T,R] = strtok(xi4s(i),':');
  xi4   = strtok(R,':');
  xi4   = str2num(xi4{1});
  
  
  xi    = [xi1 xi2 xi3 xi4];
  xis   = xis + alpha(i)* xi;

end
xis
