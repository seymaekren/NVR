function computeSvmClassification

[alpha xi1s xi2s xi3s xi4s] = textread('data_libsvm_scaled.txt','%f %s %s %s %s');
weight                      = 1.0e-03 * [-0.0445    0.0101   -0.1163   -0.0010];
rho                         = 1.00013;
fid = fopen('myClassification.txt','w');
for indexInFile = 1:length(alpha)
  combinedScore =  computeCombinedWeightedScore(indexInFile, alpha, xi1s,xi2s,xi3s,xi4s, ...
						weight, rho);
  fprintf(fid, 'alpha = %f combinedScore = %f\n', alpha(indexInFile), ...
	  combinedScore);
end
fclose(fid);
keyboard

function combinedScore =  computeCombinedWeightedScore(indexInFile, alpha, xi1s,xi2s,xi3s,xi4s, ...
						  weight, rho);
  [T,R] = strtok(xi1s(indexInFile),':');
  xi1   = strtok(R,':');
  xi1   = str2num(xi1{1});
  
  [T,R] = strtok(xi2s(indexInFile),':');
  xi2   = strtok(R,':');
  xi2   = str2num(xi2{1});
  
  
  [T,R] = strtok(xi3s(indexInFile),':');
  xi3   = strtok(R,':');
  xi3   = str2num(xi3{1});
  
  
  [T,R] = strtok(xi4s(indexInFile),':');
  xi4   = strtok(R,':');
  xi4   = str2num(xi4{1});
  
  
  xi    = [xi1 xi2 xi3 xi4];
  combinedScore = dot(xi,weight)-rho;



