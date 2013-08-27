function exists = checkIfFileExists(filename)

fid = fopen(filename,'r');

if (fid == -1)
  exists = 0;
else
  exists = 1;
end


