function [RESNUMS, ALLDISTS] = loadHSRI_Data(filename)

[RESNUMS resonanceAA_Name H_CS N_CS protonX protonY protonZ] ...
    = textread(filename,'%d %s %f %f %f %f %f');

ALLDISTS = zeros(length(RESNUMS),length(RESNUMS));
for i = 1:length(RESNUMS)
  proton1 = [protonX(i) protonY(i) protonZ(i)];
  for j = 1:length(RESNUMS)
    proton2 = [protonX(j) protonY(j) protonZ(j)];
    distance_sqr = 0;
    for k = 1:3
      distance_sqr = distance_sqr + (proton1(k)-proton2(k))* (proton1(k)-proton2(k));
    end
    ALLDISTS(i,j) = sqrt(distance_sqr);
  end
end