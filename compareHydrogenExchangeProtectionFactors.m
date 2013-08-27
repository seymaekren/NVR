name = '/mnt/sdb1/speech/users/apaydin/Workdir/InputFiles/myinput.m.1UBI';

[R T RDC1 RDC2 X Y Z H_CS NH_CS SS, HB,MOLMOL_HB HX,HY,HZ] ...
    =textread(name,'%f %s  %f %f %f %f %f %f %f %s %s %d %f %f %f');

load parsedHydrogenExchangeProtectionFactors.txt

count1 = 0; count2 = 0;
slowExchangerRates = [];
fastExchangerRates = [];

for i = 1:length(parsedHydrogenExchangeProtectionFactors)
  residueNumber      = parsedHydrogenExchangeProtectionFactors(i,1);
  exchangeRate       = parsedHydrogenExchangeProtectionFactors(i,2);
  indexInMyInputFile = find(R == residueNumber);
  slowExchangerAccordingToMyInputFile = ...
      MOLMOL_HB(indexInMyInputFile);
  fprintf(1, '%d %d %d\n', residueNumber, exchangeRate, slowExchangerAccordingToMyInputFile);
  if (slowExchangerAccordingToMyInputFile == 1)
%    figure(1)
%    plot(count1, exchangeRate,'*');

%    count1 = count1 + 1;
%    hold on
    slowExchangerRates = [slowExchangerRates exchangeRate];
  else
    fastExchangerRates = [fastExchangerRates exchangeRate];
%    figure(2)
%    plot(count2, exchangeRate,'*');
%    fprintf(1, '%d %d\n', count2, exchangeRate);
%    count2 = count2 + 1;
%    hold on
  end
end