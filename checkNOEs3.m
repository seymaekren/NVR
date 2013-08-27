%this one only takes a parsedPDB file containing H coordinates and
%NOE list containing absolute residue indices.

%load unambiguousNOEs.txt
%name = '1EZA.parsedPDB';

unambiguousNOEs = load ('1Y8B_hnhnParsedNOEs.txt');
name = '1D8C.protonatedWithAmber.parsedPDB';


[AtomType residueName resId  HX HY HZ] = textread(name,'%s %s %d %f %f %f');

ALLDISTS=zeros(length(HX));
for(i=1:size(ALLDISTS,1))
   for(j=1:size(ALLDISTS,1))
      ALLDISTS(i,j) = sqrt((HX(i)-HX(j)).^2+(HY(i)-HY(j)).^2+(HZ(i)-HZ(j)).^2);
   end   
end

NTH = 6;
for i = 1:size(unambiguousNOEs,1)
  noeResidue1 = unambiguousNOEs(i,1);
  noeResidue2 = unambiguousNOEs(i,2);
  relResId1   = find(resId == noeResidue1);
  relResId2   = find(resId == noeResidue2);
  if (ALLDISTS(relResId1,relResId2) > NTH)
    fprintf(1, 'the distance between residues %d and %d is %f\n',noeResidue1,noeResidue2,ALLDISTS(relResId1,relResId2));
  end
end
