function idx = getIndexIntoVectors(entryA, entryB, vectorA, vectorB)
idx = 0;
tol = 1e-6;
i = 1;

while(idx == 0 && i <= length(vectorA))
  if( abs( vectorA(i)-entryA ) <= tol && abs( vectorB(i)-entryB ) <= tol)
    idx = i;
  end
  i=i+1;
end
