function X = Yi2X(Y, i)
if i == 3
    X = shiftdim(Trans_Faces(Y),i+1);
elseif i == 2
    X = shiftdim(Trans_Faces(Y),i);
else
    X = shiftdim(Trans_Faces(Y),i-1);
end
end
   
