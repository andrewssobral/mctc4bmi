%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: Carla Dee Martin (carlam@math.jmu.edu)
% DATE: June 2, 2010
%
% PROGRAM: tran.m
% PURPOSE: Input a tensor, and output the tensor transpose. 
%
% VARIABLES:
% A = Input tensor
% T = Tensor transpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = tran(A)

n=size(A);
k=length(n);
a=reshape(A,prod(n),1);

t = tensortransposerecur(a,n);

T = reshape(t,[n(2),n(1),n(3:end)]);


end

%%%%%%%%%%%%%%%%%%%%%%
function t = tensortransposerecur(a,n)
  
k=length(n);

if k<=2
  A=reshape(a,n(1),n(2))'; 
  a=reshape(A,prod(n),1);
  t=a;
 

elseif k>=3
  
  % flips order
  Acell=cell(n(k),1);
  Bcell=cell(n(k),1);

  % puts into cell array
  for j = 1:n(k)
     Acell{j} = a([(j-1)*prod(n(1:k-1))+1:j*prod(n(1:k-1))]);
  end

  %flip order of cells
  y = size(Acell);
  Bcell{1} = Acell{1};
  z = 2;
  for i = y(1):-1:2
    Bcell{z} = Acell{i};
    z = z+1;
  end

  % put back into long vector
  c=[];
  for m = 1:n(k)
     c = [c;Bcell{m}];
  end
  a=c;
  % end flip order
  
  for i=1:n(k)
    v=a((i-1)*prod(n(1:k-1))+1:i*prod(n(1:k-1)));
    a((i-1)*prod(n(1:k-1))+1:i*prod(n(1:k-1))) = tensortransposerecur(v,n(1:k-1));
  end
  
end

t=a;

end




