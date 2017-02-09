function [Ztruth,Zest] = rescaleorder(Ztruth, Zest)


N = length(Ztruth);
nCom = size(Ztruth{1},2);


% standard the factors
for n=1:N
    Zest{n} = zscore(Zest{n});
    Ztruth{n} = zscore(Ztruth{n});
    
    Zest{n} = bsxfun(@minus, Zest{n}, min(Zest{n}));
    Zest{n} = bsxfun(@rdivide, Zest{n}, max(Zest{n})-min(Zest{n}));
    Zest{n} = 2*Zest{n}-1;
    
    Ztruth{n} = bsxfun(@minus, Ztruth{n}, min(Ztruth{n}));
    Ztruth{n} = bsxfun(@rdivide, Ztruth{n}, max(Ztruth{n})-min(Ztruth{n}));
    Ztruth{n} = 2*Ztruth{n}-1;
end

for n=1:N
    % rearrange the component order 
    COV = corr(Zest{n},Ztruth{n});
    [~,idx]= max(abs(COV));
    Zest{n} = Zest{n}(:,idx);
    % rearrange the component sign 
    idx = sub2ind(size(COV), idx', [1:nCom]');
    Zest{n} = Zest{n}.*repmat(sign(COV(idx))',size(Zest{n},1),1);
end



