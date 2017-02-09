

function  performance = perfscore(X_hat, X)

flag =1;

switch flag
    case 1
        err = X_hat(:) - X(:);        
        rse = sqrt(sum(err.^2)/sum(X(:).^2));        
        performance = rse;
end
