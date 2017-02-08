function rank_dec_adaptive(Xsq,n,coreNway,rank_min,rank_adj,X,X0,Y,Y0)
    % decrease the estimated rank
    sx{n} = svd(Xsq{n});
    dR = sx{n}(rank_min(n):end);
    drops = dR(1:end-1)./dR(2:end);
    [dmx,imx] = max(drops);
    rel_drp = (coreNway(n)-rank_min(n))*dmx/(sum(drops)-dmx);
    %If a large drop is found, adjust the rank estimate to imx
    if rel_drp>10
        coreNway(n) = imx+rank_min(n)-1;
        % set rank_adj(n) to 0, so only decrease the rank once
        rank_adj(n) = 0;
        [Qx,Rx] = qr(X{n},0);
        [Qy,Ry] = qr(Y{n}',0);
        [U,S,V] = svd(Rx*Ry');
        sigv = diag(S);
        X{n} = Qx*U(:,1:coreNway(n))*spdiags(sigv(1:coreNway(n)),...
            0,coreNway(n),coreNway(n)); 
        X0{n} = X{n};
        
        Yt{n} = Qy*V(:,1:coreNway(n));
        Y{n} = Yt{n}';
        Y0{n} = Y{n};
    end
end