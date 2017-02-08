function rank_inc_adaptive(M,Y,n,Yt,Y0,coreNway,rank_max,rank_inc,rank_inc_num,nstall,X0,solX,Nway,coNway)
    % increase the estimated rank
               
    [Q,R] = qr(Y{n}',0);
    for ii = 1:rank_inc(n)
        rdnx = randn(coNway(n),1);
        rdnx = rdnx-Q*(Q'*rdnx);
        rdnx = rdnx/norm(rdnx);
        Q = [Q,rdnx];
    end
    Y{n} = Q'; Yt{n} = Q; Y0{n} = Q';
    
    coreNway(n) = coreNway(n)+rank_inc(n);
    if coreNway(n) >= rank_max(n)
        rank_inc_num = rank_inc_num - 1;
    end

    if rank_inc_num == 0
        nstall = 0;
    end

    Mn = Unfold(M,Nway,n);
    X{n} = Mn*Y{n}';
    X0{n} = X{n};
    solX(n) = 0;
    
end