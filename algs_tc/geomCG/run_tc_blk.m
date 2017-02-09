function T_hat = run_tc_blk(params)
  blksize = 8;
  kblk = 0;
  sTk = size(params.T);
  T_hat = zeros(sTk); % slice3(Tk_hat); show_3dvideo(Tk_hat);
  for i = 1:blksize:sTk(1)-(blksize-1)
    for j = 1:blksize:sTk(2)-(blksize-1)
      kblk = kblk + 1;
      Tblk = params.T(i:i+blksize-1, j:j+blksize-1,:); % slice3(Tblk);
      Iblk = params.Idx(i:i+blksize-1, j:j+blksize-1,:);

      params_tc.T = Tblk;
      params_tc.Idx = Iblk;
      %Tblk_hat = zeros(size(Tblk));
      Tblk_hat = run_tc_full(params_tc); % slice3(Tblk_hat);

      T_hat(i:i+blksize-1, j:j+blksize-1,:) = Tblk_hat;
      disp(kblk);
      %break;
    end
    %break;
  end
end
