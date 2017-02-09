function T_hat = run_tc_blk(params)
  blksize = params.blksize;
  kblk = 0;
  sTk = size(params.T);
  T_hat = zeros(sTk); % slice3(Tk_hat); show_3dvideo(Tk_hat);
  disp(['frontal slice: ' num2str(sTk(1)) 'x' num2str(sTk(2))]);
  disp(['block size: ' num2str(blksize(1)) 'x' num2str(blksize(2))]);
  for i = 1:blksize(1):sTk(1)-(blksize(1)-1)
    for j = 1:blksize(2):sTk(2)-(blksize(2)-1)
      kblk = kblk + 1;
      %if(kblk == 1 || mod(kblk,10) == 0), 
        disp(['block: ' num2str(kblk) ' of ' num2str((sTk(1)*sTk(2))/(blksize(1)*blksize(2)))]);
      %end
      Tblk = params.T(i:i+blksize(1)-1, j:j+blksize(2)-1,:); % slice3(Tblk);
      Iblk = params.Idx(i:i+blksize(1)-1, j:j+blksize(2)-1,:);

      params_tc.T = Tblk;
      params_tc.Idx = Iblk;
      %Tblk_hat = zeros(size(Tblk));
      Tblk_hat = run_tc_full(params_tc); % slice3(Tblk_hat);

      T_hat(i:i+blksize(1)-1, j:j+blksize(2)-1,:) = Tblk_hat;
      %break;
    end
    %break;
  end
end
