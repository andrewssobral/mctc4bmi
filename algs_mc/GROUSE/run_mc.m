function M_hat = run_mc(params)
  M = params.M;
  Idx = params.Idx;
  maxrank = 1;
  maxCycles = 100;
  step_size = 0.1;
  [numr,numc] = size(M);
  [Usg,Vsg,err_reg] = grouse(M,Idx,numr,numc,maxrank,step_size,maxCycles);
  M_hat = Usg*Vsg';
end
