function M_hat = run_mc(params)
  M = params.M;
  Idx = params.Idx;
  M_hat = low_rank_matrix_completion_new(M,sparse(Idx));
end
