function sir = SIR(Ori,Est)

  s = Ori(:);

  t = Est(:);

  sir = 10*log10( (s'*s) / ((s-t)'*(s-t)) );


