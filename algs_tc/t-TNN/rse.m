function r = rse(Xrec, Xref)
r = -10*log10((Xrec(:) - Xref(:))'*(Xrec(:) - Xref(:))/(Xref(:)'*Xref(:)));