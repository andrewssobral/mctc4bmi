function eta = calcNewDirection( X_old, xi_old, eta_old, X, xi)
%CALCNEWDIRECTION calculates new CG direction
%   
%   ETA = CALCNEWDIRECTION( X_OLD, XI_OLD, ETA_OLD, X, XI )
%   computes the new CG direction using the Polak-Ribiere+ 
%   update formula.
%
%   See also calcGradient, transport
%

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    % transport old gradient to new X:
    xi_trans = transport( X_old, X, xi_old );

    % compute orthogonality between transported old and new gradient:
    ip_xitrans_xi = innerProduct( X, xi_trans, xi );
    ip_xi_xi      = innerProduct( X, xi, xi );

    theta = ip_xitrans_xi / ip_xi_xi;

    % compute conjugate direction
    if theta >= 0.1 
        eta = uminusFactorized(xi);
    else
        beta = max( 0, (ip_xi_xi - ip_xitrans_xi) / innerProduct( X_old, xi_old, xi_old ));
        % Fletcher-Reeves: beta = (ip_xi_xi) / innerProduct( X_old, xi_old, xi_old );

        % Calculate new new eta = -xi + beta*transport(eta_old)
        eta = addFactorized( uminusFactorized(xi), transport( X_old, X, eta_old ), beta);
    end
end

