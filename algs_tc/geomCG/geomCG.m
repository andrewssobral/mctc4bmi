function [resX, err, conv] = geomCG( A_Omega, X, A_Test, opts )
	%geomCG Tensor completion by Riemannian optimization
	%	[X, ERR, CONV] = geomCG( A_OMEGA, X0, [], OPTS) performs tensor completion by 
	%	Riemannian optimization using the algorithm described in [1]. 
	%	The given data on the sampling set is specified by the sparse tensor A_OMEGA
    %   (sptensor class of Tensor Toolbox). The initial guess on the manifold is specified
    %   by the Tucker tensor (ttensor class of Tensor Toolbox) X0. Different aspects of the
    %   algorithm can be adjusted by the OPTS struct as described below.
    %   The algorithm outputs the final iterate X (ttensor) and the relative residual 
    %   (rel. error on the sampling set Omega) ERR at each iteration step. CONV is true if 
    %   the algorithm converged and false otherwise.
	%	
	%	[X, ERR, CONV] = geomCG( A_OMEGA, X0, A_TEST, OPTS) is the same as above but 
    %   the error is not only computed on the sampling set Omega but additionally on a test
    %   set Gamma given by the sptensor A_Test. 
    %   In this case, ERR is not just a vector of size (1 x iter), but a (2 x iter) matrix. 
    %   ERR(1,:) is the rel. error on the sampling set Omega (the residual)
    %   ERR(2,:) is the rel. error on the test set Gamma
    %   This is useful to detect cases of overfitting.
    %
    %   Specify the wanted options using the OPTS argument. All options have
    %   default values (denoted in brackets). Hence, you only have to specify those you want 
    %   to change.
    %   The OPTS struct is organized as follows: 
	%       OPTS.maxiter    Maximum number of iterations to perform     [400]
	%       OPTS.verbose    Turn on iteration information (true/false)  [true]
	%       OPTS.tol        Residual tolerance stopping criterion       [1e-6]
    %       OPTS.testtol    Error on test set stopping criterion        [opts.tol]
    %                       (only used if test set A_Test specified)
    %       OPTS.gradtol    Norm of gradient stopping criterion         [10*eps]
	%       OPTS.difftol    Difference to last iterate stopping crit.   [2*eps]
	%       OPTS.reltol     Relative tolerance stopping criterion       [0.05]
	%       OPTS.allX       Output not only final iterate, but all      [false]
    %                       intermediate approximations. Results are
    %                       in the cell array X of length number of
    %                       iterations, with X{end} the final approx.
	%       OPTS.reltol     Difference to last iterate stopping crit.   [2*eps]
	%       OPTS.backtrack  Perform up to 5 backtracking steps. Tests   [false]
    %                       showed that this is almost never necessary
    %                       Hence, turned off by default.
    %
    %
	%	References: 
	%	[1] D. Kressner, M. Steinlechner, B. Vandereycken:
	%		Low-Rank Tensor Completion by Riemannian Optimization.
	%		MATHICSE Technical Report 20.2013, June 2013. Submitted to BIT Numerical Mathematics.
	%

    %   GeomCG Tensor Completion. Copyright 2013 by
    %   Michael Steinlechner
    %   Questions and contact: michael.steinlechner@epfl.ch
    %   BSD 2-clause license, see LICENSE.txt

    if ~isfield( opts, 'maxiter');  opts.maxiter = 400;       end
    if ~isfield( opts, 'verbose');  opts.verbose = 1;         end
    if ~isfield( opts, 'tol');      opts.tol = 1e-6;          end
    if ~isfield( opts, 'gradtol');  opts.gradtol = 10*eps;    end
    if ~isfield( opts, 'difftol');  opts.difftol = 2*eps;     end
    if ~isfield( opts, 'testtol');  opts.testtol = opts.tol;  end
    if ~isfield( opts, 'reltol');   opts.reltol = 1e-5;       end
    if ~isfield( opts, 'allX');     opts.allX = 0;            end
	if ~isfield( opts, 'backtrack');opts.backtrack = false;   end

    if isempty( A_Test )
        calcOnTest = false;
    else
        calcOnTest = true;
    end

    norm_A_Omega = norm(A_Omega);

    % initial error:
    err_Omega(1) = sqrt( calcFunction( A_Omega, X) ) / norm_A_Omega;
    if calcOnTest
        norm_A_Test = norm(A_Test);
        err_test(1) = sqrt( calcFunction( A_Test, X) ) / norm_A_Test;
    end

    if opts.verbose
        disp(sprintf('Error on Omega: %g', err_Omega(1)))
    end

    if opts.allX == 1
        resX{1} = X;
    end

    % first step, simple steepest descent:
    % ====================================
    Deriv = calcGradient( A_Omega, X );
    xi = calcProjection( X, Deriv );
    eta = uminusFactorized(xi);

    size(eta.Y_tilde);
    
    t = calcInitial( Deriv, X, eta ); 
    % make the step.
    Xold = X;
    X = retraction( X, eta, t );
    value_old = calcFunction( A_Omega, X);

    err_Omega(2) = sqrt( value_old ) / norm_A_Omega;
    if opts.verbose
        disp(sprintf('t: %f,   Error on Omega: %g', t, err_Omega(2)))
    end

    if calcOnTest
        err_test(2) = sqrt( calcFunction( A_Test, X) ) / norm_A_Test;
    end

    if opts.allX == 1
        resX{2} = X;
    end

    % CONJUGATE GRADIENT SCHEME
    % =========================
    backtrackiter = 0;
    for iter = 3:opts.maxiter
        if opts.verbose
            disp(sprintf('Current iteration: %i',iter))
        end
        % 1. calculate new search direction
        % 1.1. Calc euclid. derivative of obj. function
        Deriv = calcGradient( A_Omega, X );

        % 1.2. Proj. into tangent space
        xi_old = xi;
        xi = calcProjection( X, Deriv );
        
        normgradient = sqrt(innerProduct( X, xi, xi ));
        if sqrt(innerProduct( X, xi, xi )) < opts.gradtol
            disp(sprintf('CONVERGED AFTER %i STEPS. Gradient is smaller than %0.3g', iter, opts.gradtol))
            conv = 1;
            break
        end

        % 1.3 calculate conjugate direction
        eta = calcNewDirection( Xold, xi_old, eta, X, xi);

        % 2. Initial stepsize guess
        t = calcInitial( Deriv, X, eta ); 

        % 3. Possibly do backtracking
        Xold = X;
		if ~opts.backtrack
        	X = retraction( X, eta, t);
		else
       	 	X_new = retraction( X, eta, t);
	        value_new = calcFunction( A_Omega, X_new );
	        inprod_xi_eta = innerProduct( X, xi, eta);
            
	        while ( value_old - value_new ) < -1e-4 * t * inprod_xi_eta;
	            disp('BACKTRACKING!')
	            t = 0.5 * t;
	            X_new = retraction( X, eta, t);
	            value_new = calcFunction( A_Omega, X_new );
	            if backtrackiter > 4
	                disp('BACKTRACKING STOPPED AFTER 5 ITERATIONS!')
	                break
	            end
	            backtrackiter = backtrackiter + 1;
	        end
	        backtrackiter = 0;
            
	        % make the step.
	        value_old = value_new;
        	X = X_new;
		end
        value_new = calcFunction( A_Omega, X );
        err_Omega(iter) = sqrt( value_new ) / norm_A_Omega;
        diff_iter = abs( err_Omega(iter-1) - err_Omega(iter));

        if calcOnTest
            err_test(iter) = sqrt( calcFunction( A_Test, X) ) / norm_A_Test;
            if err_test(iter)  < opts.testtol 
                disp(sprintf('CONVERGED AFTER %i STEPS. Err_Test smaller than %0.3g', iter, opts.testtol))
                conv = 1;
                break
            end
        end

        if opts.verbose
            disp(sprintf('Norm of curr. grad.: %g. Diff. to last iter: %g. Rel Diff: %g', normgradient, diff_iter, diff_iter/err_Omega(iter)))
        end

        if diff_iter  < opts.difftol
            disp(sprintf('CONVERGED AFTER %i STEPS. No More Progress in err_Omega ( < %0.3g ) ', iter, opts.difftol))
            conv = 1;
            break
        end

        if diff_iter/err_Omega(iter) < opts.reltol 
            disp(sprintf('CONVERGED AFTER %i STEPS. No More Progress, diff_iter/err_Omega < %0.3g', iter, opts.reltol))
            conv = 1;
            break
        end

        if err_Omega(iter)  < opts.tol 
            disp(sprintf('CONVERGED AFTER %i STEPS. Err_Omega smaller than %0.3g', iter, opts.tol))
            conv = 1;
            break
        end

        if opts.verbose
            disp(sprintf('t: %f,  Error on Omega: %g', t, err_Omega(iter)))
        end 

        if opts.allX == 1
            resX{iter} = X;
        end

    end

    if iter == opts.maxiter
        conv = 0;
    end

    if opts.allX == 0
        resX = X;
    end

    if calcOnTest
        err = [err_Omega; err_test];
    else
        err = err_Omega;
    end

end
