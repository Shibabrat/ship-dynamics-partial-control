function func_get_saddle_eqpts


	% Parameters for Edith Terkol
    alpha0 = 0.73;
    omegaN = 0.62;
    omegaE = 0.527;
    lambda = 221.94;
    H = 4.94;

    b(1) = 0.0043;
    b(2) = 0.0225;
    c(1) = 0.384;
    c(2) = 0.1296;
    c(3) = 1.0368;
    c(4) = -4.059;
    c(5) = 2.4052;
    c = c./omegaN^2;
    b = b./omegaN^2;

    c(6) = alpha0*pi*(H/lambda);	% Non-dimensional wave height
    c(7) = omegaE/omegaN;			% Non-dimensional encounter frequency

    t = (2*pi)/c(7);
    x0 = [0.8,0];

    f = @(x)odeRHS(x,b,c,t)

	[saddleEqPts, rollVelVal, exitFlag] = fsolve(f, x0)

end
function y = odeRHS(x,bVec,cVec,t)

y = [x(2); ...
	-bVec(1)*x(2) - bVec(2)*(abs(x(2)).*x(2)) - ...
	cVec(1)*x(1) - cVec(2)*abs(x(1)).*((x(1).^3)) - ...
	cVec(3)*x(1).^3 - cVec(4)*abs(x(1)).*(x(1).^3) - ...
	cVec(5)*x(1).^5 + cVec(6)*sin(cVec(7)*t)];

end
