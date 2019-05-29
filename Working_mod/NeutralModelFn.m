function dydt = NeutralModelFn(t, y, myBeta, nu, mu, a)
% y is the vector of S and I values
% note that S = y(1), I_0 = y(2), and I_n = y(length(y))
S = y(1);

nPlus2 = length(y);  % total number of equations and dimensions of the system
n = nPlus2 - 2;
dydt = zeros(nPlus2, 1);

% first the equation for S:
dydt(1) = -myBeta * S * sum( y(2:nPlus2) );  % this is dS/dt 

% next the equation for class 0 of the influenza
dydt(2) = myBeta * S * y(2) - (nu + mu) * y(2);  % this is dI0/dt

% the last special case is the last strain, class "n" in mathematical
% notation, i.e., dIn/dt:
dydt(nPlus2) = myBeta * S * y(nPlus2) - nu * y(nPlus2) + mu * y(nPlus2 - 1);

% now dk/dt:
for k = 1:(n-1) % iterating over all strains from strain 1 to strain n - 1
    index = k + 2;  % position in vector y and in vector dydt
    Ik = y(index);
    dydt(index) = myBeta * S * Ik - ((nu + mu) * Ik) + (mu * y(index - 1));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function immunity = tau(j, k, a)
m = j + k;
immunity = exp(-a * m);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thetaVal = theta(S, Q)
thetaVal = 1 - (S/Q);
end




