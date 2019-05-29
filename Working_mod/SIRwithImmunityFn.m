function dydt = SIRwithImmunityFn(t, y, myBeta, nu, mu, a, theta)
% y is the vector of Q and i values
% suppose
Q = y(1);  % fraction of pop not currently infected
I = y(2);   % fraction of pop infected summed over all strains

addedDims = 3;  % dimensions of the system in addion to the authors' "n"
nPlus = length(y);  % total number of equations and dimensions of the system
n = nPlus - addedDims;
dydt = zeros(nPlus, 1);  % output vector pre-allocated


% first handle special cases (three of them)

% dydt(1) = dQdt, equation 3.9


% dydt(2) = dIdt, equation 3.10


% dydt(3) = di0dt, equation 3.6


% dydt(nPlus) = dindt, equation 3.8
% now dik/dt (all other classes of infection from 1 to n-1):
% equation 3.7:
for k = 1:n % iterating over all strains from strain 1 to strain n - 1
    index = k + addedDims;  % position in vector y and in vector dydt
    ik = y(index); % kth strain's frequency; focal strain for this interation of loop
    ikMinus1 = y(index - 1);  % frequency of (k-1)th strain which can mutate to be kth strain
    
    SumOverL = 0; % summation variable needed for equation 3.7 and 3.8
    for l = 0:n
        il = y(l + addedDims); % ith strain
        SumOverL = SumOverL + (1 - theta * tau(l, a)) * il;
    end
    
    dydt(index) = myBeta * Q * ik * ( (1 - theta * tau(k, a)) - SumOverL ) + (mu * ikMinus1);  % equation 3.8 and most of 3.7
    if ( k < n )
        dydt(index) = dydt(index)  - (mu * ik);  % equation 3.7
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tauValue = tau(distance, a)
tauValue = exp( -a * distance );
end

