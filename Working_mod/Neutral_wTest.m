%a = 0.12; % testing escape parameter - the higher a, the more rapid the escape
% 
% function detection = d(k, a)
% sensitivity = 0.78; % test sensitivity - from Harkins and Munson:low end of Commercial Nucleic Acid Hybridization test on 16s rRNA
% testRate = 0.0483; % https://www.cdc.gov/std/stats17/chlamydia.htm
% a = 0.12; % testing escape parameter - the higher a, the more rapid the escape
% dMax = sensitivity * testRate;
% detection = dMax * exp(-k * a);
% end

%%
% sensitivity = 0.78; % test sensitivity - from Harkins and Munson:low end of Commercial Nucleic Acid Hybridization test on 16s rRNA
% testRate = 0.0483; % https://www.cdc.gov/std/stats17/chlamydia.htm
% a = 0.12; % testing escape parameter - the higher a, the more rapid the escape
% dMax = sensitivity * testRate;

%%
function dydt = Neutral_wTest(t, y, myBeta, nu, mu, b, d, a)
    % y is the vector of S and I values
    % note that S = y(1), I_0 = y(2), and I_n = y(length(y))
    S = y(1);

    nPlus2 = length(y);  % total number of equations and dimensions of the system
    n = nPlus2 - 2;
    dydt = zeros(nPlus2, 1);

    % first the equation for S:
    dydt(1) = -myBeta * S * sum( y(2:nPlus2) )+ b;  % this is dS/dt 

    % next the equation for class 0 of the influenza
    dydt(2) = myBeta * S * y(2) - (nu + mu + d(k,a)) * y(2);  % this is dI0/dt

    % the last special case is the last strain, class "n" in mathematical
    % notation, i.e., dIn/dt:
    dydt(nPlus2) = myBeta * S * y(nPlus2) - (nu + d(k,a)) * y(nPlus2) + mu * y(nPlus2 - 1);% dIn/dt

    % now d/dt:
    for k = 1:(n-1) % iterating over all strains from strain 1 to strain n - 1
        index = k + 2;  % position in vector y and in vector dydt
        Ik = y(index);
        dydt(index) = myBeta * S * Ik - ((nu + mu + d(k,a)) * Ik) + (mu * y(index - 1));% dIk/dt
    end

end


% function detection = d(k, a)
%     detection = dMax * exp(-k * a);
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function detection = d(k, a)
sensitivity = 0.78; % test sensitivity - from Harkins and Munson:low end of Commercial Nucleic Acid Hybridization test on 16s rRNA
testRate = 0.0483; % https://www.cdc.gov/std/stats17/chlamydia.htm
a = 0.12; % testing escape parameter - the higher a, the more rapid the escape
dMax = sensitivity * testRate;
detection = dMax * exp(-k * a);
end

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function thetaVal = theta(S, Q)
% thetaVal = 1 - (S/Q);
% end




