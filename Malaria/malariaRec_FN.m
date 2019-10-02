% Next steps:
% 1.) figure out if the 'd' function is working
    % if the code is working as expected, rethink the math or parameter values
% 2.) Figure out why code crashes for b*S (rather than b*N) about half the time - depending upon max time 
% 3.) Everything else - actually get at the selection pressure and pop-gen values

function dydt = malariaRec_FN(time, y, myBeta, nu, mu, a, N, b, gamma)
    % y is the vector of S and I values
    % note that S = y(1), I_0 = y(2), and I_n = y(length(y))
    % dydt = @detection1; Probably don't have to use this unless I rewrite d as a function
        
    %% Overall function parameters for Neutral_wTest
    nPlus2 = length(y);  % total number of equations and dimensions of the system
    n = nPlus2 - 2;
    dydt = zeros(nPlus2, 1);
    S = y(1);
    R = y(nPlus2);

    %% Parameters for 'd' or detection for loop
    sensitivity = 0.78; % test sensitivity - from Harkins and Munson:low end of Commercial Nucleic Acid Hybridization test on 16s rRNA
    testRate = 0.0483; % https://www.cdc.gov/std/stats17/chlamydia.htm
    %dMax = sensitivity * testRate; % Still need to adjust testRate unless plan to model in years
    dMax = 0.5;
   
    
        %% Main equations other than Ik
        % first the equation for S:
        dydt(1) = b*N - myBeta * S * sum( y(2:(nPlus2-1)) ) + gamma * R;  % this is dS/dt 

        % next the equation for class 0 of the influenza
        det = dMax * exp(-0 * a);
        dydt(2) = myBeta * S * y(2) - (nu + mu + det) * y(2); % this is dI0/dt

        % the last special case of infecteds is the last strain, class "n" in mathematical notation, i.e., dIn/dt:
        DET = dMax * exp(-n * a);
        dydt(nPlus2-1) = myBeta * S * y(nPlus2-1) - (nu + DET) * y(nPlus2-1) + mu * y(nPlus2 - 2);% dIn/dt
        
        % "R" class: dR/dt with waning immunity:
        %d=0.3;
        %dydt(nPlus2) = ((nu + d) * sum( y(2:(nPlus2-1)) )) - gamma * R;
        

        %% now d/dt:
        for k = 1:(n-2) % iterating over all strains from strain 1 to strain n - 1
            index = k + 2;  % position in vector y and in vector dydt
            d = dMax * exp(-k * a); % Sam added this here and changed "K" to "k"
            Ik = y(index);
            dydt(index) = myBeta * S * Ik - ((nu + mu + d) * Ik) + (mu * y(index - 1)); % dIk/dt
            %dydt(nPlus2) = ((nu) * sum( y(2:(nPlus2-1)).*d )) - gamma * R; % dR/dt
            dydt(nPlus2) = (sum( y(2:(nPlus2-1)).*(d + nu) )) - gamma * R; % dR/dt

        end
   
    
end