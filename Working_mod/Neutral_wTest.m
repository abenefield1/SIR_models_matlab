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
function dydt = Neutral_wTest(time, y, myBeta, nu, mu, b, a)
    % y is the vector of S and I values
    % note that S = y(1), I_0 = y(2), and I_n = y(length(y))
    %dydt = @detection1;
    
        nPlus2 = length(y);  % total number of equations and dimensions of the system
        n = nPlus2 - 2;
        dydt = zeros(nPlus2, 1);
        S = y(1);
    
    for k = 1:(n-1)
    ind = k+2;
    d = 0.9 * exp(-ind * a);


        % first the equation for S:
        dydt(1) = -myBeta * S * sum( y(2:nPlus2) )+ b;  % this is dS/dt 

        % next the equation for class 0 of the influenza
        dydt(2) = myBeta * S * y(2) - (nu + mu + d(ind, a)) * y(2);  % this is dI0/dt
        %dydt(2) = myBeta * S * y(2) - (nu + mu + d) * y(2);

        % the last special case is the last strain, class "n" in mathematical
        % notation, i.e., dIn/dt:
        dydt(nPlus2) = myBeta * S * y(nPlus2) - (nu + d(ind, a)) * y(nPlus2) + mu * y(nPlus2 - 1);% dIn/dt
        %dydt(nPlus2) = myBeta * S * y(nPlus2) - (nu + d) * y(nPlus2) + mu * y(nPlus2 - 1);

        %%
    %     for k = 1:(n-1)
    %         index = k + 2;
    %         function detection1 = d(index, a)
    %         %detection1 = dMax * exp(-index * a);
    %         detection1 = 0.9 * exp(-index * a);
    %         end
    %     end

        %% now d/dt:
        for k = 1:(n-1) % iterating over all strains from strain 1 to strain n - 1
            index = k + 2;  % position in vector y and in vector dydt
            Ik = y(index);
            dydt(index) = myBeta * S * Ik - ((nu + mu + d(ind, a)) * Ik) + (mu * y(index - 1));% dIk/dt
            %dydt(index) = myBeta * S * Ik - ((nu + mu + d) * Ik) + (mu * y(index - 1));
        end
    end
    
end

%%

%     for k = 1:(n-1)
%         ind = k + 2;
%         d = 0.9 * exp(-ind * a);
%     end
        
        
%         function detection1 = d(index, a)
%         %detection1 = dMax * exp(-index * a);
%         detection1 = 0.9 * exp(-index * a);
%         end
    

% function detection = d(k, a)
%     %detection = dMax * exp(-k * a);
%     detection = 0.9 * exp(-k * a);
% end

% function detection = d(k, a)
%     detection = dMax * exp(-k * a);
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function detection = d(k, a)
% sensitivity = 0.78; % test sensitivity - from Harkins and Munson:low end of Commercial Nucleic Acid Hybridization test on 16s rRNA
% testRate = 0.0483; % https://www.cdc.gov/std/stats17/chlamydia.htm
% a = 0.12; % testing escape parameter - the higher a, the more rapid the escape
% dMax = sensitivity * testRate;
% detection = dMax * exp(-k * a);
% end

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function thetaVal = theta(S, Q)
% thetaVal = 1 - (S/Q);
% end




