% First order and total effect indices for a given
% model computed with Extended Fourier Amplitude
% Sensitivity Test (EFAST).
% Andrea Saltelli, Stefano Tarantola and Karen Chan.
% 1999. % "A quantitative model-independent method for global
% sensitivity analysis of model output". % Technometrics 41:39-56
clear;
close all;
%% INPUT
NR = 4; %: no. of search curves - RESAMPLING
k = 16 + 1; % # of input factors (parameters varied) + dummy parameter
NS = 513; % # of samples per search curve
wantedN=NS*k*NR; % wanted no. of sample points
 
% OUTPUT
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of k frequencies
% OMi : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group 
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output

MI = 4; %: maximum number of fourier coefficients
% that may be retained in calculating the partial
% variances without interferences between the
% assigned frequencies

%% PARAMETERS AND ODE SETTINGS (they are included in the following file)
Parameter_settings_EFAST;


% Computation of the frequency for the group
% of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1;
% if(NS*NR < 65)
%     fprintf(['Error: sample size must be >= ' ...
%     '65 per factor.\n']);
%     return;
% end


%% Pre-allocation of the output matrix Y
%% Y will save only the points of interest specified in
%% the vector time_points
Y(NS,length(time_points),length(y0),length(pmin),NR)=0;  % pre-allocation

% Loop over k parameters (input factors)
for i = 1:k % Loop over k parameters (input factors)
    try
        % Algorithm for selecting the set of frequencies
        OMci = SETFREQ(k, OMi / (2 * MI), i);   
        
        for L = 1:NR % Loop over the NR search curves
            % Setting the vector of frequencies OM
            cj = 1;
            for j = 1:k
                if j == i
                    OM(i) = OMi; % Parameter of interest
                else
                    OM(j) = OMci(cj); % Complementary group
                    cj = cj + 1;
                end
            end
            
            % Random phase shift and angular calculations
            FI = rand(1, k) * 2 * pi; % Random phase shift
            S_VEC = pi * (2 * (1:NS) - NS - 1) / NS;
            OM_VEC = OM(1:k);
            FI_MAT = FI(ones(NS, 1), 1:k)';
            ANGLE = OM_VEC' * S_VEC + FI_MAT;

            X(:, :, i, L) = 0.5 + asin(sin(ANGLE')) / pi; % Between 0 and 1
            
            % Transform distributions
            X(:, :, i, L) = parameterdist(X(:, :, i, L), pmax, pmin, 0, 1, NS, 'unif');

            % Run model evaluations
            for run_num = 1:NS
                try
                    % Perform a single model evaluation
                    [i, run_num, L] % Keeps track of [parameter, run, NR]
                    Y(run_num, :, :, i, L) = ODE_efast_wavespeed(X(:, :, i, L), run_num);
                catch ME
                    % Handle error in this specific NS run
                    warning(['Run ', num2str(run_num), ' failed for parameter ', num2str(i), ...
                             ' and search curve ', num2str(L), '. Error: ', ME.message]);
                    Y(run_num, :, :, i, L) = NaN; % Assign NaN to failed runs
                    continue; % Continue with the next run
                end
            end
        end
    catch ME
        % Handle error in the NR or parameter setup
        warning(['Parameter ', num2str(i), ' failed completely. Error: ', ME.message]);
        continue; % Continue with the next parameter
    end
end

save Model_efast.mat;
% CALCULATE Si AND STi for each resample (1,2,...,NR) [ranges]
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,1,1);
mybarplot(Si,Sti)
% Calculate Coeff. of Var. for Si and STi for Viral load (variable 4). See
% online Supplement A.5 for details.
% [CVsi CVsti]=CVmethod(Si, rangeSi,Sti,rangeSti,4)
% % T-test on Si and STi for Viral load (variable 4)
% s_HIV = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),efast_var,4,y_var_label,0.05)