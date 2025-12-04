%{
Input-Output HMM for mouse Set Shift

Here we have 4 states: explore, exploit left, exploit right, exploit light
The emission matrix (emat) depends on the input (the location of
the light--left or right)
%}
function [LLout,priorOut,transmatOut,alphaOut] = fitOreOitIOHMM_SetShift(data,verbose)

if nargin < 2
    verbose = true;
end


%%

% initialize
nChoices = 2; % number of options in the task (2 options for each environment configuration)
nStates = 4; % we need 4 states: explore, exploit right, exploit left, exploit light
nSeeds = 20; % number of times to reseed

% transition matrix:
transmat_logical = eye(nStates); % all-to-self OK
transmat_logical(1,:) = 1;
transmat_logical(:,1) = 1; % all to and from ORE

% emissions matrix depends on inputs (i.e. state of the light)
emat = zeros(4,2,2);

% hard coding emat
% light on the left
emat(:,:,1) = [.5, .5;... %explore
              1,  0;... %left
              0,  1;... %right
              1,  0];   %light
% light on the right
emat(:,:,2) = [.5, .5;... %explore
              1,  0;... %left
              0,  1;... %right
              0,  1];   %light

% initial state distribution:
initmat_logical = zeros(nStates,1); % start in ORE
initmat_logical(1) = 1;

% initmat_logical = ones(nStates,1);
 % only used if initSeeds are not provided

% tolerance and exiting stuff
MaxIter = 10^4;
MaxdLL = 10^-4;


% process the sequence into the format we want it in
[observations] = deal(cell(size(data,1),1));
        % would preallocate the state labels here
for seq = 1:size(data,1)

    obs = data{seq}; %needs to be {} for running all animals 
    choice = obs(1,:); % choices 1=left choice, 2=right choice
    light = obs(2,:); % data needs to be a cell array where first row is observed choices and 2nd row is inputs (i.e. the light location, 1=left, 2=right) and columns are trials
    
    observations{seq} = [choice' light'];
end

% For re-running one animal only:
%[observations] = deal(cell(size(data,1),1));
%obs = data;
%choice = obs(1,:);
%light = obs(2,:);
%observations = cell(1,1);
%observations{1,1} = [choice' light'];


%data.Properties.RowNames = {'choices','light'};
%observations = data;

fvalMin = -Inf;
%%

figure(); hold on;
colors = {'r','b','m','g','k','c'};

for seed = 1:nSeeds
    
    % randomly seed the starting parameters
    transmat = mk_stochastic(transmat_logical.*rand(nStates,nStates)); % doesn't depend on input

    % tie the parameters right off the bat
    transmat = tieParams(transmat); 
    
    initmat = mk_stochastic(initmat_logical.*rand(nStates,1));
    
    % initalize for EM
    LL = NaN; deltaLL = NaN; oldLL = 0; iter = 0;

    while (1)

        iter = iter + 1; LL = 0;

        % E step
        % set up sufficient statistics for the HMM
        exp_num_trans = zeros(nStates,nStates);%+10^-10;
        exp_num_visits1 = zeros(nStates,1);%+10^-10;
        exp_num_visitsT = zeros(nStates,1);%+10^-10;

        % step through ever sequence we observed
        for seq = 1:size(observations,1)

            % pull in the data for this round
            obs = observations{seq};
            choices = obs(:,1); light = obs(:,2);
            Time = length(choices);
            input = light + 1;
            choices = choices;

            % calculate the probability of each obs, given each state
            p_obs = obs_lik(choices,light,emat); 

            % now the joint probability of states and obs, given model
            %maximize = false; fwdOnly = false;
            transmat
            [~,~,gamma,currlik,xi_summed] = fwdback(initmat,transmat,p_obs,'fwd_only',0,'maximize',0); % function from HMM toolbox; removing input as the third item bc fwdback doesn't take input?
            % gamma is the posterior probabilitiy of the states
            %   p(yi | xt)
            % xi replicates the transition matrix
            %keyboard();
            % update the lik
            LL = LL+currlik;
            
            % now update the sufficient statistics
            input = input(1:end-1); % first, only keep the useful input
            
            % easy for the discrete nodes (don't need for mouse Set Shift)
            %for outcome = 1:2
                exp_num_trans(:,:) = sum(xi_summed(:,:),3); %check here, dimension
            %end % INPUT ONLY AFFECTS THE TRANSMAT - NOT THE OBSMAT
            % after learn_dhmm in bnet toolbox
            % https://github.com/the0s/DriverSim/blob/a1ff748ba21a7ac15c08501d980207b8b76ba02e/HMM_mat/learn_dhmm.m

            exp_num_visits1 = exp_num_visits1 + gamma(:,1);
            exp_num_visitsT = exp_num_visitsT + gamma(:,Time); % this is ci in the case of 1 mix

        end

        % now update the weights here
        deltaLL = oldLL-LL; % change in log likelihood
        oldLL = LL; % calculate log likelihood
        
        % give the console an update, if requested
        if verbose, fprintf(1, 'iteration %d, loglik = %f\n', iter, LL); end
        % stop now if we're not improving or we've iterated too much
        if (abs(deltaLL) < MaxdLL || iter > MaxIter); break; end

        % M step

        % first, we'll accomplish our parameter tying:
        
        % then update w/ the sufficient statistics:
        % square up our discrete nodes:
        startprob = normalise(exp_num_visits1);
%         endprob = normalise(exp_num_visitsT);
        
        initmat = startprob;

        exp_num_trans = tieParams(exp_num_trans); %check here
        
%         transmat = exp_num_trans; 

        transmat = mk_stochastic(exp_num_trans); % doesn't depend on input
        
        if round(iter/10) == (iter/10) || iter == 1
            initmat
            transmat
            emat
        end

        try
            plot(iter,LL,'.','MarkerSize',10,'Color',colors{seed});
            drawnow;
        end
    end
    
    if iter>MaxIter
        warning('fitVonMixHMM:MaxIter','This iteration did not converge.');
    elseif LL > fvalMin
        disp('model improvement!')
        
        fvalMin = LL; % set the new minimum
        
        % recalculate after final M-step (not true!)
        [loglik] = get_loglik(observations,initmat,transmat,emat)

        %each time model gets better, do viterbi here using fwdback
        %grabs the trial-by-trial state estimate 
        %maximize = 1; fwdOnly = 1;
        
        [alpha,beta,gamma,currlik,xi_summed] = fwdback(initmat,transmat,p_obs,'fwd_only',1,'maximize',1);
        
        %Convert data to cell arrays for data structure
        
        %alpha = num2cell(alpha);


        alphaOut = alpha;
        %betaOut = beta;
        %gammaOut = gamma;
        %xi_summedOut = xi_summed;
        transmatOut = transmat; 
        LLout = loglik;
        priorOut = initmat;


    end
end

if isinf(fvalMin)
    warning('fitVonMixHMM:MaxIter','Baum-Welch did not converge!');
    transmatOut = NaN(nStates,nStates); LLout = NaN;
    priotOut = initmat;
else
    disp('Baum-Welch converged just fine. Proceed with confidence!')
end

%% helpers

function mxOut = tieParams(mxIn)
% do the parameter tying where ever it's required: 
% (do these have to change?)
    
    mxOut = mxIn;
    
    for i = 1:size(mxIn,3)
    
        % first, the self-excitations are tied:
        mxOut(2:end,2:end,i) = eye(nStates-1) .* ((sum(diag(mxIn(:,:,i)))-mxIn(1,1,i))/2);

        % then the transitions into ORE
        mxOut(2:end,1,i) = deal(nanmean(mxIn(2:end,1,i)));

        % then the transitions from ORE
        mxOut(1,2:end,i) = deal(nanmean(mxIn(1,2:end,i)));
        
    end
end

function B = obs_lik(data, light, obsmat)
    % MK_DHMM_OBS_LIK  Make the observation likelihood vector for a discrete HMM.
    % B = mk_dhmm_obs_lik(data, obsmat, obsmat1)
    %
    % Inputs:
    % data(t) = y(t) = observation at time t
    % obsmat(i,o) = Pr(Y(t)=o | Q(t)=i)
    % obsmat1(i,o) = Pr(Y(1)=o | Q(1)=i). Defaults to obsmat if omitted.
    %
    % Output:
    % B(i,t) = Pr(y(t) | Q(t)=i)

    
    %keyboard();
    [Q O] = size(obsmat);
    T = length(data);
    B = zeros(Q,T);
    t = 1;
    for t=1:T
        try
      B(:,t) = obsmat(:, data(t), light(t));
        catch
            keyboard
        end
    end

end

function [loglik] = get_loglik(data,initmat,transmat,obsmat)
% LOG_LIK_DHMM Compute the log-likelihood of a dataset using a discrete HMM
% [loglik, errors] = log_lik_dhmm(data, prior, transmat, obsmat)
%
% data{m} or data(m,:) is the m'th sequence
% errors  is a list of the cases which received a loglik of -infinity
%keyboard();
[Q O] = size(obsmat);

ncases = length(data);
    
loglik = 0;
errors = [];

for m=1:ncases
    obs = data{m};
    choices = obs(:,1);
    rewards = obs(:,2); input = rewards + 1;
    T = length(obs);

    obslik = zeros(Q,T);

    for t=1:T
      obslik(:,t) = obsmat(:, choices(t));
    end
    
    fwdOnly = true; maximize = false;
    [~,~,gamma,ll,xi_summed] = fwdback(initmat, transmat, obslik, maximize, fwdOnly); 

    if ll==-inf
        errors = [errors m];
    end
    loglik = loglik + ll;
end

end

end