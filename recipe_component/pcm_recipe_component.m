function varargout = pcm_recipe_component(varargin)
%%
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% Atsushi Yokoi (2017 Feb)

%-------------------------------------------------------------------------%
% Load data: S1,M1,PMd,IPS (all left hemisphere)
%-------------------------------------------------------------------------%
load('data_component.mat');

for r=1:numel(Y)
    for s=1:numel(Y{r})
        Ghat{r}(:,:,s) = pcm_estGCrossval(Y{r}{s},partitionVec{r}{s},conditionVec{r}{s});
    end
    Ghat{r} = nanmean(Ghat{r},3);
end

%-------------------------------------------------------------------------%
% Create component G matrices
%-------------------------------------------------------------------------%
% 1st-finger model
for r=1:numel(Y)
    S = Ghat_1d{r};
    S = S/S(1,1);
    G1F{r} = [S(1,1),S(1,1),S(1,3),S(1,3),S(1,5),S(1,5);
              S(1,1),S(1,1),S(1,3),S(1,3),S(1,5),S(1,5);
              S(3,1),S(3,1),S(3,3),S(3,3),S(3,5),S(3,5);
              S(3,1),S(3,1),S(3,3),S(3,3),S(3,5),S(3,5);
              S(5,1),S(5,1),S(5,3),S(5,3),S(5,5),S(5,5);
              S(5,1),S(5,1),S(5,3),S(5,3),S(5,5),S(5,5)];
end

r = 4; % M1 (1:S1, 2:M1, 3:PMd, 4:IPS)
%-------------------------------------------------------------------------%
% Define models (now let the pcm_fitModelCrossval find the .theta0)
%-------------------------------------------------------------------------%
% 1. noiseceiling
M(1).type = 'noiseceiling'; 
M(1).numGparams = 0;
M(1).Gc = [];
M(1).name = 'noiseceiling';
% 2. null model 
M(2).type = 'fixed'; 
M(2).numGparams = 0;
M(2).Gc = zeros(6);
M(2).name = 'null';
% 3. component model 1
M(3).type = 'fixed'; 
M(3).numGparams = 0;
M(3).Gc = G1F{r}; % 1st-finger model
M(3).name = '1st-finger';
% 4. component model 2
M(4).type = 'fixed';
M(4).numGparams = 0;
M(4).Gc = Gc(1).G; % 2-finger transition model
M(4).name = '2-finger transition';
% 5. component model 3
M(5).type = 'fixed';
M(5).numGparams = 0;
M(5).Gc = Gc(5).G; % 6-finger transition model (=sequence model)
M(5).name = '6-finger transition';
% 6. component model 4
M(6).type = 'component';
M(6).numGparams = 3;
M(6).Gc(:,:,1) = M(2).Gc; % null-model
M(6).Gc(:,:,2) = M(3).Gc; % 1st-finger
%M(6).Gc(:,:,3) = M(4).Gc; % 2-finger transition
M(6).Gc(:,:,3) = M(5).Gc; % 6-finger transition (=sequence)
M(6).name = '1stF+6DT';

%-------------------------------------------------------------------------%
% Fit models 
%-------------------------------------------------------------------------%
[T,M] = pcm_fitModelCrossval(Y{r},M,partitionVec{r},conditionVec{r},...
            'isCheckDeriv',1,'runEffect','fixed');
% adjust for the number of voxels
for s=1:numel(Y{r})
    Nvox(s,1) = size(Y{r}{s},2);
end
T.likelihood        = bsxfun(@rdivide,T.likelihood,Nvox)*800;
T.likelihood_all    = bsxfun(@rdivide,T.likelihood_all,Nvox)*800;
        
%-------------------------------------------------------------------------%
% Show bar graph of crossvalidated log-likelihood
%-------------------------------------------------------------------------%
figure('unit','centimeters','position',[5,5,20,10]);
pcm_plotDeltaLogLikelihood(T,M,'Nnull',2,'Nceil',1,'normalize',0);


varargout = {T,M};