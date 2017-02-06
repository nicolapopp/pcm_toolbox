function pcm_plotDeltaLogLikelihood(T,M,varargin);
%%

Nnull       = 2;
Nceil       = 1;
modelNames  = [];
figDir      = 'horz';
normalize   = 0;
pcm_vararginoptions(varargin,{'Nnull','Nceil','modelNames','normalize'});

Nmodel = size(T.likelihood,2);
idxmodel = setdiff(1:Nmodel,[Nnull,Nceil]);

%=========================================================%
% Calculate plotting values
%=========================================================%
T.dLike     = bsxfun(@minus,T.likelihood,T.likelihood(:,Nnull)); % delta-log-likelihood over the null-model

% Re-calculate ceiling including other models, as the definition is the
% maximally possibly achievable log-likelihood.
[maxLike,idx] = max(T.likelihood_all,[],2); 

Ce.upper     = maxLike - T.likelihood(:,Nnull); % adjust for the null-model
T.nLike      = bsxfun(@rdivide,T.dLike,Ce.upper); % normalize 

% Calculate lower noiseceiling
for i=1:size(T.likelihood,1)
    Ce.lower(i,1) = T.likelihood(i,idx(i)) - T.likelihood(i,Nnull);
end
Ce.nLower    = Ce.lower ./ Ce.upper; % normalized lower ceiling
Ce.nUpper    = Ce.upper ./ Ce.upper; % normalized upper ceiling (=1)

switch normalize
    case 1
        Plotval = [T.nLike(:,idxmodel)];
        Ceilings= [1,mean(Ce.nLower)];
        label   = 'Norm \Delta log-likelihood';
        ylim    = [-0.3 1.3];
    otherwise
        Plotval = [T.dLike(:,idxmodel)];
        Ceilings= [mean(Ce.upper),mean(Ce.lower)];
        label   = '\Delta log-likelihood';
        limval  = [Plotval,Ce.lower,Ce.upper];
        meany   = mean(limval,1);
        [maxy,maxx] = max(meany);
        e = stderr(limval(:,maxx));        
        ylim = [0,maxy+1.5*e];
end
for m=1:numel(idxmodel)
    if isfield(M(idxmodel(m)),'name')
        modelnames{m} = M(idxmodel(m)).name;
    else
        modelnames{m} = num2str(m);
    end
end


%=========================================================%
% plot data
%=========================================================%
patch([0.45 8.0 8.0 0.45],... % Draw noise ceiling first
            [Ceilings(1),Ceilings(1),Ceilings(2),Ceilings(2)],...
            [0.8 0.8 0.8],'edgecolor','non'); hold on;

[xpos,ypos,e] = barplot([],Plotval,... % Plot log-likelihood
    'edgecolor',[1 1 1],...
    'facecolor',[0.5 0.5 0.5],...
    'barwidth',0.9,'gapwidth',[0.1 0.1 0.1],...
    'errorwidth',2,'capwidth',0.001);

ylabel(label)
xlabel('');
set(gca,'Xticklabel',modelnames);
set(gca,'YLim',ylim,'Xlim',[xpos(1)-0.5,xpos(end)+0.5]);
set(gca,'tickdir','out','ticklength',[0.025,0.025]);

switch (figDir)
    case 'vert'
        xlabelangle = 90;
    case 'horz'
        set(gca,'view',[90,90],'yaxislocation','left');
        xlabelangle = 0;
end

end