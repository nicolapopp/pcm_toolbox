function varargout = pcm_plotStepwiseLikelihood(Tall,Ce,M,history,varargin)
%%
% 
% 
% 
% 
% 
% 
% 
% 

normalize   = 0;
modelNames  = {};
ceilingtype = 'patch'; % control way to show noise ceiling

pcm_vararginoptions(varargin,{'normalize','modelNames',...
    'ceilingtype'});


%=========================================================%
% Choose result based on history of foward selection to get
%=========================================================%
[compLike,ceilings,models] = pcm_stepwise_getCompLike(Tall,Ce,history);
for i=1:numel(history.logBF)
    maxBF(i) = max(history.logBF{i});
end
[~,history.bestiter] = max(maxBF);
[~,ceilings,~] = pcm_stepwise_getCompLike(Tall,Ce,history);

%=========================================================%
% Plot
%=========================================================%
switch normalize
    case 1
        Plotval = bsxfun(@rdivide,compLike,ceilings(:,2));
        Ceilings= [1,mean(ceilings(:,1)./ceilings(:,2))];
        label   = 'Norm \Delta log-likelihood';
        ylim    = [-0.3 1.3];
    otherwise
        Plotval = compLike;
        Ceilings= [mean(ceilings(:,2)),mean(ceilings(:,1))];
        label   = '\Delta log-likelihood';
        limval  = [Plotval,ceilings];
        meany   = mean(limval,1);
        [maxy,maxx] = max(meany);
        e = stderr(limval(:,maxx));
        
        ylim = [0,maxy+1.5*e];
        ylim = [min(meany)*1.3,Ceilings(1)*1.3];
end

%=========================================================%
% plot data
%=========================================================%
% Draw noise ceiling and so on
switch ceilingtype
    case 'patch'        
        patch([0.05 size(Plotval,2)+1.05 size(Plotval,2)+1.05 0.05],...
            [Ceilings(1),Ceilings(1),Ceilings(2),Ceilings(2)],...
            [0.8 0.8 0.8],'edgecolor','non'); hold on;
    case 'line'
        drawline(Ceilings(1),'dir','horz','color',[0.9 0 0]);
        drawline(Ceilings(2),'dir','horz','color',[0.5 0.5 0.5]);
end

% Plot log-likelihood
[xpos,ypos,e] = barplot([],Plotval,...
    'edgecolor',[0 0 0],...
    'facecolor',[0.5 0.5 0.5],...
    'barwidth',0.8);

% set xtick label
if isempty(modelNames)
    for m=1:size(Plotval,2)-1
        modelNames{m} = sprintf('Model %d',m);
    end
    modelNames{end+1} = 'Best model combination';
end


ylabel(label)
xlabel('');
set(gca,'Xtick',xpos,'Xticklabel',modelNames);
set(gca,'YLim',ylim,'Xlim',[0 size(Plotval,2)+1],'view',[90 90]);
set(gca,'tickdir','out','ticklength',[0.025,0.025]);
title('Component log-likelihood')

drawline((xpos(end)+xpos(end-1))/2,'linestyle','--');

varargout = {Plotval};
