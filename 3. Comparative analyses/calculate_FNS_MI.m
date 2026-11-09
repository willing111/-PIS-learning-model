

clear; clc;
base=fileparts(mfilename('fullpath'));
inputFile=fullfile(base,'results','MI_MII','comparison_results.mat');
assert(isfile(inputFile),'Run calculate_MI_MII.m first.');
S=load(inputFile,'comparison');
assert(strcmp(S.comparison.SnapshotVersion,'MI_MII_GENERATED_20261109'));
preferences={ ...
    [1 1 4 5 NaN NaN;2 1 5 2 NaN NaN;2 3 2 3 5 2;3 1 3 1 NaN NaN;4 1 2 6 NaN NaN], ...
    [1 1 3 4 NaN NaN;2 1 6 1 NaN NaN;3 1 6 2 NaN NaN;4 1 5 1 NaN NaN; ...
     3 3 5 3 2 5;4 3 2 1 6 1]};
inputs=struct('Version','FNS_MI_GENERATED_20261109','Cases',{cell(2,1)});
for c=1:2
    A=S.comparison.Cases{c};
    margins=cellfun(@(x) x(29),A.MISolvedPoints);
    objectives=A.Stage1OptimalObjectives;
    calculated=-margins+A.Config.C1*cellfun(@(P) sum((P-mean(P,1)).^2,'all'),A.BaselinePNS);
    assert(max(abs(calculated-objectives))<1e-8);
    assert(isequaln(preferences{c},A.CertainPreferences));
    inputs.Cases{c}=struct('ID',sprintf('example%d',c),'Name',A.Name,'SourceIndex',c, ...
        'beta_ref',A.ReferenceData,'BaselinePNS',{A.BaselinePNS},'Delta',A.Config.Delta, ...
        'C1',A.Config.C1,'Stage1Margin',margins,'Stage1Optimum',objectives, ...
        'ReferenceMI',A.ReferenceMI,'Preferences',preferences{c},'SelectionSeed',A.MISelectionSeed);
end
out=fullfile(base,'results','FNS_MI'); if ~isfolder(out),mkdir(out);end
scales=[0 .167 .333 .5 .667 .833 1; ...
        0 .1 .3 .5 .7 .9 1; 0 .05 .25 .5 .75 .95 1];
methods={'FNS1','FNS2','FNS3','PNS'};
typeWeights=[1 1 .5 .5]; tolerance=1e-8;
WCRResults=struct('InputFile',inputFile,'SnapshotVersion',inputs.Version, ...
    'Scales',scales,'Methods',{methods},'TypeWeights',typeWeights, ...
    'SourceResultVersion',S.comparison.SnapshotVersion,'Tolerance',tolerance,'Cases',{cell(2,1)});
WCRResults.Rule=['Types 1/2: strict alternative preference/indifference; ' ...
    'types 3/4: strict pair preference/pair indifference. ' ...
    'Pair preferences require the weaker gap to be positive. ' ...
    'WCR=1-sum(typeWeight.*violations)/sum(typeWeight.*counts).'];
for c=1:2
    C=inputs.Cases{c};
    beta=C.beta_ref; [m,n,z,K]=size(beta);
    assert(m==6 && n==4 && z==7 && K==4);
    assert(max(abs(sum(beta,3)-1),[],'all')<1e-12);
    scores=zeros(m,K,4);
    for f=1:4
        for k=1:K
            if f<=3, pns=repmat(scales(f,:),n,1);
            else, pns=C.BaselinePNS{k}; end
            assert(isequal(size(pns),[n z]));
            assert(all(isfinite(pns),'all') && min(pns,[],'all')>=-tolerance ...
                && max(pns,[],'all')<=1+tolerance);
            if f==4
                assert(min(diff(pns,1,2),[],'all')>=C.Delta-tolerance);
                centered=pns-mean(pns,1);
                objective=-C.Stage1Margin(k)+C.C1*sum(centered.^2,'all');
                assert(abs(objective-C.Stage1Optimum(k))<tolerance, ...
                    'M-I PNS is not on the saved first-stage optimal face.');
            end
            for i=1:m
                for j=1:n
                    scores(i,k,f)=scores(i,k,f)+ ...
                        reshape(beta(i,j,:,k),1,z)*pns(j,:)'/n;
                end
            end
        end
    end
    assert(max(abs(scores(:,:,4)-C.ReferenceMI),[],'all')<1e-12, ...
        'Reference PNS scores differ from the aligned M-I input snapshot.');
    preferences=C.Preferences;
    caseID=C.ID;
    assert(size(preferences,2)==6 && all(ismember(preferences(:,2),1:4)));
    wcr=zeros(4,K); counts=zeros(K,4); violations=zeros(4,K,4);
    for k=1:K
        for q=1:4,counts(k,q)=sum(preferences(:,1)==k & preferences(:,2)==q);end
    end
    audit=cell(0,13);
    for f=1:4
        for p=1:size(preferences,1)
            pref=preferences(p,:);k=pref(1);q=pref(2);
            gap1=scores(pref(3),k,f)-scores(pref(4),k,f);gap2=NaN;
            if q>=3,gap2=scores(pref(5),k,f)-scores(pref(6),k,f);end
            switch q
                case 1, satisfied=gap1>tolerance;
                case 2, satisfied=abs(gap1)<=tolerance;
                case 3, satisfied=gap1-gap2>tolerance && gap2>tolerance;
                case 4, satisfied=abs(gap1-gap2)<=tolerance && gap2>tolerance;
                otherwise,error('Unknown preference type.');
            end
            violations(f,k,q)=violations(f,k,q)+double(~satisfied);
            audit(end+1,:)={methods{f},k,q,pref(3),pref(4),pref(5),pref(6), ...
                gap1,gap2,gap1-gap2,satisfied,typeWeights(q),double(~satisfied)}; %#ok<SAGROW>
        end
        for k=1:K
            denominator=sum(typeWeights.*counts(k,:));
            assert(denominator>0,'No certain preferences for DM%d.',k);
            v=reshape(violations(f,k,:),1,4);
            wcr(f,k)=1-sum(typeWeights.*v)/denominator;
        end
    end
    assert(all(abs(wcr(4,:)-1)<tolerance), ...
        'The aligned M-I point must satisfy every input certain preference.');
    summary=array2table(wcr,'VariableNames',{'DM1','DM2','DM3','DM4'});
    summary=addvars(summary,methods','Before',1,'NewVariableNames','Method');
    detail=cell2table(audit,'VariableNames',{'Method','DM','Type','i','j','x','y', ...
        'Gap_ij','Gap_xy','GapDifference','Satisfied','TypeWeight','Violation'});
    scoreRows=cell(4*K*m,4); row=0;
    for f=1:4,for k=1:K,for i=1:m
        row=row+1;scoreRows(row,:)={methods{f},k,i,scores(i,k,f)};
    end,end,end
    scoreTable=cell2table(scoreRows,'VariableNames',{'Method','DM','Alternative','CS'});
    WCRResults.Cases{c}=struct('Name',C.Name, ...
        'SourceIndex',C.SourceIndex,'Scores',scores,'Preferences',preferences, ...
        'Counts',counts,'Violations',violations,'WCR',wcr,'Summary',summary, ...
        'PreferenceAudit',detail,'ScoreTable',scoreTable,'PNS',{C.BaselinePNS}, ...
        'SelectionSeed',C.SelectionSeed);
    writetable(summary,fullfile(out,[caseID '_WCR.csv']));
    writetable(detail,fullfile(out,[caseID '_preference_audit.csv']));
    writetable(scoreTable,fullfile(out,[caseID '_reference_scores.csv']));
    fprintf('\nExample %d: reference-data WCR\n',c);disp(summary);
end
save(fullfile(out,'WCR_results.mat'),'WCRResults');
fprintf('\nSaved WCR and preference audit to: %s\n',out);
