
clear; clc;
base = fileparts(mfilename('fullpath'));
inputfile=fullfile(base,'data','comparison_data.mat');
assert(isfile(inputfile),'Run generate_comparison_data.m first.');
raw=load(inputfile,'DataCases','DataVersion');
assert(strcmp(raw.DataVersion,'COMPARISON_LDA_20261109'));
D=build_comparison_inputs(raw.DataCases);
D=prepare_direct5_inputs(D);
out = fullfile(base, 'results', 'MI_MII');
if ~isfolder(out), mkdir(out); end
sourceIndices = [1 3];
caseIDs = {'main_case','example2'};
caseNames = {'Example 1','Example 2'};
lambda = cell(2,1);

lambda{1} = [.13 .51 .26 .10; .52 .12 .25 .11; .15 .33 .20 .32; ...
             .46 .21 .16 .17; .22 .14 .29 .35]';
lambda{2} = [.16 .43 .17 .24; .31 .34 .12 .23; .21 .55 .13 .11; ...
             .13 .20 .23 .44; .27 .28 .16 .29]';
options = optimoptions('quadprog', 'Algorithm','active-set', ...
    'Display','off', 'ConstraintTolerance',1e-11, ...
    'OptimalityTolerance',1e-12, 'MaxIterations',10000);
lpoptions = optimoptions('linprog','Display','off', ...
    'ConstraintTolerance',1e-10, 'OptimalityTolerance',1e-10);
comparison = struct;
comparison.CDDefinition = 'sum_k sum_i (CSik - CSic)^2; no normalization';
comparison.Source = inputfile;
comparison.SourceSnapshotVersion = D.Experiment.Version;
comparison.OriginalSnapshotVersion = D.Experiment.OriginalVersion;
comparison.SnapshotVersion = 'MI_MII_GENERATED_20261109';
comparison.HARSeedBases = D.Experiment.HARSeedBases;
comparison.SampleSelectionRule = D.Test{1}.SampleSelectionRule;
comparison.SelectionSeeds = [NaN NaN];
save(fullfile(out,'comparison_inputs.mat'),'D','-v7');
miSeeds=[20261060 2026111503];
comparison.Cases = cell(2,1);
for c = 1:2
    T = D.Test{sourceIndices(c)};
    cfg = T.cfg;
    assert(cfg.C1 == 1 && cfg.C2 == 10 && cfg.C3 == 15);
    assert(abs(cfg.Delta-.006)<1e-12);
    assert(max(abs(sum(T.beta_ref,3)-1),[],'all')<1e-12);
    assert(max(abs(sum(T.beta_nonref,3)-1),[],'all')<1e-12);
    K=cfg.DMCount; n=cfg.Criteria; z=cfg.TermCount;
    mr=size(T.beta_ref,1); m=size(T.beta_nonref,1); r=cfg.SelectedSampleCount;
    assert(K==4 && n==4 && z==7 && r==5);
    assert(max(abs(cfg.AttributeWeights(:)-ones(n,1)/n))<1e-12);
    L=lambda{c};
    assert(all(L(:)>=.10) && max(abs(sum(L,1)-1))<1e-12);
    assert(max(abs(100*L-round(100*L)),[],'all')<1e-10);
    assert(all(arrayfun(@(k) numel(unique(round(100*L(k,:))))==5,1:K)), ...
        'Each DM must have five distinct scenario weights.');
    scoreMap=zeros(mr*K,r*K);
    nonreferenceMap=cell(K,1);
    MI=zeros(m,K); MIref=zeros(mr,K);
    baselinePNS=cell(K,1); samplePNS=cell(K,1); sampleIDs=cell(K,1);
    firstStageCheck=zeros(K,4);
    miSampleIDs=nan(K,1);
    if c==1
        rng(miSeeds(c),'twister');
        for k=1:K,miSampleIDs(k)=randi(size(T.H{k}.X,2));end
        assert(isequal(miSampleIDs(:)',[3301 1439 1148 1445]), ...
            'Main-case HAR pool differs from the retained baseline.');
    else
        assert(all(T.SpaceDimensions>0));
        randomInitials=random_feasible_initials(T.models,miSeeds(c),3,lpoptions);
    end
    originalBaselinePNS=cell(K,1);
    originalMI=zeros(m,K);originalMIref=zeros(mr,K);
    miInitialPoints=cell(K,1);miPoints=cell(K,1);
    miVerification=zeros(K,4);
    for k=1:K
        Q=T.models{k}; X=T.selected{k}.X; x0=Q.Baseline;
        B=[kron(eye(n)-ones(n)/n,eye(z)),zeros(n*z,2)];
        assert(max(abs(Q.B-B),[],'all')<1e-12, ...
            'First-stage objective differs from cross-attribute mean regularization.');
        [xr,fr,flag]=quadprog(Q.H,Q.f,Q.A,Q.b,Q.Aeq,Q.beq, ...
            Q.lb,Q.ub,x0,options);
        assert(flag>0,'First-stage verification failed, DM%d.',k);
        oldObjective=-x0(n*z+1)+cfg.C1*sum((B*x0).^2);
        assert(abs(fr-Q.Fstar)<1e-7 && abs(oldObjective-Q.Fstar)<1e-8);
        objectiveError=max(abs(-X(n*z+1,:)+cfg.C1*sum((B*X).^2,1)-Q.Fstar));
        hardError=max([Q.A*X-Q.b; Q.lb-X; X-Q.ub; ...
                       abs(Q.E*X-Q.q); zeros(1,r)],[],'all');
        assert(objectiveError<1e-8 && hardError<1e-8);
        for a=1:r
            for b=a+1:r
                assert(norm(X(1:n*z,a)-X(1:n*z,b),inf)>1e-8, ...
                    'M-II must use five distinct PNSs for every decision maker.');
            end
        end
        assert(min(X(n*z+1,:))>cfg.PositiveMarginCheckTolerance);
        firstStageCheck(k,:)=[k abs(fr-Q.Fstar) objectiveError hardError];
        Mr=make_map(T.beta_ref,k,n,z);
        Mn=make_map(T.beta_nonref,k,n,z);
        assert(max(abs(Mr*X-T.selected{k}.C),[],'all')<1e-12);
        scoreMap((k-1)*mr+(1:mr),(k-1)*r+(1:r))=Mr*X;
        nonreferenceMap{k}=Mn;
        originalBaselinePNS{k}=reshape(x0(1:n*z),z,n)';
        originalMI(:,k)=Mn*x0;originalMIref(:,k)=Mr*x0;
        if c==1
            initial=T.H{k}.X(:,miSampleIDs(k));
            [miPoint,~,miFlag]=quadprog(Q.H,Q.f,Q.A,Q.b,Q.Aeq,Q.beq, ...
                Q.lb,Q.ub,initial,options);
            assert(miFlag>0,'M-I initial-point solve failed.');
        else
            initial=randomInitials{k};
            [miPoint,~,miFlag]=quadprog(Q.H,Q.f,Q.A,Q.b,Q.Aeq,Q.beq, ...
                Q.lb,Q.ub,initial,options);
            assert(miFlag>0,'M-I first-stage solve failed.');
        end
        miObjError=abs(-miPoint(n*z+1)+cfg.C1*sum((B*miPoint).^2)-Q.Fstar);
        miPrimal=max([Q.A*miPoint-Q.b;abs(Q.E*miPoint-Q.q); ...
            Q.lb-miPoint;miPoint-Q.ub;0]);
        assert(miObjError<1e-8 && miPrimal<1e-8);
        miVerification(k,:)=[miSampleIDs(k) miObjError miPrimal norm(miPoint-initial,inf)];
        miInitialPoints{k}=initial;miPoints{k}=miPoint;
        MI(:,k)=Mn*miPoint; MIref(:,k)=Mr*miPoint;
        baselinePNS{k}=reshape(miPoint(1:n*z),z,n)';
        samplePNS{k}=X; sampleIDs{k}=T.selected{k}.IDs;
    end
    [A,b,Aeq,beq,lb,ub]=weight_constraints(T,mr,K,r);
    [feasible,~,flag]=linprog(zeros(K*r+K,1),A,b,Aeq,beq,lb,ub,lpoptions);
    assert(flag>0,'Fixed selected samples do not support the uncertain preferences.');
    P=size(L,2);
    MII=zeros(m,K,P); MIIref=zeros(mr,K,P);
    W=zeros(r,K,P); epsilon2=zeros(K,P);
    groupMI=MI*L; groupMII=zeros(m,P);
    cdMI=zeros(P,1); cdMII=zeros(P,1);
    refCDMI=zeros(P,1); refCDMII=zeros(P,1);
    verification=zeros(P,3);
    for p=1:P
        lam=L(:,p);
        G=kron(eye(K)-ones(K,1)*lam',eye(mr))*scoreMap;
        H=zeros(K*r+K);
        H(1:K*r,1:K*r)=2*cfg.C3*eye(K*r)+2*cfg.C2*(G'*G);
        H=(H+H')/2;
        f=[zeros(K*r,1);-ones(K,1)];
        [x,obj,flag,~,multipliers]=quadprog(H,f,A,b,Aeq,beq,lb,ub,feasible,options);
        assert(flag>0,'Second stage failed: %s, p%d.',caseNames{c},p);
        primal=max([A*x-b;abs(Aeq*x-beq);lb-x;x-ub;0]);
        kkt=norm(H*x+f+A'*multipliers.ineqlin+Aeq'*multipliers.eqlin ...
            -multipliers.lower+multipliers.upper,inf);
        assert(primal<1e-8 && kkt<1e-6);
        assert(all(x(K*r+1:end)>cfg.PositiveMarginCheckTolerance));
        W(:,:,p)=reshape(x(1:K*r),r,K);
        epsilon2(:,p)=x(K*r+1:end);
        MIIref(:,:,p)=reshape(scoreMap*x(1:K*r),mr,K);
        for k=1:K
            MII(:,k,p)=nonreferenceMap{k}*T.selected{k}.X*W(:,k,p);
        end
        groupMII(:,p)=MII(:,:,p)*lam;
        cdMI(p)=sum((MI-groupMI(:,p)).^2,'all');
        cdMII(p)=sum((MII(:,:,p)-groupMII(:,p)).^2,'all');
        refCDMI(p)=sum((MIref-MIref*lam).^2,'all');
        refCDMII(p)=sum((MIIref(:,:,p)-MIIref(:,:,p)*lam).^2,'all');
        verification(p,:)=[primal kkt obj];
    end
    rankMI=descending_ranks(groupMI);
    rankMII=descending_ranks(groupMII);
    rhoMI=corrcoef(rankMI); rhoMII=corrcoef(rankMII);
    assert(all(isfinite(rhoMI),'all') && all(isfinite(rhoMII),'all'));
    upper=triu(true(P),1);
    C=struct('ID',caseIDs{c},'Name',caseNames{c},'Config',cfg, ...
        'Lambda',L,'MI',MI,'MII',MII,'ReferenceMI',MIref, ...
        'ReferenceMII',MIIref,'GroupMI',groupMI,'GroupMII',groupMII, ...
        'RankMI',rankMI,'RankMII',rankMII,'RhoMI',rhoMI,'RhoMII',rhoMII, ...
        'MeanRhoMI',mean(rhoMI(upper)),'MeanRhoMII',mean(rhoMII(upper)), ...
        'CD_MI',cdMI,'CD_MII',cdMII,'ReferenceCD_MI',refCDMI, ...
        'ReferenceCD_MII',refCDMII,'Weights',W,'Epsilon2',epsilon2, ...
        'BaselinePNS',{baselinePNS},'SelectedPNS',{samplePNS}, ...
        'SampleIDs',{sampleIDs},'ReferenceData',T.beta_ref,'NonreferenceData',T.beta_nonref, ...
        'AlternativeUncertain',T.altpref, ...
        'PairUncertain',T.pairpref,'FirstStageVerification',firstStageCheck, ...
        'SecondStageVerification',verification);
    if c==2
        C.CertainPairs=T.NewCertainPairs;
        C.CertainPairPreferences=T.CertainPairPreferences;
        C.CertainPreferences=T.CertainPreferences;
        C.SpaceDimensions=T.SpaceDimensions;
        C.NonreferenceSeed=2026102301;
    else
        C.CertainConstraints='e1:a4>a5; e2:(CS2-CS3)>(CS5-CS2),a5>a2; e3:a3>a1; e4:a2>a6';
        C.NonreferenceSeed=cfg.NonreferenceSeed;
    end
    C.SelectionSeed=NaN;
    C.CertainPreferences=T.CertainPreferences;
    C.Stage1OptimalObjectives=cellfun(@(Q) Q.Fstar,T.models);
    C.HARSeedBase=cfg.HARSeedBase;
    C.SampleSelectionRule=T.SampleSelectionRule;
    C.DirectHARChecks=T.DirectHARChecks;
    C.MISelectionSeed=miSeeds(c);
    C.MISelectionProvenance='';
    C.MISelectionHARIDs=miSampleIDs;
    if c==1
       
    else
        C.MISelectionRule='Three random linear programs over the original feasible set; normalized independent -log(U) convex weights; original first-stage active-set QP; one solved point fixed across all five weight settings.';
        C.MIInitialVertexCount=3;

        C.MIInitialGenerator='twister';
        C.SampleSelectionRule=T.SampleSelectionRule;
    end
    C.MIInitialPoints=miInitialPoints;
    C.MISolvedPoints=miPoints;
    C.MIVerification=miVerification;
    C.OriginalBaselinePNS=originalBaselinePNS;
    C.OriginalMI=originalMI;
    C.OriginalReferenceMI=originalMIref;
    comparison.Cases{c}=C;
    folder=fullfile(out,C.ID);
    if ~isfolder(folder),mkdir(folder);end
    writetable(table((1:K)',repmat(miSeeds(c),K,1),miSampleIDs, ...
        miVerification(:,2),miVerification(:,3),miVerification(:,4), ...
        'VariableNames',{'DM','MISeed','HARSampleID','ObjectiveError','ConstraintViolation','MovementFromStart'}), ...
        fullfile(folder,'MI_selection_audit.csv'));
    labels=compose('p%d',(1:P)');
    writetable(table(labels,L(1,:)',L(2,:)',L(3,:)',L(4,:)',cdMI,cdMII, ...
        'VariableNames',{'Setting','Lambda1','Lambda2','Lambda3','Lambda4','CD_MI','CD_MII'}), ...
        fullfile(folder,'weights_CD.csv'));
    names=cellstr(compose('p%d',1:P));
    writetable(array2table([(1:m)' groupMI],'VariableNames',[{'Alternative'} names]),fullfile(folder,'group_scores_MI.csv'));
    writetable(array2table([(1:m)' groupMII],'VariableNames',[{'Alternative'} names]),fullfile(folder,'group_scores_MII.csv'));
    writetable(array2table([(1:m)' rankMI],'VariableNames',[{'Alternative'} names]),fullfile(folder,'ranks_MI.csv'));
    writetable(array2table([(1:m)' rankMII],'VariableNames',[{'Alternative'} names]),fullfile(folder,'ranks_MII.csv'));
    writetable(array2table(rhoMI,'VariableNames',names),fullfile(folder,'Spearman_MI.csv'));
    writetable(array2table(rhoMII,'VariableNames',names),fullfile(folder,'Spearman_MII.csv'));
    fprintf('\n%s\n',C.Name);
    disp(table(labels,cdMI,cdMII,'VariableNames',{'Setting','CD_MI','CD_MII'}));
    fprintf('Mean CD: M-I=%.9f, M-II=%.9f\n',mean(cdMI),mean(cdMII));
    fprintf('Mean Spearman (10 distinct pairs): M-I=%.9f, M-II=%.9f\n',C.MeanRhoMI,C.MeanRhoMII);
end

comparison.VerifiedMetrics=cell2mat(cellfun(@(C) [mean(C.CD_MI) mean(C.CD_MII) C.MeanRhoMI C.MeanRhoMII],comparison.Cases,'UniformOutput',false));
save(fullfile(out,'comparison_results.mat'),'comparison','-v7');
fprintf('\nSaved: %s\nNext: run calculate_FNS_MI.m or calculate_MII_MIII_MIV.m\n',fullfile(out,'comparison_results.mat'));

function M=make_map(beta,k,n,z)
    M=zeros(size(beta,1),n*z+2);
    for i=1:size(beta,1)
        for j=1:n
            M(i,(j-1)*z+(1:z))=reshape(beta(i,j,:,k),1,z)/n;
        end
    end
end

function [A,b,E,e,lb,ub]=weight_constraints(T,m,K,r)
    nv=K*r+K; A=zeros(0,nv);
    tol=T.cfg.EventSignTolerance;
    for k=1:K
        C=T.selected{k}.C;
        assert(size(C,1)==m && size(C,2)==r);
        ar=T.altpref(T.altpref(:,1)==k,:);
        pr=T.pairpref(T.pairpref(:,1)==k,:);
        assert(~isempty(ar)||~isempty(pr));
        for a=1:size(ar,1)
            d=C(ar(a,2),:)-C(ar(a,3),:);
            v=double(d>tol)-double(d< -tol);
            row=zeros(1,nv);row((k-1)*r+(1:r))=-v;row(K*r+k)=1;
            A(end+1,:)=row;
        end
        for a=1:size(pr,1)
            v=double(C(pr(a,2),:)-C(pr(a,3),:)>tol) ...
             -double(C(pr(a,4),:)-C(pr(a,5),:)>tol);
            row=zeros(1,nv);row((k-1)*r+(1:r))=-v;row(K*r+k)=1;
            A(end+1,:)=row;
        end
    end
    b=zeros(size(A,1),1);
    E=[kron(eye(K),ones(1,r)),zeros(K)];e=ones(K,1);
    lb=[zeros(K*r,1);ones(K,1)*1e-6];ub=ones(nv,1);
end

function ranks=descending_ranks(scores)
    ranks=zeros(size(scores));
    for p=1:size(scores,2)
        [values,order]=sort(scores(:,p),'descend');a=1;
        while a<=numel(values)
            b=a;
            while b<numel(values) && values(b+1)==values(a),b=b+1;end
            ranks(order(a:b),p)=(a+b)/2;
            a=b+1;
        end
    end
end

function D = prepare_direct5_inputs(D)
seedBases = [2026362000 2032691700 2026111800 2026113068; ...
             2026362000 2032691700 2026111800 2026111900];
indices = [1 3];
for c = 1:2
    T = D.Test{indices(c)};
    cfg = T.cfg;
    cfg.HARChains = 5;
    cfg.HARBurnIn = 2000;
    cfg.HARThin = 10;
    cfg.HARPerChain = 1;
    cfg.HARSeedBase = seedBases(c,:);
    cfg.SelectedSampleCount = 5;
    if isfield(cfg,'SelectionSeed')
        cfg = rmfield(cfg,'SelectionSeed');
    end
    assert(cfg.C1==1 && cfg.C2==10 && cfg.C3==15 && abs(cfg.Delta-.006)<1e-12);
    assert(cfg.DMCount==4 && cfg.Criteria==4 && cfg.TermCount==7);
    selected = cell(cfg.DMCount,1);
    dimensions = zeros(cfg.DMCount,1);
    checks = zeros(cfg.DMCount,3);
    for k = 1:cfg.DMCount
        Q = T.models{k};
        [X,chainIDs,dim] = optimal_face_har(Q,cfg,k);
        assert(isequal(size(X),[numel(Q.Baseline) 5]));
        objectiveError = max(abs(-X(29,:)+cfg.C1*sum((Q.B*X).^2,1)-Q.Fstar));
        hardError = max([Q.A*X-Q.b;Q.lb-X;X-Q.ub;zeros(1,5)],[],'all');
        equalityError = max(abs(Q.E*X-Q.q),[],'all');
        assert(max([objectiveError hardError equalityError])<cfg.FeasibilityTolerance);
        assert(min(Q.D*X,[],'all')>=cfg.Delta-cfg.FeasibilityTolerance);
        assert(min(X(29,:))>cfg.PositiveMarginCheckTolerance);
        selected{k} = struct('IDs',(1:5)','X',X,'C',Q.C*X,'ChainIDs',chainIDs(:));
        dimensions(k) = dim;
        checks(k,:) = [objectiveError hardError equalityError];
    end
    T.cfg = cfg;
    T.selected = selected;
    T.DirectHARSpaceDimensions = dimensions;
    T.SpaceDimensions = dimensions;
    T.DirectHARChecks = checks;
    T.SampleSelectionRule = 'Five chains; seed=HARSeedBase(k)+chain; burn-in 2000; retain iteration 2010; no pool subsampling.';
    D.Test{indices(c)} = T;
end
D.Experiment.OriginalVersion = D.Experiment.Version;
D.Experiment.Version = 'COMPARISON_MODELS_20261109';
D.Experiment.HARSeedBases = seedBases;
end
function [X, chainids, dim] = optimal_face_har(Q, cfg, k, legacyOrigin)
    if nargin<4,legacyOrigin=false;end
    lp = optimoptions('linprog', 'Display', 'off', ...
        'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
    B = Q.Reduction;
    yy = zeros(size(B.A, 2), double(legacyOrigin));
    universal = false(size(B.b));
    for r = 1:numel(B.b)
        [y, mn, flag] = linprog(B.A(r, :)', B.A, B.b, [], [], [], [], lp);
        assert(flag > 0, '最优集合约束分析失败。');
        universal(r) = B.b(r) - mn < 2e-9;
        yy(:, end + 1) = y;
    end

    origin = mean(yy, 2);
    F = B.A(universal, :);
    g = B.b(universal);
    if isempty(F)
        Z = eye(size(B.A, 2));
    else
        [U, D, V] = svd(F, 'econ');
        values = diag(D);
        rk = sum(values > 1e-9);
        if rk > 0
            eq = V(:, 1:rk)';
            rhs = (U(:, 1:rk)' * g) ./ values(1:rk);
            origin = origin + eq' * (rhs - eq * origin);
            Z = null(eq, 1e-9);
        else
            Z = eye(size(B.A, 2));
        end
    end

    dim = size(Z, 2);
    if dim == 0
        X = B.Origin + B.Z * origin;
        chainids = 1;
        return;
    end
    G = B.A * Z;
    h = B.b - B.A * origin;
    norms = vecnorm(G, 2, 2);
    nz = norms > 1e-10;
    assert(all(h(~nz) > -1e-8));
    G = G(nz, :);
    h = h(nz);

    [cc, ~, flag] = linprog([zeros(dim, 1); -1], ...
        [G, vecnorm(G, 2, 2)], h, [], [], [-inf(dim, 1); 0], [], lp);
    assert(flag > 0 && cc(end) > 1e-10, '最优集合内点求解失败。');
    center = cc(1:dim);
    slack = h - G * center;
    HH = G' * (G ./ slack.^2);
    [V, D] = eig((HH + HH') / 2);
    scale = V * diag(1 ./ sqrt(diag(D))) * V';
    GG = G * scale;
    hh = h - G * center;
    count = cfg.HARChains * cfg.HARPerChain;
    X = zeros(numel(Q.Baseline), count);
    chainids = repelem(1:cfg.HARChains, cfg.HARPerChain);
    for chain = 1:cfg.HARChains
        rng(cfg.HARSeedBase(k) + chain, 'twister');
        z = zeros(dim, 1);
        saved = 0;
        for it = 1:(cfg.HARBurnIn + cfg.HARThin * cfg.HARPerChain)
            direction = randn(dim, 1);
            direction = direction / norm(direction);
            slopes = GG * direction;
            s = hh - GG * z;
            lo = max(s(slopes < -1e-14) ./ slopes(slopes < -1e-14));
            hi = min(s(slopes > 1e-14) ./ slopes(slopes > 1e-14));
            assert(isfinite(lo) && isfinite(hi) && lo <= hi + 1e-10);
            z = z + (lo + rand * (hi - lo)) * direction;

            if it > cfg.HARBurnIn && mod(it - cfg.HARBurnIn, cfg.HARThin) == 0
                saved = saved + 1;
                ix = (chain - 1) * cfg.HARPerChain + saved;
                X(:, ix) = B.Origin + B.Z * (origin + Z * (center + scale * z));
            end
        end
    end
end

function D=build_comparison_inputs(DataCases)
assert(numel(DataCases)==2);
seeds=[42 2026;101 2026102301];
seedBases=[2026362000 2032691700 2026111800 2026113068; ...
    2026362000 2032691700 2026111800 2026111900];
preferences={ ...
    [1 1 4 5 NaN NaN;2 1 5 2 NaN NaN;2 3 2 3 5 2;3 1 3 1 NaN NaN;4 1 2 6 NaN NaN], ...
    [1 1 3 4 NaN NaN;2 1 6 1 NaN NaN;3 1 6 2 NaN NaN;4 1 5 1 NaN NaN; ...
     3 3 5 3 2 5;4 3 2 1 6 1]};
alternativeUncertain={ [3 5 6],zeros(0,3) };
pairUncertain={ [1 5 1 1 4;2 3 6 6 5;4 1 2 2 3], ...
    [1 1 2 1 3;2 6 1 1 6;3 1 3 5 1;4 3 1 3 5] };
D.Test=cell(3,1);
indices=[1 3];
for c=1:2
    R=DataCases{c};cfg=R.GeneratorSettings;
    assert(cfg.ReferenceSeed==seeds(c,1) && cfg.NonreferenceSeed==seeds(c,2));
    cfg.Delta=.006;cfg.C1=1;cfg.C2=10;cfg.C3=15;cfg.MinimumWeight=0;
    cfg.DMWeights=ones(4,1)/4;cfg.AttributeWeights=ones(4,1)/4;
    cfg.HARChains=5;cfg.HARBurnIn=2000;cfg.HARThin=10;cfg.HARPerChain=1;
    cfg.HARSeedBase=seedBases(c,:);cfg.SelectedSampleCount=5;
    cfg.FeasibilityTolerance=1e-8;cfg.PositiveMarginCheckTolerance=1e-7;
    cfg.EventSignTolerance=1e-12;
    assert(isequal(size(R.beta_ref),[6 4 7 4]) && isequal(size(R.beta_nonref),[6 4 7 4]));
    T=struct('Name',R.Name,'cfg',cfg,'beta_ref',R.beta_ref,'beta_nonref',R.beta_nonref, ...
        'CertainPreferences',preferences{c},'altpref',alternativeUncertain{c},'pairpref',pairUncertain{c});
    T.models=build_stage1_models(T,c);
    T.H=cell(4,1);
    if c==1
        poolCfg=cfg;poolCfg.HARChains=4;poolCfg.HARPerChain=1250;
        poolCfg.HARSeedBase=[2026361600 2032691700 2026111800 2026111900];
        for k=1:4
            [X,~,dimension]=optimal_face_har(T.models{k},poolCfg,k,true);
            T.H{k}=struct('X',X,'C',T.models{k}.C*X,'Dimension',dimension);
        end
    else
        T.NewCertainPairs=[3 4;6 1;6 2;5 1];
        T.CertainPairPreferences=[3 5 3 2 5;4 2 1 6 1];
    end
    D.Test{indices(c)}=T;
end
D.Experiment=struct('Version','COMPARISON_LDA_20261109');
end

function models=build_stage1_models(T,caseIndex)
cfg=T.cfg;beta_ref=T.beta_ref;m=cfg.ReferenceCount;n=cfg.Criteria;z=cfg.TermCount;
np=n*z;N=np+2;ie=np+1;id=np+2;e=zeros(1,N);e(ie)=1;ed=zeros(1,N);ed(id)=1;
lp=optimoptions('linprog','Display','off','ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);
qp=optimoptions('quadprog','Display','off','Algorithm','active-set', ...
    'ConstraintTolerance',1e-11,'OptimalityTolerance',1e-12,'MaxIterations',10000);
models=cell(cfg.DMCount,1);
for k=1:cfg.DMCount
    D=zeros(n*(z-1),N);
    for j=1:n
        for t=1:z-1
            r=(j-1)*(z-1)+t;D(r,(j-1)*z+t)=-1;D(r,(j-1)*z+t+1)=1;
        end
    end
    U=zeros(m*n,N);C=zeros(m,N);
    for i=1:m
        for j=1:n
            r=(i-1)*n+j;U(r,(j-1)*z+(1:z))=reshape(beta_ref(i,j,:,k),1,z);
            C(i,:)=C(i,:)+U(r,:)/n;
        end
    end
    pref=T.CertainPreferences(T.CertainPreferences(:,1)==k,:);
    G=zeros(0,N);extra=zeros(0,N);
    for a=1:size(pref,1)
        left=C(pref(a,3),:)-C(pref(a,4),:);
        if pref(a,2)==1
            if caseIndex==1 && any(pref(:,2)==3 & pref(:,5)==pref(a,3) & pref(:,6)==pref(a,4))
                continue;
            end
            G(end+1,:)=left;
        else
            assert(pref(a,2)==3);
            right=C(pref(a,5),:)-C(pref(a,6),:);
            if caseIndex==1,G=[G;left-right;right];
            else,extra=[extra;left-right;right];end
        end
    end
    B=[kron(eye(n)-ones(n)/n,eye(z)),zeros(np,2)];H=2*cfg.C1*(B'*B);f=-e';
    A=[-D+repmat(ed,size(D,1),1);U;-U;C;-C;-G+repmat(e,size(G,1),1); ...
        -extra+repmat(e,size(extra,1),1)];
    b=[zeros(size(D,1),1);ones(m*n,1);zeros(m*n,1);ones(m,1);zeros(m,1); ...
        zeros(size(G,1)+size(extra,1),1)];
    marginLower=1e-6;if caseIndex==1,marginLower=-inf;end
    lb=[zeros(np,1);marginLower;cfg.Delta];ub=[ones(np,1);inf;1];Aeq=ed;beq=cfg.Delta;
    [x0,~,flag]=linprog(f,A,b,Aeq,beq,lb,ub,lp);assert(flag>0);
    [x,~,flag,~,mu]=quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,qp);assert(flag>0);
    Fstar=-x(ie)+cfg.C1*sum((B*x).^2);
    primal=max([A*x-b;abs(Aeq*x-beq);lb-x;x-ub;0]);
    kkt=norm(H*x+f+A'*mu.ineqlin+Aeq'*mu.eqlin-mu.lower+mu.upper,inf);
    assert(primal<cfg.FeasibilityTolerance && kkt<1e-7 && x(ie)>cfg.PositiveMarginCheckTolerance);
    E=[Aeq;B;e];q=[beq;B*x;x(ie)];Z=null(E,1e-10);I=eye(N);
    il=isfinite(lb);iu=isfinite(ub);allA=[A;-I(il,:);I(iu,:)];allb=[b;-lb(il);ub(iu)];
    reducedA=allA*Z;reducedb=allb-allA*x;norms=vecnorm(reducedA,2,2);keep=norms>1e-10;
    assert(min(reducedb)>-cfg.FeasibilityTolerance);reducedb(abs(reducedb)<1e-10)=0;
    reduction=struct('Origin',x,'Z',Z,'A',reducedA(keep,:)./norms(keep),'b',reducedb(keep)./norms(keep));
    models{k}=struct('DM',k,'A',A,'b',b,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub, ...
        'B',B,'H',H,'f',f,'D',D,'U',U,'C',C,'G',G,'Baseline',x,'Fstar',Fstar,'E',E,'q',q,'Reduction',reduction);
end
end

function initials=random_feasible_initials(models,seed,vertexCount,lp)
rng(seed,'twister');initials=cell(numel(models),1);
for k=1:numel(models)
    Q=models{k};N=numel(Q.Baseline);V=zeros(N,vertexCount);
    for v=1:vertexCount
        direction=[randn(N-2,1);0;0];
        [V(:,v),~,flag]=linprog(direction,Q.A,Q.b,Q.Aeq,Q.beq,Q.lb,Q.ub,lp);
        assert(flag>0,'Random feasible initial-point construction failed.');
    end
    w=-log(rand(vertexCount,1));w=w/sum(w);x0=V*w;
    violation=max([Q.A*x0-Q.b;abs(Q.Aeq*x0-Q.beq);Q.lb-x0;x0-Q.ub;0]);
    assert(violation<1e-8);initials{k}=x0;
end
end
