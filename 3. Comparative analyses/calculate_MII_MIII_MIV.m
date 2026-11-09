

clear;clc;
base=fileparts(mfilename('fullpath'));addpath(base);
miiFile=fullfile(base,'results','MI_MII','comparison_results.mat');
inputFile=fullfile(base,'results','MI_MII','comparison_inputs.mat');
assert(isfile(miiFile) && isfile(inputFile),'请先运行新版calculate_MI_MII.m。');
S=load(miiFile,'comparison');
P=load(inputFile,'D');D=P.D;
settings=struct;
settings.ClassMeaning='1 low, 2 middle, 3 high; reference alternatives only';
settings.ClassificationRule='Fixed exploratory reference classification shared by Li and Xiao.';
settings.Li=struct('MinGap',1e-6,'Alpha',.5,'Beta',.5);
settings.XiaoParameters=[1 .01 .001];
settings.ClassesExample1=[2 1 1 3 2 3;1 2 1 2 3 3;2 2 3 3 1 1;3 2 3 1 2 1];
settings.ClassesExample2=[2 1 3 2 1 3;1 2 3 3 1 2;2 3 1 1 2 3;1 3 1 2 2 3];
assert(strcmp(D.Experiment.Version,'COMPARISON_MODELS_20261109'));
assert(strcmp(S.comparison.SourceSnapshotVersion,D.Experiment.Version));
assert(strcmp(S.comparison.SnapshotVersion,'MI_MII_GENERATED_20261109'), ...
    '请先运行新版calculate_MI_MII.m，生成匹配的M-II结果。');
settings.Version='MII_MIII_MIV_GENERATED_20261109';
settings.HARSeedBases=S.comparison.HARSeedBases;
out=fullfile(base,'results','MII_MIII_MIV');if ~isfolder(out),mkdir(out);end
comparison=struct('Version',settings.Version,'MethodNames',{{'M-II','M-III','M-IV'}}, ...
 'CDDefinition','sum_k sum_i (CSik-CSic)^2, no normalization', ...
 'SpearmanDefinition','Mean of 10 off-diagonal upper-triangle correlations between five rankings', ...
 'Classification',settings,'SourceExperiment',D.Experiment, ...
 'SourceResult',miiFile,'SourceResultVersion',S.comparison.SnapshotVersion,'Cases',{cell(2,1)});
indices=[1 3];
for c=1:2
 T=D.Test{indices(c)};prior=S.comparison.Cases{c};L=prior.Lambda;
 [mr,n,z,K]=size(T.beta_ref);mn=size(T.beta_nonref,1);
 assert(isequal(T.cfg,prior.Config));
 assert(isequal(T.beta_ref,prior.ReferenceData) && isequal(T.beta_nonref,prior.NonreferenceData));
 assert(isequal(T.altpref,prior.AlternativeUncertain) && isequal(T.pairpref,prior.PairUncertain));
 assert(all(L(:)>=.10) && max(abs(sum(L,1)-1))<1e-12);
 assert(max(abs(100*L-round(100*L)),[],'all')<1e-10);
 assert(all(arrayfun(@(k) numel(unique(round(100*L(k,:))))==5,1:K)));
 for k=1:K
  assert(isequal(T.selected{k}.IDs(:),prior.SampleIDs{k}(:)));
  assert(max(abs(T.selected{k}.X-prior.SelectedPNS{k}),[],'all')<1e-12);
 end
 classes=settings.(sprintf('ClassesExample%d',c));
 assert(isequal(size(classes),[4 6]));
 assert(all(sum(classes==1,2)==2) && all(sum(classes==2,2)==2) && all(sum(classes==3,2)==2));
 if c==1,strictPairs={ [4 5],[2 3;5 2],[3 1],[2 6] };
 else,strictPairs=arrayfun(@(k) T.NewCertainPairs(k,:),1:K,'UniformOutput',false);end
 for k=1:K
  pairs=strictPairs{k};assert(all(classes(k,pairs(:,1))>=classes(k,pairs(:,2))));
 end
 fprintf('Calculating %s\n',prior.Name);
 [mii,miiAudit]=solve_mii(T,L);
 assert(max(abs(mii-prior.MII),[],'all')<2e-6,'M-II differs from aligned experiment.');
 li=cell(K,1);liNon=zeros(mn,K);liRef=zeros(mr,K);
 for k=1:K
  li{k}=solve_Li_P2(T.beta_ref(:,:,:,k),classes(k,:), ...
      MinGap=settings.Li.MinGap,Alpha=settings.Li.Alpha,Beta=settings.Li.Beta);
  liRef(:,k)=li{k}.CS;
  liNon(:,k)=score_scale(T.beta_nonref(:,:,:,k),li{k}.PNS);
 end
 xiao=solve_Xiao_transfer(T.beta_ref,T.beta_nonref,classes,L);
 predictions={mii,repmat(liNon,1,1,5),xiao.PersonalScores};
 metrics=cell(3,1);
 for method=1:3,metrics{method}=evaluate(predictions{method},L);end
  assert(max(abs(metrics{1}.CD-prior.CD_MII))<2e-6,'两组比较的M-II CD不一致。');
 assert(max(abs(metrics{1}.GroupScores-prior.GroupMII),[],'all')<2e-6, ...
     '两组比较的M-II群体评分不一致。');
 assert(max(abs(metrics{1}.Rho-prior.RhoMII),[],'all')<1e-10, ...
     '两组比较的M-II排序相关矩阵不一致。');
 C=struct('ID',prior.ID,'Name',prior.Name,'Lambda',L,'Classes',classes, ...
  'SelectionSeed',NaN,'HARSeedBase',T.cfg.HARSeedBase, ...
  'SampleSelectionRule',T.SampleSelectionRule,'SelectedSampleIDs',{prior.SampleIDs}, ...
  'CertainAlternativePairs',{strictPairs},'AlternativeUncertain',T.altpref, ...
  'PairUncertain',T.pairpref,'LiModels',{li},'LiReferenceScores',liRef, ...
  'Xiao',xiao,'MIIAudit',miiAudit,'Metrics',{metrics});
 if c==2,C.CertainPairPreferences=T.CertainPairPreferences;end
 comparison.Cases{c}=C;
 folder=fullfile(out,prior.ID);if ~isfolder(folder),mkdir(folder);end
 names={'p1','p2','p3','p4','p5'};
 writetable(array2table([(1:K)' classes],'VariableNames', ...
  {'DM','a1','a2','a3','a4','a5','a6'}),fullfile(folder,'classes_Li_Xiao.csv'));
 writetable(array2table([(1:5)' L' metrics{1}.CD metrics{2}.CD metrics{3}.CD], ...
  'VariableNames',{'Setting','Lambda1','Lambda2','Lambda3','Lambda4','CD_MII','CD_MIII','CD_MIV'}),fullfile(folder,'weights_CD.csv'));
 for method=1:3
  M=metrics{method};tag=comparison.MethodNames{method};
  writetable(array2table(M.Rho,'VariableNames',names),fullfile(folder,['Spearman_' tag '.csv']));
  writetable(array2table([(1:mn)' M.GroupScores],'VariableNames',[{'Alternative'} names]),fullfile(folder,['group_scores_' tag '.csv']));
  writetable(array2table([(1:mn)' M.Ranks],'VariableNames',[{'Alternative'} names]),fullfile(folder,['ranks_' tag '.csv']));
  writematrix(reshape(M.PersonalScores,mn,[]),fullfile(folder,['personal_scores_' tag '.csv']));
  fprintf('%s meanCD=%.9f meanSpearman=%.9f\n',tag,M.MeanCD,M.MeanRho);
 end
end
save(fullfile(out,'comparison_results.mat'),'comparison','-v7');
fprintf('Saved numerical results. Next: plot_MII_MIII_MIV.m\n');

function cs=score_scale(beta,pns)
[m,n,z]=size(beta);cs=zeros(m,1);
for j=1:n,cs=cs+reshape(beta(:,j,:),m,z)*pns(j,:)'/n;end
end
function M=evaluate(personal,L)
[m,~,P]=size(personal);group=zeros(m,P);ranks=zeros(m,P);CD=zeros(P,1);
for p=1:P
 group(:,p)=personal(:,:,p)*L(:,p);
 CD(p)=sum((personal(:,:,p)-group(:,p)).^2,'all');
 ranks(:,p)=tiedrank(-group(:,p));
end
rho=corrcoef(ranks);assert(all(isfinite(rho),'all'));
M=struct('PersonalScores',personal,'GroupScores',group,'Ranks',ranks, ...
 'CD',CD,'MeanCD',mean(CD),'Rho',rho,'MeanRho',mean(rho(triu(true(P),1)))) ;
end

function [personal,audit]=solve_mii(T,L)
cfg=T.cfg;[mr,n,z,K]=size(T.beta_ref);r=cfg.SelectedSampleCount;
mn=size(T.beta_nonref,1);maps=cell(K,1);scoreMap=zeros(mr*K,K*r);
for k=1:K
 Q=T.models{k};X=T.selected{k}.X;
 B=[kron(eye(n)-ones(n)/n,eye(z)),zeros(n*z,2)];
 assert(max(abs(Q.B-B),[],'all')<1e-12);
 assert(max(abs(-X(n*z+1,:)+cfg.C1*sum((B*X).^2,1)-Q.Fstar))<1e-8);
 assert(max([Q.A*X-Q.b;abs(Q.E*X-Q.q);Q.lb-X;X-Q.ub;zeros(1,r)],[],'all')<1e-8);
 scoreMap((k-1)*mr+(1:mr),(k-1)*r+(1:r))=T.selected{k}.C;
 maps{k}=zeros(mn,n*z+2);
 for j=1:n,maps{k}(:,(j-1)*z+(1:z))=reshape(T.beta_nonref(:,j,:,k),mn,z)/n;end
end
[A,b,E,e,lb,ub]=weight_constraints(T,mr,K,r);
lp=optimoptions('linprog','Display','off','ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);
[x0,~,flag]=linprog(zeros(K*r+K,1),A,b,E,e,lb,ub,lp);assert(flag>0);
opt=optimoptions('quadprog','Algorithm','active-set','Display','off', ...
 'ConstraintTolerance',1e-11,'OptimalityTolerance',1e-12,'MaxIterations',10000);
P=size(L,2);personal=zeros(mn,K,P);audit=struct;W=zeros(r,K,P);epsilon=zeros(K,P);residual=zeros(P,2);
for p=1:P
 G=kron(eye(K)-ones(K,1)*L(:,p)',eye(mr))*scoreMap;
 H=zeros(K*r+K);H(1:K*r,1:K*r)=2*cfg.C3*eye(K*r)+2*cfg.C2*(G'*G);H=(H+H')/2;
 f=[zeros(K*r,1);-ones(K,1)];
 [x,~,flag,~,mu]=quadprog(H,f,A,b,E,e,lb,ub,x0,opt);assert(flag>0);
 primal=max([A*x-b;abs(E*x-e);lb-x;x-ub;0]);
 kkt=norm(H*x+f+A'*mu.ineqlin+E'*mu.eqlin-mu.lower+mu.upper,inf);
 assert(primal<1e-8 && kkt<1e-6);
 W(:,:,p)=reshape(x(1:K*r),r,K);epsilon(:,p)=x(K*r+1:end);residual(p,:)=[primal kkt];
 for k=1:K,personal(:,k,p)=maps{k}*T.selected{k}.X*W(:,k,p);end
end
assert(all(epsilon(:)>cfg.PositiveMarginCheckTolerance));
audit.Weights=W;audit.Epsilon2=epsilon;audit.PrimalKKT=residual;
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


function result=solve_Li_P2(beta,classID,options)
arguments
    beta double
    classID double
    options.MinGap (1,1) double = 1e-6
    options.Alpha (1,1) double = .5
    options.Beta (1,1) double = .5
end
[m,n,z]=size(beta);q=z-1;eta=options.MinGap;
classID=classID(:)';Hc=max(classID);
assert(Hc==3 && numel(classID)==m,'This implementation supports three classes.');
assert(all(ismember(1:3,classID)) && eta>0 && eta*q<1);
assert(abs(options.Alpha+options.Beta-1)<1e-12 && min([options.Alpha options.Beta])>=0);
assert(max(abs(sum(beta,3)-1),[],'all')<1e-10);
nv=n*q; V=zeros(m,nv);u0=zeros(m,1);scale=1-q*eta;
for j=1:n
    for t=1:q
        tail=sum(beta(:,j,t+1:z),3)/n;
        V(:,(j-1)*q+t)=scale*tail;
        u0=u0+eta*tail;
    end
end
G=zeros(2,m);Q=zeros(m);
for h=1:3
    ix=find(classID==h);nh=numel(ix);
    Q(ix,ix)=2*nh*(eye(nh)-ones(nh)/nh);
    if h<3
        iy=find(classID==h+1);G(h,iy)=1/numel(iy);G(h,ix)=-1/nh;
    end
end
L=G*V;g0=G*u0;eq=kron(eye(n),ones(1,q));
qp=optimoptions('quadprog','Display','off','Algorithm','interior-point-convex', ...
    'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-11,'MaxIterations',2000);
lp=optimoptions('linprog','Display','off','ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);
[xm,~,flag]=linprog([zeros(nv,1);-1],[-L ones(2,1)],g0, ...
    [eq zeros(n,1)],ones(n,1),zeros(nv+1,1),ones(nv+1,1),lp);
assert(flag>0,'Li classification infeasible or margin LP failed.');dMax=xm(end);
HQ=2*V'*Q*V;HQ=(HQ+HQ')/2;fQ=2*V'*Q*u0;
[xmin,~,flag]=quadprog(HQ,fQ,-L,g0,eq,ones(n,1),zeros(nv,1),ones(nv,1),[],qp);
assert(flag>0,'Li minimum-dispersion QP failed.');umin=u0+V*xmin;Dmin=umin'*Q*umin;
[Dmax,vertexInfo]=maximum_dispersion_vertices(V,u0,G,Q,n,q);
assert(Dmax>=Dmin-1e-8 && dMax>1e-10 && Dmax-Dmin>1e-10, ...
    'Degenerate objective range; original membership normalization is undefined.');
cd=options.Alpha/dMax;cg=options.Beta/(Dmax-Dmin);
H=zeros(nv+1);H(1:nv,1:nv)=cg*HQ;
f=[cg*fQ;-cd];A=[-L ones(2,1)];Aeq=[eq zeros(n,1)];
[x,~,flag,~,mult]=quadprog(H,f,A,g0,Aeq,ones(n,1),zeros(nv+1,1),ones(nv+1,1),[],qp);
assert(flag>0,'Li normalized P2 QP failed.');
primal=max([A*x-g0;abs(Aeq*x-1);-x;x-1;0]);
kkt=norm(H*x+f+A'*mult.ineqlin+Aeq'*mult.eqlin-mult.lower+mult.upper,inf);
assert(primal<1e-7 && kkt<1e-5);
eps=eta+scale*reshape(x(1:nv),q,n)';pns=[zeros(n,1) cumsum(eps,2)];
cs=u0+V*x(1:nv);D=cs'*Q*cs;
assert(max(abs(pns(:,end)-1))<1e-9 && min(diff(pns,1,2),[],'all')>=eta-1e-9);
assert(D>=Dmin-1e-7 && D<=Dmax+1e-7);
membership=[x(end)/dMax (Dmax-D)/(Dmax-Dmin)];
result=struct('PNS',pns,'CS',cs,'ClassID',classID,'Margin',x(end), ...
    'Dispersion',D,'MarginRange',[0 dMax],'DispersionRange',[Dmin Dmax], ...
    'EffectiveCoefficients',[cd cg],'Membership',membership, ...
    'Objective',options.Alpha*membership(1)+options.Beta*membership(2), ...
    'MinGap',eta,'Alpha',options.Alpha,'Beta',options.Beta, ...
    'PrimalError',primal,'KKTError',kkt,'VertexAudit',vertexInfo);
end

function [maximum,audit]=maximum_dispersion_vertices(V,u0,G,Q,n,q)
maximum=-Inf;bestScore=[];counts=zeros(1,4);tol=1e-9;
sets=choices(q,n);U=base_scores(V,u0,1:n,sets,q);consider(U,1);
pairs=nchoosek(1:q,2);triples=nchoosek(1:q,3);
for j=1:n
    others=setdiff(1:n,j);otherChoices=choices(q,numel(others));
    B=base_scores(V,u0,others,otherChoices,q);
    for a=1:size(pairs,1)
        ids=(j-1)*q+pairs(a,:);base=B+V(:,ids(1));direction=V(:,ids(2))-V(:,ids(1));
        gaps=G*base;dg=G*direction;
        for active=1:2
            if abs(dg(active))<1e-13,continue;end
            theta=-gaps(active,:)/dg(active);ok=theta>=-tol & theta<=1+tol;
            consider(base(:,ok)+direction*theta(ok),2);
        end
    end
    for a=1:size(triples,1)
        ids=(j-1)*q+triples(a,:);base=B+V(:,ids(1));directions=V(:,ids(2:3))-V(:,ids(1));
        A=G*directions;if abs(det(A))<1e-13,continue;end
        theta=-(A\(G*base));ok=all(theta>=-tol,1) & sum(theta,1)<=1+tol;
        consider(base(:,ok)+directions*theta(:,ok),3);
    end
end
attributes=nchoosek(1:n,2);
for jj=1:size(attributes,1)
    j=attributes(jj,1);k=attributes(jj,2);others=setdiff(1:n,[j k]);
    B=base_scores(V,u0,others,choices(q,numel(others)),q);
    for a=1:size(pairs,1)
        ia=(j-1)*q+pairs(a,:);
        for b=1:size(pairs,1)
            ib=(k-1)*q+pairs(b,:);base=B+V(:,ia(1))+V(:,ib(1));
            directions=[V(:,ia(2))-V(:,ia(1)) V(:,ib(2))-V(:,ib(1))];
            A=G*directions;if abs(det(A))<1e-13,continue;end
            theta=-(A\(G*base));ok=all(theta>=-tol & theta<=1+tol,1);
            consider(base(:,ok)+directions*theta(:,ok),4);
        end
    end
end
assert(isfinite(maximum));
audit=struct('Method','Exact vertex enumeration of faces of dimension 0/1/2', ...
    'FeasibleCandidateCounts',counts,'MaximizingReferenceScores',bestScore, ...
    'MaximumClassGapViolation',max([0;-G*bestScore]));
    function consider(U,kind)
        if isempty(U),return;end
        good=all(G*U>=-tol,1);U=U(:,good);counts(kind)=counts(kind)+size(U,2);
        if isempty(U),return;end
        values=sum(U.*(Q*U),1);[value,idx]=max(values);
        if value>maximum,maximum=value;bestScore=U(:,idx);end
    end
end

function U=base_scores(V,u0,attributes,sets,q)
U=repmat(u0,1,size(sets,1));
for a=1:numel(attributes),U=U+V(:,(attributes(a)-1)*q+sets(:,a));end
end

function sets=choices(q,n)
if n==0,sets=zeros(1,0);return;end
c=cell(1,n);[c{:}]=ndgrid(1:q);sets=zeros(numel(c{1}),n);
for a=1:n,sets(:,a)=c{a}(:);end
end

function result=solve_Xiao_transfer(br,bn,classes,L)
[m,n,z,K]=size(br);r=5;
assert(isequal(size(classes),[K m]) && size(L,1)==K);
opt=optimoptions('quadprog','Display','off','Algorithm','interior-point-convex', ...
    'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-11,'MaxIterations',2000);
rawCandidates=cell(K,1);fitCandidates=cell(K,1);
rawScore=cell(K,1);fitRefScore=cell(K,1);fitNonScore=cell(K,1);
modelAudit=zeros(K,r,4);fitSSE=0;fitDenominator=0;
fitRMSE=zeros(K,r);fitClassViolations=zeros(K,r);fitScoreError=zeros(K,r);
 for k=1:K
  beta=br(:,:,:,k);
  nv=m*n*z;U=zeros(m,nv);
  for i=1:m,for j=1:n,for t=1:z,U(i,sub2ind([m n z],i,j,t))=beta(i,j,t)/n;end,end,end
  center=class_center(classes(k,:));
  H=zeros(nv+1);H(1:nv,1:nv)=2*(U'*center*U+.01*eye(nv));
  A=zeros(2,nv+1);
  for h=1:2,d=mean(U(classes(k,:)==h+1,:),1)-mean(U(classes(k,:)==h,:),1);A(h,:)=[-d 1];end
  lb=zeros(nv+1,1);ub=[double(beta(:)>0);Inf];
  accepted=zeros(nv,r);acceptedScores=zeros(m,r);Qfit=zeros(n,z,r);SN=zeros(size(bn,1),r);SR=zeros(m,r);
  f0=[zeros(nv,1);-1];num=0;
  for trial=1:80
   rng(1000*k+trial,'twister');noise=randn(m,n,z);
   f=f0+[.01*.01*noise(:);0];
   [x,~,flag,~,mult]=quadprog(H,f,A,zeros(2,1),[],[],lb,ub,[],opt);
   if flag<=0,continue;end
   if num>0 && any(vecnorm(accepted(:,1:num)-x(1:nv),2,1)<1e-4),continue;end
   primal=max([A*x;lb-x;x-ub;0]);
   kkt=norm(H*x+f+A'*mult.ineqlin-mult.lower+mult.upper,inf);
   assert(primal<1e-7 && kkt<1e-5);
   num=num+1;accepted(:,num)=x(1:nv);acceptedScores(:,num)=U*x(1:nv);
   modelAudit(k,num,:)=[trial primal kkt .5*x'*H*x+f0'*x];
   raw=reshape(x(1:nv),m,n,z);fitted=zeros(n,z);err=0;den=0;
   for j=1:n
    b=squeeze(beta(:,j,:));pp=squeeze(raw(:,j,:));wt=sum(b,1);target=sum(b.*pp,1);
    observed=find(wt>1e-12);l=numel(observed);
    Am=diff(eye(l),1,1)*-1;
    [qo,~,fitflag]=quadprog(2*diag(wt(observed)),-2*target(observed)', ...
       Am,zeros(max(0,l-1),1),[],[],zeros(l,1),ones(l,1),[],opt);
    assert(fitflag>0);
    if l==1,qv=repmat(qo,1,z);else,qv=interp1(observed,qo,1:z,'linear');qv(1:observed(1)-1)=qo(1);qv(observed(end)+1:z)=qo(end);end
    assert(all(diff(qv)>=-1e-8));fitted(j,:)=qv;
    err=err+sum(b.*(pp-qv).^2,'all');den=den+sum(b,'all');
   end
   Qfit(:,:,num)=fitted;SN(:,num)=xiao_score_scale(bn(:,:,:,k),fitted);SR(:,num)=xiao_score_scale(beta,fitted);
   fitRMSE(k,num)=sqrt(err/den);fitSSE=fitSSE+err;fitDenominator=fitDenominator+den;
   fitScoreError(k,num)=max(abs(SR(:,num)-acceptedScores(:,num)));
   for h=1:2,fitClassViolations(k,num)=fitClassViolations(k,num)+double(mean(SR(classes(k,:)==h+1,num))<mean(SR(classes(k,:)==h,num))-1e-8);end
   if num==r,break;end
  end
  assert(num==r,'Could not obtain 5 Xiao candidates.');
  rawCandidates{k}=accepted;fitCandidates{k}=Qfit;rawScore{k}=acceptedScores;fitRefScore{k}=SR;fitNonScore{k}=SN;
 end
 scoreMap=blkdiag(rawScore{:});eq=kron(eye(K),ones(1,r));
 weights=zeros(r,K,5);MIV=zeros(size(bn,1),K,5);MIVref=zeros(m,K,5);MIVrefFit=zeros(m,K,5);
 weightAudit=zeros(5,2);
 for p=1:5
  lam=L(:,p);G=kron(eye(K)-ones(K,1)*lam',eye(m))*scoreMap;
  H=2*(G'*G+.001*eye(K*r));H=(H+H')/2;
  [x,~,flag,~,mult]=quadprog(H,zeros(K*r,1),[],[],eq,ones(K,1),zeros(K*r,1),ones(K*r,1),[],opt);
  assert(flag>0);primal=max([abs(eq*x-1);-x;x-1;0]);kkt=norm(H*x+eq'*mult.eqlin-mult.lower+mult.upper,inf);
  assert(primal<1e-7 && kkt<1e-5);weightAudit(p,:)=[primal kkt];
  weights(:,:,p)=reshape(x,r,K);
  for k=1:K
   MIV(:,k,p)=fitNonScore{k}*weights(:,k,p);MIVref(:,k,p)=rawScore{k}*weights(:,k,p);MIVrefFit(:,k,p)=fitRefScore{k}*weights(:,k,p);
  end
 end
O=struct('RawCandidates',{rawCandidates},'FittedCandidates',{fitCandidates}, ...
 'Weights',weights);
nv=n*z;
opt.MaxIterations=3000;
A=kron(eye(n),-diff(eye(z),1,1));b=zeros(n*(z-1),1);
projectionOpt=optimoptions('quadprog','Display','off','Algorithm','active-set', ...
 'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-11,'MaxIterations',10000);
 fitted=cell(K,1);nonScores=cell(K,1);refScores=cell(K,1);
 audit=zeros(K,r,7);oldErrors=zeros(m,K,r);newErrors=zeros(m,K,r);
 oldRef=cell(K,1);rawRef=cell(K,1);projectionDistances=zeros(K,r);
 fprintf('\nXiao reference-score transfer fit\n');
 for k=1:K
  M=score_map(br(:,:,:,k));Mn=score_map(bn(:,:,:,k));
  H=2*(M'*M);H=(H+H')/2;Qfit=zeros(n,z,r);SN=zeros(size(bn,1),r);SR=zeros(m,r);
  rawReference=zeros(m,r);oldReference=zeros(m,r);
  for v=1:r
   raw=reshape(O.RawCandidates{k}(:,v),m,n,z);
   target=squeeze(sum(sum(br(:,:,:,k).*raw,3),2))/n;
   q0=reshape(O.FittedCandidates{k}(:,:,v)',nv,1);
   [qmin,~,flag,~,mult]=quadprog(H,-2*M'*target,A,b,[],[],zeros(nv,1),ones(nv,1),[],opt);
   assert(flag>0,'Reference score LS failed.');
   firstError=max([A*qmin;-qmin;qmin-1;0]);
   firstKKT=norm(H*qmin-2*M'*target+A'*mult.ineqlin-mult.lower+mult.upper,inf);
   assert(firstError<1e-7 && firstKKT<1e-5);
   optimalPrediction=M*qmin;
   [q,~,flag,~,mult]=quadprog(2*eye(nv),-2*q0,A,b,M,optimalPrediction,zeros(nv,1),ones(nv,1),qmin,projectionOpt);
   assert(flag>0,'Closest optimal score-fit projection failed.');
   primal=max([A*q;abs(M*q-optimalPrediction);-q;q-1;0]);
   kkt=norm(2*q-2*q0+A'*mult.ineqlin+M'*mult.eqlin-mult.lower+mult.upper,inf);
   sse=sum((M*q-target).^2);firstSSE=sum((optimalPrediction-target).^2);
   oldSSE=sum((M*q0-target).^2);
   assert(primal<1e-7 && kkt<1e-5 && abs(sse-firstSSE)<1e-8 && sse<=oldSSE+1e-8);
   Qfit(:,:,v)=reshape(q,z,n)';SN(:,v)=Mn*q;SR(:,v)=M*q;
   rawReference(:,v)=target;oldReference(:,v)=M*q0;
   oldErrors(:,k,v)=M*q0-target;newErrors(:,k,v)=M*q-target;
   projectionDistances(k,v)=norm(q-q0);
   audit(k,v,:)=[firstError firstKKT primal kkt firstSSE sse oldSSE];
   fprintf('DM%d candidate%d reference RMSE %.9g -> %.9g, optimal-face residual %.3g\n', ...
     k,v,sqrt(oldSSE/m),sqrt(sse/m),norm(M*q-optimalPrediction,inf));
  end
  fitted{k}=Qfit;nonScores{k}=SN;refScores{k}=SR;oldRef{k}=oldReference;rawRef{k}=rawReference;
 end
 MIV=zeros(size(bn,1),K,5);refFinal=zeros(m,K,5);refOriginal=zeros(m,K,5);refOldFit=zeros(m,K,5);
 group=zeros(size(bn,1),5);ranks=zeros(size(group));CD=zeros(5,1);
 originalRefCD=zeros(5,1);oldFitRefCD=zeros(5,1);newFitRefCD=zeros(5,1);
 for p=1:5
  for k=1:K
   w=O.Weights(:,k,p);MIV(:,k,p)=nonScores{k}*w;
   refFinal(:,k,p)=refScores{k}*w;refOriginal(:,k,p)=rawRef{k}*w;refOldFit(:,k,p)=oldRef{k}*w;
  end
   group(:,p)=MIV(:,:,p)*L(:,p);ranks(:,p)=tiedrank(-group(:,p));
  CD(p)=sum((MIV(:,:,p)-group(:,p)).^2,'all');
  originalRefCD(p)=sum((refOriginal(:,:,p)-refOriginal(:,:,p)*L(:,p)).^2,'all');
  oldFitRefCD(p)=sum((refOldFit(:,:,p)-refOldFit(:,:,p)*L(:,p)).^2,'all');
  newFitRefCD(p)=sum((refFinal(:,:,p)-refFinal(:,:,p)*L(:,p)).^2,'all');
 end
result=struct('PersonalScores',MIV,'ReferenceScores',refFinal, ...
 'ReferenceScoresOriginal',refOriginal,'RawCandidates',{rawCandidates}, ...
 'PriorFittedCandidates',{fitCandidates},'FittedCandidates',{fitted}, ...
 'Weights',weights,'Classes',classes,'Parameters',[1 .01 .001], ...
 'NoiseStrength',.01,'CandidateSeedRule','1000*k+trial, twister', ...
 'ModelAudit',modelAudit,'WeightAudit',weightAudit,'ScoreFitAudit',audit, ...
 'ReferenceCandidateRMSE',[sqrt(mean(oldErrors.^2,'all')) sqrt(mean(newErrors.^2,'all'))], ...
 'ReferenceFinalRMSE',[sqrt(mean((refOldFit-refOriginal).^2,'all')) sqrt(mean((refFinal-refOriginal).^2,'all'))], ...
 'ReferenceCDOriginal',originalRefCD,'ReferenceCDFitted',newFitRefCD, ...
 'ProjectionDistances',projectionDistances, ...
 'TransferRule','');
end

function M=class_center(cid)
m=numel(cid);M=zeros(m);for h=1:3,ix=find(cid==h);M(ix,ix)=eye(numel(ix))-ones(numel(ix))/numel(ix);end
end

function cs=xiao_score_scale(beta,pns)
[m,n,z]=size(beta);cs=zeros(m,1);for i=1:m,for j=1:n,cs(i)=cs(i)+reshape(beta(i,j,:),1,z)*pns(j,:)'/n;end,end
end

function M=score_map(beta)
[m,n,z]=size(beta);M=zeros(m,n*z);for i=1:m,for j=1:n,M(i,(j-1)*z+(1:z))=reshape(beta(i,j,:),1,z)/n;end,end
end
