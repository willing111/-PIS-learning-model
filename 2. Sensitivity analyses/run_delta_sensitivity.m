function Summary = run_delta_sensitivity()
% 依次计算C1=1、C1=2的完整两阶段delta敏感性分析。
% 两组结果分别保存，保持现有绘图程序的读取位置。
c1Values = [1 2];
caseSummaries = cell(numel(c1Values),1);
for index = 1:numel(c1Values)
    caseSummaries{index} = run_delta_case(c1Values(index));
end
Summary = vertcat(caseSummaries{:});
end

function Summary = run_delta_case(C1)
validateattributes(C1, {'numeric'}, {'scalar','>=',1,'<=',2});
base = fileparts(mfilename('fullpath'));
caseRoot = fullfile(fileparts(base), '1.Hypothetical case');
assert(isfolder(caseRoot),'找不到模型源代码目录。');
if ~isfile(fullfile(caseRoot,'results','data.mat'))
    run_case_script(fullfile(caseRoot,'generate_data.m'));
end
input = load(fullfile(caseRoot,'results','data.mat'));
input.cfg.HARChains = 5;
input.cfg.HARBurnIn = 2000;
input.cfg.HARThin = 10;
input.cfg.HARPerChain = 1;
input.cfg.HARSeedBase = [2026362000 2032691700 2026111800 2026113068];
input.cfg.SelectedSampleCount = 5;
input.cfg.C1 = C1;
input.cfg.C2 = 10;
input.cfg.C3 = 15;
epsilonLowerBound = 1e-6;
assert(epsilonLowerBound > input.cfg.PositiveMarginCheckTolerance);
deltaValues = [0.002; 0.004; 0.006; 0.008; 0.010];
assert(input.cfg.DMCount == 4 && input.cfg.NonreferenceCount == 6);

% 两组C1的计算结果分别保存，R1和MeanR1分别表示各DM指标及其均值。
out = fullfile(base,'results',sprintf('delta_C1_%g',C1));
if ~isfolder(out), mkdir(out); end
Summary = struct('C1',C1,'C2',10,'C3',15,'Delta',deltaValues, ...
    'EpsilonLowerBound',epsilonLowerBound,'BaseConfig',input.cfg, ...
    'AlternativeIDs',string(input.nonreference_ids(:)), ...
    'R1',nan(5,4),'Epsilon1',nan(5,4),'Epsilon2',nan(5,4), ...
    'CD',nan(5,1),'WeightSquareSum',nan(5,1),'CSic',nan(5,6),'Ranking',nan(5,6),'Success',false(5,1), ...
    'CompletedStage',repmat("pending",5,1),'Message',strings(5,1));
scriptNames = {'run_stage1.m','run_har.m','run_stage2.m'};
for h = 1:numel(deltaValues)
    delta = deltaValues(h);
    caseDir = fullfile(out,'runs',['delta_' strrep(sprintf('%.3f',delta),'.','p')]);
    if ~isfolder(fullfile(caseDir,'results')), mkdir(fullfile(caseDir,'results')); end
    fprintf('\n=== C1=%g, delta=%.3f, C2=10, C3=15 ===\n',C1,delta);
    attemptedStage = "prepare";
    try
        % 每个delta独立复制主流程，主算例的脚本和results均不会被覆盖。
        for j = 1:numel(scriptNames)
            text = fileread(fullfile(caseRoot,scriptNames{j}));
            if j == 1
                expression = 'lb\s*=\s*\[zeros\(np,\s*1\);\s*[^;]+;\s*delta\];';
                assert(numel(regexp(text,expression,'match')) == 1, ...
                    '无法定位第一阶段epsilon1下界，请检查主代码。');
                text = regexprep(text,expression,'lb = [zeros(np, 1); 1e-6; delta];');
            elseif j == 3
                expression = 'lb\s*=\s*\[cfg\.MinimumWeight\s*\*\s*ones\(20,\s*1\);[^\]]+\];';
                assert(numel(regexp(text,expression,'match')) == 1, ...
                    '无法定位第二阶段epsilon2下界，请检查主代码。');
                text = regexprep(text,expression, ...
                    'lb = [cfg.MinimumWeight * ones(20, 1); 1e-6 * ones(4, 1)];');
                text = regexprep(text,'cfg\.C2\s*=\s*[^;]+;','cfg.C2 = 10;');
                text = regexprep(text,'cfg\.C3\s*=\s*[^;]+;','cfg.C3 = 15;');
                functionStart = 'function R = solve_reference_stage2';
                assert(contains(text,functionStart),'无法定位第二阶段结果保存位置。');
                capture = sprintf(['save(fullfile(out,''sensitivity_stage2.mat''), ' ...
                    '''group_score'',''rank'',''R'',''-v7'');\n\n']);
                text = strrep(text,functionStart,[capture functionStart]);
            end
            write_script(fullfile(caseDir,scriptNames{j}),text);
        end
        data = input;
        data.cfg.Delta = delta;
        save(fullfile(caseDir,'results','data.mat'),'-struct','data','-v7');

        attemptedStage = "stage1";
        run_case_script(fullfile(caseDir,'run_stage1.m'));
        first = load(fullfile(caseDir,'results','stage1.mat'));
        for k = 1:4
            Q = first.models{k};
            Summary.R1(h,k) = sum((Q.B * Q.Baseline).^2);
            Summary.Epsilon1(h,k) = Q.Baseline(data.cfg.Criteria * data.cfg.TermCount + 1);
        end
        Summary.CompletedStage(h) = "stage1";

        attemptedStage = "HAR";
        run_case_script(fullfile(caseDir,'run_har.m'));
        Summary.CompletedStage(h) = "HAR";

        attemptedStage = "stage2";
        run_case_script(fullfile(caseDir,'run_stage2.m'));
        second = load(fullfile(caseDir,'results','sensitivity_stage2.mat'));
        Summary.CSic(h,:) = second.group_score(:)';
        Summary.Ranking(h,:) = second.rank(:)';
        Summary.Epsilon2(h,:) = second.R.eps2(:)';
        referenceScores = second.R.ReferenceScores;
        groupReference = referenceScores * data.cfg.DMWeights;
        Summary.CD(h,1) = sum((referenceScores-groupReference).^2,'all');
        Summary.WeightSquareSum(h,1) = sum(second.R.w.^2,'all');
        assert(isequal(sort(Summary.Ranking(h,:)),1:6));
        assert(all(Summary.Epsilon1(h,:) > data.cfg.PositiveMarginCheckTolerance));
        assert(all(Summary.Epsilon2(h,:) > data.cfg.PositiveMarginCheckTolerance));
        Summary.Success(h) = true;
        Summary.CompletedStage(h) = "complete";
    catch ME
        Summary.Message(h) = attemptedStage + ": " + string(ME.message);
        warning('DeltaSensitivity:CaseFailed','C1=%g, delta=%.3f: %s', ...
            C1,delta,char(Summary.Message(h)));
    end
    % 每组结束即保存状态；失败组不会被当成有效排名绘图。
    save(fullfile(out,'delta_stage1_summary.mat'),'Summary','-v7');
end
Summary = summarize_metrics(Summary,out,"Delta",'delta_stage1_summary');
save(fullfile(out,'delta_stage1_summary.mat'),'Summary','-v7');
write_summary_tables(Summary,out);
fprintf('\nC1=%g: %d/%d cases completed.\n',C1,nnz(Summary.Success),numel(deltaValues));
disp(array2table([deltaValues Summary.CSic], ...
    'VariableNames',[{'Delta'},cellstr(Summary.AlternativeIDs')]));
fprintf('Ranking: 1 = best; columns retain the original alternative order.\n');
disp(array2table([deltaValues Summary.Ranking], ...
    'VariableNames',[{'Delta'},cellstr(Summary.AlternativeIDs')]));
end

function run_case_script(scriptPath)
% 主脚本中的clear仅清理本辅助函数，不影响外层delta循环。
run(scriptPath);
end

function write_script(path,text)
fid = fopen(path,'w','n','UTF-8');
assert(fid >= 0,'无法写入实验脚本：%s',path);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',text);
end

function write_summary_tables(Summary,out)
file = fullfile(out,'delta_stage1_summary.xlsx');
dmNames = {'Delta','DM1','DM2','DM3','DM4'};
altNames = [{'Delta'},cellstr(Summary.AlternativeIDs')];
fields = {'R1','Epsilon1','Epsilon2','CSic','Ranking'};
for j = 1:numel(fields)
    field = fields{j};
    if j <= 3, names = dmNames; else, names = altNames; end
    T = array2table([Summary.Delta Summary.(field)],'VariableNames',names);
    writetable(T,file,'Sheet',field,'WriteMode','overwritesheet','UseExcel',false);
end
Status = table(Summary.Delta,Summary.Success,Summary.CompletedStage,Summary.Message, ...
    'VariableNames',{'Delta','Success','CompletedStage','Message'});
writetable(Status,file,'Sheet','Status','WriteMode','overwritesheet','UseExcel',false);
end

function S = summarize_metrics(S,out,parameter,stem)
% MeanR1 is the arithmetic mean over four decision makers, without C1 scaling.
S.MeanR1 = mean(S.R1,2);
S.MeanEpsilon1 = mean(S.Epsilon1,2);
n = numel(S.(parameter));
count = size(S.CSic,3);
S.MeanEpsilon2 = reshape(mean(S.Epsilon2,2),n,count);
S.MeanCD = S.CD / 4;
S.MeanR2 = S.WeightSquareSum / 4;
S.MeanCSic = reshape(mean(S.CSic,2),n,count);
S.MeanResults = cell(count,1);
for e = 1:count
    if parameter == "Delta"
        c1 = repmat(S.C1,n,1);
        c2 = repmat(S.C2,n,1);
        c3 = repmat(S.C3,n,1);
        delta = S.Delta(:);
        r1 = S.MeanR1;
        ep1 = S.MeanEpsilon1;
        suffix = 'All';
    elseif parameter == "C1"
        c1 = S.C1(:);
        c2 = repmat(S.Environments(e).C2,n,1);
        c3 = repmat(S.Environments(e).C3,n,1);
        delta = repmat(S.Delta,n,1);
        r1 = S.MeanR1;
        ep1 = S.MeanEpsilon1;
        suffix = S.Environments(e).Name;
    else
        c1 = repmat(S.Environments(e).C1,n,1);
        delta = repmat(S.Delta,n,1);
        if parameter == "C2"
            c2 = S.C2(:);
            c3 = repmat(S.Environments(e).C3,n,1);
        else
            c2 = repmat(S.Environments(e).C2,n,1);
            c3 = S.C3(:);
        end
        r1 = repmat(S.MeanR1(e),n,1);
        ep1 = repmat(S.MeanEpsilon1(e),n,1);
        suffix = S.Environments(e).Name;
    end
    ok = S.Success(:,e);
    ep2 = S.MeanEpsilon2(:,e);
    cd = S.MeanCD(:,e);
    r2 = S.MeanR2(:,e);
    cs = S.MeanCSic(:,e);
    scores = S.CSic(:,:,e);
    ep2(~ok) = NaN; cd(~ok) = NaN; r2(~ok) = NaN;
    cs(~ok) = NaN; scores(~ok,:) = NaN;
    S.MeanEpsilon2(:,e) = ep2;
    S.MeanCD(:,e) = cd;
    S.MeanR2(:,e) = r2;
    S.MeanCSic(:,e) = cs;
    T = table(delta,c1,c2,c3,ok,ep1,r1,ep2,cd,r2,cs, ...
        'VariableNames',{'Delta','C1','C2','C3','Success','MeanEpsilon1', ...
        'MeanR1','MeanEpsilon2','MeanCD','MeanR2','MeanCSic'});
    for i = 1:numel(S.AlternativeIDs)
        T.(['CSic_' char(S.AlternativeIDs(i))]) = scores(:,i);
    end
    S.MeanResults{e} = T;
    fprintf('\nEnvironment %s: mean indicators and CSic in original alternative order.\n',suffix);
    fprintf('MeanR1 denotes R1-bar; means over decision makers use divisor 4.\n');
    disp(T);
    writetable(T,fullfile(out,[char(stem) '.xlsx']),'Sheet',['Means_' suffix], ...
        'WriteMode','overwritesheet','UseExcel',false);
end
end
