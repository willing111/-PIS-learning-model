function Summary = run_C3_sensitivity(sourceRoot)
% C3敏感性分析：各环境完整运行阶段一/HAR，五个C3共用该环境的固定样本。
% 不绘图；输出原始共识偏差CD、非参考方案CSic/排序及epsilon2。
base = fileparts(mfilename('fullpath'));
if nargin < 1, sourceRoot = fullfile(fileparts(base), '1.Hypothetical case'); end
assert(isfolder(sourceRoot),'找不到模型源代码目录。');
if ~isfile(fullfile(sourceRoot,'results','data.mat'))
    isolated_run(fullfile(sourceRoot,'generate_data.m'));
end
input = load(fullfile(sourceRoot,'results','data.mat'));
input.cfg.HARChains = 5;
input.cfg.HARBurnIn = 2000;
input.cfg.HARThin = 10;
input.cfg.HARPerChain = 1;
input.cfg.HARSeedBase = [2026362000 2032691700 2026111800 2026113068];
input.cfg.SelectedSampleCount = 5;
c3Values = [5;10;15;20;25];
delta = .006; margin = 1e-6;
environments = struct('Name',{'A','B'},'C1',{1,2},'C2',{10,15});
assert(margin > input.cfg.PositiveMarginCheckTolerance);
assert(input.cfg.DMCount == 4 && input.cfg.NonreferenceCount == 6);
out = fullfile(base,'results','C3');
if ~isfolder(out), mkdir(out); end
Summary = struct('C3',c3Values,'Delta',delta,'EpsilonLowerBound',margin, ...
    'Environments',environments,'AlternativeIDs',string(input.nonreference_ids(:)), ...
    'SourceConfig',input.cfg,'R1',nan(2,4),'Epsilon1',nan(2,4), ...
    'Epsilon2',nan(5,4,2),'CSic',nan(5,6,2),'Ranking',nan(5,6,2), ...
    'CD',nan(5,2),'Objective',nan(5,2),'WeightSquareSum',nan(5,2), ...
    'Stage1OK',false(2,1),'HAROK',false(2,1),'Success',false(5,2), ...
    'Message',strings(5,2));
for e = 1:2
    caseDir = fullfile(out,'runs',['environment_' environments(e).Name]);
    if ~isfolder(fullfile(caseDir,'results')), mkdir(fullfile(caseDir,'results')); end
    fprintf('\n=== Environment %s: C1=%g, C2=%g, delta=%g ===\n', ...
        environments(e).Name,environments(e).C1,environments(e).C2,delta);
    attempted = "prepare";
    try
        data = input;
        data.cfg.C1 = environments(e).C1; data.cfg.Delta = delta;
        data.cfg.C3 = 15; data.cfg.C2 = environments(e).C2;
        save(fullfile(caseDir,'results','data.mat'),'-struct','data','-v7');
        first = fileread(fullfile(sourceRoot,'run_stage1.m'));
        expression = 'lb\s*=\s*\[zeros\(np,\s*1\);\s*[^;]+;\s*delta\];';
        assert(isscalar(regexp(first,expression,'match')),'无法定位eps1下界。');
        first = regexprep(first,expression,'lb = [zeros(np, 1); 1e-6; delta];');
        write_text(fullfile(caseDir,'run_stage1.m'),first);
        copyfile(fullfile(sourceRoot,'run_har.m'),fullfile(caseDir,'run_har.m'));
        attempted = "stage1";
        captured = evalc('isolated_run(fullfile(caseDir,''run_stage1.m''));');
        write_text(fullfile(caseDir,'stage1.log'),captured);
        st = load(fullfile(caseDir,'results','stage1.mat'));
        for k = 1:4
            Q = st.models{k};
            Summary.R1(e,k) = sum((Q.B*Q.Baseline).^2);
            Summary.Epsilon1(e,k) = Q.Baseline(29);
        end
        assert(all(Summary.Epsilon1(e,:) > input.cfg.PositiveMarginCheckTolerance));
        assert(all(Summary.Epsilon1(e,:) >= margin-input.cfg.FeasibilityTolerance));
        Summary.Stage1OK(e) = true;
        attempted = "HAR";
        captured = evalc('isolated_run(fullfile(caseDir,''run_har.m''));');
        write_text(fullfile(caseDir,'HAR.log'),captured);
        Summary.HAROK(e) = true;
        fprintf('  Stage1 and HAR PASS; eps1 = %s\n',mat2str(Summary.Epsilon1(e,:),9));
    catch ME
        Summary.Message(:,e) = attempted + ": " + string(ME.message);
        warning('C3Sensitivity:PreparationFailed','Environment %s: %s',environments(e).Name,ME.message);
        save(fullfile(out,'C3_sensitivity_summary.mat'),'Summary','-v7');
        continue;
    end
    for h = 1:numel(c3Values)
        runDir = fullfile(caseDir,['C3_' strrep(sprintf('%g',c3Values(h)),'.','p')]);
        if ~isfolder(fullfile(runDir,'results')), mkdir(fullfile(runDir,'results')); end
        try
            for name = {'data.mat','stage1.mat','har.mat'}
                copyfile(fullfile(caseDir,'results',name{1}),fullfile(runDir,'results',name{1}));
            end
            second = fileread(fullfile(sourceRoot,'run_stage2.m'));
            second = regexprep(second,'cfg\.C3\s*=\s*[^;]+;',sprintf('cfg.C3 = %g;',c3Values(h)));
            second = regexprep(second,'cfg\.C2\s*=\s*[^;]+;',sprintf('cfg.C2 = %g;',environments(e).C2));
            expression = 'lb\s*=\s*\[cfg\.MinimumWeight\s*\*\s*ones\(20,\s*1\);[^\]]+\];';
            assert(isscalar(regexp(second,expression,'match')),'无法定位eps2下界。');
            second = regexprep(second,expression, ...
                'lb = [cfg.MinimumWeight * ones(20, 1); 1e-6 * ones(4, 1)];');
            functionStart = 'function R = solve_reference_stage2';
            assert(contains(second,functionStart),'无法定位第二阶段结果保存位置。');
            capture = sprintf('save(fullfile(out,''sensitivity_stage2.mat''),''R'',''group_score'',''rank'',''-v7'');\n\n');
            second = strrep(second,functionStart,[capture functionStart]);
            write_text(fullfile(runDir,'run_stage2.m'),second);
            captured = evalc('isolated_run(fullfile(runDir,''run_stage2.m''));');
            write_text(fullfile(runDir,'stage2.log'),captured);
            result = load(fullfile(runDir,'results','sensitivity_stage2.mat'));
            assert(all(result.R.eps2 > input.cfg.PositiveMarginCheckTolerance));
            assert(all(result.R.eps2 >= margin-input.cfg.FeasibilityTolerance));
            Summary.CSic(h,:,e) = result.group_score(:)';
            Summary.Ranking(h,:,e) = result.rank(:)';
            Summary.Epsilon2(h,:,e) = result.R.eps2(:)';
            referenceScores = result.R.ReferenceScores;
            groupReference = referenceScores * input.cfg.DMWeights;
            Summary.CD(h,e) = sum((referenceScores-groupReference).^2,'all');
            Summary.WeightSquareSum(h,e) = sum(result.R.w.^2,'all');
            Summary.Objective(h,e) = result.R.Objective;
            direct = -sum(result.R.eps2) + environments(e).C2*Summary.CD(h,e) ...
                + c3Values(h)*Summary.WeightSquareSum(h,e);
            assert(abs(direct-result.R.Objective)<1e-8);
            assert(isequal(sort(Summary.Ranking(h,:,e)),1:6));
            Summary.Success(h,e) = true;
            fprintf('  C3=%g PASS; eps2=%s; CD=%.10g\n',c3Values(h), ...
                mat2str(result.R.eps2(:)',9),Summary.CD(h,e));
        catch ME
            Summary.Message(h,e) = "stage2: " + string(ME.message);
            warning('C3Sensitivity:SecondStageFailed','Environment %s C3=%g: %s', ...
                environments(e).Name,c3Values(h),ME.message);
        end
        save(fullfile(out,'C3_sensitivity_summary.mat'),'Summary','-v7');
    end
end
Summary.Epsilon2Total = reshape(sum(Summary.Epsilon2,2),numel(c3Values),2);
Summary.MinEpsilon2 = reshape(min(Summary.Epsilon2,[],2),numel(c3Values),2);
save(fullfile(out,'C3_sensitivity_summary.mat'),'Summary','-v7');
Summary = summarize_metrics(Summary,out,"C3",'C3_sensitivity_summary');
save(fullfile(out,'C3_sensitivity_summary.mat'),'Summary','-v7');
write_tables(Summary,out);
fprintf('\nCompleted: %d/%d combinations. No figures generated.\n',nnz(Summary.Success),numel(Summary.Success));
end

function isolated_run(scriptPath)
run(scriptPath);
end

function write_text(path,text)
fid = fopen(path,'w','n','UTF-8');
assert(fid>=0,'无法写入文件：%s',path);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',text);
end

function write_tables(S,out)
file = fullfile(out,'C3_sensitivity_summary.xlsx');
dmNames = {'C3','DM1','DM2','DM3','DM4'};
altNames = [{'C3'},cellstr(S.AlternativeIDs')];
for e = 1:2
    suffix = S.Environments(e).Name;
    for field = {'CSic','Ranking','Epsilon2'}
        if strcmp(field{1},'Epsilon2'), names=dmNames; else, names=altNames; end
        T = array2table([S.C3 S.(field{1})(:,:,e)],'VariableNames',names);
        writetable(T,file,'Sheet',[field{1} '_' suffix],'WriteMode','overwritesheet','UseExcel',false);
    end
    E2Table = table(S.C3,S.Epsilon2Total(:,e),S.CD(:,e),S.WeightSquareSum(:,e),S.MinEpsilon2(:,e), ...
        'VariableNames',{'C3','E2','CD','Rw','MinimumIndividualEpsilon2'});
    writetable(E2Table,file,'Sheet',['E2_' suffix],'WriteMode','overwritesheet','UseExcel',false);
    T = table(S.C3,repmat(S.Stage1OK(e),5,1),repmat(S.HAROK(e),5,1), ...
        S.Success(:,e),S.CD(:,e),S.Objective(:,e),S.Message(:,e), ...
        'VariableNames',{'C3','Stage1OK','HAROK','Stage2OK','CD','Objective','Message'});
    writetable(T,file,'Sheet',['Status_' suffix],'WriteMode','overwritesheet','UseExcel',false);
end
for field = {'R1','Epsilon1'}
    T = array2table([(1:2)' S.(field{1})],'VariableNames',{'Environment','DM1','DM2','DM3','DM4'});
    writetable(T,file,'Sheet',field{1},'WriteMode','overwritesheet','UseExcel',false);
end
writetable(struct2table(S.Environments),file,'Sheet','Environments','WriteMode','overwritesheet','UseExcel',false);
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
