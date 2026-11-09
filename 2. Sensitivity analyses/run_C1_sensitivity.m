function Summary = run_C1_sensitivity()
% C1敏感性分析：同一个C1先运行一次阶段一和HAR，两组环境共用样本。
base = fileparts(mfilename('fullpath'));
sourceRoot = fullfile(fileparts(base), '1.Hypothetical case');
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
c1Values = [0.6; 1.0; 1.4; 1.8; 2.2];
delta = 0.006;
margin = 1e-6;
environments = struct('Name',{'A','B'},'C2',{10,15},'C3',{15,20});
assert(margin > input.cfg.PositiveMarginCheckTolerance);
assert(input.cfg.DMCount == 4 && input.cfg.NonreferenceCount == 6);
out = fullfile(base,'results','C1');
if ~isfolder(out), mkdir(out); end
Summary = struct('C1',c1Values,'Delta',delta,'EpsilonLowerBound',margin, ...
    'Environments',environments,'AlternativeIDs',string(input.nonreference_ids(:)), ...
    'SourceConfig',input.cfg,'R1',nan(5,4),'Epsilon1',nan(5,4), ...
    'Epsilon2',nan(5,4,2),'CD',nan(5,2),'WeightSquareSum',nan(5,2), ...
    'CSic',nan(5,6,2),'Ranking',nan(5,6,2), ...
    'Stage1OK',false(5,1),'HAROK',false(5,1),'Success',false(5,2), ...
    'Message',strings(5,2));
for h = 1:numel(c1Values)
    caseDir = fullfile(out,'runs',['C1_' strrep(sprintf('%.1f',c1Values(h)),'.','p')]);
    if ~isfolder(fullfile(caseDir,'results')), mkdir(fullfile(caseDir,'results')); end
    fprintf('\n=== C1=%.1f, delta=%.3f ===\n',c1Values(h),delta);
    attempted = "prepare";
    try
        data = input;
        data.cfg.C1 = c1Values(h); data.cfg.Delta = delta;
        data.cfg.C2 = environments(1).C2; data.cfg.C3 = environments(1).C3;
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
            Summary.R1(h,k) = sum((Q.B * Q.Baseline).^2);
            Summary.Epsilon1(h,k) = Q.Baseline(29);
        end
        Summary.Stage1OK(h) = true;
        attempted = "HAR";
        captured = evalc('isolated_run(fullfile(caseDir,''run_har.m''));');
        write_text(fullfile(caseDir,'HAR.log'),captured);
        Summary.HAROK(h) = true;
    catch ME
        Summary.Message(h,:) = attempted + ": " + string(ME.message);
        warning('C1Sensitivity:FirstStageFailed','C1=%.1f: %s',c1Values(h),char(Summary.Message(h,1)));
        save(fullfile(out,'C1_sensitivity_summary.mat'),'Summary','-v7');
        continue;
    end
    for e = 1:2
        envDir = fullfile(caseDir,['environment_' environments(e).Name]);
        if ~isfolder(fullfile(envDir,'results')), mkdir(fullfile(envDir,'results')); end
        try
            % 两组环境复制并共用同一份HAR样本。
            for name = {'data.mat','stage1.mat','har.mat'}
                copyfile(fullfile(caseDir,'results',name{1}),fullfile(envDir,'results',name{1}));
            end
            second = fileread(fullfile(sourceRoot,'run_stage2.m'));
            second = regexprep(second,'cfg\.C2\s*=\s*[^;]+;',sprintf('cfg.C2 = %g;',environments(e).C2));
            second = regexprep(second,'cfg\.C3\s*=\s*[^;]+;',sprintf('cfg.C3 = %g;',environments(e).C3));
            expression = 'lb\s*=\s*\[cfg\.MinimumWeight\s*\*\s*ones\(20,\s*1\);[^\]]+\];';
            assert(isscalar(regexp(second,expression,'match')),'无法定位eps2下界。');
            second = regexprep(second,expression, ...
                'lb = [cfg.MinimumWeight * ones(20, 1); 1e-6 * ones(4, 1)];');
            functionStart = 'function R = solve_reference_stage2';
            assert(contains(second,functionStart),'无法定位第二阶段结果保存位置。');
            capture = sprintf('save(fullfile(out,''sensitivity_stage2.mat''),''R'',''group_score'',''rank'',''-v7'');\n\n');
            second = strrep(second,functionStart,[capture functionStart]);
            write_text(fullfile(envDir,'run_stage2.m'),second);
            captured = evalc('isolated_run(fullfile(envDir,''run_stage2.m''));');
            write_text(fullfile(envDir,'stage2.log'),captured);
            result = load(fullfile(envDir,'results','sensitivity_stage2.mat'));
            Summary.CSic(h,:,e) = result.group_score(:)';
            Summary.Ranking(h,:,e) = result.rank(:)';
            Summary.Epsilon2(h,:,e) = result.R.eps2(:)';
            referenceScores = result.R.ReferenceScores;
            groupReference = referenceScores * data.cfg.DMWeights;
            Summary.CD(h,e) = sum((referenceScores-groupReference).^2,'all');
            Summary.WeightSquareSum(h,e) = sum(result.R.w.^2,'all');
            assert(isequal(sort(Summary.Ranking(h,:,e)),1:6));
            Summary.Success(h,e) = true;
            fprintf('  Environment %s: C2=%g, C3=%g, PASS; min eps2=%.8g\n', ...
                environments(e).Name,environments(e).C2,environments(e).C3,min(result.R.eps2));
        catch ME
            Summary.Message(h,e) = "stage2: " + string(ME.message);
            warning('C1Sensitivity:SecondStageFailed','C1=%.1f, environment %s: %s', ...
                c1Values(h),environments(e).Name,ME.message);
        end
        save(fullfile(out,'C1_sensitivity_summary.mat'),'Summary','-v7');
    end
end
Summary = summarize_metrics(Summary,out,"C1",'C1_sensitivity_summary');
save(fullfile(out,'C1_sensitivity_summary.mat'),'Summary','-v7');
write_tables(Summary,out);
fprintf('\nCompleted: %d/%d C1-environment combinations.\n',nnz(Summary.Success),numel(Summary.Success));
for e = 1:2
    fprintf('Environment %s CSic (original alternative order):\n',environments(e).Name);
    disp(array2table([c1Values Summary.CSic(:,:,e)], ...
        'VariableNames',[{'C1'},cellstr(Summary.AlternativeIDs')]));
end

end

function isolated_run(scriptPath)
run(scriptPath);
end

function write_text(path,text)
fid = fopen(path,'w','n','UTF-8');
assert(fid >= 0,'无法写入文件：%s',path);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',text);
end

function write_tables(Summary,out)
file = fullfile(out,'C1_sensitivity_summary.xlsx');
dmNames = {'C1','DM1','DM2','DM3','DM4'};
altNames = [{'C1'},cellstr(Summary.AlternativeIDs')];
for field = {'R1','Epsilon1'}
    T = array2table([Summary.C1 Summary.(field{1})],'VariableNames',dmNames);
    writetable(T,file,'Sheet',field{1},'WriteMode','overwritesheet','UseExcel',false);
end
for e = 1:2
    suffix = Summary.Environments(e).Name;
    for field = {'CSic','Ranking','Epsilon2'}
        if strcmp(field{1},'Epsilon2'), names = dmNames; else, names = altNames; end
        T = array2table([Summary.C1 Summary.(field{1})(:,:,e)],'VariableNames',names);
        writetable(T,file,'Sheet',[field{1} '_' suffix],'WriteMode','overwritesheet','UseExcel',false);
    end
end
Status = table(Summary.C1,Summary.Stage1OK,Summary.HAROK,Summary.Success(:,1), ...
    Summary.Success(:,2),Summary.Message(:,1),Summary.Message(:,2), ...
    'VariableNames',{'C1','Stage1OK','HAROK','EnvironmentA_OK','EnvironmentB_OK','MessageA','MessageB'});
writetable(Status,file,'Sheet','Status','WriteMode','overwritesheet','UseExcel',false);
Parameters = struct2table(Summary.Environments);
writetable(Parameters,file,'Sheet','Environments','WriteMode','overwritesheet','UseExcel',false);
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
