clear;
clc;

base = fileparts(mfilename('fullpath'));
out = fullfile(base, 'results');
D = load(fullfile(out, 'data.mat'));
S = load(fullfile(out, 'stage1.mat'));
V = load(fullfile(out, 'har.mat'));
assert(isequal(D.cfg, S.cfg) && isequal(D.cfg, V.cfg));
assert(isequal(D.beta_ref, S.beta_ref) && isequal(D.beta_ref, V.beta_ref), ...
    '各阶段的参考方案数据不一致。');

cfg = D.cfg;
cfg.MinimumWeight = 0;
cfg.C2 = 10;
cfg.C3 = 15;
selected = V.selected;
m = cfg.NonreferenceCount;
n = cfg.Criteria;
T = cfg.DMCount;
r = cfg.SelectedSampleCount;
assert(r == 5, '样本数量应为5。');

altpref = [3 5 6];

pairpref = [1 5 1 1 4;
            2 3 6 6 5;
            4 1 2 2 3];

for k = 1:T
    C = zeros(cfg.ReferenceCount, r);
    for i = 1:cfg.ReferenceCount
        for j = 1:n
            p = selected{k}.X((j - 1) * cfg.TermCount + (1:cfg.TermCount), :);
            C(i, :) = C(i, :) + reshape(D.beta_ref(i, j, :, k), 1, []) * p / n;
        end
    end
    assert(max(abs(C - selected{k}.C), [], 'all') < 1e-12);
    selected{k}.C = C;

    Q = S.models{k};
    X = selected{k}.X;
    assert(max(abs(-X(29, :) + cfg.C1 * sum((Q.B * X).^2, 1) - Q.Fstar)) < 1e-8);
    assert(Q.Baseline(29) > cfg.PositiveMarginCheckTolerance);
end

R = solve_reference_stage2(selected, altpref, pairpref, cfg);
assert(all(R.eps2 > cfg.PositiveMarginCheckTolerance), ...
    'epsilon2未满足正值要求。');

UNR = zeros(m, n, T);
CSNR = zeros(m, T);
averagePNS = zeros(n, cfg.TermCount, T);
conversion_rows = cell(0, 4);
mapping_error = 0;

for k = 1:T
    for j = 1:n
        p = selected{k}.X((j - 1) * cfg.TermCount + (1:cfg.TermCount), :);
        averagePNS(j, :, k) = (p * R.w(:, k))';
        for i = 1:m
            sample_u = reshape(D.beta_nonref(i, j, :, k), 1, []) * p;
            UNR(i, j, k) = sample_u * R.w(:, k);
            via_average = reshape(D.beta_nonref(i, j, :, k), 1, []) ...
                * averagePNS(j, :, k)';
            mapping_error = max(mapping_error, abs(via_average - UNR(i, j, k)));
            conversion_rows(end + 1, :) = {k, D.nonreference_ids(i), ...
                j, UNR(i, j, k)};
        end
    end
    CSNR(:, k) = UNR(:, :, k) * cfg.AttributeWeights;
end

assert(mapping_error < 1e-12);
assert(min(UNR, [], 'all') >= -cfg.FeasibilityTolerance);
assert(max(UNR, [], 'all') <= 1 + cfg.FeasibilityTolerance);

group_score = CSNR * cfg.DMWeights;
[~, order] = sort(group_score, 'descend');
rank = zeros(m, 1);
rank(order) = 1:m;

Weights = array2table([(1:T)' R.w' R.eps2], 'VariableNames', ...
    {'DM', 'w1', 'w2', 'w3', 'w4', 'w5', 'Epsilon2'});
ReferenceScores = table(D.reference_ids, ...
    R.ReferenceScores(:, 1), R.ReferenceScores(:, 2), ...
    R.ReferenceScores(:, 3), R.ReferenceScores(:, 4), ...
    'VariableNames', {'Alternative', 'DM1', 'DM2', 'DM3', 'DM4'});
Conversion = cell2table(conversion_rows, 'VariableNames', ...
    {'DM', 'Alternative', 'Criterion', 'NumericalAssessment'});
Scores = table(D.nonreference_ids, group_score, rank, ...
    'VariableNames', {'Alternative', 'CScic', 'Rank'});
Scores = sortrows(Scores, 'Rank');

file = fullfile(out, '第二阶段结果.xlsx');
writetable(Weights, file, 'Sheet', 'Weights', ...
    'WriteMode', 'overwritesheet', 'UseExcel', false);
writetable(ReferenceScores, file, 'Sheet', 'Reference_Scores', ...
    'WriteMode', 'overwritesheet', 'UseExcel', false);
writetable(Conversion, file, 'Sheet', 'Nonreference_U', ...
    'WriteMode', 'overwritesheet', 'UseExcel', false);
writetable(Scores, file, 'Sheet', 'Nonreference_Scores', ...
    'WriteMode', 'overwritesheet', 'UseExcel', false);

fprintf('PNS权重及epsilon2：\n');
disp(Weights);
fprintf('参考方案的个人综合得分（DM1-DM4）：\n');
disp(ReferenceScores);
fprintf('非参考方案的群体综合得分CS_i^c及排名：\n');
fprintf('%-15s %8s %6s\n', 'Alternative', 'CScic', 'Rank');
for i = 1:height(Scores)
    fprintf('%-15s %8.4f %6d\n', char(string(Scores.Alternative(i))), ...
        Scores.CScic(i), Scores.Rank(i));
end
fprintf('第二阶段目标函数值：%.10f\n', R.Objective);

function R = solve_reference_stage2(selected, altpref, pairpref, cfg)
    C2 = cfg.C2;
    C3 = cfg.C3;
    dm_weights = cfg.DMWeights;
    tol = cfg.EventSignTolerance;

    scoreMap = zeros(24, 20);
    for k = 1:4
        scoreMap((k - 1) * 6 + (1:6), (k - 1) * 5 + (1:5)) = selected{k}.C;
    end
    consensusMap = kron(eye(4) - ones(4, 1) * dm_weights', eye(6)) * scoreMap;

    H = diag([2 * C3 * ones(20, 1); zeros(4, 1)]);
    H(1:20, 1:20) = H(1:20, 1:20) + 2 * C2 * (consensusMap' * consensusMap);
    H = (H + H') / 2;
    f = [zeros(20, 1); -ones(4, 1)];

    Aeq = [kron(eye(4), ones(1, 5)), zeros(4)];
    beq = ones(4, 1);
   lb = [cfg.MinimumWeight * ones(20, 1); 1e-6 * ones(4, 1)];
    ub = ones(24, 1);

    A = zeros(0, 24);
    for k = 1:4
        C = selected{k}.C;
        ar = altpref(altpref(:, 1) == k, :);
        if ~isempty(ar)
            a = ar(2:3);
            d = C(a(1), :) - C(a(2), :);
            left = double(d > tol);
            right = double(d < -tol);
            row = zeros(1, 24);
            row((k - 1) * 5 + (1:5)) = -(left - right);
            row(20 + k) = 1;
            A(end + 1, :) = row;
        end

        pr = pairpref(pairpref(:, 1) == k, :);
        if ~isempty(pr)
            left = double(C(pr(2), :) - C(pr(3), :) > tol);
            right = double(C(pr(4), :) - C(pr(5), :) > tol);
            row = zeros(1, 24);
            row((k - 1) * 5 + (1:5)) = -(left - right);
            row(20 + k) = 1;
            A(end + 1, :) = row;
        end
        assert(~isempty(ar) || ~isempty(pr), '决策者的偏好输入为空。');
    end
    b = zeros(size(A, 1), 1);

    lp = optimoptions('linprog', 'Display', 'off', ...
        'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
    qp = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'active-set', ...
        'ConstraintTolerance', 1e-11, 'OptimalityTolerance', 1e-12, ...
        'MaxIterations', 10000);
    [x0, ~, lp_flag] = linprog(zeros(24, 1), A, b, Aeq, beq, lb, ub, lp);
    assert(lp_flag > 0, '第二阶段约束不可行。');
    [x, objective, flag, ~, lambda] = quadprog(H, f, A, b, Aeq, beq, lb, ub, x0, qp);
    assert(flag > 0, '第二阶段求解失败。');

    violation = max([A * x - b; abs(Aeq * x - beq); lb - x; x - ub; 0]);
    kkt = norm(H * x + f + A' * lambda.ineqlin + Aeq' * lambda.eqlin ...
        - lambda.lower + lambda.upper, inf);
    assert(violation < cfg.FeasibilityTolerance && kkt < 1e-7);

    w = reshape(x(1:20), 5, 4);
    epsilon2 = x(21:24);
    scores = zeros(6, 4);
    for k = 1:4
        scores(:, k) = selected{k}.C * w(:, k);
    end
    group = scores * dm_weights;
    direct_objective = -sum(epsilon2) ...
        + C2 * sum((scores - group).^2, 'all') + C3 * sum(w.^2, 'all');
    assert(abs(direct_objective - objective) < 1e-9);

    R = struct('w', w, 'eps2', epsilon2, 'Objective', objective, ...
        'ReferenceScores', scores);
end
