clear;
clc;
base = fileparts(mfilename('fullpath'));
out = fullfile(base, 'results');
load(fullfile(out, 'data.mat'), 'cfg', 'beta_ref');
cfg.HARSeedBase = [2026362000 2032691700 2026111800 2026113068];

m = cfg.ReferenceCount;
n = cfg.Criteria;
T = cfg.DMCount;
z = cfg.TermCount;
C1 = cfg.C1;
delta = cfg.Delta;
lp = optimoptions('linprog', 'Display', 'off', ...
    'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
qp = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'active-set', ...
    'ConstraintTolerance', 1e-11, 'OptimalityTolerance', 1e-12, ...
    'MaxIterations', 10000);

models = cell(T, 1);
summary = zeros(T, 3);
for k = 1:T
    np = n * z;
    N = np + 2;
    ie = np + 1;
    id = np + 2;
    e = zeros(1, N);
    e(ie) = 1;
    ed = zeros(1, N);
    ed(id) = 1;
    D = zeros(n * (z - 1), N);
    for j = 1:n
        for t = 1:z - 1
            r = (j - 1) * (z - 1) + t;
            D(r, (j - 1) * z + t) = -1;
            D(r, (j - 1) * z + t + 1) = 1;
        end
    end
    U = zeros(m * n, N);
    C = zeros(m, N);
    for i = 1:m
        for j = 1:n
            r = (i - 1) * n + j;
            U(r, (j - 1) * z + (1:z)) = reshape(beta_ref(i, j, :, k), 1, z);
            C(i, :) = C(i, :) + U(r, :) / n;
        end
    end
    switch k
        case 1
            G = C(4, :) - C(5, :);
        case 2
            left = C(2, :) - C(3, :);
            right = C(5, :) - C(2, :);
            G = [left - right; right];
        case 3
            G = C(3, :) - C(1, :);
        case 4
            G = C(2, :) - C(6, :);
    end
    B = [kron(eye(n) - ones(n) / n, eye(z)), zeros(np, 2)];
    H = 2 * C1 * (B' * B);
    f = -e';
    A = [-D + repmat(ed, size(D, 1), 1); ...
         U; -U; C; -C; ...
         -G + repmat(e, size(G, 1), 1)];
    b = [zeros(size(D, 1), 1); ...
         ones(m * n, 1); zeros(m * n, 1); ...
         ones(m, 1); zeros(m, 1); zeros(size(G, 1), 1)];
    lb = [zeros(np, 1); 1e-6; delta];
    ub = [ones(np, 1); inf; 1];
    Aeq = ed;
    beq = delta;
    [x0, ~, lp_flag] = linprog(f, A, b, Aeq, beq, lb, ub, lp);
    assert(lp_flag > 0, 'DM%d：约束不可行或目标无界。', k);
    [x, ~, flag, ~, lambda] = quadprog(H, f, A, b, Aeq, beq, lb, ub, x0, qp);
    assert(flag > 0, 'DM%d：第一阶段求解失败。', k);

    Fstar = -x(ie) + C1 * sum((B * x).^2);
    violation = max([A * x - b; abs(Aeq * x - beq); lb - x; x - ub; 0]);
    kkt = norm(H * x + f + A' * lambda.ineqlin + Aeq' * lambda.eqlin ...
        - lambda.lower + lambda.upper, inf);
    assert(violation < cfg.FeasibilityTolerance && kkt < 1e-7);
    assert(x(ie) > cfg.PositiveMarginCheckTolerance, ...
        'DM%d：epsilon1未满足正值要求。', k);
    E = [Aeq; B; e];
    q = [beq; B * x; x(ie)];
    Z = null(E, 1e-10);
    I = eye(N);
    il = isfinite(lb);
    iu = isfinite(ub);
    allA = [A; -I(il, :); I(iu, :)];
    allb = [b; -lb(il); ub(iu)];
    reducedA = allA * Z;
    reducedb = allb - allA * x;
    norms = vecnorm(reducedA, 2, 2);
    keep = norms > 1e-10;
    assert(min(reducedb) > -cfg.FeasibilityTolerance);
    reducedb(abs(reducedb) < 1e-10) = 0;

    reduction = struct('Origin', x, 'Z', Z, ...
        'A', reducedA(keep, :) ./ norms(keep), ...
        'b', reducedb(keep) ./ norms(keep));
    models{k} = struct('DM', k, 'A', A, 'b', b, ...
        'Aeq', Aeq, 'beq', beq, 'lb', lb, 'ub', ub, ...
        'B', B, 'H', H, 'f', f, 'D', D, 'U', U, 'C', C, 'G', G, ...
        'Baseline', x, 'Fstar', Fstar, 'E', E, 'q', q, 'Reduction', reduction);
    summary(k, :) = [k Fstar x(ie)];
end
Stage1Summary = array2table(summary, ...
    'VariableNames', {'DM', 'OptimalObjective', 'Epsilon1'});
writetable(Stage1Summary, fullfile(out, '第一阶段结果.xlsx'), ...
    'Sheet', 'Summary', 'WriteMode', 'overwritesheet', 'UseExcel', false);
save(fullfile(out, 'stage1.mat'), 'cfg', 'beta_ref', 'models', '-v7');
save(fullfile(out, 'data.mat'), 'cfg', '-append');
disp(Stage1Summary);
