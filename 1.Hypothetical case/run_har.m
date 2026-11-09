clear;
clc;
base = fileparts(mfilename('fullpath'));
out = fullfile(base, 'results');
S = load(fullfile(out, 'stage1.mat'));
D = load(fullfile(out, 'data.mat'));
assert(isequal(S.cfg, D.cfg) && isequal(S.beta_ref, D.beta_ref), ...
    '输入数据与第一阶段数据不一致。');

cfg = D.cfg;
beta_ref = D.beta_ref;
models = S.models;
T = cfg.DMCount;
r = cfg.SelectedSampleCount;
assert(r == 5, '样本数量应为5。');
assert(cfg.HARChains * cfg.HARPerChain == r, 'HAR输出数量与样本数量不一致。');
H = cell(T, 1);
for k = 1:T
    Q = models{k};
    [X, ~, dim] = optimal_face_har(Q, cfg, k);
    assert(size(X, 2) == cfg.HARChains * cfg.HARPerChain);

    objective = -X(29, :) + cfg.C1 * sum((Q.B * X).^2, 1);
    objective_error = max(abs(objective - Q.Fstar));
    hard_error = max([Q.A * X - Q.b; Q.lb - X; X - Q.ub; ...
        zeros(1, size(X, 2))], [], 'all');
    equality_error = max(abs(Q.E * X - Q.q), [], 'all');

    assert(objective_error < cfg.FeasibilityTolerance);
    assert(hard_error < cfg.FeasibilityTolerance);
    assert(equality_error < cfg.FeasibilityTolerance);
    assert(min(Q.D * X, [], 'all') >= cfg.Delta - cfg.FeasibilityTolerance);
    assert(min(X(29, :)) > cfg.PositiveMarginCheckTolerance);

    H{k} = struct('X', X, 'C', Q.C * X);
    fprintf('DM%d：最优空间维数%d，HAR样本%d。\n', k, dim, size(X, 2));
end

selected = cell(T, 1);
for k = 1:T
    ids = 1:r;
    selected{k} = struct('IDs', ids(:), ...
        'X', H{k}.X(:, ids), 'C', H{k}.C(:, ids));
end
rows = [];
for k = 1:T
    Q = models{k};
    X = selected{k}.X;
    loss = -X(29, :) + cfg.C1 * sum((Q.B * X).^2, 1) - Q.Fstar;
    assert(max(abs(loss)) < cfg.FeasibilityTolerance);

    for v = 1:r
        for j = 1:cfg.Criteria
            rows(end + 1, :) = [k v selected{k}.IDs(v) j ...
                X((j - 1) * cfg.TermCount + (1:cfg.TermCount), v)'];
        end
    end
end

SampledPNS = array2table(rows, 'VariableNames', ...
    {'DM', 'Sample', 'HARSampleID', 'Criterion', ...
     'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6'});
writetable(SampledPNS, fullfile(out, '抽样PNS.xlsx'), ...
    'Sheet', 'SampledPNS', 'WriteMode', 'overwritesheet', 'UseExcel', false);
save(fullfile(out, 'har.mat'), 'cfg', 'beta_ref', 'selected', '-v7');

fprintf('每位决策者的%d组PNS已保存。\n', r);
function [X, chainids, dim] = optimal_face_har(Q, cfg, k)
    lp = optimoptions('linprog', 'Display', 'off', ...
        'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
    B = Q.Reduction;
    yy = zeros(size(B.A, 2), 0);
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
