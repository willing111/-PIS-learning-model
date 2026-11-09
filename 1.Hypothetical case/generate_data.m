clear;
clc;
base = fileparts(mfilename('fullpath'));
out = fullfile(base, 'results');
if ~isfolder(out)
    mkdir(out);
end
cfg = struct( ...
    'ReferenceCount', 6, 'NonreferenceCount', 6, ...
    'Criteria', 4, 'DMCount', 4, 'TermCount', 7, ...
    'ReferenceSeed', 42, 'NonreferenceSeed', 2026, ...
    'Generator', 'twister', ...
    'SupportProbabilities', [0.15 0.50 0.35], 'SupportSigma', 1.2, ...
    'SupportSizes', [2 3 4], 'TwoTermBase', 0.35, 'TwoTermWidth', 0.30, ...
    'ThreeTermBase', [0.25 0.50 0.25], ...
    'FourTermBase', [0.15 0.35 0.35 0.15], ...
    'JitterAmplitude', 0.03, 'RoundingDigits', 4, ...
    'FlipTwoTermProbability', 0.5, ...
    'Delta', 0.006, 'C1', 1, 'C2', 10, 'C3', 15, ...
    'MinimumWeight', 0, ...
    'DMWeights', ones(4, 1) / 4, ...
    'AttributeWeights', ones(4, 1) / 4, ...
    'HARChains', 5, 'HARBurnIn', 2000, ...
    'HARThin', 10, 'HARPerChain', 1, ...
    'HARSeedBase', [2026362000 2032691700 2026111800 2026113068], ...
    'SelectedSampleCount', 5, ...
    'FeasibilityTolerance', 1e-8, ...
    'PositiveMarginCheckTolerance', 1e-7, 'EventSignTolerance', 1e-12);
beta_ref = generate_beta(cfg.ReferenceSeed, cfg.ReferenceCount, cfg);
beta_nonref = generate_beta(cfg.NonreferenceSeed, cfg.NonreferenceCount, cfg);

reference_ids = compose('aR%d', (1:cfg.ReferenceCount)');
nonreference_ids = compose('b%d', (1:cfg.NonreferenceCount)');

assert(min(beta_ref, [], 'all') >= 0);
assert(min(beta_nonref, [], 'all') >= 0);
assert(max(abs(sum(beta_ref, 3) - 1), [], 'all') < 1e-12);
assert(max(abs(sum(beta_nonref, 3) - 1), [], 'all') < 1e-12);
ReferenceLDA = lda_table(beta_ref, reference_ids, cfg);
NonreferenceLDA = lda_table(beta_nonref, nonreference_ids, cfg);

file = fullfile(out, '输入数据.xlsx');
writetable(ReferenceLDA, file, 'Sheet', 'Reference_LDA', ...
    'WriteMode', 'overwritesheet', 'UseExcel', false);
writetable(NonreferenceLDA, file, 'Sheet', 'Nonreference_LDA', ...
    'WriteMode', 'overwritesheet', 'UseExcel', false);

save(fullfile(out, 'data.mat'), 'cfg', 'beta_ref', 'beta_nonref', ...
    'reference_ids', 'nonreference_ids', '-v7');

fprintf('语言分布评估已保存：%s\n', file);
function beta = generate_beta(seed, n_alternatives, cfg)
    rng(seed, 'twister');
    n_criteria = cfg.Criteria;
    n_dms = cfg.DMCount;
    z = cfg.TermCount;
    beta = zeros(n_alternatives, n_criteria, z, n_dms);

    p2 = cfg.SupportProbabilities(1);
    p3 = cfg.SupportProbabilities(2);

    for i = 1:n_alternatives
        for j = 1:n_criteria
            for k = 1:n_dms
                r = rand;
                if r < p2
                    s = 2;
                elseif r < p2 + p3
                    s = 3;
                else
                    s = 4;
                end

                start_candidates = 1:(z - s + 1);
                mid = (z + 1) / 2;
                center_start = round(mid - (s - 1) / 2);
                sigma = cfg.SupportSigma;
                w = exp(-((start_candidates - center_start).^2) / (2 * sigma^2));
                w = w / sum(w);
                start = randsample(start_candidates, 1, true, w);
                idx = start:(start + s - 1);

                switch s
                    case 2
                        u = cfg.TwoTermBase + cfg.TwoTermWidth * rand;
                        temp = [u; 1 - u];
                        if rand < cfg.FlipTwoTermProbability
                            temp = flipud(temp);
                        end
                    case 3
                        base = cfg.ThreeTermBase(:);
                        temp = base + cfg.JitterAmplitude * rand(3, 1);
                    case 4
                        base = cfg.FourTermBase(:);
                        temp = base + cfg.JitterAmplitude * rand(4, 1);
                end

                temp = temp / sum(temp);
                temp = round(temp, cfg.RoundingDigits);
                if sum(temp) == 0
                    temp = ones(s, 1) / s;
                end
                temp = temp / sum(temp);
                temp = round(temp, cfg.RoundingDigits);
                difference = 1 - sum(temp);
                temp(end) = temp(end) + difference;
                beta(i, j, idx, k) = temp;
            end
        end
    end
end

function T = lda_table(beta, ids, cfg)
    rows = zeros(numel(ids) * cfg.Criteria * cfg.DMCount, 10);
    labels = strings(size(rows, 1), 1);
    r = 0;

    for k = 1:cfg.DMCount
        for i = 1:numel(ids)
            for j = 1:cfg.Criteria
                r = r + 1;
                rows(r, :) = [k i j reshape(beta(i, j, :, k), 1, [])];
                labels(r) = ids(i);
            end
        end
    end

    T = array2table(rows, 'VariableNames', ...
        {'DM', 'Alternative', 'Criterion', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6'});
    T = addvars(T, labels, 'After', 'Alternative', ...
        'NewVariableNames', 'AlternativeID');
end
