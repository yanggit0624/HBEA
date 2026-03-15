classdef HBEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Hyperbolic Busemann Evolutionary Algorithm (HBEA)

    methods
        function main(Algorithm, Problem)
            %% Initialization
            Population = Problem.Initialization();
            [Population, FrontNo, Iso, BInfo] = EnvironmentalSelection(Population, Problem.N);

            %% Optimization loop
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population, FrontNo, Iso, BInfo);
                Offspring  = OperatorGA(Problem, Population(MatingPool));
                [Population, FrontNo, Iso, BInfo] = EnvironmentalSelection([Population, Offspring], Problem.N);
            end
        end
    end
end

%% ================================================================
%  Mating selection
%% ================================================================
function MatingPool = MatingSelection(Population, FrontNo, Iso, BInfo)
    EPS = 1e-12;
    N   = length(Population);
    Delta = BInfo.Delta;
    C     = BInfo.C;

    % Parent 1: Tournament selection by rank and isolation
    Parent1 = TournamentSelection(2, N, FrontNo, -Iso);

    % Parent 2: Complementary pairing in Delta-space
    nrm = max(sqrt(sum(Delta.^2, 2)), EPS);
    U   = Delta ./ nrm;
    U1  = U(Parent1, :);

    % Cosine similarity matrix
    S = U1 * U';
    S = max(-1, min(1, S));

    % Exclude self-pairing
    idx_lin = sub2ind([N, N], (1:N)', Parent1(:));
    S(idx_lin) = +inf;

    % Top-K most complementary candidates
    K = max(2, ceil(sqrt(N)));
    [~, sortedIdx] = sort(S, 2, 'ascend');
    TopK = sortedIdx(:, 1:K);

    % Select the candidate with the smallest C
    CTopK = reshape(C(TopK), N, K);
    [~, bestInK] = min(CTopK, [], 2);
    Parent2 = TopK(sub2ind([N, K], (1:N)', bestInK))';

    % Generate mating pool
    MatingPool = zeros(1, 2*N);
    MatingPool(1:2:end) = Parent1;
    MatingPool(2:2:end) = Parent2;
end

%% ================================================================
%  Environmental selection
%% ================================================================
function [Population, FrontNo, Iso, BInfo] = EnvironmentalSelection(Population, N)
    EPS = 1e-12;
    PopObj = Population.objs;
    PopCon = Population.cons;
    [NQ, M] = size(PopObj);

    %% Stage 1: Constrained non-dominated sorting
    [FrontNoAll, MaxFNo] = NDSort(PopObj, PopCon, N);
    Next   = FrontNoAll < MaxFNo;
    S0Idx  = find(Next);
    Last   = find(FrontNoAll == MaxFNo);
    nLast  = numel(Last);
    Remain = N - sum(Next);

    %% Auxiliary representation
    [XiAll, ~]       = EmbedToPoincareBall(PopObj, PopCon);
    BAll             = BusemannFingerprint(XiAll);
    [CAll, DeltaAll] = BusemannDecompose(BAll);

    %% Degenerate cases
    if Remain <= 0 || nLast == 0
        if Remain <= 0
            SelIdx = S0Idx;
            if numel(SelIdx) > N, SelIdx = SelIdx(1:N); end
        else
            Next(Last) = true;
            SelIdx = find(Next);
            if numel(SelIdx) > N, SelIdx = SelIdx(1:N); end
        end
        [Population, FrontNo, Iso, BInfo] = PackOutput(Population, FrontNoAll, SelIdx, DeltaAll, CAll);
        return;
    end

    %% Stage 2: Hyperbolic convergence gap for the last front
    nF1 = sum(FrontNoAll == 1);
    rho = nF1 / max(NQ, 1);

    GapLast = ComputeHyperbolicGap(XiAll, S0Idx, Last);
    gMin = min(GapLast); gMax = max(GapLast);
    GapNorm = (GapLast - gMin) / (gMax - gMin + EPS);

    %% Prepare last-front data
    DeltaLast = DeltaAll(Last, :);
    CLast     = CAll(Last);

    cMin = min(CLast); cMax = max(CLast);
    CNorm = (CLast - cMin) / (cMax - cMin + EPS);

    % Initialize diversity distance
    if ~isempty(S0Idx)
        DistToS0 = pdist2(DeltaLast, DeltaAll(S0Idx, :), 'euclidean');
        deltaInit = min(DistToS0, [], 2);
    else
        deltaInit = inf(nLast, 1);
    end

    %% Extreme preservation
    extLocal = [];
    for m = 1:M
        [~, pos] = min(PopObj(Last, m));
        extLocal(end+1) = pos; %#ok<AGROW>
    end
    extLocal = unique(extLocal(:));

    selected = false(nLast, 1);
    if numel(extLocal) > Remain
        key = GapNorm(extLocal) + CNorm(extLocal);
        [~, ord] = sort(key, 'ascend');
        extLocal = extLocal(ord(1:Remain));
    end
    selected(extLocal) = true;

    delta = deltaInit;
    for t = 1:numel(extLocal)
        bestIdx = extLocal(t);
        remaining = find(~selected);
        if ~isempty(remaining)
            diff = DeltaLast(remaining,:) - DeltaLast(bestIdx,:);
            delta(remaining) = min(delta(remaining), sqrt(sum(diff.^2, 2)));
        end
    end

    %% Stage 3: Greedy filling in Delta-space
    wMin = 1 / max(M, 1);
    wC   = max(0, 1 - rho);
    Remain2 = max(0, Remain - sum(selected));

    for k = 1:Remain2
        candIdx = find(~selected);
        if isempty(candIdx), break; end

        wk = max(rho * (Remain2 - k + 1) / max(Remain2, 1), wMin);

        deltaCand = delta(candIdx);
        gapCand   = GapNorm(candIdx);
        cCand     = CNorm(candIdx);

        if all(isinf(deltaCand)), deltaCand = zeros(size(deltaCand)); end
        dMin = min(deltaCand); dMax = max(deltaCand);
        deltaNorm = (deltaCand - dMin) / (dMax - dMin + EPS);

        Score = deltaNorm - wk * gapCand - wC * cCand;
        [~, bestLocal] = max(Score);
        bestIdx = candIdx(bestLocal);
        selected(bestIdx) = true;

        remaining = find(~selected);
        if ~isempty(remaining)
            diff = DeltaLast(remaining,:) - DeltaLast(bestIdx,:);
            delta(remaining) = min(delta(remaining), sqrt(sum(diff.^2, 2)));
        end
    end

    %% Merge selected indices
    Next(Last(selected)) = true;
    SelIdx = find(Next);
    if numel(SelIdx) > N, SelIdx = SelIdx(1:N); end

    [Population, FrontNo, Iso, BInfo] = PackOutput(Population, FrontNoAll, SelIdx, DeltaAll, CAll);
end

%% ================================================================
%  Output packaging
%% ================================================================
function [Population, FrontNo, Iso, BInfo] = PackOutput(Population, FrontNoAll, SelIdx, DeltaAll, CAll)
    Population  = Population(SelIdx);
    FrontNo     = reshape(FrontNoAll(SelIdx), 1, []);
    DeltaSel    = DeltaAll(SelIdx, :);
    CSel        = CAll(SelIdx);

    Iso         = ComputeAngularIsolation(DeltaSel);
    BInfo.Delta = DeltaSel;
    BInfo.C     = CSel;
end

%% ================================================================
%  Poincaré-ball embedding
%% ================================================================
function [Xi, Fbar] = EmbedToPoincareBall(PopObj, PopCon)
    EPS = 1e-12;
    [NP, M] = size(PopObj);

    if nargin >= 2 && ~isempty(PopCon)
        Feasible = all(PopCon <= 0, 2);
        if ~any(Feasible), Feasible = true(NP,1); end
    else
        Feasible = true(NP, 1);
    end

    zIdeal = min(PopObj(Feasible,:), [], 1);
    zNadir = max(PopObj(Feasible,:), [], 1);
    range  = zNadir - zIdeal + EPS;
    Fbar   = max((PopObj - zIdeal) ./ range, 0);

    % Radius
    s = sum(Fbar, 2);

    % Direction on simplex
    d = Fbar ./ (s + EPS);
    zmask = s < EPS;
    if any(zmask), d(zmask,:) = 1/M; end

    % Simplex to sphere
    u  = sqrt(max(d, 0));
    u  = u ./ (sqrt(sum(u.^2, 2)) + EPS);

    % Radial mapping
    sf = s(Feasible);
    if isempty(sf), sf = s; end
    gamma = atanh(0.6) / (median(sf) + EPS);
    
    Xi = tanh(gamma * s) .* u;

    % Safeguard boundary
    normXi = sqrt(sum(Xi.^2, 2));
    mask   = normXi >= 1 - 1e-10;
    if any(mask)
        Xi(mask,:) = Xi(mask,:) .* ((1-1e-10) ./ (normXi(mask) + EPS));
    end
end

%% ================================================================
%  Busemann fingerprint
%% ================================================================
function B = BusemannFingerprint(Xi)
    EPS     = 1e-12;
    normXi2 = sum(Xi.^2, 2);
    denom   = max(1 - normXi2, EPS);

    dist2 = max(normXi2 + 1 - 2*Xi, EPS);
    B     = log(dist2 ./ denom);
end

%% ================================================================
%  Busemann decomposition
%% ================================================================
function [C, Delta] = BusemannDecompose(B)
    C     = mean(B, 2);
    Delta = B - C;
end

%% ================================================================
%  Hyperbolic convergence gap
%% ================================================================
function Gap = ComputeHyperbolicGap(XiAll, S0Idx, LastIdx)
    EPS = 1e-12;
    XiL = XiAll(LastIdx, :);

    if isempty(S0Idx)
        normXi = min(sqrt(sum(XiL.^2, 2)), 1-EPS);
        Gap    = 2 * atanh(normXi);
    else
        D   = HyperbolicDistMatrix(XiL, XiAll(S0Idx, :));
        Gap = min(D, [], 2);
    end
    Gap = reshape(Gap, [], 1);
end

%% ================================================================
%  Hyperbolic distance matrix in the Poincaré ball
%% ================================================================
function D = HyperbolicDistMatrix(X, Y)
    EPS   = 1e-12;
    nx    = sum(X.^2, 2);
    ny    = sum(Y.^2, 2)';

    dist2 = max(nx + ny - 2*(X*Y'), 0);
    denom = max((1-nx) .* (1-ny), EPS);

    Z = max(1 + 2*dist2./denom, 1+EPS);
    D = log(Z + sqrt(Z.*Z - 1));
    D(~isfinite(D)) = 0;
end

%% ================================================================
%  Angular isolation in Delta-space
%% ================================================================
function Iso = ComputeAngularIsolation(Delta)
    EPS = 1e-12;
    NP  = size(Delta, 1);
    if NP <= 1
        Iso = inf(1, NP);
        return;
    end

    nrm = sqrt(sum(Delta.^2, 2));
    zeroMask = nrm < EPS;
    nrm(zeroMask) = 1;
    U = Delta ./ nrm;

    S = max(-1, min(1, U * U'));
    S(1:NP+1:end) = 1;

    A = real(acos(S));
    A(1:NP+1:end) = inf;

    isoCol = min(A, [], 2);
    isoCol(zeroMask) = 0;
    Iso = isoCol';
end