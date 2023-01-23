%
%clone from https://github.com/kepecsLab/CellBase/
%bbisstim.m
%-------------------------------------------------------------------------
function [p_value Idiff] = isikldist(spt, dt, win, res)

    spt_baseline = spt;
    % Trial number and epoch length
    [tno tl] = size(spt_baseline); %baseline

    % Number of bins for ISI histograms
    nmbn = round(res / dt);

    % Pre-stimulus time window to consider for null hypothesis
    st = abs(win(1)) / dt; % number of pre-stim values in 'spt'

    % ISI histogram - baseline
    edges = 0:nmbn + 1;
    nm = floor(st / nmbn);
    lsi = zeros(tno, nm); % ISI's
    slsi = zeros(tno, nm); % sorted ISI's
    hlsi = zeros(nmbn + 1, nm); % ISI hist.; +1: zero when no spike in the segment
    nhlsi = zeros(nmbn + 1, nm); % normalized ISI histogram
    next = 1;

    for t = 1:nmbn:st

        for k = 1:tno
            cspt = spt_baseline(k, t:t + nmbn - 1);
            pki = find(cspt, 1, 'first');

            if ~isempty(pki)
                lsi(k, next) = pki;
            else
                lsi(k, next) = 0;
            end

        end

        slsi(:, next) = sort(lsi(:, next));
        hst = hist(slsi(:, next), edges);
        hlsi(:, next) = hst(1:end - 1);
        nhlsi(:, next) = hlsi(:, next) / sum(hlsi(:, next));
        next = next + 1;
    end
    
    % ISI histogram - test
    spt_test = spt;
    tno_test = size(spt_test, 1);
    lsi_tt = nan(tno_test, 1);

    for k = 1:tno_test
        cspt = spt_test(k, st + 1:st + nmbn);
        pki = find(cspt, 1, 'first');

        if ~isempty(pki)
            lsi_tt(k, 1) = pki;
        else
            lsi_tt(k, 1) = 0;
        end

    end

    slsi_tt = sort(lsi_tt(:, 1));
    hst = hist(slsi_tt, edges);
    hlsi(:, next) = hst(1:end - 1);
    nhlsi(:, next) = hlsi(:, next) / sum(hlsi(:, next));

   
    %figure      % plot ISIs
    %imagesc(lsi)
    %figure      % plot sorted ISIs
    %imagesc(slsi)
    %figure      % plot ISI histograms
    % imagesc(hlsi(2:end,:))
    % figure, hold on
    % plot(cumsum(mean(nhlsi(:,1:end-1),2)))
    %plot(cumsum(nhlsi(:,end)))

   
    % Symmetric KL-divergence and JS-divergence
    kn = st / nmbn + 1;
    jsd = nan(kn, kn); % pairwise modified JS-divergence (which is a metric!)
   
    for k1 = 1:kn
        D1 = nhlsi(:, k1);

        for k2 = k1 + 1:kn
            D2 = nhlsi(:, k2);
            jsd(k1, k2) = sqrt(JSdiv(D1, D2) * 2);
        end

    end

    % figure    % plot KL-distance
    % imagesc(kld)

    % Calculate p-value and information difference
    [p_value Idiff] = makep(jsd, kn);
    % keyboard
   
return;
    % -------------------------------------------------------------------------
    function [p_value Idiff] = makep(kld, kn)
        % Calculates p value from distance matrix.

        pnhk = kld(1:kn - 1, 1:kn - 1);
        nullhypkld = pnhk(~isnan(pnhk)); % nullhypothesis
        testkld = median(kld(1:kn - 1, kn)); % value to test
        sno = length(nullhypkld(:)); % sample size for nullhyp. distribution
        p_value = length(find(nullhypkld >= testkld)) / sno;
        Idiff = testkld - median(nullhypkld);
return;
%%%%%%%%%%%%%%%%%%
function D = JSdiv(P,Q)
    %JSDIV   Jensen-Shannon divergence.
    %   D = JSDIV(P,Q) calculates the Jensen-Shannon divergence of the two 
    %   input distributions.
    %
    %   See also KLDIST, FDIV, HDISC, BDIST, CHISQUAREDIV, VARDIST and 
    %   HARMONICMEAN.
    
    % Input argument check
    error(nargchk(2,2,nargin))
    if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
        error('Input arguments must be probability distributions.')
    end
    if ~isequal(size(P),size(Q))
        error('Input distributions must be of the same size.')
    end
    
    % JS-divergence
    M = (P + Q) / 2;
    D1 = KLdist(P,M);
    D2 = KLdist(Q,M);
    D = (D1 + D2) / 2;
    return;
    %%%%%%%%%%
    function D = KLdist(P,Q)
        %KLDIST   Kullbach-Leibler distance.
        %   D = KLDIST(P,Q) calculates the Kullbach-Leibler distance (information
        %   divergence) of the two input distributions.
        %
        %   See also JSDIV, FDIV, HDISC, BDIST, CHISQUAREDIV, VARDIST and 
        %   HARMONICMEAN.
        
        % Input argument check
        error(nargchk(2,2,nargin))
        if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
            error('Input arguments must be probability distributions.')
        end
        if ~isequal(size(P),size(Q))
            error('Input distributions must be of the same size.')
        end
        
        % KL-distance
        P2 = P(P.*Q>0);     % restrict to the common support
        Q2 = Q(P.*Q>0);
        P2 = P2 / sum(P2);  % renormalize
        Q2 = Q2 / sum(Q2);
        
        D = sum(P2.*log(P2./Q2));
        
        % Alternative way of computation:
        % HPQ = -sum(P2.*log(Q2));      % cross-entropy
        % HP = -sum(P2.*log(P2));       % entropy
        % D = HPQ - HP;
    