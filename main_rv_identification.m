%% ========================================================================
%  IDENTIFICAZIONE ARMA DELLA VOLATILITA' REALIZZATA — versione 2
%
%  Modifiche rispetto alla v1:
%   (A) Validazione dei residui: auto-correlation test + Ljung-Box
%   (B) Intervalli di confidenza sui parametri -> errorbar su a1
%   (C) Modello pooled con medie di periodo rimosse (merge multi-esperimento)
%   (D) PRESS (LOO esatta) e GCV al posto/accanto all'hold-out 70/30
%   (E) Controllo automatico delle quasi-cancellazioni polo-zero
%   (F) Correzione del bias di campione finito su a1
%   (G) Dispersione cross-market dei poli per periodo
%   (H) RPEM con forgetting factor: a1(t) continuo nel tempo
%   (I) Fix minori: parsing date robusto, legenda 'data1', RV normalizzata
%% ========================================================================

clear all
close all
clc

index_files = {'^tsx_d.csv', '^dax_d.csv', '^nkx_d.csv', ...
               '^kospi_d.csv', '^ndx_d.csv', '^fmib_d.csv'};
index_names = {'TSX', 'DAX', 'Nikkei225', 'KOSPI', 'Nasdaq100', 'FTSE-MIB'};
n_idx       = length(index_files);

period_names  = {'2000-2006 (Dot-com recovery)', ...
                 '2007-2010 (Financial crisis)', ...
                 '2011-2017 (Bull market)', ...
                 '2018-2026 (COVID + rate hikes)'};
period_short  = {'P1: Dot-com','P2: Crisis','P3: Bull','P4: COVID'};
period_bounds = {'2006-12', '2010-12', '2017-12', ''};
n_per         = 4;

Ts        = 1;
na_range  = 1:6;
nc_range  = 0:4;
h_cv      = 1;
epsilon   = 1e-6;

% -- opzioni nuove ---------------------------------------------------------
L_acf        = 20;      % lag massimo per i test sui residui
alpha_bound  = 2.58;    % 99% -> +/- 2.58/sqrt(N)  (par. 8.6.1)
pz_tol       = 0.10;    % soglia per segnalare quasi-cancellazioni polo-zero
lambda_rls   = 0.97;    % forgetting factor RPEM -> memoria ~ 1/(1-l) = 33 mesi
S0_rls       = 1e-2;    % inizializzazione della matrice di accumulo
normalize_RV = true;    % normalizza la RV per il n. di sedute del mese

%% ── 1. STRUTTURE DI MEMORIZZAZIONE ──────────────────────────────────────
all_log_vol = cell(n_idx, 1);
all_labels  = cell(n_idx, 1);
all_months  = cell(n_idx, 1);
all_periods = cell(n_idx, 1);

na_bic_tab    = zeros(n_idx, n_per);
nc_bic_tab    = zeros(n_idx, n_per);
na_cv_tab     = zeros(n_idx, n_per);
nc_cv_tab     = zeros(n_idx, n_per);
na_gcv_tab    = zeros(n_idx, n_per);   % (D)
nc_gcv_tab    = zeros(n_idx, n_per);   % (D)
na_press_tab  = zeros(n_idx, n_per);   % (D)
noisevar_tab  = zeros(n_idx, n_per);
models_period = cell(n_idx, n_per);

a1_tab        = zeros(n_idx, n_per);   % (B)
a1_std_tab    = NaN(n_idx, n_per);     % (B)
a1_corr_tab   = zeros(n_idx, n_per);   % (F)
Nseg_tab      = zeros(n_idx, n_per);

acf_out_tab   = NaN(n_idx, n_per);     % (A) frazione di lag fuori banda
lb_p_tab      = NaN(n_idx, n_per);     % (A) p-value Ljung-Box
pz_min_tab    = NaN(n_idx, n_per);     % (E) distanza minima polo-zero

% pooled: (a) media globale  (b) medie di periodo
models_pooled     = cell(n_idx, 1);
noisevar_pooled   = zeros(n_idx, 1);
a1_pooled         = zeros(n_idx, 1);
a1_pooled_std     = NaN(n_idx, 1);
acf_out_pooled    = NaN(n_idx, 1);
lb_p_pooled       = NaN(n_idx, 1);

models_pooled_dm   = cell(n_idx, 1);   % (C)
noisevar_pooled_dm = zeros(n_idx, 1);  % (C)
a1_pooled_dm       = zeros(n_idx, 1);  % (C)
a1_pooled_dm_std   = NaN(n_idx, 1);    % (C)

theta_rls = cell(n_idx, 1);            % (H)

%% ── 2. CARICAMENTO E PREPROCESSING ──────────────────────────────────────
for ix = 1:n_idx

    fprintf('\n########################################\n');
    fprintf('  LOADING INDEX: %s (%s)\n', index_names{ix}, index_files{ix});
    fprintf('########################################\n');

    data_raw     = readtable(index_files{ix});
    close_prices = data_raw.Close;

    % (I) parsing robusto: readtable puo' aver gia' convertito la colonna
    dates = data_raw.Date;
    if isdatetime(dates)
        dates_dt = dates;
    elseif isnumeric(dates)
        dates_dt = datetime(dates, 'ConvertFrom', 'datenum');
    else
        dates_dt = datetime(string(dates), 'InputFormat', 'yyyy-MM-dd');
    end

    % rendimenti log giornalieri (%)
    daily_returns = diff(log(close_prices)) * 100;
    dates_returns = dates_dt(2:end);

    % volatilita' realizzata mensile
    year_month    = dateshift(dates_returns, 'start', 'month');
    unique_months = unique(year_month);
    M_all         = length(unique_months);

    volatility   = zeros(M_all, 1);
    n_days       = zeros(M_all, 1);
    month_labels = strings(M_all, 1);
    for i = 1:M_all
        sel             = (year_month == unique_months(i));
        n_days(i)       = sum(sel);
        volatility(i)   = sqrt(sum(daily_returns(sel).^2));
        month_labels(i) = string(unique_months(i), 'yyyy-MM');
    end

    % (I) normalizzazione per il numero di sedute: rimuove la componente
    %     deterministica 0.5*log(n_days) che ha struttura stagionale
    if normalize_RV
        volatility = volatility ./ sqrt(n_days);
    end

    volatility    = volatility(1:end-1);
    month_labels  = month_labels(1:end-1);
    months_dt     = unique_months(1:end-1);
    M             = length(volatility);

    log_vol = log(volatility + epsilon);

    all_log_vol{ix} = log_vol;
    all_labels{ix}  = month_labels;
    all_months{ix}  = months_dt;

    idx1_end = find(month_labels == period_bounds{1});
    idx2_end = find(month_labels == period_bounds{2});
    idx3_end = find(month_labels == period_bounds{3});
    idx4_end = M;

    if isempty(idx1_end) || isempty(idx2_end) || isempty(idx3_end)
        error('Index %s non copre l''intero range 2000-2026.', index_names{ix});
    end

    all_periods{ix} = {1:idx1_end, idx1_end+1:idx2_end, ...
                       idx2_end+1:idx3_end, idx3_end+1:idx4_end};

    fprintf('Total months: %d | Periods: [%d %d %d %d]\n', M, ...
            idx1_end, idx2_end-idx1_end, idx3_end-idx2_end, idx4_end-idx3_end);
end

%% ── 3. MODELLI PERIOD-SPECIFIC ──────────────────────────────────────────
for ix = 1:n_idx
    fprintf('\n============================================\n');
    fprintf(' INDEX: %s — Period-specific model fitting\n', index_names{ix});
    fprintf('============================================\n');

    log_vol = all_log_vol{ix};
    periods = all_periods{ix};

    for p = 1:n_per
        seg    = log_vol(periods{p});
        n_seg  = length(seg);
        mu_seg = mean(seg);
        seg_d  = iddata(seg - mu_seg, [], Ts);

        Nseg_tab(ix,p) = n_seg;

        % ---- selezione dell'ordine: quattro criteri ---------------------
        [na_bic, nc_bic]   = select_order_bic(seg_d, na_range, nc_range);
        [na_gcv, nc_gcv]   = select_order_gcv(seg_d, na_range, nc_range);
        na_press           = select_order_press(seg - mu_seg, na_range);

        n_tr   = floor(0.7 * n_seg);
        seg_tr = iddata(seg(1:n_tr) - mean(seg(1:n_tr)), [], Ts);
        seg_va = iddata(seg(n_tr+1:end) - mean(seg(1:n_tr)), [], Ts);
        [na_cv, nc_cv] = select_order_cv(seg_tr, seg_va, na_range, nc_range, h_cv);

        % ---- modello finale (ordine BIC) --------------------------------
        M_p = armax(seg_d, [na_bic nc_bic]);

        % (B) deviazione standard di a1
        a1 = -M_p.A(2);
        s1 = NaN;
        try
            [~, dpv] = getpvec(M_p);
            if ~isempty(dpv) && ~isnan(dpv(1))
                s1 = dpv(1);
            end
        catch
        end

        % (F) correzione bias campione finito (Marriott-Pope / Kendall)
        a1_c = a1 + (1 + 3*a1)/n_seg;

        % (A) test sui residui
        E  = resid(M_p, seg_d);
        e  = E.OutputData;
        pM = na_bic + nc_bic;
        [frac_out, ~, pval] = whiteness_test(e, L_acf, pM);

        % (E) quasi-cancellazioni polo-zero
        dmin = min_pole_zero_dist(M_p);

        % ---- memorizzazione ---------------------------------------------
        na_bic_tab(ix,p)   = na_bic;    nc_bic_tab(ix,p)  = nc_bic;
        na_cv_tab(ix,p)    = na_cv;     nc_cv_tab(ix,p)   = nc_cv;
        na_gcv_tab(ix,p)   = na_gcv;    nc_gcv_tab(ix,p)  = nc_gcv;
        na_press_tab(ix,p) = na_press;
        noisevar_tab(ix,p) = M_p.NoiseVariance;
        models_period{ix,p}= M_p;
        a1_tab(ix,p)       = a1;
        a1_std_tab(ix,p)   = s1;
        a1_corr_tab(ix,p)  = a1_c;
        acf_out_tab(ix,p)  = frac_out;
        lb_p_tab(ix,p)     = pval;
        pz_min_tab(ix,p)   = dmin;

        fprintf(['  %-32s | BIC=(%d,%d) CV=(%d,%d) GCV=(%d,%d) PRESS=(%d,0)\n' ...
                 '  %-32s | a1=%.3f+/-%.3f (corr %.3f) | NoiseVar=%.5f | ACF out=%4.1f%% | LB p=%.3f\n'], ...
                period_names{p}, na_bic, nc_bic, na_cv, nc_cv, na_gcv, nc_gcv, na_press, ...
                '', a1, s1, a1_c, M_p.NoiseVariance, 100*frac_out, pval);
    end
end

%% ── 4. MODELLI POOLED: media globale vs medie di periodo ────────────────
fprintf('\n============================================\n');
fprintf(' CROSS-PERIOD (POOLED) MODELS\n');
fprintf('============================================\n');

for ix = 1:n_idx
    log_vol = all_log_vol{ix};
    periods = all_periods{ix};

    % (a) versione originale: UNA sola media globale
    data_full = iddata(log_vol - mean(log_vol), [], Ts);
    [na_p, nc_p] = select_order_bic(data_full, na_range, nc_range);
    Mp  = armax(data_full, [na_p nc_p]);

    E   = resid(Mp, data_full);
    [fo, ~, pv_lb] = whiteness_test(E.OutputData, L_acf, na_p + nc_p);

    models_pooled{ix}   = Mp;
    noisevar_pooled(ix) = Mp.NoiseVariance;
    a1_pooled(ix)       = -Mp.A(2);
    acf_out_pooled(ix)  = fo;
    lb_p_pooled(ix)     = pv_lb;
    try
        [~, dpv] = getpvec(Mp);
        a1_pooled_std(ix) = dpv(1);
    catch
    end

    % (C) versione corretta: coefficienti comuni, MEDIE PER PERIODO.
    %     merge() crea un iddata multi-esperimento: armax stima un unico
    %     theta senza mai costruire regressori a cavallo dei confini.
    segs = cell(1, n_per);
    for p = 1:n_per
        s       = log_vol(periods{p});
        segs{p} = iddata(s - mean(s), [], Ts);
    end
    data_dm = merge(segs{:});
    Mdm     = armax(data_dm, [na_p nc_p]);   % stesso ordine, confronto equo

    models_pooled_dm{ix}   = Mdm;
    noisevar_pooled_dm(ix) = Mdm.NoiseVariance;
    a1_pooled_dm(ix)       = -Mdm.A(2);
    try
        [~, dpv] = getpvec(Mdm);
        a1_pooled_dm_std(ix) = dpv(1);
    catch
    end

    fprintf(['  %-11s | order=(%d,%d) | a1 pooled=%.3f  a1 pooled-dm=%.3f  (delta %+0.3f)\n' ...
             '  %-11s | NoiseVar=%.5f -> %.5f | ACF out=%4.1f%% | LB p=%.3f\n'], ...
            index_names{ix}, na_p, nc_p, a1_pooled(ix), a1_pooled_dm(ix), ...
            a1_pooled_dm(ix) - a1_pooled(ix), ...
            '', noisevar_pooled(ix), noisevar_pooled_dm(ix), 100*fo, pv_lb);
end

%% ── 5. HEATMAP DEI COEFFICIENTI AR ──────────────────────────────────────
na_max   = max(na_bic_tab(:));
coeff_AR = NaN(n_idx, n_per, na_max);
for ix = 1:n_idx
    for p = 1:n_per
        a = -models_period{ix,p}.A(2:end);
        coeff_AR(ix,p,1:length(a)) = a;
    end
end

figure('Name', 'AR Coefficient Heatmaps by Lag');
for k = 1:na_max
    subplot(1, na_max, k);
    Ck = squeeze(coeff_AR(:,:,k));
    imagesc(Ck, 'AlphaData', ~isnan(Ck));
    colorbar; colormap(parula);
    set(gca, 'XTick', 1:n_per, 'XTickLabel', {'P1','P2','P3','P4'});
    set(gca, 'YTick', 1:n_idx, 'YTickLabel', index_names);
    title(['a_' num2str(k) ' coefficient']);
    xtickangle(45);
end
sgtitle('AR Coefficients Across Indexes and Periods');

%% ── 6. VARIANZA DEL RUMORE ──────────────────────────────────────────────
figure('Name', 'Noise Variance Comparison');
bar(noisevar_tab');
set(gca, 'XTickLabel', period_short); xtickangle(20);
ylabel('Noise Variance (log-vol space)');
title('Innovation Noise Variance by Period and Index');
legend(index_names, 'Location', 'northoutside', 'Orientation', 'horizontal');
grid on;

%% ── 7. MAPPE POLO-ZERO + DISPERSIONE CROSS-MARKET ───────────────────────
colors    = lines(n_idx);
pole_disp = zeros(n_per, 1);

for p = 1:n_per
    figure('Name', ['Pole-Zero Map - ' period_names{p}]);
    hold on;
    theta = linspace(0, 2*pi, 200);
    % (I) fix: il cerchio unitario non deve comparire in legenda
    plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

    pr = zeros(n_idx, 1);
    for ix = 1:n_idx
        M_p    = models_period{ix,p};
        poles  = roots(M_p.A);
        zeros_ = roots(M_p.C);
        pr(ix) = max(abs(poles));
        plot(real(poles), imag(poles), 'x', 'Color', colors(ix,:), ...
             'MarkerSize', 10, 'LineWidth', 2, ...
             'DisplayName', [index_names{ix} ' poles']);
        if ~isempty(zeros_)
            plot(real(zeros_), imag(zeros_), 'o', 'Color', colors(ix,:), ...
                 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    pole_disp(p) = std(pr);

    axis equal; grid on;
    xlabel('Real'); ylabel('Imaginary');
    title(sprintf('Poles (x) and Zeros (o): %s   [std raggio = %.3f]', ...
                  period_names{p}, pole_disp(p)));
    legend('show', 'Location', 'eastoutside');
    hold off;
end

%% ── 8. a1 PERIOD-SPECIFIC vs POOLED, CON BARRE D'ERRORE ─────────────────
figure('Name', 'Period-specific vs Pooled AR(1) Coefficient');
hold on;
for ix = 1:n_idx
    e1 = a1_std_tab(ix,:);
    e1(isnan(e1)) = 0;
    errorbar(1:n_per, a1_tab(ix,:), 1.96*e1, '-o', 'Color', colors(ix,:), ...
             'LineWidth', 1.5, 'CapSize', 8, ...
             'DisplayName', [index_names{ix} ' (period-specific)']);
    yline(a1_pooled(ix), '--', 'Color', colors(ix,:), 'HandleVisibility', 'off');
    yline(a1_pooled_dm(ix), ':',  'Color', colors(ix,:), 'HandleVisibility', 'off');
end
xticks(1:n_per); xticklabels(period_short); xlim([0.7 n_per+0.3]);
ylabel('a_1 coefficient');
title({'AR(1): period-specific (o, IC 95%) vs pooled', ...
       'tratteggio = pooled media globale, punteggiato = pooled medie di periodo'});
legend('show', 'Location', 'eastoutside'); grid on; hold off;

% versione con correzione del bias
figure('Name', 'a1: stima grezza vs corretta per bias');
plot(1:n_per, mean(a1_tab,1), '-o', 'LineWidth', 2, 'DisplayName', 'a_1 stimato');
hold on;
plot(1:n_per, mean(a1_corr_tab,1), '-s', 'LineWidth', 2, ...
     'DisplayName', 'a_1 corretto (bias campione finito)');
xticks(1:n_per); xticklabels(period_short);
ylabel('a_1 medio sui 6 indici'); grid on; legend('show'); hold off;

%% ── 9. RPEM: a1(t) COME FUNZIONE CONTINUA DEL TEMPO ─────────────────────
% Regressore aumentato phi(t) = [y(t-1); 1] cosi' il livello locale viene
% tracciato insieme al coefficiente, senza detrending a priori.
figure('Name', 'RPEM — a1(t) con forgetting factor');
hold on;
burn = ceil(3/(1 - lambda_rls));

for ix = 1:n_idx
    y  = all_log_vol{ix};
    th = rls_forget(y, lambda_rls, S0_rls);
    theta_rls{ix} = th;

    tt = all_months{ix};
    plot(tt(burn:end), th(burn:end,1), 'Color', colors(ix,:), ...
         'LineWidth', 1.3, 'DisplayName', index_names{ix});
end

% confini di periodo scelti a priori: il grafico dice se sono nei punti giusti
per_dates = datetime({'2006-12-01','2010-12-01','2017-12-01'}, ...
                     'InputFormat','yyyy-MM-dd');
for k = 1:numel(per_dates)
    xline(per_dates(k), 'k--', 'HandleVisibility', 'off');
end
ylabel('a_1(t)'); ylim([0 1]); grid on;
title(sprintf('RPEM, \\lambda = %.2f (memoria \\approx %.0f mesi)', ...
              lambda_rls, 1/(1-lambda_rls)));
legend('show', 'Location', 'eastoutside'); hold off;

%% ── 10. TABELLE RIASSUNTIVE ─────────────────────────────────────────────
fprintf('\n\n========== ORDINI SELEZIONATI: BIC / CV / GCV / PRESS ==========\n');
fprintf('%-11s', 'Index');
for p = 1:n_per, fprintf('  %-22s', period_short{p}); end
fprintf('\n');
for ix = 1:n_idx
    fprintf('%-11s', index_names{ix});
    for p = 1:n_per
        fprintf('  (%d,%d)/(%d,%d)/(%d,%d)/%d   ', ...
            na_bic_tab(ix,p), nc_bic_tab(ix,p), na_cv_tab(ix,p), nc_cv_tab(ix,p), ...
            na_gcv_tab(ix,p), nc_gcv_tab(ix,p), na_press_tab(ix,p));
    end
    fprintf('\n');
end

fprintf('\n========== VALIDAZIONE RESIDUI (99%%, L=%d lag) ==========\n', L_acf);
fprintf('%-11s %-8s', 'Index', 'model');
for p = 1:n_per, fprintf('  %-16s', period_short{p}); end
fprintf('  %-16s\n', 'POOLED');
for ix = 1:n_idx
    fprintf('%-11s %-8s', index_names{ix}, 'ACF/LB');
    for p = 1:n_per
        fprintf('  %4.1f%% p=%.3f %s', 100*acf_out_tab(ix,p), lb_p_tab(ix,p), ...
                verdict(lb_p_tab(ix,p)));
    end
    fprintf('  %4.1f%% p=%.3f %s\n', 100*acf_out_pooled(ix), lb_p_pooled(ix), ...
            verdict(lb_p_pooled(ix)));
end
fprintf('\n[OK] = residui compatibili con rumore bianco (p > 0.05)\n');
fprintf('[NO] = struttura residua non modellata -> ordine insufficiente\n');

fprintf('\n========== POOLED: media globale vs medie di periodo ==========\n');
fprintf('%-11s %8s %8s %8s %10s %10s\n', 'Index', 'a1 glob', 'a1 dm', 'delta', 'NV glob', 'NV dm');
for ix = 1:n_idx
    fprintf('%-11s %8.3f %8.3f %+8.3f %10.5f %10.5f\n', index_names{ix}, ...
        a1_pooled(ix), a1_pooled_dm(ix), a1_pooled_dm(ix)-a1_pooled(ix), ...
        noisevar_pooled(ix), noisevar_pooled_dm(ix));
end
fprintf('\nSe delta < 0: parte della persistenza del pooled era SPURIA,\n');
fprintf('generata dai salti di livello tra regimi (cfr. par. 8.4).\n');

fprintf('\n========== CANCELLAZIONI POLO-ZERO (soglia %.2f) ==========\n', pz_tol);
n_flag = 0;
for ix = 1:n_idx
    for p = 1:n_per
        if ~isnan(pz_min_tab(ix,p)) && pz_min_tab(ix,p) < pz_tol
            fprintf('  ATTENZIONE %-11s %-6s : dist min = %.4f -> sovra-parametrizzazione\n', ...
                    index_names{ix}, period_short{p}, pz_min_tab(ix,p));
            n_flag = n_flag + 1;
        end
    end
end
if n_flag == 0
    fprintf('  Nessuna quasi-cancellazione rilevata.\n');
end

fprintf('\n========== CONFRONTO PREDITTIVO (par. 8.6.2, G ~ s2(1+p/N)) ==========\n');
fprintf('%-11s %12s %12s %10s\n', 'Index', 'G pooled', 'G per-period', 'guadagno');
for ix = 1:n_idx
    Ntot   = length(all_log_vol{ix});
    p_pool = na_pooled(models_pooled{ix}) + 1;                 % coeff + 1 media
    p_per  = sum(na_bic_tab(ix,:) + nc_bic_tab(ix,:)) + n_per; % coeff + 4 medie
    G_pool = noisevar_pooled(ix) * (1 + p_pool/Ntot);
    G_per  = mean(noisevar_tab(ix,:)) * (1 + p_per/Ntot);
    fprintf('%-11s %12.5f %12.5f %9.1f%%\n', index_names{ix}, G_pool, G_per, ...
            100*(G_pool - G_per)/G_pool);
end

fprintf('\n========== DISPERSIONE CROSS-MARKET DEI POLI ==========\n');
for p = 1:n_per
    fprintf('  %-14s : std raggio poli = %.4f\n', period_short{p}, pole_disp(p));
end

fprintf('\n========== COSTANTI DI TEMPO (mesi) ==========\n');
fprintf('%-11s', 'Index');
for p = 1:n_per, fprintf('  %-10s', period_short{p}); end
fprintf('\n');
for ix = 1:n_idx
    fprintf('%-11s', index_names{ix});
    for p = 1:n_per
        fprintf('  %10.2f', -1/log(abs(a1_tab(ix,p))));
    end
    fprintf('\n');
end


%% ════════════════════════════════════════════════════════════════════════
%  FUNZIONI LOCALI
%% ════════════════════════════════════════════════════════════════════════

function [na_best, nc_best, bic_table] = select_order_bic(data, na_range, nc_range)
% Ordine che minimizza il BIC (par. 8.6.4).
    bic_table = NaN(length(na_range), length(nc_range));
    best = Inf; na_best = na_range(1); nc_best = nc_range(1);
    for i = 1:length(na_range)
        for j = 1:length(nc_range)
            try
                M = armax(data, [na_range(i) nc_range(j)]);
                b = bic(M);
                bic_table(i,j) = b;
                if b < best
                    best = b; na_best = na_range(i); nc_best = nc_range(j);
                end
            catch
            end
        end
    end
end

function [na_best, nc_best, gcv_table] = select_order_gcv(data, na_range, nc_range)
% Generalized cross-validation (par. 8.6.5): J = V_N / (1 - p/N)^2.
% Approssima la LOO senza risolvere N problemi, applicabile anche a nc>0.
    N = size(data.OutputData, 1);
    gcv_table = NaN(length(na_range), length(nc_range));
    best = Inf; na_best = na_range(1); nc_best = nc_range(1);
    for i = 1:length(na_range)
        for j = 1:length(nc_range)
            p = na_range(i) + nc_range(j);
            if p >= N, continue; end
            try
                M = armax(data, [na_range(i) nc_range(j)]);
                J = M.NoiseVariance / (1 - p/N)^2;
                gcv_table(i,j) = J;
                if J < best
                    best = J; na_best = na_range(i); nc_best = nc_range(j);
                end
            catch
            end
        end
    end
end

function [na_best, press_vec] = select_order_press(y, na_range)
% LOO ESATTA via identita' di Woodbury (par. 8.6.5), valida solo per
% modelli AR puri (nc=0), che sono regressioni lineari:
%     J_PRESS = sum_t  eps(t)^2 / (1 - h_t)^2,   H = Phi(Phi'Phi)^-1 Phi'
% Costa quanto una singola stima e usa TUTTI i dati: decisivo su segmenti
% corti, dove l'hold-out 70/30 valida su pochissimi punti.
    y = y(:);
    N = length(y);
    press_vec = NaN(length(na_range), 1);
    best = Inf; na_best = na_range(1);
    for i = 1:length(na_range)
        na = na_range(i);
        if na >= N - 5, continue; end
        rows = na+1 : N;
        Phi  = zeros(length(rows), na);
        for k = 1:na
            Phi(:,k) = y(rows - k);
        end
        yv = y(rows);
        if rank(Phi) < na, continue; end
        th  = Phi \ yv;
        eps = yv - Phi*th;
        % h_t senza costruire H per intero
        R  = Phi' * Phi;
        h  = sum((Phi / R) .* Phi, 2);
        h  = min(h, 1 - 1e-9);
        J  = sum(eps.^2 ./ (1 - h).^2);
        press_vec(i) = J;
        if J < best
            best = J; na_best = na;
        end
    end
end

function [na_best, nc_best, fit_table] = select_order_cv(data_tr, data_va, na_range, nc_range, h)
% Hold-out cross-validation, mantenuta come termine di paragone.
    fit_table = NaN(length(na_range), length(nc_range));
    best = -Inf; na_best = na_range(1); nc_best = nc_range(1);
    opt = compareOptions('InitialCondition', 'z');
    for i = 1:length(na_range)
        for j = 1:length(nc_range)
            try
                M = armax(data_tr, [na_range(i) nc_range(j)]);
                [~, f] = compare(data_va, M, h, opt);
                fit_table(i,j) = f;
                if f > best
                    best = f; na_best = na_range(i); nc_best = nc_range(j);
                end
            catch
            end
        end
    end
end

function [frac_out, Q, pval, r_norm] = whiteness_test(e, L, p_par)
% Auto-correlation test (par. 8.6.1) + statistica portmanteau di Ljung-Box.
%   frac_out : frazione di lag s=1..L con |r(s)/r(0)| oltre 2.58/sqrt(N)
%   Q, pval  : Ljung-Box e p-value (chi2 con L-p gradi di liberta')
    e = e(:); e = e - mean(e);
    N = numel(e);
    L = min(L, N - 2);
    r = zeros(L+1, 1);
    for s = 0:L
        r(s+1) = sum(e(1:N-s) .* e(1+s:N)) / N;
    end
    r_norm = r / r(1);

    bound    = 2.58 / sqrt(N);
    frac_out = mean(abs(r_norm(2:end)) > bound);

    s_idx = (1:L)';
    Q     = N * (N + 2) * sum(r_norm(2:end).^2 ./ (N - s_idx));

    k = max(L - p_par, 1);                       % gradi di liberta'
    % approssimazione di Wilson-Hilferty a chi2 (evita chi2cdf)
    z    = ((Q/k)^(1/3) - (1 - 2/(9*k))) / sqrt(2/(9*k));
    pval = 0.5 * erfc(z / sqrt(2));
end

function d = min_pole_zero_dist(M)
% Distanza minima polo-zero (par. 8.6.3). NaN se il modello non ha zeri.
    P = roots(M.A);
    Z = roots(M.C);
    if isempty(Z) || isempty(P)
        d = NaN; return;
    end
    D = abs(P(:) - Z(:).');
    d = min(D(:));
end

function theta_hist = rls_forget(y, lambda, S0)
% RPEM / RLS con forgetting factor (par. 7.2), AR(1) con livello:
%     phi(t) = [y(t-1); 1],   theta = [a1; c]
%     S_t     = lambda*S_{t-1} + phi phi'
%     theta_t = theta_{t-1} + S_t^-1 phi [y(t) - phi' theta_{t-1}]
    y = y(:);
    N = numel(y);
    theta = [0; mean(y(1:min(12,N)))];
    S     = S0 * eye(2);
    theta_hist = NaN(N, 2);
    for t = 2:N
        phi   = [y(t-1); 1];
        S     = lambda * S + phi * phi';
        theta = theta + S \ (phi * (y(t) - phi' * theta));
        theta_hist(t,:) = theta';
    end
    theta_hist(1,:) = theta_hist(2,:);
end

function n = na_pooled(M)
% Numero di coefficienti AR di un modello idpoly.
    n = length(M.A) - 1;
end

function s = verdict(p)
    if isnan(p)
        s = '[??]';
    elseif p > 0.05
        s = '[OK]';
    else
        s = '[NO]';
    end
end