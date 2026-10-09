---
title: "simplePHENOTYPES — Equation-to-Code Audit Map"
subtitle: "Every theory-derived equation, its primary reference, and where it is implemented (function, file, line)"
author: "Independent dual-model audit (Claude Fable + OpenAI Codex), assembled from per-group reports"
date: "2026-09-29 audit at `c511c6f`; V2 line references re-checked against `20cd2c9` (2026-10-08); V1 references as remapped to `e666a2c` (2026-10-04)"
toc: true
toc-depth: 2
---

# How to read this document

Each entry gives the **equation as implemented**, the **variables**, the **primary reference** (a page or equation number is quoted only when the auditor verified it; otherwise it is marked *page unverified*), the **implementing function**, and the **exact `file:line` range**. V1 ranges are at `e666a2c` (remapped by git line mapping from the audit commit `c511c6f`; see Appendix C). V2 ranges, and the V2 rows of Appendix B, are at `20cd2c9` (2026-10-08): every V2 *Location* was re-mapped from the commit it was actually written against (the audit commit for `:NNN` shorthands, `e666a2c` or `4959077` for file-qualified ranges) and checked against the code at `20cd2c9` (content of the mapped lines, plus an automated check that the first R range of each entry falls inside its named function). Line hints inside *Notes* fields were updated only where an entry was revised (selection, OCS/usefulness, crossing, transcriptome); other *Notes* hints may still refer to `c511c6f`. Appendix A is the original automated check at the audit commit.

**Revision 2026-10-08 (equations changed in code after the audit).** The following entries were rewritten to match the current source; their page-verification status is unchanged: V2-crossing-3 (default meiosis is the two-pathway gamma interference model, $\nu = 2.6$, $p = 0$, since 2.0.0.9003, DECISION-047; Poisson is `interference = "poisson"`), with Poisson-only notes on V2-crossing-4/5/9 and V2-rust-core-1/13; V2-transcriptome-12 (realized cis fraction $v_{cis}/\mathrm{Var}(G)$); V2-transcriptome-19 ($\hat\kappa$ now subtracts the GREML-implied genetic trans share); V2-transcriptome-20 (`observe_counts` carries `@references`); V2-transcriptome-C2 (TX-F2 note); V2-grammar-25 and V2-transcriptome-25 (polynomial hash of the layer type in the sub-seed); V2-grammar-10 and V2-effects-arch-18 (`residual_mode`); V2-effects-arch-20 (`partner = "random"`); V2-selection-15 (multinomial bulk); V2-ocs-usefulness-marker-8 (`sample_parents()` allocates by default).

Conventions: dosage is coded −1/0/1 unless stated; `p` is the allele frequency of the +1 allele; `q = 1 − p`. "Own design" means the equation is a package design decision (documented in `docs/DECISIONS.md`) rather than a literature formula.

Sections V2 and V1 are separate: V1 is the frozen legacy `create_phenotypes()` engine (DECISION-003/008), V2 is the composable grammar, crossing, selection, prediction and transcriptome layers.


# Version 2 (grammar, crossing, selection, prediction, transcriptome, I/O)


## Grammar: variance budget, layers, realization, orthogonal model, vQTL, seeds

*Source report:* `v2-grammar.md`  
*Entries:* 28


**V2-grammar-1. Layer scaling to prop**


$\tilde c = \dfrac{c - \bar c}{\mathrm{sd}(c)}\sqrt{\pi_\ell}$


- *Variables:* c raw component, $\pi_\ell$ layer prop

- *Reference:* Package convention (SPEC §2)

- *Function:* `.genetic_matrix`

- *Location:* R/grammar_realize.R:139-161 (146)

- *Notes:* Sample sd (n−1). Layers summed; no orthogonalization.


**V2-grammar-2. Additive component**


$c_A = \sum_j e_j x_{ij}$, $x\in\{-1,0,1\}$


- *Variables:* e geometric/custom effects

- *Reference:* Fernandes & Lipka 2020 (v1 convention)

- *Function:* `.component_raw`

- *Location:* R/grammar_realize.R:498-509 (500), centred 432


**V2-grammar-3. Orthogonal genotypic value**


$g_{ij} = a_j x_{ij} + d_j\,\mathbb 1[x_{ij}=0]$ i.e. $-a/+d/+a$


- *Variables:* a, d per locus

- *Reference:* Falconer & Mackay 1996 (genotypic values; page unverified)

- *Function:* `.component_raw`

- *Location:* R/grammar_realize.R:501-507 (506)

- *Notes:* DECISION-020


**V2-grammar-4. Dominance component**


$c_D=\sum_j e_j\,\mathbb 1[x_{ij}=0]$


- *Variables:* het indicator

- *Reference:* Package convention

- *Function:* `.component_raw`

- *Location:* R/grammar_realize.R:510

- *Notes:* Not Fisher's D; see GRAM-F1


**V2-grammar-5. Epistasis unit column**


$z_{i}=\prod_{k}\big(t_{ik}-\bar t_k\big)$, $t_k$ = dosage ("a") or het indicator ("d")


- *Variables:* interaction positions k

- *Reference:* Package convention (centred products; cf. Cockerham/NOIA only under HWE+LE)

- *Function:* `.epi_unit_column`

- *Location:* R/grammar_realize.R:579-589 (586, 588)

- *Notes:* Single source for realization and pleio normalizer


**V2-grammar-6. Epistasis component**


$c_E=\sum_p e_p z_{ip}$


- *Variables:* e per set

- *Reference:* as above

- *Function:* `.component_raw`

- *Location:* R/grammar_realize.R:511-529 (524)


**V2-grammar-7. Geometric effect series**


$e_k = b^{k}$, $k=1..n$, default $b=0.5$


- *Variables:* b base

- *Reference:* Fernandes & Lipka 2020 (v1 `geometric`)

- *Function:* `.effect_series`

- *Location:* R/effects_series.R:23-72 (64)

- *Notes:* length-n vector used verbatim


**V2-grammar-8. Repulsion phase**


$e_k \leftarrow e_k(-1)^{k+1}$


- *Reference:* Package convention

- *Function:* `.apply_phase`

- *Location:* R/grammar_layers.R:1482-1490 (1487)

- *Notes:* positional, not LD-derived (documented)


**V2-grammar-9. Residual budget**


$\sigma^2_e = \max(0,\,1-\sum_\ell \pi_\ell - \sum_v \pi_v)$


- *Variables:* mean-effect and vqtl props

- *Reference:* SPEC §2

- *Function:* `.realize_phenotype`

- *Location:* R/grammar_realize.R:55-59 (59)


**V2-grammar-10. Residual draw**


$e \sim N(0,\sigma_e^2)$ then $e\leftarrow \dfrac{e-\bar e}{\mathrm{sd}(e)}\sigma_e$


- *Reference:* Package convention

- *Function:* `.draw_residual`

- *Location:* R/effects_series.R:85-102 (99)

- *Notes:* exact sample variance; same standardized vector per (seed,trait,rep). Default `residual_mode = "fixed"`; `residual_mode = "random"` (added after the audit) returns the raw $N(0,\sigma_e^2)$ draw (`:94`)


**V2-grammar-11. Phenotype**


$y_{it} = G_{it} + T_{it} + e_{it} + \mu_t$


- *Variables:* G marker genetic, T transcriptome, μ mean

- *Reference:* —

- *Function:* `.realize_phenotype`

- *Location:* R/grammar_realize.R:96


**V2-grammar-12. Requested budget identity**


$\sum_{A,D,E}\pi = h^2$; $\sum_{\text{all}}\pi\le 1$


- *Reference:* SPEC §2

- *Function:* `.resolve_prop`, `.add_layer`, `.check_h2_complete`

- *Location:* R/grammar_layers.R:1085-1098, 1102-1115; R/grammar_realize.R:1153-1166

- *Notes:* tol 1e-8


**V2-grammar-13. One-call split**


$\pi_c = h^2/\vert{}\text{model}\vert{}$


- *Variables:* components A,D,E

- *Reference:* SPEC §4.1

- *Function:* `.build_one_call`

- *Location:* R/grammar_simulate_phenotype.R:767-780 (769)

- *Notes:* "AE": n_pairs = n_qtn


**V2-grammar-14. Realized H²**


$\hat H^2_t = \frac{1}{R}\sum_r \dfrac{\mathrm{Var}(G_{\cdot t r})}{\mathrm{Var}(y_{\cdot t r})}$


- *Variables:* G = marker + genetic-mediated Tx

- *Reference:* SPEC §2 (V2)

- *Function:* `.realized_h2`

- *Location:* R/grammar_realize.R:1078-1106 (716-720)

- *Notes:* vqtl in denominator only


**V2-grammar-15. Average effect of substitution**


$\alpha_j = a_j + d_j(q_j-p_j) = a_j + d_j(1-2p_j)$


- *Variables:* p freq of +1 allele

- *Reference:* Falconer 1985 *Genet. Res.* 46:337-347; Falconer & Mackay 1996 (page unverified); Fisher 1918

- *Function:* `.avg_effect`

- *Location:* R/grammar_realize.R:182-184 (183)

- *Notes:* verified numerically


**V2-grammar-16. Breeding value**


$A_i=\sum_j \alpha_j (x_{ij}-2p_j)$, $x\in\{0,1,2\}$


- *Reference:* Falconer & Mackay 1996; Lynch & Walsh 1998

- *Function:* `.breeding_value_matrix`

- *Location:* R/grammar_realize.R:282-303 (227-235)

- *Notes:* a, d rescaled by $\sqrt{\pi}/\mathrm{sd}$ (273); DECISION-019


**V2-grammar-17. Orthogonal split**


$\mathrm{Var}(g)=\mathrm{Var}(A)+\mathrm{Var}(D)+2\mathrm{Cov}(A,D)$; rows $\tfrac{\mathrm{Var}(A)}{\mathrm{Var}(g)},\tfrac{\mathrm{Var}(D)}{\mathrm{Var}(g)},1-\text{add}-\text{dom}$


- *Variables:* D = g − A

- *Reference:* Fisher 1918; Falconer & Mackay 1996

- *Function:* `.orthogonal_var_split`, `.variance_budget`

- *Location:* R/grammar_realize.R:778-794 (646-650), 536-551

- *Notes:* Cov = 0 under random mating (verified)


**V2-grammar-18. HWE theory (check only)**


$V_A=\sum 2p_jq_j\alpha_j^2$, $V_D=\sum (2p_jq_jd_j)^2$


- *Reference:* Falconer & Mackay 1996 (page unverified)

- *Function:* — (audit check)

- *Location:* —

- *Notes:* matches 0.8536 vs 0.8539


**V2-grammar-19. vQTL loading**


$L_i=\sum_v \sqrt{\pi_v}\,\dfrac{s_{iv}-\bar s_v}{\mathrm{sd}(s_v)}$, $s_v = X_v e_v$


- *Reference:* Package choice (DGLM-style log link; Rönnegård & Valdar 2011, Murphy et al. 2022 cited for context only)

- *Function:* `.apply_vqtl`

- *Location:* R/grammar_realize.R:605-618 (616)


**V2-grammar-20. vQTL residual**


$\log \mathrm{Var}(E_v\mid g)=\text{const}+L_i$; $E_{vi}=z_i e^{L_i/2}$, then centred and scaled to $\sqrt{\pi_v}$


- *Variables:* z standardized normal

- *Reference:* as above

- *Function:* `.apply_vqtl`

- *Location:* R/grammar_realize.R:627-641

- *Notes:* verified numerically


**V2-grammar-21. Transcriptome score**


$c_T=\sum_g w_g z_g$, $z_g=(E_g-\bar E_g)/\mathrm{sd}(E_g)$, $w\leftarrow w/\max\vert{}w\vert{}$; genetic part same denominator


- *Reference:* SPEC-transcriptome (not this group)

- *Function:* `.tx_raw`, `.transcriptome_matrix`

- *Location:* R/grammar_realize.R:454-476 (470, 473, 475), 308-338


**V2-grammar-22. Mediation split**


$\tfrac{\mathrm{Var}(T_g)}{V_P},\tfrac{\mathrm{Var}(T_e)}{V_P},\tfrac{2\mathrm{Cov}(T_g,T_e)}{V_P}$


- *Reference:* —

- *Function:* `.mediation_budget`

- *Location:* R/grammar_realize.R:738-750 (604-606)

- *Notes:* stale in complex (F3)


**V2-grammar-23. Complex combination**


$G^{(c)}_t=\dfrac{\sum_m G^{(m)}_t-\overline{\cdot}}{\mathrm{sd}}\sqrt{h^2_t}$; $e\sim$ residual with $1-h^2_t$


- *Reference:* SPEC §4.3

- *Function:* `complex_phenotypes`

- *Location:* R/grammar_complex.R:122-134 (123, 127), 96-107 (99)

- *Notes:* inputs weighted by sd ($\propto$$\surd$prop)


**V2-grammar-24. MAF**


$p_j=\overline{(x_j+1)/2}$, $\mathrm{MAF}=\min(p,1-p)$


- *Reference:* —

- *Function:* `.marker_maf_ref` (via `.marker_stats_ref`)

- *Location:* R/grammar_simulate_phenotype.R:1022-1052 (1040, 1044)

- *Notes:* candidate = MAF>0 (F2)


**V2-grammar-25. Sub-seed**


$b \leftarrow (257\,b + u_k) \bmod 2147483629$ over the UTF-8 codes $u_k$ of the layer type ($b_0 = 0$); $s_\ell = (1009\,s + 7919\,b + 104729\,\text{occ}) \bmod (2^{31}-1)$


- *Reference:* Package scheme (SPEC §6)

- *Function:* `.layer_seed`

- *Location:* R/grammar_simulate_phenotype.R:1079-1094 (1083, 1092)

- *Notes:* no collisions found (audit, with the former $\sum \mathrm{utf8}$ type term); the type term is now an order-sensitive polynomial hash (`:1083-1086`)


**V2-grammar-26. Per-QTN variance**


$v_j = \dfrac{\pi_\ell}{\mathrm{sd}(c)^2}\,\dfrac{\mathrm{Var}(\text{col}_j)}{V_P^{real}}$


- *Reference:* Package convention

- *Function:* `.qtn_var`

- *Location:* R/io_write.R:1288-1312 (217-218)

- *Notes:* marginal; sums $\neq$ prop under LD/shared loci


**V2-grammar-27. Fixed-scale phenotype**


$y=g+e$, $e\sim N(0,\sigma^2_e)$, $\sigma^2_e=\mathrm{Var}(g_{ref})\dfrac{1-h^2}{h^2}$


- *Reference:* Falconer & Mackay 1996 ($h^2=V_A/(V_A+V_E)$)

- *Function:* `phenotype_value`

- *Location:* R/cross_population.R:1004-1104 (1061)

- *Notes:* true normal, not standardized (DECISION-021)


**V2-grammar-28. Genetic correlation (plot)**


$r=\mathrm{cor}(G_{\cdot1},G_{\cdot2})$


- *Reference:* —

- *Function:* `.plot_cor`

- *Location:* R/grammar_plot.R:144-162 (150)

- *Notes:* rep 1 only


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V2-grammar-C1. PleioArch shared covariance**


$\Sigma_{ii}=\pi_i V_i,\quad \Sigma_{ij}=r_{ij}\sqrt{V_iV_j}$


- *Variables:* $\pi_i$ shared (pleiotropic) fraction, $V_i$ layer variance (= prop), $r_{ij}$ target `cor`

- *Reference:* Prado et al., in preparation (PleioArch); page unverified; reference code absent from checkout

- *Function:* `.pleio_draw`, `.pleio_nonadditive_draw`

- *Location:* R/effects_pleioarch.R:65-93 (sigma at 78-80), 345-358 (356-357)

- *Notes:* UNVERIFIABLE against primary source; internal form explicit. DECISION-007/013


**V2-grammar-C2. PleioArch effect allocation (non-additive)**


$w_u\sim\mathrm{MVN}(0,\Sigma/m),\quad b_u = w_u/\mathrm{sd}(z_u)$ for the $m$ informative shared units; trait-specific units $\sim N(0,(1-\pi_t)V_t)$


- *Variables:* $z_u$ realized unit design column, $m$ shared units with sd > 0

- *Reference:* Prado et al., in preparation (additive form); non-additive extension is package design; page unverified

- *Function:* `.pleio_unit_effects`, `.draw_mvnorm`

- *Location:* R/effects_pleioarch.R:583-690 (535-537), 780-800

- *Notes:* Constant (hetless) units get effect 0 and are excluded from the allocation count


**V2-grammar-C3. PleioArch attainability**


$M\succeq 0,\ M_{ii}=\pi_i,\ M_{ij}=r_{ij}$; two traits: $r^2\le\pi_1\pi_2$


- *Variables:* variance-free feasibility matrix $M$

- *Reference:* Prado et al., in preparation (two-trait constraint); PSD generalisation is algebraic

- *Function:* `.check_pleio_feasible`

- *Location:* R/effects_pleioarch.R:870-917

- *Notes:* Multi-trait check on the variance-free matrix avoids tolerance masking


**V2-grammar-C4. Total pleiotropic correlation target**


$r^{tot}_{ij}=\dfrac{\sum_c r_{c,ij}\sqrt{V_{c,i}V_{c,j}}}{\sqrt{\sum_c V_{c,i}\,\sum_c V_{c,j}}}$


- *Variables:* $c$ = mean-effect layers (A, D, E), $V_{c,t}$ per-layer per-trait prop

- *Reference:* Package derivation (covariance addition + Cauchy–Schwarz); no external source

- *Function:* `.pleio_total_cor_check`

- *Location:* R/effects_pleioarch.R:724-788 (roxygen 543-574)

- *Notes:* Equals $r_{ij}$ iff per-trait variance profiles are proportional or $r_{ij}=0$; warns when off by > 1 %


**V2-grammar-C5. Additive MAF effect scaling**


$s_j = 1/\sqrt{2\,\mathrm{MAF}_j(1-\mathrm{MAF}_j)}$


- *Variables:* marker MAF; $s_j\leftarrow 1$ when MAF = 0 or non-finite

- *Reference:* Claimed PleioArch `scaleQTNEffects.R`; unavailable, page unverified

- *Function:* `.pleio_draw`

- *Location:* R/effects_pleioarch.R:201-217 (144-149)

- *Notes:* P3: implements the documented factor; source parity UNVERIFIABLE


**V2-grammar-C6. LD window (two-trait linked distinct loci)**


$r^2_{jk}=\mathrm{Cor}(X_j,X_k)^2\in[r^2_{\min},r^2_{\max}]$, same chromosome, causal loci disjoint


- *Variables:* dosage columns $X$; `ld_type` direct/indirect

- *Reference:* Package v2 implementation; Fernandes & Lipka 2020 *BMC Bioinformatics* 21:491 (v1 LD concept; page unverified for the v2 algorithm)

- *Function:* `.draw_qtn_ld`

- *Location:* R/arch_ld.R:62-225 (window 69-95, sampling 105-175)

- *Notes:* DECISION-014; indirect path also checks the causal pair's own $r^2$


## Effects and architectures: PleioArch pleiotropy, effect series, independent and LD architectures

*Source report:* `v2-effects-arch.md`  
*Entries:* 23


**V2-effects-arch-1. Pleiotropic covariance**


$\Sigma_{ii} = \pi_i V_i,\quad \Sigma_{ij} = \rho_{ij}\sqrt{V_i V_j}$


- *Variables:* $\pi_i$ shared share, $V_i$ = layer `prop` of trait i (phenotypic variance 1), $\rho_{ij}$ = `cor`

- *Reference:* `simulateEffects.R:34, 38-39` (Prado et al., in prep.; page n/a — bundled code)

- *Function:* `.pleio_draw`, `.pleio_nonadditive_draw`

- *Location:* `effects_pleioarch.R:90-91`, `:504-505`

- *Notes:* F1 bit-exact vs reference; generalised to n traits (DECISION-013)


**V2-effects-arch-2. Trait-specific variance**


$V^{spec}_i = (1-\pi_i) V_i$, i.i.d. per trait


- *Variables:* as above

- *Reference:* `simulateEffects.R:45-46`

- *Function:* `.pleio_draw`, `.pleio_unit_effects`

- *Location:* `:197-199`, `:685`

- *Notes:* zero-variance classes still drawn (EFF-F4)


**V2-effects-arch-3. Major/minor split**


$\Sigma_{maj} = \Sigma\,\phi,\ \Sigma_{min} = \Sigma(1-\phi)$; per-unit $\Sigma_{maj}/n_{maj}$, $\Sigma_{min}/n_{min}$


- *Variables:* $\phi$ = `prop_var_major`

- *Reference:* `simulateEffects.R:49-60, 69`

- *Function:* `.pleio_draw`, `.draw_mvnorm`

- *Location:* `:195-196`, `:934`

- *Notes:* additive layer only (`:63-68`)


**V2-effects-arch-4. MVN draw**


$E = Z\,\Sigma_{per}^{1/2},\ \Sigma_{per}^{1/2} = U\,\mathrm{diag}(\sqrt{\lambda_+})\,U^{\top}$


- *Variables:* $Z \sim N(0, I)$ $n\times n_t$

- *Reference:* reference uses `chol` (`:73`); symmetric root is standard linear algebra (page unverified)

- *Function:* `.draw_mvnorm`

- *Location:* `:929-939`

- *Notes:* same law, different realizations; handles singular Σ


**V2-effects-arch-5. Univariate draw**


$e_k \sim N\!\big(0,\ V^{spec}_i / n_{spec}\big)$


- *Reference:* `simulateEffects.R:78-83`

- *Function:* `.draw_univariate`

- *Location:* `:944-949`


**V2-effects-arch-6. Shared-unit count**


$n_{pleio} = \mathrm{round}\big(\bar\pi\, n\big)$, $n_{spec} = n - n_{pleio}$


- *Variables:* $\bar\pi$ = mean of `pi`

- *Reference:* package design (reference takes `pleioSnps` as input)

- *Function:* `.pleio_partition`

- *Location:* `:291-292, 320`

- *Notes:* R half-to-even (EFF-F9)


**V2-effects-arch-7. Attainability (2 traits)**


$\rho^2 \le \pi_1 \pi_2$


- *Reference:* `simulateEffects.R:17`

- *Function:* `.check_pleio_feasible`

- *Location:* `:876-890`

- *Notes:* tol `8\,\epsilon\max(\cdot)`


**V2-effects-arch-8. Attainability (n traits)**


$M \succeq 0,\ M_{ii} = \pi_i,\ M_{ij} = \rho_{ij}$; $\Sigma = D^{1/2} M D^{1/2}$


- *Variables:* $D = \mathrm{diag}(V)$

- *Reference:* congruence preserves PSD (standard; page unverified)

- *Function:* `.check_pleio_feasible`

- *Location:* `:898-915`

- *Notes:* E6 catches min eig -0.2


**V2-effects-arch-9. Allele -> genotype scaling**


$e^{*}_k = e_k / \sqrt{2\,\mathrm{MAF}_k(1-\mathrm{MAF}_k)}$


- *Variables:* MAF from $p = \overline{(g+1)/2}$, $\mathrm{MAF} = \min(p, 1-p)$

- *Reference:* `scaleQTNEffects.R:26-37`; $\mathrm{Var}(\mathrm{Bin}(2,p)) = 2p(1-p)$ (elementary)

- *Function:* `.pleio_draw`; `.marker_maf_ref` (via `.marker_stats_ref`)

- *Location:* `:201-217`; `grammar_simulate_phenotype.R:1040-1044`

- *Notes:* HWE sd of -1/0/1 dosage


**V2-effects-arch-10. Non-additive normaliser**


$e^{*}_u = e_u / \widehat{\mathrm{sd}}(z_u)$, $z_u$ = het indicator or $\prod_k (x_k - \bar x_k)$


- *Variables:* sample sd (n-1)

- *Reference:* DECISION-023 (own design, `:332`)

- *Function:* `.pleio_unit_effects`; `.epi_unit_column`

- *Location:* `:601-604, 686`; `grammar_realize.R:579-589`; dominance design `:518` vs `grammar_realize.R:510`

- *Notes:* constant columns -> effect 0, excluded from allocation `:454-457`


**V2-effects-arch-11. Effective per-component target**


$\rho^{eff}_{ij} = \Sigma_{ij}\,\mathbb{1}_{sh} / \sqrt{v^{eff}_i v^{eff}_j}$, $v^{eff}_i = \pi_i V_i \mathbb{1}_{sh} + (1-\pi_i)V_i \mathbb{1}_{spec,i}$


- *Variables:* indicators = live units exist

- *Reference:* DECISION-023 round 16

- *Function:* `.pleio_unit_effects`

- *Location:* `:659-664`


**V2-effects-arch-12. Total-correlation target**


$\rho^{tot}_{ij} = \dfrac{\sum_c \rho^{(c)}_{ij}\sqrt{V_{ci}V_{cj}}}{\sqrt{\sum_c V_{ci}\sum_c V_{cj}}} \le \rho_{ij}$ (Cauchy–Schwarz)


- *Variables:* c over A/D/E layers

- *Reference:* DECISION-023 §3 (own derivation)

- *Function:* `.pleio_total_cor_check`

- *Location:* `:754-759`

- *Notes:* warns if `abs(tot - cor) > 0.01*abs(cor)` `:616`; E12 gives 0.400 for (.4,.1)/(.1,.4) at 0.5


**V2-effects-arch-13. Layer rescale (context)**


$c_t \leftarrow (c_t - \bar c_t)\,\sqrt{\mathrm{prop}_t}/\widehat{\mathrm{sd}}(c_t)$


- *Reference:* SPEC §2

- *Function:* `.genetic_matrix`

- *Location:* `grammar_realize.R:144-146`

- *Notes:* makes realized r a random ratio (P1 wording)


**V2-effects-arch-14. Few-unit attenuation**


$\mathbb{E}[r]/\rho \approx 0.66, 0.82, 0.92, 0.96, 0.98, 0.99$ at $K = 1,2,5,10,20,60$ ($\rho = .5$, $\pi = 1$, orthogonal designs)


- *Variables:* K shared units

- *Reference:* numerical (F2, 40 000 draws); no closed form cited

- *Function:* doc claim

- *Location:* `:26-28`; `grammar_simulate_phenotype.R:318-320`

- *Notes:* verified


**V2-effects-arch-15. Complete-LD ensemble mean**


$\mathbb{E}[r] = 2\arcsin(\rho)/\pi$ from $P(XY>0) = \tfrac12 + \arcsin(\rho)/\pi$


- *Variables:* bivariate normal orthant

- *Reference:* Sheppard's orthant formula (page unverified)

- *Function:* doc claim

- *Location:* `:28-30`; DECISIONS.md 023 Scope

- *Notes:* 1/3 at ρ = .5; E7 realized ±1


**V2-effects-arch-16. Geometric effect series**


$a_k = b^{k},\ k = 1..n$ (default $b = 0.5$)


- *Reference:* v1 `sim_method = "geometric"` (Fernandes & Lipka 2020; page unverified)

- *Function:* `.effect_series`

- *Location:* `effects_series.R:55, 64`

- *Notes:* custom vector used verbatim `:32-38`


**V2-effects-arch-17. Repulsion phase**


$a_k \leftarrow a_k(-1)^{k+1}$, same for every trait


- *Reference:* package design

- *Function:* `.apply_phase`

- *Location:* `grammar_layers.R:1482-1490`

- *Notes:* applied to PleioArch effects too (EFF-F11)


**V2-effects-arch-18. Residual draw**


$e \sim N(0, \sigma^2_e)$, then $e \leftarrow (e-\bar e)\sigma_e/\widehat{\mathrm{sd}}(e)$, $\sigma^2_e = 1 - \sum \mathrm{prop}$


- *Reference:* SPEC §2

- *Function:* `.draw_residual`

- *Location:* `effects_series.R:85-102`

- *Notes:* exact-variance residual under the default `residual_mode = "fixed"`; `"random"` skips the rescaling (`:94`)


**V2-effects-arch-19. LD measure**


$r^2 = \mathrm{cor}(g_a, g_b)^2$ on -1/0/1 dosage (composite)


- *Variables:* same-chromosome candidates

- *Reference:* standard Pearson (genotypic/“composite” r²; page unverified)

- *Function:* `window_partners`, `r2_pair`

- *Location:* `arch_ld.R:111, 124-126`

- *Notes:* inclusive window `:82`


**V2-effects-arch-20. Direct pair choice**


$t_2 = \arg\max_{j \in W(t_1)} r^2_{t_1 j}$ (default `partner = "strongest"`); $t_2 \sim U\{W(t_1)\}$ for `partner = "random"`


- *Variables:* $W$ = in-window unused partners

- *Reference:* package design (undocumented)

- *Function:* `.draw_qtn_ld`

- *Location:* `arch_ld.R:149-160` (`:151-152`); `partner` `:73-74`

- *Notes:* EFF-F2


**V2-effects-arch-21. Indirect flanks**


$\mathrm{pos}(u) < \mathrm{pos}(c) < \mathrm{pos}(d)$, $r^2_{uc}, r^2_{dc}, r^2_{ud} \in [r^2_{min}, r^2_{max}]$, strongest-to-cause first


- *Variables:* c = hidden cause

- *Reference:* DECISION-014

- *Function:* `.draw_qtn_ld`

- *Location:* `arch_ld.R:161-197`

- *Notes:* reported r² = $r^2_{ud}$ `:170-174` (EFF-F7)


**V2-effects-arch-22. Layer sub-seed**


$s = (1009\,\mathrm{seed} + 7919 \sum \mathrm{utf8}(\mathrm{type}) + 104729\,\mathrm{occ}) \bmod (2^{31}-1)$


- *Reference:* package design

- *Function:* `.layer_seed`

- *Location:* `grammar_simulate_phenotype.R:1083-1093`

- *Notes:* permutation-invariant (EFF-F1)


**V2-effects-arch-23. Gabriel blocks tag**


keep $\arg\max_{m \in B}\mathrm{MAF}_m$ per block $B$; MAF floor 0.05; window $\lfloor 1000\,kb\,(1+\epsilon)\rfloor$ bp


- *Reference:* Gabriel et al. 2002 Science 296:2225-2229 (page unverified); PLINK 1.9 constants (unverified)

- *Function:* `.gabriel_blocks`

- *Location:* `qc_ld_methods.R:24-29, 38-43`

- *Notes:* UNVERIFIABLE parity


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V2-effects-arch-C1. Zero-variance correlation guard**


$V_i = 0 \Rightarrow \mathrm{Cor}(g_i, g_j) = 0/0$ undefined


- *Variables:* $V_i$ = `vg[i]` (layer `prop`), $R_{ij}$ = requested `cor`

- *Reference:* Definition of correlation; no external source

- *Function:* `.pleio_check_zero_var()`

- *Location:* `R/effects_pleioarch.R:246-275`

- *Notes:* nonzero requested pair with a zero-variance trait -> error (`:186`); zero requested pair -> warning, reported NA (`:193`)


**V2-effects-arch-C2. Shared/specific unit count (clamped)**


$n_{sh} = \min\!\big(n, \max(0, \mathrm{round}(n\,\bar\pi))\big),\ n_{sp} = n - n_{sh}$


- *Variables:* $n$ = `n_qtn`/`n_pairs`; $\bar\pi$ = `mean(pi)`

- *Reference:* Package generalisation (reference takes `pleioSnps` as input)

- *Function:* `.pleio_partition()`

- *Location:* `R/effects_pleioarch.R:288-329` (clamp `:291-292`)

- *Notes:* refines Fable's row (which omits the clamp); errors when a requested class rounds to zero (`:241-256`); R half-to-even (EFF-F9)


**V2-effects-arch-C3. One shared unit**


$g_i = a_i z \Rightarrow \lvert\mathrm{Cor}(g_i, g_j)\rvert = 1$ when neither trait has specific variance; otherwise one noisy draw around `cor`


- *Variables:* $z$ = the single shared design column; $a_i$ = trait-$i$ effect

- *Reference:* Algebraic consequence; no external source

- *Function:* `.pleio_single_unit_consequence()`

- *Location:* `R/effects_pleioarch.R:412-441`

- *Notes:* pairs with $\lvert R_{ij}\rvert \ge 1$ or zero variance skipped (`:277-278`); feeds the warning at `:235-241`


**V2-effects-arch-C4. Correlation input expansion**


scalar $r \mapsto R,\ R_{ii} = 1,\ R_{ij} = r$; matrix must be $n_t \times n_t$, finite, in $[-1, 1]$, symmetric, unit diagonal


- *Variables:* `cor`, `n_traits`

- *Reference:* Correlation-matrix definition; no external source

- *Function:* `.pleio_cor_matrix()`

- *Location:* `R/effects_pleioarch.R:957-986`

- *Notes:* PSD is checked later by `.check_pleio_feasible()` after the diagonal is replaced by $\pi_i$


**V2-effects-arch-C5. Pleiotropic-share argument mapping**


$\pi = (\pi_T, \pi_S)$ from `pi_target`/`pi_secondary`, or scalar/length-$n_t$ `pi`; all in $[0, 1]$


- *Variables:* `pi`, `pi_target`, `pi_secondary`

- *Reference:* Bundled PleioArch interface (`simulateEffects.R:10`: `piT`, `piS`); Prado et al. in prep., page unverified

- *Function:* `.pleio_pi_vector()`

- *Location:* `R/effects_pleioarch.R:829-860`

- *Notes:* spellings mutually exclusive (`:684-686`); two-trait spelling rejected for $n_t > 2$ (`:689-691`)


**V2-effects-arch-C6. Architecture-specific QTN sampling**


independent: per-trait SRSWOR from $C$; pleiotropy: one shared SRSWOR replicated to every trait; ld: delegated to `.draw_qtn_ld()`


- *Variables:* $C$ = candidate markers

- *Reference:* Fernandes & Lipka 2020 BMC Bioinf 21:491 (page unverified); package v2 design

- *Function:* `.draw_qtn()`

- *Location:* `R/arch_independent.R:21-45`

- *Notes:* sub-seed set and restored (`:27-31`); the pleiotropy branch `:34-37` is reached only by `vqtl(same_as_add = FALSE)` (EFF-F8) — the additive PleioArch path uses `.pleio_draw()`


**V2-effects-arch-C7. Distinct-chromosome allocation**


chromosome $k \mapsto 1 + ((k-1) \bmod T)$ (first-occurrence order), then SRSWOR within the assigned pool


- *Variables:* $T$ = `n_traits`

- *Reference:* Package design; no external source

- *Function:* `.draw_qtn_distinct_chr()`

- *Location:* `R/arch_independent.R:56-75`

- *Notes:* errors if fewer than $T$ chromosomes (`:60-63`) or a pool is short (`:68-72`)


**V2-effects-arch-C8. Epistatic-set sampling**


draw $n_{pairs} \times \mathrm{interaction}$ markers SRSWOR, reshape `byrow`


- *Variables:* candidate markers; `interaction`

- *Reference:* Package design; no external source

- *Function:* `.draw_qtn_pairs()`

- *Location:* `R/arch_independent.R:81-103`

- *Notes:* pleiotropy branch `:98-101` unreachable (EFF-F8); PleioArch non-additive uses `.pleio_units()`


**V2-effects-arch-C9. Candidate-marker rule**


$C = \{ j : \mathrm{isfinite}(\mathrm{MAF}_j) \wedge \mathrm{MAF}_j > 0 \}$


- *Variables:* `sim$maf`

- *Reference:* Package design

- *Function:* `.candidate_markers()`

- *Location:* `R/arch_independent.R:114-134`

- *Notes:* makes the `s == 0` guard at `effects_pleioarch.R:205` dead


**V2-effects-arch-C10. Layer occurrence index**


$o_t = \sum_l \mathbb{1}(\mathrm{type}_l = t)$


- *Variables:* prior layers, requested type

- *Reference:* Package RNG design (grammar seed contract)

- *Function:* `.type_occurrence()`

- *Location:* `R/arch_independent.R:139-141`

- *Notes:* feeds `.layer_seed()`; type-local, so inserting a different layer type does not perturb an additive draw


## Crossing: populations, meiosis draws, maps, pedigree, mating designs, crossbreeding

*Source report:* `v2-crossing.md`  
*Entries:* 24


**V2-crossing-1. Heterozygote phasing**


$\text{cis}_j = \mathbf 1(x_j \ge 0),\ \text{trans}_j = \mathbf 1(x_j > 0)$


- *Variables:* $x_j \in \{-1,0,1\}$ dosage

- *Reference:* package convention (documented assumption)

- *Function:* `as_population`

- *Location:* `R/cross_population.R:176-177`

- *Notes:* every het $\rightarrow$ allele 1 on cis (coupling)


**V2-crossing-2. Dosage from strands**


$x_j = \text{cis}_j + \text{trans}_j - 1$


- *Variables:* strands $\in$ {0,1}

- *Reference:* isqg `genotype_num` (Genetics.cpp:363)

- *Function:* `dosages`

- *Location:* `R/cross_population.R:623`; `genome.rs:273-278`

- *Notes:* integer output


**V2-crossing-3. Crossover count**


Default (`interference = NULL`, since 2.0.0.9003): two-pathway gamma model with $\nu = 2.6$, $p = 0$. Chiasmata on the bivalent at 2 per Morgan; a share $p$ is a Poisson process (drawn at the gamete level, $n \sim \text{Poisson}(pL)$, positions $U(0,L)$); the rest is a stationary renewal process with gaps $\sim \text{Gamma}(\nu,\ 2\nu(1-p))$, first point $U\cdot\text{Gamma}(\nu+1,\ 2\nu(1-p))$, each chiasma kept in the gamete with probability $1/2$. `interference = "poisson"` (the default before 2.0.0.9003, isqg stream): $n_x \sim \text{Poisson}(L)$. $L = \text{last cM}/100$ in both


- *Variables:* L in Morgans

- *Reference:* Poisson branch: Karlin & Liberman 1978 PNAS 75:6332–6336 (as cited by isqg; not re-verified); isqg Genetics.cpp:84. Gamma model: McPeek & Speed 1995 *Genetics* 139(2):1031–1044; two-pathway extension Housworth & Stahl 2003 *Am. J. Hum. Genet.* 73(1):188–197 (as cited in `man/cross.Rd`; pages not verified here)

- *Function:* `.draw_meiosis`, `.draw_meiosis_interference`, `.check_interference`

- *Location:* `R/cross_mating.R:53-56` (dispatch), Poisson branch `:77-78`; gamma model `:215-275`; default `.INTERFERENCE_DEFAULT` `:117` (DECISION-047); cM$\rightarrow$M `:339,341`

- *Notes:* last position, not span (CROSS-F5), in both models; seeded meiosis results changed with the new default (NEWS, 2.0.0.9003)


**V2-crossing-4. Chiasma positions**


$u_1,\dots,u_{n_x} \overset{iid}{\sim} U(0,L)$, sorted


- *Reference:* isqg Genetics.cpp:93-95

- *Function:* `.draw_meiosis`

- *Location:* `R/cross_mating.R:79-87`

- *Notes:* not drawn when $n_x=0$; Poisson branch only (`interference = "poisson"`); default gamma model in V2-crossing-3


**V2-crossing-5. Strand choice**


$f \sim \text{Bernoulli}(1/2)$, always drawn


- *Reference:* isqg Genetics.cpp:74

- *Function:* `.draw_meiosis`

- *Location:* `R/cross_mating.R:90`

- *Notes:* unconditional; the gamma branch also draws one flip per slot (`:272`)


**V2-crossing-6. Ancestry mask**


$m_j = f \oplus \bigoplus_{k} \mathbf 1(j \ge b_k),\ b_k = \#\{p_i \le u_k\}$


- *Variables:* $p_i$ marker positions

- *Reference:* isqg Genetics.cpp:65-69, 74-75

- *Function:* `chromosome_mask`, `breaks_at`

- *Location:* `src/rust/src/meiosis.rs:49-51, 63-72`

- *Notes:* `<=` (upper_bound)


**V2-crossing-7. Gamete**


$g_j = m_j\,\text{cis}_j + (1-m_j)\,\text{trans}_j$


- *Reference:* isqg Genetics.cpp:356

- *Function:* `recombine`

- *Location:* `meiosis.rs:224-237`


**V2-crossing-8. Cross / self / DH**


cross: (cis,trans) = (g^{(1)}_{p_1}, g^{(2)}_{p_2}); self: same with $p_1=p_2$, 2 independent meioses; DH: (g, g)


- *Reference:* isqg Mating.cpp:91-95, 112-116, 127-130; Genetics.cpp:521-541

- *Function:* `mate_haplotypes`

- *Location:* `meiosis.rs:350-371`; events per progeny `R/cross_mating.R:345`

- *Notes:* progeny-major, parent-1 then parent-2


**V2-crossing-9. Recombination fraction (implied, verified)**


$r = \tfrac12\left(1-e^{-2d}\right)$


- *Variables:* d Morgans

- *Reference:* Haldane 1919 J. Genet. 8:299–309 (page unverified; standard)

- *Function:* — (property of the Poisson process; holds for `interference = "poisson"`, or $\nu = 1$ / $p = 1$, not for the default gamma model, whose $r(d) = [1 - P_0(d)]/2$ is given in `man/cross.Rd`)

- *Location:* verified e4.R/e5.R; `tests/testthat/test-cross.R:179-202` (file runs under `simplePHENOTYPES.interference = "poisson"`, `:13`)

- *Notes:* holds for any map origin


**V2-crossing-10. Heterozygosity under selfing**


$H_t = H_0 / 2^t$


- *Reference:* Falconer & Mackay 1996 ch. 5 (page unverified)

- *Function:* `selfcross` (doc)

- *Location:* `R/cross_mating.R:627-629`; verified e1.R E3


**V2-crossing-11. Synthetic map rate**


$r(p) = 1 - s\exp\!\left(-\tfrac12\left(\tfrac{p-c}{wL}\right)^2\right)$


- *Variables:* s suppression, w width, c centromere, L bp span

- *Reference:* phenomenological (package); scale $\approx$ Bauer et al. 2013 GB 14:R103 (UNVERIFIABLE numbers)

- *Function:* `synthetic_map`

- *Location:* `R/cross_map.R:157`


**V2-crossing-12. Synthetic map integration**


$\text{cM}_i = \dfrac{\sum_{k<i}(p_{k+1}-p_k)\,\tfrac{r_k+r_{k+1}}2}{\sum_{k}(\cdot)}\cdot \text{len}_{cM}$; $\text{len}_{cM} = \text{cM/Mb}\cdot L/10^6$ if `total_cm` NULL


- *Reference:* package

- *Function:* `synthetic_map`

- *Location:* `R/cross_map.R:168-173`, `:133`

- *Notes:* trapezoid; monotone


**V2-crossing-13. Additive value (fixed scale)**


$A_i = \sum_j x_{ij}\,\beta_j$


- *Reference:* Falconer & Mackay 1996 (a on −1/0/1 scale; page unverified)

- *Function:* `additive_value`

- *Location:* `R/cross_population.R:746`

- *Notes:* no centring


**V2-crossing-14. Genotypic value**


$G_i = \sum_j [a_j x_{ij} + d_j\,\mathbf 1(x_{ij}=0)]$


- *Reference:* Falconer & Mackay 1996 ch. 7 (page unverified)

- *Function:* `genotypic_value`

- *Location:* `R/cross_population.R:821`


**V2-crossing-15. Residual from h²**


$\sigma_e^2 = \text{Var}(g_{\text{ref}})\,(1-h^2)/h^2$


- *Reference:* $h^2 = V_A/(V_A+V_E)$, F&M 1996 (page unverified)

- *Function:* `phenotype_value`

- *Location:* `R/cross_population.R:1061`


**V2-crossing-16. Breed composition**


$F_{i,b} = \tfrac12(F_{m(i),b} + F_{f(i),b})$, founders $F = \mathbf 1(\text{pool}=b)$


- *Reference:* SPEC-block3b §7 (package)

- *Function:* `breed_composition`

- *Location:* `R/cross_breed.R:37-43`

- *Notes:* expected, not realized


**V2-crossing-17. Realized heterosis**


$H = \bar G_{\text{cross}} - \sum_b \bar F_b\,\bar G_b$


- *Variables:* $\bar F_b$ mean composition

- *Reference:* mid-parent heterosis, F&M 1996 ch. 14 (page unverified)

- *Function:* `heterosis`

- *Location:* `R/cross_breed.R:135-139, 152`

- *Notes:* composition-weighted


**V2-crossing-18. Expected F1 heterosis**


$H_{F1} = \sum_j d_j\,[h_{AB,j} - \tfrac12(h_{A,j}+h_{B,j})]$, $h_{AB}=p_A(1-p_B)+(1-p_A)p_B$; HWE: $\sum_j d_j (p_{Aj}-p_{Bj})^2$; inbred lines: $\sum_j d_j(p_A+p_B-2p_Ap_B)$


- *Variables:* $p$ = allele-1 gamete freq

- *Reference:* F&M 1996 $H_{F1}=\sum d y^2$ (page unverified); general form = package derivation

- *Function:* `heterosis` + `.expected_cross_means`

- *Location:* `R/cross_breed.R:140-150`; `R/select_combining.R:244-255`

- *Notes:* verified exactly (E12)


**V2-crossing-19. Per-locus cross mean**


$E = a(g_i+g_k-1) + d(g_i+g_k-2g_ig_k)$


- *Reference:* expansion of $p_{AA}a + p_{Aa}d - p_{aa}a$

- *Function:* `.expected_cross_means`

- *Location:* `R/select_combining.R:249-254`


**V2-crossing-20. Retention fractions (doc only)**


F2, BC: 1/2 H_F1; 2-breed rotation: 2/3; n-breed: $(2^n-2)/(2^n-1)$


- *Reference:* SPEC §7 derivation (breed-origin dominance model), **requires within-breed HWE**

- *Function:* roxygen

- *Location:* `R/cross_breed.R:68-72`; `docs/SPEC-block3b.md` §7

- *Notes:* CROSS-F1


**V2-crossing-21. Rotation sire sequence**


sire$_g$ = breed $((g-1) \bmod n_b)+1$


- *Reference:* —

- *Function:* `crossbreed`

- *Location:* `R/cross_breed.R:302-305`


**V2-crossing-22. Pedigree key**


key = FNV-1a-128(canonical(design, keys$_{p_1}$, keys$_{p_2}$, RNG state before/after, counts, flips, chiasma bits, n)) ‖ "_i"


- *Reference:* FNV spec (128-bit params)

- *Function:* `.mating_pedigree`, `.stable_key`

- *Location:* `R/cross_pedigree.R:25-37` (`.stable_key`), `:42-59` (`.key_part`), `:205-207` (mating key); `hash.rs:12-23`

- *Notes:* founder key `:97-138` (`.founder_pedigree`, batched since the post-audit perf change)


**V2-crossing-23. Generation**


$\text{gen}_i = 1 + \max(\text{gen}_{m}, \text{gen}_{f})$


- *Reference:* —

- *Function:* `.mating_pedigree`

- *Location:* `R/cross_pedigree.R:208-209`


**V2-crossing-24. A-matrix (other group, consistency only)**


$A_{ij} = \tfrac12(A_{i,m(j)} + A_{i,f(j)})$, $A_{jj} = 1 + \tfrac12 A_{m(j)f(j)}$, DH: 2


- *Reference:* Wright 1922; Emik & Terrill 1949 (pages as cited in roxygen, not re-verified)

- *Function:* `a_matrix`

- *Location:* `R/select_blup.R:69-90`

- *Notes:* consistent with this group's pedigree (E10)


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V2-crossing-C1. Average-effect vs genotypic-value distinction (doc)**


$\alpha_j = a_j + d_j(1-2p_j),\quad A_i = \sum_j \alpha_j (x_{ij} - 2p_j)$


- *Variables:* $p_j$ counted-allele frequency, $x_{ij}$ gene count

- *Reference:* Falconer & Mackay 1996 (secondary; page unverified); Fisher 1918 average effect (not cited in roxygen)

- *Function:* roxygen of `genotypic_value` (implemented in `select_ind(on = "bv")`, other group)

- *Location:* `R/cross_population.R:774-780` (verified: α text at 776-777)

- *Notes:* Doc correctly states `additive_value(alpha)` is only ranking-equivalent (differs by an additive constant)


**V2-crossing-C2. Full-sib family identity**


$\text{fam}(i) = \{k_{m(i)}, k_{f(i)}\}$ (unordered key pair); maternal/paternal half-sib = $k_m$ / $k_f$; selfed = $k_m$ iff design $\in$ {self, dh}


- *Variables:* stable parent keys $k$

- *Reference:* package bookkeeping (no external source needed)

- *Function:* `families`

- *Location:* `R/cross_pedigree.R:307-336` (verified: `families <-` at 307, closes at 336)

- *Notes:* `ifelse(m <= f, paste(m,f), paste(f,m))` merges reciprocals; labels made unique


**V2-crossing-C3. Mating designs (counts)**


factorial $\vert{}M\vert{}\times\vert{}F\vert{}$ (minus selfs); diallel $n(n-1)$ (+$n$ with selfs); half-diallel $\binom{n}{2}$ (+$n$); nested: $m$ distinct mothers per father via bipartite matching; random: mother uniform, then admissible father uniform


- *Variables:* parent keys, `progeny_per_cross`, `allow_self`

- *Reference:* Griffing 1956 Aust. J. Biol. Sci. 9:463–493 for the diallel methods (standard citation, pages not verified here; not cited in roxygen); NC designs I/II uncited

- *Function:* `mating_design`

- *Location:* `R/cross_mate.R:360-486` (verified: `mating_design <-` at 229; file is 350 lines)

- *Notes:* Fable's rubric line covers the counts (e2.R E13) but Fable's equation map omitted the row


**V2-crossing-C4. Within-family relationship validation (test)**


expected BV correlations: S1 sibs $2/3$ ($A_{ij}=1$, $A_{ii}=1.5$), full sibs and DH $1/2$; $\text{ICC} = \dfrac{MS_B - MS_W}{MS_B + (k-1)MS_W}$


- *Variables:* $k$ balanced family size, $MS_B$, $MS_W$

- *Reference:* Wright 1922 (relationship; page unverified); one-way ANOVA ICC estimator (Fisher 1925 / standard; uncited)

- *Function:* test helper `.icc`

- *Location:* `tests/testthat/test-family-relationship.R:1-58` (verified: 58 lines, `.icc` at 26, targets at 52-57)

- *Notes:* validates the crossing engine's realized relationships, not S3 index weights


## Rust core: count-location meiosis, genome masks, numericalization (isqg parity)

*Source report:* `v2-rust-core.md`  
*Entries:* 16


**V2-rust-core-1. Crossover count (count-location)**


$n_x \sim \mathrm{Poisson}(L)$ (`interference = "poisson"` only; the default gamma model is V2-crossing-3)


- *Variables:* $L$ = last map position of the chromosome, Morgans

- *Reference:* Karlin & Liberman 1978, PNAS 75(12):6332–6336 (as cited in isqg `Genetics.cpp:79`; page unverified by me)

- *Function:* `.draw_meiosis`

- *Location:* R `R/cross_mating.R:77-78`; isqg `Genetics.cpp:44,84`; Rust: none (input `counts`)

- *Notes:* $L$ is NOT the span (last − first); RUST-F10 docstring


**V2-rust-core-2. Chiasma locations**


$x_i \overset{iid}{\sim} U(0, L),\; i=1..n_x$, sorted


- *Variables:* $x_i$

- *Reference:* Karlin & Liberman 1978 (as above)

- *Function:* `.draw_meiosis`

- *Location:* R `:79-87`; isqg `Genetics.cpp:86-97`

- *Notes:* R `runif` is on the open interval; Rust never sorts (`meiosis.rs:60-62`)


**V2-rust-core-3. Strand raffle**


$f \sim \mathrm{Bernoulli}(1/2)$, drawn unconditionally


- *Variables:* $f$

- *Reference:* isqg convention (`Genetics.cpp:74`)

- *Function:* `.draw_meiosis` / `chromosome_mask`

- *Location:* R `:90`; Rust `meiosis.rs:68-70`; isqg `:74-75`

- *Notes:* skipping it when $n_x=0$ desynchronises the stream


**V2-rust-core-4. Breakpoint rank**


$b(x) = \#\{j : p_j \le x\}$


- *Variables:* $p_j$ ascending map positions

- *Reference:* `std::upper_bound` semantics (isqg)

- *Function:* `breaks_at`

- *Location:* Rust `meiosis.rs:49-51`; isqg `Genetics.cpp:67`

- *Notes:* 0-based; $b = n$ $\Rightarrow$ no-op; marker exactly at $x$ stays upstream


**V2-rust-core-5. Ancestry mask**


$m_j = f \oplus \bigoplus_{i} \mathbb{1}[\,j \ge b(x_i)\,]$


- *Variables:* $m_j \in\{0,1\}$

- *Reference:* isqg `Genetics.cpp:54-77`

- *Function:* `chromosome_mask`, `Bits::toggle_from`, `flip_all`

- *Location:* Rust `meiosis.rs:63-72`, `genome.rs:69-87`; isqg `:58-70` (XOR loop), `:74-75` (flip)

- *Notes:* duplicates cancel; order irrelevant


**V2-rust-core-6. Gamete assembly**


$g_j = m_j\, c_j + (1-m_j)\, t_j$


- *Variables:* $c,t$ = parental cis/trans strands

- *Reference:* isqg `DNA::recombination`

- *Function:* `recombine`

- *Location:* Rust `meiosis.rs:224-237`; isqg `Genetics.cpp:356`


**V2-rust-core-7. Cross**


progeny $i$: $(\text{cis},\text{trans}) = (g^{(P_1)}_{2i},\, g^{(P_2)}_{2i+1})$


- *Variables:* event index progeny-major

- *Reference:* isqg `Mating.cpp:92-94`

- *Function:* `mate_haplotypes`

- *Location:* Rust `meiosis.rs:359-369`; R `.mate` `cross_mating.R:345,453,460-467`

- *Notes:* cis from parent 1 (tested by `cross_swapped`)


**V2-rust-core-8. Self**


as cross with $P_2 = P_1$, two independent meioses


- *Reference:* isqg `Mating.cpp:113-115`

- *Function:* `selfcross` $\rightarrow$ `.mate(parent, parent)`

- *Location:* R `cross_mating.R:647-652`; Rust same path


**V2-rust-core-9. Doubled haploid**


$(\text{cis},\text{trans}) = (g^{(P)}_{i}, g^{(P)}_{i})$


- *Variables:* one event per progeny

- *Reference:* isqg `Mating.cpp:128-129` $\rightarrow$ `duplication` `Genetics.cpp:532-541`

- *Function:* `mate_haplotypes` (Dh)

- *Location:* Rust `meiosis.rs:263-269,351-358`; R `cross_mating.R:686-691`


**V2-rust-core-10. Genotype projection**


$G_j = \begin{cases}1 & c_j \wedge t_j\\ 0 & c_j \oplus t_j\\ -1 & \neg c_j \wedge \neg t_j\end{cases}$


- *Reference:* isqg `DNA::genotype_num` `Genetics.cpp:358-369`

- *Function:* `write_genotype`, `dosages`

- *Location:* Rust `genome.rs:273-284`; R `cross_population.R:623` (`cis + trans - 1`)

- *Notes:* lossy (phase)


**V2-rust-core-11. Founder phasing**


$c_j = \mathbb{1}[d_j \ge 0],\; t_j = \mathbb{1}[d_j > 0]$


- *Variables:* $d_j\in\{-1,0,1\}$ dosage

- *Reference:* isqg `founder("Aa")` `Genetics.cpp:307-309`

- *Function:* `as_population`

- *Location:* R `cross_population.R:176-177`

- *Notes:* all hets coupling-phased (RUST-F11)


**V2-rust-core-12. Map units**


$p_j = \mathrm{cM}_j / 100$


- *Reference:* definition

- *Function:* `.mate`

- *Location:* R `cross_mating.R:339,341`


**V2-rust-core-13. Haldane map function (consequence, tested only)**


$r(d) = \tfrac{1}{2}\left(1 - e^{-2d}\right)$


- *Variables:* $d$ Morgans between two markers

- *Reference:* Haldane 1919, J. Genet. 8:299–309 (page unverified)

- *Function:* —

- *Location:* test `tests/testthat/test-cross.R:179-202`; auditor `T2–T4`

- *Notes:* follows from a rate-1 Poisson process with no interference (`interference = "poisson"`; not the default gamma model)


**V2-rust-core-14. Numericalization — orientation**


allele-2 major $\iff n_{2} > n_{0}$ ($\iff 2n_2+n_1 > 2n_0+n_1$)


- *Variables:* $n_k$ = count of raw dosage $k$

- *Reference:* package contract (`io_as_numeric.R:8-12`)

- *Function:* `compute_flip`

- *Location:* R `io_detect_format.R:362-364`; reference mode `:354-359`

- *Notes:* strict `>`: tie $\rightarrow$ allele-1 (RUST-F9)


**V2-rust-core-15. Numericalization — codes**


`-101`: (major, het, minor) = (1,0,−1); `012`: (2,1,0); Dom: non-het$\rightarrow$minor; Left: het$\rightarrow$minor; Right: het$\rightarrow$major; impute Middle/Minor/Major $\rightarrow$ het/minor/major class


- *Reference:* package contract; v1 `numericalization.R` (Dom/Left/Right lines)

- *Function:* `numericalize_core`

- *Location:* Rust `numeric.rs:76-88,135-180`

- *Notes:* model applied after imputation (fixes a v1 gap, `numeric.rs:127-134`)


**V2-rust-core-16. Pedigree key hash**


$h_0 = \text{offset}_{128};\; h_{k+1} = (h_k \oplus b_k)\cdot p \bmod 2^{128}$, $p = 2^{88}+2^{8}+\mathrm{0x3b}$


- *Variables:* $b_k$ UTF-8 bytes

- *Reference:* Fowler–Noll–Vo FNV-1a (IETF draft-eastlake-fnv; section unverified)

- *Function:* `fnv1a_128`, `.stable_key`

- *Location:* Rust `hash.rs:12-23`; R `cross_pedigree.R:22-59`

- *Notes:* verified vs Python (RUST-F15)


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V2-rust-core-C1. Public ingestion dispatch**


$\texttt{as\_numeric}(x) = \texttt{format\_conversion}(x, \texttt{to = "numeric"})$ after path/object disambiguation


- *Variables:* $x$: bare character scalar path, or in-memory genotype object (character matrix, data.frame, …)

- *Reference:* Package API contract; no external theory source

- *Function:* `as_numeric()`

- *Location:* `R/io_as_numeric.R:192-207` (verified: `if (is.character(x) && is.null(dim(x)))` at :198, `format_conversion(...)` at :206)

- *Notes:* Deterministic; character matrices deliberately not treated as paths (`:101-104`); F4 lives downstream in `io_detect_format.R:286-297`


**V2-rust-core-C2. Genetic-model post-transform (split out of Fable's single "codes" row)**


Dom: $z = c_H$ if het else $c_m$; Left: het $\to c_m$; Right: het $\to c_M$; Add: unchanged


- *Variables:* $z$: output code; $c_M, c_H, c_m$: coding-specific major/het/minor values

- *Reference:* Package-defined simulation coding (v1 `numericalization.R` Dom/Left/Right lines); no external theory source

- *Function:* `numericalize_core()`

- *Location:* `src/rust/src/numeric.rs:149-180` (verified: `match model` at :156, `"Dom"` :157, `"Left"` :164, `"Right"` :171, default :178)

- *Notes:* Applied *after* imputation (`:64-71`); Codex's audit-only 64-case grid passed; unknown `model` silently = Add (F7/RC-04)


## Selection engine: truncation, indices (Smith-Hazel, Lush, QGSI), culling, schemes

*Source report:* `v2-selection.md`  
*Entries:* 18


**V2-selection-1. Truncation count from prop**


$k=\max(1,\operatorname{round}(pN))$


- *Variables:* p prop, N candidates

- *Reference:* — (implementation)

- *Function:* `.resolve_keep`

- *Location:* `R/select_ind.R:574`

- *Notes:* R half-to-even


**V2-selection-2. Selection intensity (large N)**


$i(p)=\dfrac{\varphi(\Phi^{-1}(1-p))}{p}$


- *Variables:* p = k/N

- *Reference:* Falconer & Mackay 1996 (page unverified)

- *Function:* `.resolve_keep`

- *Location:* `R/select_ind.R:585–587`

- *Notes:* count with nearest i; infinite-N form


**V2-selection-3. Selection differential**


$S=\bar{y}_{sel}-\bar{y}$


- *Variables:* y ranking criterion (sign-flipped for low)

- *Reference:* Falconer & Mackay 1996 (page unverified)

- *Function:* `select_ind`

- *Location:* `R/select_ind.R:460, 482`

- *Notes:* natural sign restored for low


**V2-selection-4. Realized intensity**


$i=S/\hat\sigma_{y}$, $\hat\sigma$ with n−1


- *Reference:* Falconer & Mackay 1996 (page unverified)

- *Function:* `select_ind`, `.select_culling`

- *Location:* `R/select_ind.R:461–478, 1089–1090`

- *Notes:* 0 when sd = 0


**V2-selection-5. Response (documented, tested here)**


$R=i\,h^2\sigma_P=h^2 S$


- *Variables:* h² narrow-sense

- *Reference:* Falconer & Mackay 1996 (page unverified)

- *Function:* roxygen only

- *Location:* `R/select_ind.R:16–21`

- *Notes:* verified numerically (exp4)


**V2-selection-6. Average effect**


$\alpha_j=a_j+d_j(q_j-p_j)=a_j+d_j(1-2p_j)$


- *Variables:* a, d per-locus (realized scale), p counted-allele freq

- *Reference:* Fisher 1918 52:399–433; Falconer & Mackay 1996 (page unverified)

- *Function:* `.avg_effect`

- *Location:* `R/grammar_realize.R:182–184`

- *Notes:* DECISION-019


**V2-selection-7. Breeding value**


$A_i=\sum_j \alpha_j (x_{ij}-2p_j)$


- *Variables:* x gene content 0/1/2

- *Reference:* as above

- *Function:* `.breeding_value_matrix`

- *Location:* `R/grammar_realize.R:236–305` (loop 282–303)

- *Notes:* refused under epistasis 250–263


**V2-selection-8. Smith–Hazel index**


$b=P^{-1}Ga,\ I=b'y$


- *Variables:* P = sample Cov(y), G = sample Cov(A), a economic

- *Reference:* Smith 1936 7:240–250; Hazel 1943 28:476–490

- *Function:* `.index_score`, `.index_weights`

- *Location:* `R/select_ind.R:769–782; 799–842`

- *Notes:* SVD pseudo-inverse if P singular; G $\neq$ Cov(y,A) off HWE (F4)


**V2-selection-9. QGSI**


$\hat I=w'\hat\gamma+\hat\gamma'W\hat\gamma$, $W\leftarrow (W+W')/2$


- *Variables:* γ̂ = true BV per trait

- *Reference:* Ceron-Rojas et al. 2026 (art. no. unverified)

- *Function:* `.quadratic_index_score`

- *Location:* `R/select_ind.R:737, 744–747`

- *Notes:* true BV, not GEBV


**V2-selection-10. Lush combined index**


$b=V^{-1}c$; $b_1=\dfrac{h^2(1-r)}{1-rh^2}$, $b_2=\dfrac{h^2 n r(1-h^2)}{[1+(n-1)rh^2](1-rh^2)}$, $I=b_1(y-\bar y)+b_2(\bar y_f-\bar y)$


- *Variables:* $V=v_P\begin{pmatrix}1&m\\m&m\end{pmatrix}$, $m=\frac{1+(n-1)t}{n}$, $c=h^2v_P(1,\frac{1+(n-1)r}{n})'$, $t=rh^2$

- *Reference:* Lush 1947 81:241–261, 362–379 (verified); Hazel 1943 (deviations)

- *Function:* `.combined_score`

- *Location:* `R/select_ind.R:867, 872–877, 889–902`

- *Notes:* family of one / t = 1 $\rightarrow$ h²·dev (890–892)


**V2-selection-11. Within-family allocation**


largest-remainder of $k\,n_f/N$ capped at $n_f$


- *Reference:* — (implementation)

- *Function:* `.sel_within_family`

- *Location:* `R/select_ind.R:939–949`

- *Notes:* tie by label order


**V2-selection-12. Among-family**


keep whole families by $\bar y_f$ until $\ge k$


- *Reference:* —

- *Function:* `.sel_among_family`

- *Location:* `R/select_ind.R:968–978`

- *Notes:* overshoots k


**V2-selection-13. Independent culling**


keep $\bigcap_t \text{top}_{\lceil c_t N\rfloor}(y_t)$; sequential: nested


- *Variables:* c_t per-trait proportions

- *Reference:* Hazel & Lush 1942 33(11):393–399 (verified)

- *Function:* `.select_culling`

- *Location:* `R/select_ind.R:1067–1081`

- *Notes:* DECISION-028


**V2-selection-14. Index : culling : tandem**


$\sqrt{T}\,i(p) : T\,i(p^{1/T}) : i(p)$


- *Variables:* T traits, p overall

- *Reference:* Hazel & Lush 1942 + package derivation

- *Function:* roxygen; test-culling.R:71–89

- *Location:* `R/select_ind.R:157–164`

- *Notes:* 1 : 0.907 : 0.707 reproduced


**V2-selection-15. Bulk pool**


$(n_1,\dots,n_N)\sim\text{Multinomial}(n, 1/N)$ seeds per current plant; each contributing plant selfed into $n_k$ seeds


- *Reference:* Bernardo 2020 (page unverified)

- *Function:* `bulk`

- *Location:* `R/select_schemes.R:156–169`; `.bulk_counts` `:470–472`

- *Notes:* multinomial since the post-audit fix of F3 (was $n_{each}=\max(2,\lceil n/N\rceil)$ then n drawn without replacement)


**V2-selection-16. Pedigree family size**


$n_{each}=\max(1,\lceil \text{pop\_size}/n_{sel}\rceil)$, trim at random


- *Reference:* Bernardo 2020 (page unverified)

- *Function:* `pedigree`

- *Location:* `R/select_schemes.R:304–307`


**V2-selection-17. Tandem schedule**


$t_g=\text{trait}[((g-1)\bmod L)+1]$


- *Reference:* DECISION-028

- *Function:* `pedigree`, `recurrent_selection`

- *Location:* `R/select_schemes.R:291, 398`


**V2-selection-18. Intermating**


each cross: pair ~ `sample.int(np, 2)` (no self), with replacement across crosses


- *Reference:* —

- *Function:* `.intermate`

- *Location:* `R/select_schemes.R:542–545`


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V2-selection-C1. General index response (context, not implemented)**


$\Delta H = i\,\operatorname{Cov}(H, I)/\sigma_I$


- *Variables:* $H = a'A$ aggregate breeding value; $I$ selection criterion; $i$ intensity

- *Reference:* Smith 1936 *Ann. Eugen.* 7:240–250 (Codex asserts pp. 241–242 Eq. 4; not re-verified here $\rightarrow$ page unverified); Hazel 1943 *Genetics* 28:476–490 (page unverified)

- *Function:* roxygen only

- *Location:* `R/select_ind.R:16-21`

- *Notes:* The general form; reduces to $i h^2 \sigma_P$ (and to $G = \operatorname{Cov}(A)$ in Smith–Hazel) only when $\operatorname{Cov}(A, P) = V_A$, i.e. Cov(A, D) = 0 under HWE. Source of F4/O4.


**V2-selection-C2. Singular Smith–Hazel extension**


$b = P^{+} G a$, $P^{+}$ = Moore–Penrose via SVD, rank tol $= \max(\dim P)\,\varepsilon\, d_1$


- *Variables:* $P^{+}$ pseudo-inverse; $d_1$ largest singular value

- *Reference:* Package numerical extension; no genetic primary source

- *Function:* `.index_weights`

- *Location:* `R/select_ind.R:799-830`

- *Notes:* Minimum-norm solution; warns on rank deficiency; exact inverse when full rank.


**V2-selection-C3. Selfed / DH family BV correlation**


$r_{S1} = \dfrac{2(1+F)}{3+F}$, $r_{DH} = \dfrac{1+F}{2}$


- *Variables:* $F$ parental inbreeding; $r = A_{ij}/\sqrt{A_{ii}A_{jj}}$

- *Reference:* Package derivation from the tabular $A$ (reconciler check: S1 sibs $A_{ij} = A_{PP} = 1+F$, $A_{ii} = 1 + (1+F)/2$; DH $A_{ij} = 1+F$, $A_{ii} = 2$ — algebra consistent); no primary page $\rightarrow$ UNVERIFIABLE attribution

- *Function:* `select_ind()` docs (`family_relationship`)

- *Location:* `R/select_ind.R:94-98`

- *Notes:* Reduces to 2/3 and 1/2 at $F = 0$.


**V2-selection-C4. Sequential culling**


$S_0 = \{1..N\};\ S_t = \operatorname{Top}_{\max(1,\operatorname{round}(c_t \vert{}S_{t-1}\vert{}))}(z_t; S_{t-1})$


- *Variables:* $c_t$ kept fraction per trait; $S_t$ survivors

- *Reference:* Hazel & Lush 1942 *J. Hered.* 33:393–399 (page unverified)

- *Function:* `.select_culling`

- *Location:* `R/select_ind.R:1067-1074`

- *Notes:* Trait order matters; `kept_at` records $\vert{}S_t\vert{}$. Fable folded this into its culling row.


**V2-selection-C5. Population pooling**


$C = [C_1\ C_2 \cdots],\ T = [T_1\ T_2 \cdots]$ subject to one shared map; ids `make.unique(sep = "_")`; pedigrees unioned by key


- *Variables:* $C, T$ cis/trans haplotype matrices

- *Reference:* Package orchestration; no primary source

- *Function:* `c.Population`

- *Location:* `R/select_schemes.R:25-68`

- *Notes:* Map equality via `identical(snp, chr)` + `all.equal(pos, cm)` (O7).


**V2-selection-C6. Single seed descent**


$P_{g+1} = \{\operatorname{self}(i, 1) : i \in P_g\}$


- *Variables:* one selfed seed per line per generation

- *Reference:* Bernardo 2020 (secondary; page unverified)

- *Function:* `single_seed_descent`

- *Location:* `R/select_schemes.R:102-118`

- *Notes:* Line count preserved; no selection.


**V2-selection-C7. Recurrent-selection cycle**


$P_{c+1} = \bigcup_{j=1}^{n_x} \operatorname{cross}(u_j, v_j, n_p)$, $(u_j, v_j) \sim$ `sample.int(np, 2)` (no selfing)


- *Variables:* $n_x$ crosses; $n_p$ progeny per cross

- *Reference:* Bernardo 2020; Falconer & Mackay 1996 (secondary; page unverified)

- *Function:* `recurrent_selection`, `.intermate`

- *Location:* `R/select_schemes.R:387-410; 354-370`

- *Notes:* Fable's "Intermating" row covers only the pair draw (361-364); this row adds the cycle.


**V2-selection-C8. Scheme RNG threading**


one `set.seed(seed)` per wrapper, then every `selfcross()`/`cross()`/`sample.int()` consumes the ambient stream (`seed = NULL`)


- *Variables:* `seed` scheme seed

- *Reference:* DECISION-006/015 (project policy; not an equation)

- *Function:* all four wrappers, `.self_each`, `.intermate`

- *Location:* `R/select_schemes.R:111, 155, 283, 384; 521; 544`

- *Notes:* Valid seeds reproduce (both auditors); `.Random.seed` not restored (F13); seed not validated (O9).


## Genomic relationship, optimum contribution, cross usefulness, MAS, MABC

*Source report:* `v2-ocs-usefulness-marker.md`  
*Entries:* 19


**V2-ocs-usefulness-marker-1. Genomic relationship**


$G = \dfrac{ZZ'}{2\sum_j p_j(1-p_j)},\; Z = M - 2p$


- *Variables:* $M$ ind$\times$marker gene content 0/1/2; $p_j$ allele freq (sample or `base_freq`)

- *Reference:* VanRaden 2008 *J Dairy Sci* 91:4414–4423, method 1 (paper verified; page unverified)

- *Function:* `g_matrix`

- *Location:* `R/select_ocs.R:71–92`

- *Notes:* monomorphic $p_j\in\{0,1\}$ dropped `:83–89`; ridge blend `:93–95`


**V2-ocs-usefulness-marker-2. Genomic inbreeding**


$F_i = G_{ii} - 1$


- *Reference:* VanRaden 2008 (page unverified)

- *Function:* `g_matrix` (doc)

- *Location:* `R/select_ocs.R:16, 42–43`

- *Notes:* relative to the base used


**V2-ocs-usefulness-marker-3. OCS objective**


$\max_{c\ge0,\;1'c=1}\; c'g - \tfrac{\lambda}{2}c'Gc$


- *Variables:* $c$ contributions; $g$ merit; $G$ from `g_matrix`

- *Reference:* Meuwissen 1997 *J Anim Sci* 75:934–940 (verified; page unverified)

- *Function:* `optimum_contribution`, `.frank_wolfe`

- *Location:* `R/select_ocs.R:105–108, 601–624`

- *Notes:* `direction="low"` negates $g$ `:264`


**V2-ocs-usefulness-marker-4. Group coancestry**


$\bar f = \tfrac12 c'Gc$


- *Reference:* Meuwissen 1997 (page unverified)

- *Function:* `optimum_contribution`, `.tune_lambda`

- *Location:* `R/select_ocs.R:319, 666`


**V2-ocs-usefulness-marker-5. FW gradient / duality gap**


$\nabla = g - \lambda Gc;\; \text{gap} = \max_k \nabla_k - \nabla'c$


- *Reference:* Lacoste-Julien & Jaggi 2015 NIPS 28 (verified)

- *Function:* `.frank_wolfe`

- *Location:* `R/select_ocs.R:601–607`

- *Notes:* tol on gap `:604`


**V2-ocs-usefulness-marker-6. Away-step and line search**


$d_{fw}=e_s-c,\; d_{aw}=c-e_a,\; \gamma^* = \min\!\big(\gamma_{max}, \tfrac{\nabla'd}{\lambda d'Gd}\big),\; \gamma_{max}^{aw}=\tfrac{c_a}{1-c_a}$


- *Reference:* Lacoste-Julien & Jaggi 2015 (verified)

- *Function:* `.frank_wolfe`

- *Location:* `R/select_ocs.R:609–625`

- *Notes:* exact maximiser of concave quadratic along $d$


**V2-ocs-usefulness-marker-7. Penalty tuning**


bisection on $\lambda$ for $\tfrac12 c^*(\lambda)'Gc^*(\lambda) = \bar f_{target}$


- *Reference:* (monotone in $\lambda$: exchange argument)

- *Function:* `.tune_lambda`

- *Location:* `R/select_ocs.R:652–718`

- *Notes:* doubling to bracket `:697–701`


**V2-ocs-usefulness-marker-8. Parent sampling**


default `method = "allocate"`: $n_i=\lfloor nc_i\rfloor$ plus one slot each for the $n-\sum_i\lfloor nc_i\rfloor$ largest remainders (slot order shuffled); `method = "multinomial"`: $n_i \sim \text{Multinomial}(n, c)$


- *Reference:* package definition (documented)

- *Function:* `sample_parents`

- *Location:* `R/select_ocs.R:429–441`; `.allocate_slots` `:473–494`

- *Notes:* allocation is the default since the post-audit fix of OCS-F2 (owner decision D5, `docs/DECISIONS.md`); multinomial only on request


**V2-ocs-usefulness-marker-9. Usefulness**


$U = \mu + i\,\sigma$ (or $\mu - i\sigma$ for `"low"`)


- *Variables:* $\mu,\sigma$ realized family mean/sd of additive value

- *Reference:* Zhong & Jannink 2007 *Genetics* 177:567–576; Lehermeier et al. 2017 *Genetics* 207:1651–1661 (verified; pages unverified)

- *Function:* `cross_usefulness`

- *Location:* `R/select_usefulness.R:142–144`

- *Notes:* family from crossing engine `:140`, `:297–315`


**V2-ocs-usefulness-marker-10. Selection intensity**


$i(p) = \dfrac{\varphi(\Phi^{-1}(1-p))}{p}$


- *Variables:* $p$ = `select_top`

- *Reference:* standard normal-theory (Falconer & Mackay 1996, author-year)

- *Function:* `.intensity_from_p`

- *Location:* `R/select_usefulness.R:161–163`

- *Notes:* infinite-population value (1.755 at p=0.1)


**V2-ocs-usefulness-marker-11. Fixed additive score**


$g_i = \sum_j x_{ij}\,e_j,\; e_j = e^{raw}_j\sqrt{\pi_t}/\text{sd}(\text{comp})$


- *Variables:* $x\in\{-1,0,1\}$; $\pi_t$ layer `prop`; sd over template

- *Reference:* package definition (mirrors `.genetic_matrix`)

- *Function:* `.additive_model`, `.additive_gv`

- *Location:* `R/select_usefulness.R:204–209, 228`

- *Notes:* orthogonal layer uses $\alpha$ `:194–198`


**V2-ocs-usefulness-marker-12. Average effect**


$\alpha = a + d(1-2p) = a + d(q-p)$


- *Reference:* Falconer & Mackay 1996 (author-year)

- *Function:* `.avg_effect`

- *Location:* `R/grammar_realize.R:182–184`

- *Notes:* used by `.additive_model:197`, `marker_select` docs


**V2-ocs-usefulness-marker-13. MAS feasibility**


$z_{mi} = x_{mi}\cdot f_m;\; \text{met}_{mi} = [z_{mi} \ge \tau_m],\; \tau=0$ carrier, $1$ homozygote; feasible $\iff \sum_m \text{met}_{mi} \ge k_{min}$


- *Variables:* $f_m\in\{-1,+1\}$ favourable homozygote

- *Reference:* DECISION-029 (package definition)

- *Function:* `marker_select`

- *Location:* `R/select_marker.R:116–119`

- *Notes:* ranking `:167`


**V2-ocs-usefulness-marker-14. MARS/oracle index**


$I_i = \sum_j w_j x_{ij}$


- *Variables:* $w_j = a_j$ or $\alpha_j$

- *Reference:* Lande & Thompson 1990 *Genetics* 124:743–756 (verified; page unverified)

- *Function:* `additive_value` (documented in `marker_select`)

- *Location:* `R/cross_population.R:735–747`; doc `R/select_marker.R:20–33`


**V2-ocs-usefulness-marker-15. Recurrent-allele score**


$s_{mi} = \dfrac{x_{mi} r_m + 1}{2}$


- *Variables:* $r_m = \pm1$ recurrent homozygote dosage

- *Reference:* package definition

- *Function:* `.mabc_founders`

- *Location:* `R/select_mabc.R:311`

- *Notes:* informative iff $\vert{}r_m\vert{}=1,\; d_m=-r_m$ `:309`


**V2-ocs-usefulness-marker-16. Recovery**


$R_i = \dfrac{\sum_{m\in B} w_m s_{mi}}{\sum_{m\in B} w_m}$


- *Variables:* $B$ background set; $w$ equal / interval / user

- *Reference:* package definition (marker-observed, not IBD)

- *Function:* `.mabc_recovery`

- *Location:* `R/select_mabc.R:464`


**V2-ocs-usefulness-marker-17. Expected recovery**


$R_t = 1 - 2^{-(t+1)}$, $R_t = \tfrac12 + \tfrac12 R_{t-1}$


- *Variables:* $t$ backcrosses after F1

- *Reference:* Mendelian; attributed to Frisch & Melchinger 2005 *Genetics* 170:909–917 (paper verified; location unverified)

- *Function:* doc only

- *Location:* `R/select_mabc.R:25–31, 227–229`

- *Notes:* MC: 0.756/0.870/0.934


**V2-ocs-usefulness-marker-18. Interval weights**


$w_m = \big[\tfrac{p_m+p_{m+1}}{2} - \tfrac{p_{m-1}+p_m}{2}\big] - \text{overlap with excluded intervals}$, end cells to chromosome ends


- *Variables:* $p$ cM positions

- *Reference:* package definition

- *Function:* `.mabc_weights`

- *Location:* `R/select_mabc.R:520–529`

- *Notes:* per-chromosome sum = mapped − excluded length


**V2-ocs-usefulness-marker-19. Foreground / recombinant / background staging**


lexicographic: feasible $\rightarrow$ #flanks recurrent-homozygous $\rightarrow$ $R_i$ $\rightarrow$ seeded tie


- *Reference:* Frisch & Melchinger 2001 *Crop Sci* 41:1485–1494 (verified; page unverified)

- *Function:* `mabc_select`

- *Location:* `R/select_mabc.R:170–177, 194`


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V2-ocs-usefulness-marker-C1. Ridge blend of G**


$G^{*} = (1-r)\,G + r\,I$


- *Variables:* $r$ = `ridge` $\in[0,1]$

- *Reference:* package extension (no primary source claimed)

- *Function:* `g_matrix`

- *Location:* `R/select_ocs.R:25-26, 44-50, 93-95`

- *Notes:* Identity regularisation for a PD matrix; not VanRaden's pedigree-$A$ blend. Fable mentioned it only in a Notes cell.


**V2-ocs-usefulness-marker-C2. MABC foreground feasibility**


carrier: $s_{ji} < 1\ \forall j\in T$; donor homozygote: $s_{ji} = 0\ \forall j\in T$, with $s_{ji} = (x_{ji} r_j + 1)/2$


- *Variables:* $T$ target markers; $r_j=\pm1$ recurrent homozygote code

- *Reference:* Frisch & Melchinger 2001 *Crop Sci* 41:1485-1494 (foreground stage; p. 1485 per Codex, page unverified by reconciler); rule itself is package/DECISION-029-style definition

- *Function:* `mabc_select`

- *Location:* `R/select_mabc.R:168-174`

- *Notes:* Under the documented shared-allele-coding precondition. Fable covered this in finding MABC-P1 but not as a map row.


## Prediction: pedigree A matrix, BLUP/GBLUP, combining ability, progeny testing

*Source report:* `v2-prediction.md`  
*Entries:* 19


**V2-prediction-1. Tabular relationship (off-diagonal)**


$A_{kj} = \tfrac12\,(A_{m(k),j} + A_{f(k),j})$, $j < k$


- *Variables:* m(k), f(k) parents of k (processed parents-first)

- *Reference:* Emik & Terrill 1949 *J. Hered.* 40:51–55 (page unverified); Wright 1922 *Am. Nat.* 56:330–338 (page unverified)

- *Function:* `a_matrix`

- *Location:* `R/select_blup.R:73-79`

- *Notes:* unknown parent contributes 0; DH: m = f so row = A_{P,·}


**V2-prediction-2. Inbreeding / diagonal**


$A_{kk} = 1 + F_k,\; F_k = \tfrac12 A_{m(k)f(k)}$; founder $A_{kk}=1$; DH $A_{kk}=2$


- *Variables:* —

- *Reference:* as above

- *Function:* `a_matrix`

- *Location:* `R/select_blup.R:80-88`

- *Notes:* selfing: $F_t = 1-(1/2)^t$ from a non-inbred founder (verified t = 1,2,3); founders forced F = 0 (PRED-F1)


**V2-prediction-3. Variance ratio**


$\lambda = \sigma^2_e/\sigma^2_A$; with `h2`: $\sigma^2_A = h^2\sigma^2_P,\ \sigma^2_e = (1-h^2)\sigma^2_P$, $\sigma^2_P = \mathrm{var}(\text{ref or pheno})$


- *Variables:* —

- *Reference:* Henderson 1975 *Biometrics* 31:423–447 (page unverified)

- *Function:* `predict_ebv`, `.blup_variances`

- *Location:* `R/select_blup.R:365`, `:676-681`

- *Notes:* not inverted (verified)


**V2-prediction-4. GLS mean**


$\hat\mu = (\mathbf 1'V^{-1}\mathbf 1)^{-1}\mathbf 1'V^{-1}y$, $V = K_{rr} + \lambda I$


- *Variables:* r = phenotyped rows; V on the σ²_A scale

- *Reference:* Henderson 1975 (page unverified)

- *Function:* `predict_ebv`

- *Location:* `R/select_blup.R:376-386`

- *Notes:* Cholesky of V is the PD check


**V2-prediction-5. BLUP of u**


$\hat u = K_{\cdot r}\,V^{-1}(y - \mathbf 1\hat\mu)$ $\equiv$ MME $\begin{bmatrix}X'X & X'Z\\ Z'X & Z'Z+\lambda K^{-1}\end{bmatrix}\begin{bmatrix}\hat\mu\\ \hat u\end{bmatrix}=\begin{bmatrix}X'y\\ Z'y\end{bmatrix}$


- *Variables:* Z selects one record per individual

- *Reference:* Henderson 1975 (page unverified)

- *Function:* `predict_ebv`

- *Location:* `R/select_blup.R:387-388`

- *Notes:* no K⁻¹ needed; equals MME to 1e-10 (executed) and RR-BLUP for singular G


**V2-prediction-6. PEV / reliability**


$\mathrm{PEV}_i/\sigma^2_A = K_{ii} - K_{i r}\,P\,K_{r i}$, $P = V^{-1} - V^{-1}\mathbf 1(\mathbf 1'V^{-1}\mathbf 1)^{-1}\mathbf 1'V^{-1}$; $\ r^2_i = 1 - \mathrm{PEV}_i/(K_{ii}\sigma^2_A)$


- *Variables:* K_ii = 1 + F_i

- *Reference:* Henderson 1975 (page unverified)

- *Function:* `predict_ebv`

- *Location:* `R/select_blup.R:395-398`

- *Notes:* K = I: r² = h²(1 − 1/n) (test line 24); equals 1 − C^{uu}_ii λ/K_ii


**V2-prediction-7. GBLUP marker back-solve**


$\hat u_m = Z_r'\alpha \,/\, [2\sum_j p_j(1-p_j)]$, $\alpha = V^{-1}(y-\mathbf 1\hat\mu)$; $Z\hat u_m = \widehat{\mathrm{GEBV}}$


- *Variables:* Z = M − 2p (VanRaden method 1)

- *Reference:* VanRaden 2008 *J. Dairy Sci.* 91:4414–4423 (page unverified)

- *Function:* `.gblup_marker_effects`

- *Location:* `R/select_blup.R:695-703`

- *Notes:* equals RR-BLUP (executed)


**V2-prediction-8. Accuracy / bias**


$\mathrm{acc} = \mathrm{cor}(\hat u, u)$; slope $= \mathrm{Cov}(u,\hat u)/\mathrm{Var}(\hat u)$


- *Variables:* —

- *Reference:* package (BLUP property E[u | û] = û)

- *Function:* `prediction_accuracy`

- *Location:* `R/select_blup.R:752-761`

- *Notes:* —


**V2-prediction-9. Expected cross mean (one locus)**


$E[G_{ik}] = g_ig_k a + [g_i(1-g_k)+(1-g_i)g_k]\,d - (1-g_i)(1-g_k)\,a = -2d\,g_ig_k + (g_i+g_k)(a+d) - a$, $g = x/2$


- *Variables:* x $\in$ {0,1,2} gene content; genotypic values −a, d, +a

- *Reference:* package derivation from Falconer & Mackay 1996 genotypic values (page unverified)

- *Function:* `.expected_cross_means`

- *Location:* `R/select_combining.R:244-255`

- *Notes:* summed over loci; linearity $\Rightarrow$ linkage-free


**V2-prediction-10. Tester-referenced average effect**


$\alpha_T = a + d(1 - 2p_T) = a + d(q_T - p_T)$; $\mathrm{GCA}_i \propto \tfrac12\sum_j \alpha_{T,j}(x_{ij} - 2p_{C,j})$


- *Variables:* p_T tester gamete freq

- *Reference:* Falconer & Mackay 1996 (page unverified)

- *Function:* `combining_ability` (doc)

- *Location:* `R/select_combining.R:30-37`

- *Notes:* = 1/2 DECISION-019 BV when p_T = p (executed exact)


**V2-prediction-11. Topcross / factorial GCA, SCA**


$g_i = \bar Y_{i\cdot} - \mu,\; h_k = \bar Y_{\cdot k} - \mu,\; s_{ik} = Y_{ik} - \mu - g_i - h_k$


- *Variables:* μ grand mean of Y

- *Reference:* Sprague & Tatum 1942 *J. Am. Soc. Agron.* 34:923–932 (page unverified)

- *Function:* `.decompose_ca`

- *Location:* `R/select_combining.R:321-324`

- *Notes:* balanced two-way LS


**V2-prediction-12. Diallel GCA (method 4)**


$g_i = (m_i - \mu)\dfrac{p-1}{p-2} = \dfrac{pY_{i\cdot} - 2Y_{\cdot\cdot}}{p(p-2)}$


- *Variables:* m_i mean of i's p−1 crosses; μ mean of p(p−1)/2 crosses

- *Reference:* Griffing 1956 *Aust. J. Biol. Sci.* 9:463–493 (page unverified)

- *Function:* `.decompose_ca`

- *Location:* `R/select_combining.R:312-316`

- *Notes:* algebraically identical to Griffing method 4 (executed)


**V2-prediction-13. Diallel SCA (method 4)**


$s_{ij} = Y_{ij} - \mu - g_i - g_j = Y_{ij} - \dfrac{Y_{i\cdot}+Y_{j\cdot}}{p-2} + \dfrac{2Y_{\cdot\cdot}}{(p-1)(p-2)}$


- *Variables:* —

- *Reference:* Griffing 1956 (page unverified)

- *Function:* `.decompose_ca`

- *Location:* `R/select_combining.R:317-318`

- *Notes:* Σ_j s_ij = 0 (executed)


**V2-prediction-14. Broad-sense residual for simulated crosses**


$\sigma^2_e = \mathrm{Var}(G_{\text{ref}})\,(1-h^2)/h^2$


- *Variables:* G = A + D of progeny (default ref)

- *Reference:* Falconer & Mackay 1996 H² (page unverified)

- *Function:* `phenotype_value` via `.simulate_cross_means`

- *Location:* `R/cross_population.R:1061`; `R/select_combining.R:290-293`

- *Notes:* —


**V2-prediction-15. Realized-scale template effects**


$a_j^{\ast} = a_j\sqrt{\pi_\ell}/\mathrm{sd}(\text{raw component}_\ell)$ (same for d)


- *Variables:* π_ℓ layer prop

- *Reference:* package (SPEC §2)

- *Function:* `.layer_scaled_effects` $\leftarrow$ `template_effects`

- *Location:* `R/grammar_realize.R:364`; `R/select_combining.R:400-406`

- *Notes:* reproduces `genetic_values()` up to a constant (executed)


**V2-prediction-16. Progeny-test expected mean**


$E[\bar y_i] = \tfrac12\sum_j \alpha_{M,j}(x_{ij} - 2p_{M,j}) + c$, $\alpha_M = a + d(q_M - p_M)$


- *Variables:* p_M mates' frequency

- *Reference:* package derivation; Falconer & Mackay 1996 α (page unverified)

- *Function:* `progeny_test` (doc)

- *Location:* `R/select_progeny.R:17-35`

- *Notes:* slope 0.503 on 200 progeny (executed)


**V2-prediction-17. Progeny-test accuracy**


$r = \dfrac{V_A/2}{\sqrt{V_A\,V_P[1+(n-1)h^2/4]/n}} = \sqrt{\dfrac{n h^2}{4 + (n-1)h^2}}$; crossover $n > \dfrac{4-h^2}{1-h^2}$


- *Variables:* n half-sib records

- *Reference:* package derivation (classical; no page asserted)

- *Function:* `progeny_test` (doc)

- *Location:* `R/select_progeny.R:37-53`

- *Notes:* validated by simulation in `test-progeny.R:52-70`


**V2-prediction-18. Mate sampling**


mates for parent j: `sample.int` without replacement from distinct mates $\neq$ j


- *Variables:* —

- *Reference:* package

- *Function:* `progeny_test`

- *Location:* `R/select_progeny.R:121-135`

- *Notes:* half-sib by construction


**V2-prediction-19. Family mean**


$\bar y_i = \mathrm{mean}\{y_k : \text{mother}(k) = i\}$


- *Variables:* —

- *Reference:* package

- *Function:* `progeny_test`

- *Location:* `R/select_progeny.R:149-150`

- *Notes:* —


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V2-prediction-C1. Known-variance mixed model (the model BLUP solves)**


$y = \mathbf{1}\mu + Zu + e,\quad u \sim N(0, K\sigma^2_A),\quad e \sim N(0, I\sigma^2_e)$


- *Variables:* $K$: $A$ or $G$; $Z$: one-record-per-individual incidence; $\mu$: single fixed mean

- *Reference:* Henderson 1975 *Biometrics* 31:423–447 (page unverified)

- *Function:* `predict_ebv`

- *Location:* `R/select_blup.R:174-181` (roxygen), `:300-321` (inputs)

- *Notes:* Observable phenotypes only; intercept-only fixed effects; variance components known, no REML (DECISION-030). Fable's map starts at the GLS solution and never states the model.


**V2-prediction-C2. Supplied-covariance validation (correlation-scale PSD test)**


$C_{ij} = K_{ij}/\sqrt{K_{ii}K_{jj}}$ over rows with $K_{ii} > 0$; accept iff $\lambda_{\min}(C) \ge -\,n\,\epsilon\,\max(1, \lvert\lambda(C)\rvert)$ and every zero-variance row is all-zero


- *Variables:* $\epsilon$: `.Machine$double.eps`; $n$: number of eigenvalues

- *Reference:* Numerical-linear-algebra tolerance; no genetics source

- *Function:* `.check_relationship`

- *Location:* `R/select_blup.R:587-650` (scaling `:630-632`, eigen `:640`, tolerance `:644-645`)

- *Notes:* PRED-01: when $C$ is non-finite the code sets `ev <- -Inf`, which makes the tolerance $\infty$ and the test vacuous (CONFIRMED). Also `(K + t(K))/2` at `:265` can overflow for entries near `.Machine$double.xmax`.


**V2-prediction-C3. VanRaden method-1 genomic relationship (consumed by GBLUP)**


$G = \dfrac{ZZ'}{2\sum_j p_j(1-p_j)},\quad Z = M - 2p$, markers with $p_j \in \{0,1\}$ dropped; optional $G^* = (1-\text{ridge})\,G + \text{ridge}\,I$


- *Variables:* $M$: individuals $\times$ markers gene content 0/1/2; $p$: `base_freq` or current frequencies

- *Reference:* VanRaden 2008 *J. Dairy Sci.* 91:4414–4423 (page unverified)

- *Function:* `predict_ebv` $\rightarrow$ `g_matrix`

- *Location:* `R/select_blup.R:337-343` (Codex wrote 168-173); `R/select_ocs.R:44-98` (Codex wrote 44-92; core at `:83, 90-92`)

- *Notes:* Fable covered this only in the O1 checklist, not in its map. Inbred founders give $G_{ii} \approx 2$ while `a_matrix()` gives 1 (PRED-F1).


**V2-prediction-C4. Progeny-test RNG draw order (reproducibility contract)**


`set.seed(seed)` $\rightarrow$ for each parent $j$: `sample.int` of $n$ distinct mates $\ne j$ $\rightarrow$ `mate()` meioses in plan order $\rightarrow$ residual draws


- *Variables:* one R stream; `seed = NULL` uses the ambient stream

- *Reference:* Package design, DECISION-012/027

- *Function:* `progeny_test`

- *Location:* `R/select_progeny.R:113-114` (seed), `:132` (mates), `:138` (meioses), `:144-145` (residual)

- *Notes:* Fable's "Mate sampling" row covers only the first step. Ambient stream is not restored afterwards (PRED-F4).


## Transcriptome: expression model, eQTL, mediation, NB counts, GREML mimic

*Source report:* `v2-transcriptome.md`  
*Entries:* 25


**V2-transcriptome-1. Expression model**


$E_{gi} = \mu_g + G_{gi} + R_{gi}$


- *Variables:* mu: gene location (0 unless mimic); G genetic; R residual

- *Reference:* Falconer & Mackay 1996 (additive model, definition of h2; page unverified)

- *Function:* `simulate_transcriptome`

- *Location:* `R/transcriptome_simulate.R:616-618`

- *Notes:* G, R drawn independently; mimic affine at :602-614


**V2-transcriptome-2. Reference-centered dosage**


$Z_{ij} = x_{ij} - \bar x_j$


- *Variables:* x in {-1,0,1}; marker mean on the reference panel

- *Reference:* own design (DECISION-020/021 fixed reference)

- *Function:* `simulate_transcriptome`

- *Location:* `:308-310`; predict `:845-847`

- *Notes:* mean, not 2p; equals x - 2p + 1 shift


**V2-transcriptome-3. cis score**


$c_g = \sum_{j \in W(g)} \beta_{gj} Z_j$, $\beta \sim N(0,1)$, $k = 1 + \mathrm{Bin}(2, 0.25)$


- *Variables:* W(g): same chr, $\lvert pos - tss\rvert$ <= cis_window, MAF >= 0.05

- *Reference:* GTEx 2020 Science 369:1318-1330 (multiple cis-eQTL; verified)

- *Function:* `simulate_transcriptome`

- *Location:* `:501-511`

- *Notes:* window in physical bp


**V2-transcriptome-4. trans factor**


$f_q = \sum_{k \in H(q)} \gamma_{qk} Z_k$, $\vert{}H(q)\vert{} = 1 + \mathrm{Bern}(0.5)$, loading 1


- *Variables:* hubs outside all module genes' cis windows

- *Reference:* Albert et al. 2018 eLife 7:e35471 (trans hotspots; verified)

- *Function:* `simulate_transcriptome`

- *Location:* `:451-471`, `:512`

- *Notes:* factor inert if no distant marker


**V2-transcriptome-5. Standardization**


$\tilde v = (v - \bar v)/s_v$, $s_v$ with $n-1$


- *Variables:* any component

- *Reference:* --

- *Function:* `z1`

- *Location:* `:492-496`

- *Notes:* returns NULL if s < 1e-9


**V2-transcriptome-6. Additive genetic score**


$G^0 = \sqrt{\omega}\,\tilde c + \sqrt{1-\omega}\,\tilde t$; $G = \sqrt{h^2}\, G^0 / s_{G^0}$


- *Variables:* omega cis fraction (marginal)

- *Reference:* own design; SPEC section 4

- *Function:* `simulate_transcriptome`

- *Location:* `:548-578`

- *Notes:* joint scaling; cancellation fallback :555-558


**V2-transcriptome-7. Epistatic blend**


$\tilde e = \widetilde{\sum_p \beta_p (Z_{j_p} Z_{k_p} - \overline{Z_j Z_k})}$; $G^0 = \sqrt{1-\epsilon}\,\tilde{G}^{0}_{ct} + \sqrt{\epsilon}\,\tilde e$


- *Variables:* epsilon epistatic fraction, 1-2 pairs

- *Reference:* own design (a x a centered product)

- *Function:* `simulate_transcriptome`

- *Location:* `:522-572`

- *Notes:* prod_mean stored for predict


**V2-transcriptome-8. Residual**


$R^0 = \sqrt{\kappa}\,\tilde u_{q(g)} + \sqrt{1-\kappa}\,\tilde\varepsilon_g$; $R = \sqrt{1-h^2}\, R^0/s_{R^0}$


- *Variables:* u_q ~ N(0,1) shared per module; kappa residual module fraction

- *Reference:* own design

- *Function:* `simulate_transcriptome`

- *Location:* `:585-595`

- *Notes:* co-expression independent of h2


**V2-transcriptome-9. Realized heritability**


$h^2_{real} = \mathrm{Var}(G)/\mathrm{Var}(E) = h^2/(1 + 2\mathrm{Cov}(G,R))$


- *Variables:* gr_cov = 2 Cov(G,R)

- *Reference:* Falconer & Mackay 1996 (definition); V2

- *Function:* `simulate_transcriptome`, `predict`

- *Location:* `:619-628`; `:943-948`

- *Notes:* unbounded above (TX-F1); the bounded allocation `h2_allocated` = Var(G)/(Var(G)+Var(R)) is reported separately (`:627`, `:949`)


**V2-transcriptome-10. Effective coefficients**


$s_c = \sqrt{h^2(1-\epsilon)\omega}/(s_{G^0_{ct}} s_{G^0})$, $s_t$, $s_e = \sqrt{h^2\epsilon}/s_{G^0}$; cis effect $= s_c\beta_j/s_c^{raw}$; trans_scale $= s_t/s_t^{raw}$


- *Variables:* unit-scale truth

- *Reference:* own design

- *Function:* `simulate_transcriptome`

- *Location:* `:636-646`, `:666`, `:673`, `:684`

- *Notes:* reconstruct G exactly (verified 1e-8)


**V2-transcriptome-11. Genetic budget**


$\mathrm{Var}(G) = v_{cis} + v_{trans} + v_{epi} + 2\mathrm{Cov}(c,t) + 2\mathrm{Cov}(c,e) + 2\mathrm{Cov}(t,e)$


- *Variables:* --

- *Reference:* cf. DECISION-020 `add_dom_cov` pattern

- *Function:* `simulate_transcriptome`

- *Location:* `:651-658`

- *Notes:* closes to 1e-9


**V2-transcriptome-12. Realized cis fraction**


$\omega_{real} = v_{cis}/\mathrm{Var}(G)$, $\mathrm{Var}(G)$ the realized genetic variance (all covariance terms of V2-transcriptome-11 included)


- *Variables:* --

- *Reference:* SPEC section 3 ("marginal")

- *Function:* `simulate_transcriptome`

- *Location:* `:659-661`; predict `:950`

- *Notes:* reported as `cis_fraction_realized`; no longer identically omega on the reference (TX-F2 addressed): it equals omega only when the epistatic and covariance shares are zero


**V2-transcriptome-13. Cross-population genetic value**


$G^{new}_g = Z^{new}_c \beta_g + \text{trans\_scale}_g Z^{new}_H \gamma + \sum_p e_p(Z_jZ_k - \overline{Z_jZ_k})$; $Z^{new} = x^{new} - \bar x^{ref}$


- *Variables:* fixed reference means

- *Reference:* own design (fixed-scale principle)

- *Function:* `predict.transcriptome_sim`

- *Location:* `:829-848`, `:860-886`

- *Notes:* scaled by stored gene_scale :921-927


**V2-transcriptome-14. Mimic affine**


$E = \mu_g + \frac{\sqrt{V_g}}{s_{G+R}}\,[(G-\bar G) + (R - \bar R)]$


- *Variables:* mu_g, V_g from user matrix

- *Reference:* own design

- *Function:* `simulate_transcriptome`

- *Location:* `:602-616`

- *Notes:* exact moments (verified 1e-10)


**V2-transcriptome-15. GRM**


$K = \frac{ZZ'}{m} \Big/ \overline{\mathrm{diag}}$


- *Variables:* Z reference-centered dosages

- *Reference:* VanRaden 2008 J. Dairy Sci. 91:4414-4423 method 1 up to a constant (page listed in THEORY_REVIEW as "verify"; not re-verified here)

- *Function:* `.tx_grm`

- *Location:* `R/transcriptome_mimic.R:24-30`

- *Notes:* ratio to VanRaden G = 0.502 (constant)


**V2-transcriptome-16. REML objective**


$-\ell_R(\delta) = \tfrac12\big[(n-1)\log \mathrm{RSS}(\delta) + \sum_i \log(\xi_i + \delta) + \log \sum_i \omega_i^2/(\xi_i+\delta)\big]$, $\mathrm{RSS} = \sum w_i\eta_i^2 - (\sum w_i\omega_i\eta_i)^2/\sum w_i\omega_i^2$


- *Variables:* K = U diag(xi) U'; eta = U'y; omega = U'1; delta = sigma_e^2/sigma_g^2

- *Reference:* Kang et al. 2008 Genetics 178(3):1709-1723 (verified); Patterson & Thompson 1971 Biometrika 58:545-554 (page unverified)

- *Function:* `.greml_h2`

- *Location:* `R/transcriptome_mimic.R:76-84`

- *Notes:* matches rrBLUP 8e-6


**V2-transcriptome-17. GREML heritability**


$h^2 = \sigma_g^2/(\sigma_g^2+\sigma_e^2) = 1/(1+\hat\delta)$


- *Variables:* mean diag(K) = 1

- *Reference:* Yang et al. 2010 Nat. Genet. 42:565-569 (page unverified)

- *Function:* `.greml_h2`

- *Location:* `R/transcriptome_mimic.R:94-97`

- *Notes:* boundaries 0/1 by comparison at log delta = +/-20


**V2-transcriptome-18. Factor count**


$Q = \#\{\lambda_i(\tfrac1T E_s'E_s) > (1+\sqrt{n/T})^2\}$, capped to [1, min(50, n-2)]


- *Variables:* E_s standardized genes x ind

- *Reference:* Marchenko & Pastur 1967 Mat. Sb. 72:507-536 (page unverified)

- *Function:* `.tx_estimate_factors`

- *Location:* `R/transcriptome_mimic.R:111-119`

- *Notes:* floor 1 on noise (verified)


**V2-transcriptome-19. kappa proxy**


$\hat\kappa = \mathrm{clamp}_{[0,1]}\!\Big[\big(r_w - \overline{\sqrt{h^2(1-\omega)}}^{\,2}\big)\big/\overline{\sqrt{1-h^2}}^{\,2}\Big]$, $r_w = \sum_{i\le Q} s_i/(T-Q)$, $s_i = \ell_i - 1$ from $\lambda_i = \ell_i\big(1 + y/(\ell_i - 1)\big)$, $y = T/n$; $s_i = 0$ if $\lambda_i \le (1+\sqrt y)^2$


- *Variables:* $\lambda_i$ leading $Q$ eigenvalues of the gene-gene correlation matrix; $h^2$ per-gene GREML estimates; $\omega$ the assumed cis fraction (default 0.25)

- *Reference:* own estimator; the spike inversion is the Baik–Ben Arous–Péché spiked-model relation, named in the roxygen without a page-level source

- *Function:* `.tx_estimate_kappa`

- *Location:* `R/transcriptome_mimic.R:149-172` (roxygen `:120-148`); called at `R/transcriptome_simulate.R:395`

- *Notes:* subtracts the GREML-implied genetic trans share (TX-F3 addressed); with $h^2 = 0$ no correction is applied


**V2-transcriptome-20. NB observation**


$\log\mu_{gi} = \log L_i + \alpha_g + \sigma_g z_{gi}$; $Y_{gi} \sim \mathrm{NB}(\mu_{gi}, \phi_g)$, $\mathrm{Var} = \mu + \phi\mu^2$, `size = 1/phi`, $\phi = 0 \Rightarrow$ Poisson


- *Variables:* z per-gene standardized latent E

- *Reference:* McCarthy, Chen & Smyth 2012 NAR 40:4288-4297; Robinson & Smyth 2008 Biostatistics 9:321-332 (pages unverified); cited in the `@references` of `observe_counts` (`R/transcriptome_counts.R:42`, `man/observe_counts.Rd`)

- *Function:* `observe_counts`

- *Location:* `R/transcriptome_counts.R:75-80`, `:98-102`

- *Notes:* moments verified within 1.3%


**V2-transcriptome-21. Layer score**


$T_x = c\sum_g w_g z_g$, $w = \text{slope}/\max\vert{}\text{slope}\vert{}$, $z_g = (E_g - \bar E_g)/s_{E_g}$, $c = \sqrt{\text{prop}}/s_{\sum w z}$


- *Variables:* causal genes sparse

- *Reference:* own design (SPEC section 3; mirrors `additive()` scaling)

- *Function:* `.tx_raw`, `.transcriptome_matrix`

- *Location:* `R/grammar_realize.R:454-476`, `:414-440`

- *Notes:* Var(Tx) = prop exactly


**V2-transcriptome-22. Mediation split**


$T_{x,g} = c\sum_g w_g (G_g - \bar G_g)/s_{E_g}$, $T_{x,e} = T_x - T_{x,g}$; shares $\mathrm{Var}(T_{x,g})/V_P$, $\mathrm{Var}(T_{x,e})/V_P$, $2\mathrm{Cov}/V_P$


- *Variables:* same denominator sd(E_g) so z_total = z_gen + z_env

- *Reference:* own design; MESC estimand analogue (benchmarks README)

- *Function:* `.tx_raw(which="genetic")`, `.mediation_budget`, `mediation_split`

- *Location:* `R/grammar_realize.R:458-476`, `:730-761`; `R/io_write.R:1267-1270`

- *Notes:* sums to realized share (1e-8)


**V2-transcriptome-23. Realized H2**


$H^2 = \mathrm{Var}(\text{Gen} + T_{x,g})/V_P$


- *Variables:* Gen marker layers

- *Reference:* V2

- *Function:* `.genetic_value_matrix`, `.realized_h2`

- *Location:* `R/grammar_realize.R:394-400`, `:1078-1106`

- *Notes:* real source: Tx_g = 0


**V2-transcriptome-24. Per-gene share**


$\text{var\_explained}_g = (c\,w_g)^2 \cdot 1 / V_P$


- *Variables:* --

- *Reference:* mirrors `.qtn_var`

- *Function:* `.tx_qtn_var`

- *Location:* `R/io_write.R:1328-1345`

- *Notes:* marginal; sum 0.36 vs prop 0.4 (co-expression cov)


**V2-transcriptome-25. Layer sub-seed**


$b \leftarrow (257\,b + u_k) \bmod 2147483629$ over the UTF-8 codes $u_k$ of the layer type ($b_0 = 0$); $s = (1009\,\text{seed} + 7919\,b + 104729\,\text{occ}) \bmod (2^{31}-1)$


- *Variables:* --

- *Reference:* R2 rule

- *Function:* `.layer_seed`, `transcriptome`

- *Location:* `R/grammar_simulate_phenotype.R:1079-1092`; `R/transcriptome_layer.R:257-276`

- *Notes:* order invariance verified


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V2-transcriptome-C1. Single-component LMM (model statement)**


$y = \mathbf{1}\mu + u + e,\ u \sim N(0, \sigma_g^2 K),\ e \sim N(0, \sigma_e^2 I);\ h^2 = \sigma_g^2/(\sigma_g^2+\sigma_e^2) = 1/(1+\delta)$


- *Variables:* $K$: GRM with mean diagonal 1; $\delta = \sigma_e^2/\sigma_g^2$

- *Reference:* Kang et al. 2008 *Genetics* 178(3):1709-1723 (EMMA; page/eq unverified)

- *Function:* `.greml_h2`

- *Location:* `R/transcriptome_mimic.R:32-53` (roxygen), fit `:58-97`

- *Notes:* Fable's map has the profile objective and $h^2$ rows but not the model statement; no identifiability guard (O4)


**V2-transcriptome-C2. Calibration metrics**


$\mathrm{bias} = \overline{\hat\theta - \theta},\ \mathrm{RMSE} = \sqrt{\overline{(\hat\theta-\theta)^2}}$


- *Variables:* $\theta$ target, $\hat\theta$ realized (h2, cis fraction)

- *Reference:* descriptive statistics; no source needed

- *Function:* benchmark 01

- *Location:* `benchmarks/01_h2_calibration.R:61-62` (h2), `:93-94` (cis)

- *Notes:* Cis-panel RMSE was 4e-16 by construction under the old $v_{cis}/(v_{cis}+v_{trans})$ definition (TX-F2); the benchmark now scores `cis_fraction_realized` = $v_{cis}/\mathrm{Var}(G)$ (`:83`), not re-run here


**V2-transcriptome-C3. Marginal eQTL scan**


$r_j = \dfrac{(D_j-\bar D_j)^\top(y-\bar y)}{\lVert D_j-\bar D_j\rVert\,\lVert y-\bar y\rVert}$; rank by $\lvert r_j\rvert$; best true rank $= \min_{j\in\text{true}} \mathrm{rank}_j$


- *Variables:* $D_j$ dosage of candidate $j$; $y$ gene expression

- *Reference:* Pearson correlation (monotone in $\lvert t\rvert$ for one predictor); page unverified

- *Function:* `scan_one`

- *Location:* `benchmarks/02_eqtl_recovery.R:61-80` ($r$ at `:75`, best rank `:80`)

- *Notes:* Printed chance comparator `:112-113` is $K/(D+1)$; exact is $1-\binom{D}{K}/\binom{D+m}{K}$ (O7)


**V2-transcriptome-C4. Co-expression rank test**


one-sided Wilcoxon rank-sum: $H_1$: within-module $\lvert r\rvert$ stochastically $>$ between-module $\lvert r\rvert$; flag if $p<0.05$


- *Variables:* pairwise gene-gene correlations of a genotype-free transcriptome

- *Reference:* Wilcoxon 1945 *Biometrics Bull.* 1(6):80-83 (page unverified)

- *Function:* benchmark 03

- *Location:* `benchmarks/03_coexpression_fp_control.R:55-70`

- *Notes:* Explicitly framed as the *naive* inference; "FPR" is relative to the stated truth $h^2 = 0$, not to the Wilcoxon null (O6, INFO)


**V2-transcriptome-C5. Benchmark-04 mediation target**


claimed $m_g \approx p\,h^2_*$; exact $m_g = \mathrm{Var}\!\big(c\sum_g w_g (G_g-\bar G_g)/s_{E_g}\big)/V_P$ (reduces to $p\sum w_g^2 h_g^2/\sum w_g^2 \cdot (p/V_P\text{-adj})$ only for independent genes)


- *Variables:* $p$ = `prop`; $w_g$ slopes; $s_{E_g}$ sd of expression

- *Reference:* none cited; algebra

- *Function:* benchmark 04

- *Location:* `benchmarks/04_mediation_recovery.R:9-12, 89-90, 114-116`

- *Notes:* Generator run: actual 0.2771, $p\,\overline{h^2}$ 0.2560, exact 0.2771 (G)


**V2-transcriptome-C6. TWAS Pearson test**


$t = r\sqrt{(n-2)/(1-r^2)}$, $p = 2F_{t_{n-2}}(-\lvert t\rvert)$


- *Variables:* $r$ gene-phenotype correlation

- *Reference:* standard Pearson-correlation $t$ test; page unverified

- *Function:* `cor_p`

- *Location:* `benchmarks/05_twas_power.R:103-112` ($t$ at `:110`, $p$ at `:111`)

- *Notes:* Equation correct; single-realization "power" (O10)


**V2-transcriptome-C7. Bonferroni threshold**


$\alpha_{\text{gene}} = 0.05/T$


- *Variables:* $T$ genes tested

- *Reference:* Bonferroni; page unverified

- *Function:* benchmark 05

- *Location:* `benchmarks/05_twas_power.R:114, 146-147`

- *Notes:* Threshold correct


## Input/output and QC: MAF, LD r2, PLINK pruning and blocks

*Source report:* `v2-io-formats.md`  
*Entries:* 14


**V2-io-formats-1. Allele frequency / MAF**


$p=\dfrac{\sum_i d_i}{2\,n_{\text{called}}},\quad \text{MAF}=\min(p,1-p)$


- *Variables:* $d_i\in\{0,1,2\}$ dosage of the allele coded 2 (from `code_as`), $n_{\text{called}}$ = non-NA calls

- *Reference:* PLINK 1.9 `--maf` semantics (docs; Purcell 2007 for the tool)

- *Function:* `filter_geno`

- *Location:* `R/qc_filter_geno.R:290-297`

- *Notes:* `n_called == 0 -> MAF 0`; verified 4/10 on `0 0 0 2 2 NA NA NA`


**V2-io-formats-2. Monomorphic / MAF filters**


keep iff $\text{MAF}>0$; $\text{MAF}\ge m_{\text{above}}$; $\text{MAF}\le m_{\text{below}}$


- *Variables:* —

- *Reference:* roxygen `@param` (`>=`, `<=`)

- *Function:* `filter_geno`

- *Location:* `:311`, `:314`, `:318`

- *Notes:* boundaries executed exact (exp3 F1)


**V2-io-formats-3. Heterozygosity**


$n_{\text{het}}=\#\{i: d_i=1\}$; include iff $n_{\text{het}}>0$, remove iff $n_{\text{het}}=0$


- *Variables:* —

- *Reference:* roxygen

- *Function:* `filter_geno`

- *Location:* `:301`, `:322`, `:324`

- *Notes:* count, not rate; no rate threshold exists


**V2-io-formats-4. Composite (genotypic) $r^2$**


$r^2=\dfrac{(n\sum xy-\sum x\sum y)^2}{(n\sum x^2-(\sum x)^2)(n\sum y^2-(\sum y)^2)}$


- *Variables:* $x,y$ dosages of two markers over $n$ individuals (integer sums)

- *Reference:* PLINK 1.9 docs `--indep-pairwise` ("correlations between genotype allele counts"); `plink_ld.c` `ld_prune` (`cov12*cov12*dxx`)

- *Function:* `.ld_sweep`

- *Location:* `:556-565`

- *Notes:* = squared Pearson correlation; NA case falls back to `cor(use="pairwise.complete.obs")` `:402`


**V2-io-formats-5. Pruning decision**


prune iff $r^2>\theta(1+\varepsilon)$, $\varepsilon=2^{-44}$; drop $i$ iff $\text{MAF}_i<(1-\varepsilon)\text{MAF}_j$ else drop $j$


- *Variables:* $\theta$ user threshold

- *Reference:* `plink_ld.c` `ld_prune` (fetched; "remove marker with lower MAF")

- *Function:* `.ld_sweep`

- *Location:* `:532`, `:663`, `:665-668`

- *Notes:* tie $\rightarrow$ later marker pruned (earlier kept); executed exp3 F7


**V2-io-formats-6. Window / step**


window of $w$ markers (or all markers within $w\cdot 1000$ bp of the start), slide by $s$, per chromosome in position order


- *Variables:* $w,s$

- *Reference:* PLINK 1.9 docs

- *Function:* `.ld_prune`, `.ld_sweep`

- *Location:* `:492-506`, `:605-617` (`win_at`), `:737-782` (slide)

- *Notes:* `step` rounded to integer $\geq$ 1; `window` fractional accepted (exp3 F7)


**V2-io-formats-7. Haplotypic $r^2$**


$r^2=\dfrac{(f_{11}-f_{1\cdot}f_{\cdot 1})^2}{f_{1\cdot}f_{2\cdot}f_{\cdot1}f_{\cdot2}}=\dfrac{D^2}{p_1q_1p_2q_2}$


- *Variables:* $f_{11}$ ML haplotype frequency (cubic root with max log-likelihood), marginals $f_{1\cdot},f_{\cdot1}$

- *Reference:* Hill & Robertson 1968 *TAG* 38:226–231 (page of eq. not verified); PLINK `em_phase_hethet` port

- *Function:* `.plink_hap_rsq`, `.plink_em_hethet`, `.plink_cubic_roots`

- *Location:* `:796-809` (r² at `:808`), `:822-879`, `:902-943`

- *Notes:* independent EM agrees to 7 digits (exp3 F8)


**V2-io-formats-8. VIF**


$\text{VIF}_k=[R^{-1}]_{kk}$; prune $\arg\max_k \text{VIF}_k$ while $>v$


- *Variables:* $R$ correlation matrix of the window's surviving markers

- *Reference:* PLINK 1.9 docs `--indep`; Purcell 2007 (VIF pruner)

- *Function:* `.ld_sweep` (vif branch)

- *Location:* `:683-736` (`solve` `:697`,`:721`; test `:725`)

- *Notes:* singular ($\text{rcond}<10^{-14}$) handling `:362`,`:435-438`,`:533-557` UNVERIFIABLE vs source


**V2-io-formats-9. Gabriel block MAF floor / tag**


keep markers with $\text{MAF}\ge 0.05(1-\varepsilon)$; tag = $\arg\max$ MAF in block


- *Variables:* —

- *Reference:* Gabriel 2002; PLINK `--blocks` (Haploview)

- *Function:* `.gabriel_blocks`

- *Location:* `R/qc_ld_methods.R:25`, `:40`

- *Notes:* classification thresholds (D' CI percentiles `:814-815`, informative fraction 0.95 `:906`) in `.plink_blocks_classify` `R/qc_filter_geno.R:978-1053` and `.plink_blocks_chrom` `:901-1015`


**V2-io-formats-10. Major-allele flip (frequency)**


flip iff $\#\{d=2\}>\#\{d=0\}$


- *Variables:* homozygote counts (hets cancel)

- *Reference:* package convention (v1)

- *Function:* `compute_flip`

- *Location:* `R/io_detect_format.R:362-364`

- *Notes:* tie $\rightarrow$ no flip (root of IO-F2)


**V2-io-formats-11. Reference flip**


flip iff $\text{allele}_1\ne\text{ref}$


- *Variables:* per-marker labels

- *Reference:* package convention

- *Function:* `compute_flip`

- *Location:* `:359`

- *Notes:* only hapmap/table may use it (`R/io_format_conversion.R:214-219`)


**V2-io-formats-12. Numeric coding**


$-101$: $(0,\text{het},2)\mapsto(+1,0,-1)$ unflipped, $(−1,0,+1)$ flipped; $012$: $(2,1,0)$ / $(0,1,2)$; Dom: non-het $\rightarrow$ minor; Left: het $\rightarrow$ minor; Right: het $\rightarrow$ major


- *Variables:* raw dosage, flip

- *Reference:* v1 `create_phenotypes` transforms

- *Function:* `numericalize_core`

- *Location:* `src/rust/src/numeric.rs:76-76`, `:138-145`, `:156-178`

- *Notes:* imputation `:47-52`, `:72-76` before the model transform


**V2-io-formats-13. HWE**


—


- *Variables:* —

- *Reference:* —

- *Function:* —

- *Location:* —

- *Notes:* **not implemented anywhere** (no exact/χ² test)


**V2-io-formats-14. Missing-rate filter**


—


- *Variables:* —

- *Reference:* —

- *Function:* —

- *Location:* —

- *Notes:* **not implemented** (`--geno`/`--mind` absent)


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V2-io-formats-C1. Raw genotype parsing and coding dispatch**


call $\to g\in\{0,1,2,\mathrm{NA}\}$ (0 = hom allele-1, 1 = het, 2 = hom allele-2); then major/het/minor $=(1,0,-1)$ or $(2,1,0)$ after optional flip


- *Variables:* $g$ raw dosage; flip from `compute_flip`

- *Reference:* Package-defined convention; no external source

- *Function:* `parse_hapmap_chars_to_raw()`, `.apply_coding()`

- *Location:* `R/io_detect_format.R:223-326`; `R/io_read_formats.R:18-194`

- *Notes:* Codex wrote `195-276`; function closes at `:277`. Contract is violated by the SNPRelate readers (IO-F2): they deliver raw 2 = hom allele-1.


**V2-io-formats-C2. Exact two-locus ML cubic**


$z^3+az^2+bz+c=0$, $a=\tfrac12(f_{11}+f_{22}-f_{12}-f_{21}-3h)$, $b=\tfrac12\big(f_{11}f_{22}+f_{12}f_{21}+h(f_{12}+f_{21}-f_{11}-f_{22}+h)\big)$, $c=-\tfrac12\,h\,f_{11}f_{22}$


- *Variables:* $f_{ab}$ known-haplotype shares; $h$ = double-het count / $2N$ (`hhs`); $z$ coupling increment

- *Reference:* Gaunt, Rodríguez & Day 2007 *BMC Bioinformatics* 8:428 (page unverified; **not in `@references`**); PLINK 1.9 `plink_ld.c` `em_phase_hethet` (no page)

- *Function:* `.plink_em_hethet()`

- *Location:* `R/qc_filter_geno.R:822-879`

- *Notes:* Codex wrote `646-713`; verified coefficients at `:674-677`. Feasible roots clipped to $[0,h]$ (`:679-686`) and compared by log-likelihood — exact root selection, not iterative EM.


**V2-io-formats-C3. Two-locus haplotype log-likelihood**


$\ell=c_{DH}\log(f_{11}f_{22}+f_{12}f_{21})+\sum_{ab}c_{ab}\log f_{ab}$, with $f_{11},f_{22}\mathrel{+}=z$, $f_{12},f_{21}\mathrel{+}=h-z$


- *Variables:* $c_{ab}$ unambiguous haplotype counts; $c_{DH}$ double-heterozygotes

- *Reference:* Gaunt et al. 2007 (page unverified); PLINK `plink_ld.c`

- *Function:* `.plink_calc_lnlike()`

- *Location:* `R/qc_filter_geno.R:884-894`

- *Notes:* Codex wrote `716-728`. Zero-count terms skipped (`:724-727`).


**V2-io-formats-C4. Cubic real roots**


$Q=(a^2-3b)/9$, $R=(2a^3-9ab+27c)/54$; if $R^2<Q^3$ three real roots $-2\sqrt{Q}\cos\!\big(\tfrac{\theta+2k\pi}{3}\big)-a/3$, $\theta=\arccos(R/Q^{3/2})$


- *Variables:* monic-cubic coefficients $a,b,c$

- *Reference:* Classical trigonometric Cardano (no page applicable); PLINK 1.9 `plink_common.c` `cubic_real_roots`

- *Function:* `.plink_cubic_roots()`

- *Location:* `R/qc_filter_geno.R:902-943`

- *Notes:* Codex wrote `731-777`. PLINK constants `EPSILON`, $\pi$ literal, and operation order retained (`:738-749`).


**V2-io-formats-C5. D' likelihood surface and class codes**


$f_{11}(q)=f_{1\cdot}f_{\cdot1}+q\cdot\mathrm{denom}$, $q=0..100$ (D' percentile), other cells from marginals; $\ell(q)$ as above; classify by thresholds `recomb_highci=89`, `strong_highci=97`, `strong_lowci=72`, `strong_lowci_outer=71`


- *Variables:* $q$ percentile of $D'$; `denom` $=D_{\max}/100$ after sign orientation

- *Reference:* Gabriel et al. 2002 *Science* 296:2225-2229 (method page unverified); PLINK `plink_ld.c` `haploview_blocks_classify`

- *Function:* `.plink_calc_lnlike_quantile()`, `.plink_blocks_classify()`

- *Location:* `R/qc_filter_geno.R:952-966`, `:978-1053`

- *Notes:* Codex wrote `780-887`. Thresholds verified at `:814-815`; 0.95 mass rule at `:874`.


**V2-io-formats-C6. Gabriel/Haploview block acceptance**


informative fraction $=0.95+\varepsilon_{ish}$; two-marker span threshold $1+\lfloor 3\cdot\mathrm{frac}\rfloor$, three-marker $\lfloor 6\cdot\mathrm{frac}\rfloor$; strong-pair fraction of informative pairs $>0.95$; greedy largest-span non-overlapping blocks


- *Variables:* strong / recombination pair counts; span in bp

- *Reference:* Gabriel et al. 2002 (method page unverified); PLINK `plink_ld.c` `haploview_blocks`

- *Function:* `.plink_blocks_chrom()`, `.gabriel_blocks()`

- *Location:* `R/qc_filter_geno.R:1066-1180`; `R/qc_ld_methods.R:23-46`

- *Notes:* Codex wrote `890-1014` / `15-37`. `inform_frac` at `:906-908`. Tag = highest-MAF member is a package convenience (`qc_ld_methods.R:40`).


**V2-io-formats-C7. Per-QTN marginal realized variance share**


$k^2=\pi_t/s_{\text{raw}}^2$; $v_j=k^2\,\mathrm{Var}(c_j)/V_P$


- *Variables:* $\pi_t$ layer prop for trait $t$; $c_j=e_jG_j$ (additive), $e_j\,I(G_j=0)$ (dominance), $a_jG_j+d_jI(G_j=0)$ (orthogonal); $V_P$ realized $\mathrm{var}(y)$

- *Reference:* Package convention, `docs/SPEC.md` §2; no external source

- *Function:* `.qtn_var()`, `qtn_table()`

- *Location:* `R/io_write.R:1288-1312`; `:1400-1403`

- *Notes:* Marginal shares; with LD among causal loci they need not sum to $\pi_t$ (roxygen `:188-191`). NA for vqtl/epistasis (`:196-198`).


**V2-io-formats-C8. Per-gene transcriptome marginal share**


$k=\sqrt{\pi_t}/s_{\text{comp}}$; $w_g=\text{slope}_g/\max_g\vert{}\text{slope}_g\vert{}$; $v_g=(k\,w_g)^2\,\mathrm{Var}(z_g)/V_P$ with $\mathrm{Var}(z_g)\in\{0,1\}$


- *Variables:* $z_g$ standardized expression; constant gene $\rightarrow$ 0

- *Reference:* Package convention, `docs/SPEC.md`; no external source

- *Function:* `.tx_qtn_var()`

- *Location:* `R/io_write.R:1328-1345`

- *Notes:* Codex wrote `235-251`; closes at `:252`. Co-expression covariance not attributed to single genes.


**V2-io-formats-C9. PLINK fixture allele serialization**


$-1\mapsto A_1A_1$, $0\mapsto A_1A_2$, $1\mapsto A_2A_2$


- *Variables:* package dosage; `allele` column labels

- *Reference:* Package fixture; PLINK TPED format (manual page not cited)

- *Function:* fixture generator

- *Location:* `data-raw/plink_parity_fixtures.R:26-34`

- *Notes:* Codex wrote `19-35` (includes comments). r², VIF, MAF invariant to a consistent A1/A2 swap. NA dosage would serialize as `NA NA` (IO-F15).


# Version 1 (frozen legacy create_phenotypes engine)


## V1 legacy: create_phenotypes core, residual/h2 machinery, linkage architecture, user QTNs

*Source report:* `v1-core-linkage.md`  
*Entries:* 19


**V1-core-linkage-1. Additive genetic value**


$g^{A}_i=\sum_k a_k x_{ik}$, $x\in\{-1,0,1\}$


- *Variables:* $a_k$ effect of QTN k; $x_{ik}$ dosage code

- *Reference:* Fernandes & Lipka 2020 (additive model; page unverified); F&M genotypic value scale $-a,d,+a$

- *Function:* `genetic_effect`

- *Location:* `legacy_genetic_effect.R:109-112`

- *Notes:* Coding centred at the heterozygote; not Fisher's average effect


**V1-core-linkage-2. Dominance value**


$g^{D}_i=\sum_k d_k\,\mathbf 1[x_{ik}=0]$


- *Variables:* $d_k$ dominance effect

- *Reference:* idem (page unverified)

- *Function:* `genetic_effect`

- *Location:* `:61-75`

- *Notes:* Het indicator; if a QTN has no het, $V_{D,k}=0$ silently


**V1-core-linkage-3. Epistatic value**


$g^{E}_i=\sum_k e_k\prod_{m\in S_k}x_{im}$


- *Variables:* $S_k$ set of `epi_interaction` markers

- *Reference:* idem, A x A only (page unverified)

- *Function:* `genetic_effect`

- *Location:* `:77-89`

- *Notes:* Product of -1/0/1 codes (sign flips)


**V1-core-linkage-4. Centering**


$g_i \leftarrow g_i-\bar g$


- *Reference:* -

- *Function:* `genetic_effect`

- *Location:* `:93`

- *Notes:* Hence phenotype mean = `mean`


**V1-core-linkage-5. Degree of dominance**


$d_k=\delta\,a_k$


- *Variables:* $\delta$ = `degree_of_dom`

- *Reference:* F&M partial/complete/over-dominance ($d/a$)

- *Function:* `check_in`

- *Location:* `legacy_check_in.R:274`

- *Notes:* Only with `same_add_dom_QTN`


**V1-core-linkage-6. Geometric series**


$a_k=a^{k}$, $k=1..n$; with major QTN: $(a_{big},a^{1},\dots,a^{n-1})$


- *Reference:* Fernandes & Lipka 2020 (geometric series; page unverified)

- *Function:* `check_in`

- *Location:* `:786-788, 796-798, 810-812, 823-825`

- *Notes:* Also applied to full-length custom vectors when mixed (F3)


**V1-core-linkage-7. Residual variance**


$\sigma_e^2=\frac{V_G(1-h^2)}{h^2}$, $V_G=\frac1{n-1}\sum(g_i-\bar g)^2$


- *Variables:* $h^2$ requested; $V_G$ realized sample variance

- *Reference:* $h^2=V_G/V_P$, $V_P=V_G+V_E$ (F&M; page unverified)

- *Function:* `phenotypes`

- *Location:* `legacy_Phenotypes.R:216, 650, 1065, 1511`

- *Notes:* Per rep for vary_QTN; n-1 denominator; no clamp on $h^2$ (F8)


**V1-core-linkage-8. Single-trait phenotype**


$y_{ij}=g_i+e_{ij}+\mu$, $e_{ij}\sim N(0,\sigma_e^2)$


- *Variables:* j = rep

- *Reference:* —

- *Function:* `phenotypes`

- *Location:* `:669-674, 1517-1522`

- *Notes:* `rnorm` inversion; seed per rep (F2)


**V1-core-linkage-9. Multi-trait residual**


$\Sigma_e=D^{1/2}RD^{1/2}$, $D=\mathrm{diag}(\sigma^2_{e,t})$, $R=$`cor_res`


- *Reference:* —

- *Function:* `phenotypes`

- *Location:* `:238-245, 1080-1086`

- *Notes:* `mvtnorm::rmvnorm` (eigen method); $R$ unvalidated


**V1-core-linkage-10. Null trait (h2=0)**


$y_{ij}\sim N(\mu,1)$; multi: $N(\mu,R)$


- *Reference:* —

- *Function:* `phenotypes`

- *Location:* `:533-534, 86-92`

- *Notes:* Unit variance, not $V_G$-scaled


**V1-core-linkage-11. Realized heritability**


$\hat h^2_{tj}=V_{G,t}/\widehat{\mathrm{Var}}(y_{tj})$, averaged over j


- *Reference:* —

- *Function:* `phenotypes`

- *Location:* `:248-251, 336, 675-676, 876`

- *Notes:* Diagnostic only


**V1-core-linkage-12. Per-QTN PVE**


$\mathrm{PVE}_{kj}=\mathrm{Var}(a_k x_k)/\widehat{\mathrm{Var}}(y_j)$


- *Reference:* -

- *Function:* `phenotypes`

- *Location:* `:260-268, 685-693`

- *Notes:* Ignores LD covariance; rep-1 numerator under vary_QTN multi-trait (F12)


**V1-core-linkage-13. Residual seed**


$s_{j}=\lfloor (s_0+j)\cdot\mathrm{round}(10h^2_1)\rfloor$; null: $s_0+j$


- *Variables:* $s_0$ master seed

- *Reference:* - (roxygen `:161-167` states a different formula)

- *Function:* `phenotypes`

- *Location:* `:234, 665, 1076, 1513; 82, 529, 924, 1374`

- *Notes:* F2


**V1-core-linkage-14. LD-QTN seed**


$s=(s_0\cdot a)+z$ (add), $+\,$`rep` (dom), $+\,x$ (direct resample)


- *Variables:* $a$ attempt, $z$ rep, $x$ QTN index

- *Reference:* -

- *Function:* `qtn_linkage`

- *Location:* `legacy_QTN_linkage.R:110, 383, 645, 971, 1220, 1243, 1455, 1479`

- *Notes:* Same seed re-set before each inner resample -> deterministic walk (F4)


**V1-core-linkage-15. LD measure**


$\ell=\lvert\mathrm{LD}_{method}(x_j,x_i)\rvert$, accept iff $\ell_{min}\le\ell\le\ell_{max}$


- *Variables:* method in {composite, r, dprime, corr}

- *Reference:* SNPRelate `snpgdsLDpair` docs (primary source for the composite measure not cited in package: UNVERIFIABLE)

- *Function:* `qtn_linkage`

- *Location:* `:150-151, 176-177, 1025-1028, 1270-1273`

- *Notes:* Correlation scale, not $r^2$ (INFO-3); `LD_between_QTNs` unsigned-ness not applied (`:198-199`)


**V1-core-linkage-16. Indirect selection**


QTN$_1$ = first $i>j$ with $\ell(j,i)\le\ell_{max}$; QTN$_2$ = first $i'<j$ likewise; require both $\ge\ell_{min}$


- *Variables:* j = intermediate ("cause") marker

- *Reference:* Fernandes & Lipka 2020 "indirect" (page unverified)

- *Function:* `qtn_linkage`

- *Location:* `:125-212`

- *Notes:* Index-based walk, chromosome-agnostic; trait 1 <- higher index (`sup`)


**V1-core-linkage-17. Direct selection**


QTN$_2$ = whichever of $j\pm1$ has $\ell$ in window, nearest $\ell_{max}$; else resample $j$


- *Reference:* idem "direct"

- *Function:* `qtn_linkage`

- *Location:* `:1229-1308` (fixed), `:978-1062, 1465-1546` (unfixed)

- *Notes:* Fixed branch adjacent-only (INFO-4); unfixed branches F5


**V1-core-linkage-18. MAF**


$p=\frac{\sum_i x_i+n}{2n}$, $\mathrm{MAF}=\min(p,1-p)$


- *Variables:* n taxa

- *Reference:* standard allele frequency from dosage

- *Function:* `constraint`, `qtn_linkage`, `qtn_from_user`

- *Location:* `legacy_constraint.R:47-50`; `QTN_linkage.R:314-317`; `qtn_from_user.R:237-240`

- *Notes:* Strict bounds in constraint (F23)


**V1-core-linkage-19. Genetic correlation (given `cor`)**


$G_s = \mathrm{scale}(G)$; $W=L_g^{-1}G_s^\top$; $T=L_cW$ with $L_gL_g^\top=\mathrm{cov}(G_s)$, $L_cL_c^\top=$`cor`; $T_t\leftarrow T_t\,s_t+\bar g_t$


- *Reference:* Cholesky whitening/colouring (standard); PleioArch not used in v1

- *Function:* `base_line_multi_traits` (not in group)

- *Location:* `legacy_Base_line_multi_traits.R:167-183`

- *Notes:* Realized sample correlation == `cor` exactly (INFO-2); upper triangle only (F18)


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V1-core-linkage-C1. Heterozygote eligibility indicator**


$I_{het}(x)=\mathbf 1[x=0]$


- *Variables:* $x\in\{-1,0,1\}$; `0` = Aa under the documented coding

- *Reference:* Package coding contract (`create_phenotypes.R:100-102`, roxygen "aa = -1, Aa = 0 and AA = 1"); external page UNVERIFIABLE

- *Function:* `constraint()` (correct); `create_phenotypes()` guard (wrong code)

- *Location:* `legacy_constraint.R:28-42`; `legacy_create_phenotypes.R:100-102, 848-935`

- *Notes:* `constraint()` tests `== 0`; the guard tests `b == 1` (V1C-F11 / Codex F06). Verified: `hets <- NULL` at `:848`, `if (any(!unlist(hets)))` at `:922`, block closes at `:935`.


**V1-core-linkage-C2. Reported residual "sample correlation"**


$R_{reported}=\mathrm{cov2cor}\!\left(\tfrac1r\sum_{z=1}^{r}\Sigma_{E,z}\right)$, $\Sigma_{E,z}=D^{1/2}RD^{1/2}$


- *Variables:* $r$ replicates; $\Sigma_{E,z}$ the *target* covariance, not realized residuals

- *Reference:* Implementation-defined; the "sample" label conflicts with the paper's stated diagnostic (page unverified)

- *Function:* `phenotypes()`

- *Location:* `legacy_Phenotypes.R:483-491, 1335-1343` (verified: `matrix(0,...)` at `:484/:1336`, `cov2cor` at `:488/:1340`)

- *Notes:* Equals `cor_res` exactly whenever `cor_res` is a valid correlation matrix (V1C-F10 / Codex F10). Not an estimate.


**V1-core-linkage-C3. QTN anchor sampling without replacement**


$J\sim\mathrm{Sample}(\mathcal E,m,\text{replace}=\mathrm{FALSE})$ after `set.seed(seed*s+z)` (add) / `set.seed(seed*s+z+rep)` (dom)


- *Variables:* $\mathcal E$ = `index` (constraint-filtered rows or all rows); $m$ = QTN count; $s$ outer attempt, $z$ rep

- *Reference:* Implementation-defined; Fernandes & Lipka 2020 (random QTN selection; page unverified)

- *Function:* `qtn_linkage()`

- *Location:* `legacy_QTN_linkage.R:109-113, 382-386, 644-648, 970-974, 1219-1223, 1454-1458` (verified: `while (s <= 10 & border)` at `:109,382,644,970,1219,1454`; `sample(index, ...)` at `:113,386,648,974,1223,1458`)

- *Notes:* Uses R's `sample()` under `RNGversion("3.5.1")` ("Rounding"). Constraints apply to anchors only; walked partners are unfiltered (Codex F14 / Fable F23).


**V1-core-linkage-C4. User marker resolution and order**


$I=(\mathrm{match}(m_1,S),\dots,\mathrm{match}(m_k,S))$ realised as `genotypes[genotypes$snp %in% i, ]` then `a[i, ]` (named-row indexing in requested order)


- *Variables:* $m$ requested IDs; $S$ genotype IDs

- *Reference:* Implementation-defined API algorithm

- *Function:* `qtn_from_user()`

- *Location:* `legacy_qtn_from_user.R:1358-1388` (add), `:1447-1460` (dom), `:1517-1530` (epi), `:1589-1602` (var) (verified: `selected_snps <- genotypes$snp %in% i` at `:1361, 1450, 1520, 1592`)

- *Notes:* Only `sum(selected_snps) == 0` is rejected; a partially missing list yields an all-`NA` row named `Chr_NA_NA` (Codex F08; public outcome is a cryptic error, R5). Single-trait assembly wraps the result one level too deep (V1C-F1).


**V1-core-linkage-C5. Epistatic interaction diagnostic expansion**


$n_{groups}=n_{markers}/k$; `QTN = rep(1:n_groups, each = k)`; effect column `rep(effect, each = n_groups)` (should be `each = k`)


- *Variables:* $k$ = `epi_interaction`

- *Reference:* Implementation-defined API grouping

- *Function:* `qtn_from_user()`

- *Location:* `legacy_qtn_from_user.R:1582-1605` (verified: `e_len <- lengths(QTN_list$epi)/epi_interaction` at `:1532`, `each = epi_interaction` at `:1541`, `each = e_len[i]` at `:1544`)

- *Notes:* Correct only when $n_{groups}=k$ (V1C-F13 / Codex F09).


**V1-core-linkage-C6. *(reconciler-added)* Direct-LD acceptance sentinel**


Loop runs while $(\ell_{sup}\notin[\ell_{min},\ell_{max}])\wedge(\ell_{inf}\notin[\ell_{min},\ell_{max}])$ with initial $\ell_{sup}=\ell_{inf}=1$; on exit partner $=i-1$ or $i_2+1$


- *Variables:* sentinel value 1

- *Reference:* -

- *Function:* `qtn_linkage()`

- *Location:* `legacy_QTN_linkage.R:1269-1269` (fixed add branch); same pattern `:123-125` (indirect add), `:976-983`, `:1467-1473`

- *Notes:* When $\ell_{max}=1$ the sentinel is "in window", the loop never executes and the partner is the anchor itself (RX-1, R6e); indirect errors "object 'i' not found" (R6f).


## V1 legacy: pleiotropic and partially pleiotropic QTNs, genetic effects, vQTL, baselines

*Source report:* `v1-pleio-effects.md`  
*Entries:* 18


**V1-pleio-effects-1. Additive genetic value**


$g^{A}_i=\sum_{k=1}^{n_A} x_{ik}\,a_k$


- *Variables:* $x_{ik}\in\{-1,0,1\}$ (+1 = major-allele homozygote), $a_k$ additive effect

- *Reference:* Fernandes & Lipka 2020, Implementation (page unverified); Falconer & Mackay 1996 ch. 7 (page unverified)

- *Function:* `genetic_effect`

- *Location:* `R/legacy_genetic_effect.R:111–57`

- *Notes:* Not Fisher's average-effect decomposition; no allele-frequency scaling


**V1-pleio-effects-2. Geometric effect series**


$a_k=a^{k},\ k=1,\dots,n$; with big QTN $(a_{\text{big}},a^{1},\dots,a^{n-1})$


- *Variables:* $a$ = `add_effect`, $n$ = QTN count

- *Reference:* Fernandes & Lipka 2020 ("geometric series of effect sizes", page unverified)

- *Function:* `check_in`

- *Location:* `R/legacy_check_in.R:905–837`

- *Notes:* Length-mismatch recycling (F3); mixed-type conversion (F4); applies equally to `dom_`, `epi_`, `var_effect`


**V1-pleio-effects-3. Dominance genetic value**


$g^{D}_i=\sum_{k=1}^{n_D}\mathbb{1}[x_{ik}=0]\,d_k$


- *Variables:* $d_k$ dominance effect; `same_add_dom_QTN`: $d_k=a_k\cdot\text{degree\_of\_dom}$

- *Reference:* paper Table 1 (UNVERIFIABLE); code

- *Function:* `genetic_effect`

- *Location:* `R/legacy_genetic_effect.R:120–77`; `R/legacy_check_in.R:274`

- *Notes:* Heterozygote indicator only; hetless locus contributes 0 with `var_dom=0`


**V1-pleio-effects-4. Epistatic genetic value**


$g^{E}_i=\sum_{m=1}^{n_E} e_m\prod_{l\in S_m} x_{il}$, $\vert{}S_m\vert{}=$ `epi_interaction`


- *Variables:* $e_m$ interaction effect

- *Reference:* paper "additive x additive epistatic" (page unverified)

- *Function:* `genetic_effect`

- *Location:* `R/legacy_genetic_effect.R:151–91`

- *Notes:* **Uncentered** dosage product; any heterozygote zeroes the term; $k$-way allowed; `epi_interaction=1` errors


**V1-pleio-effects-5. Total genetic value**


$g_i=(g^{A}_i+g^{D}_i+g^{E}_i)-\overline{g}$


- *Variables:* —

- *Reference:* —

- *Function:* `genetic_effect`

- *Location:* `R/legacy_genetic_effect.R:186–96`

- *Notes:* Mean-centred; `mean` added later


**V1-pleio-effects-6. Component variances**


$V_A=\widehat{\mathrm{Var}}(g^{A})$ etc. (sample, $n-1$); per-QTN $\mathrm{Var}(x_{\cdot k}a_k)$


- *Variables:* —

- *Reference:* —

- *Function:* `genetic_effect`

- *Location:* `R/legacy_genetic_effect.R:115,59,69,76,87,90`

- *Notes:* Per-QTN variances ignore LD covariance; Σ $\neq$ $V_A$


**V1-pleio-effects-7. MAF of a QTN**


$p=\tfrac12\left(\tfrac{\sum_i x_i+n}{n}\right)$, $\mathrm{MAF}=\min(p,1-p)$


- *Variables:* $n$ = individuals

- *Reference:* —

- *Function:* `qtn_pleiotropic`, `qtn_partially_pleiotropic`, `constraint`

- *Location:* `R/legacy_QTN_pleiotropic.R:205–165`; `…partially….R:165–169`; `R/legacy_constraint.R:47–43`

- *Notes:* Rounded to 4 dp in files; NA-unsafe


**V1-pleio-effects-8. QTN sampling seeds (pleiotropic)**


add: $s+i$; dom: $s+i+r$; var: $2s+i$; epi: $2s+i+r$


- *Variables:* $s$ = seed, $i$ = rep index, $r$ = `rep` (1 unless `vary_QTN`)

- *Reference:* DECISION‑008 (legacy seed math)

- *Function:* `qtn_pleiotropic`

- *Location:* `R/legacy_QTN_pleiotropic.R:174,318,408,485`

- *Notes:* Collisions (F8); sets not disjoint


**V1-pleio-effects-9. QTN sampling seeds (partially)**


pleio add $s+j$; spec add $s+i+j$; pleio dom $s+j+r$; spec dom $s+i+j+r$; pleio epi $2s+j$; spec epi $2s+i+j$


- *Variables:* $j$ = rep, $i$ = trait

- *Reference:* —

- *Function:* `qtn_partially_pleiotropic`

- *Location:* `…partially….R:125,153,446,475,621,637`

- *Notes:* Pools: `setdiff` (disjoint) except F6


**V1-pleio-effects-10. Het re-sample loop**


re-draw while $\neg\exists (i,k): x_{ik}=0$, $\leq$10 tries, same seed, pool $\setminus$ rejected


- *Variables:* —

- *Reference:* —

- *Function:* both QTN samplers

- *Location:* `R/legacy_QTN_pleiotropic.R:180–145,324–333`; `…partially….R:131–140,160–170,452–461,481–493`

- *Notes:* "any het in the whole set", not per QTN


**V1-pleio-effects-11. Whitening**


$G_s=\mathrm{scale}(G)$, $C_g=\mathrm{cov}(G_s)$, $W=G_s\,L_g^{-\top}$, $L_gL_g^{\top}=C_g$


- *Variables:* $G$ = $n\times t$ genetic values

- *Reference:* Fernandes & Lipka 2020, Implementation (whitening/colouring paragraph; page unverified); Kessy et al. 2018 (page unverified)

- *Function:* `base_line_multi_traits`

- *Location:* `R/legacy_Base_line_multi_traits.R:151–109,178–184`

- *Notes:* $C_g$ passed through `make_pd`


**V1-pleio-effects-12. Colouring & rescale**


$T=W\,L_c^{\top}$, $L_cL_c^{\top}=\text{cor}$; $\tilde g_{\cdot i}=T_{\cdot i}\,\mathrm{sd}(G_{\cdot i})+\overline{G_{\cdot i}}$


- *Variables:* `cor` user matrix ($t\times t$ only)

- *Reference:* same

- *Function:* `base_line_multi_traits`

- *Location:* `…multi_traits.R:160–120,185–194`

- *Notes:* Exact sample correlation; trait $\ge2$ effects/PVE not effective (F2); variance preserved iff $\mathrm{diag}(\text{cor})=1$


**V1-pleio-effects-13. PD repair**


$\tau=\max(0,\ 2n\,\varepsilon\max\vert{}\lambda\vert{}-\lambda)$, $m'=\mathrm{round}(m+V\mathrm{diag}(\tau)V^{\top},\,2)$


- *Variables:* $\lambda,V$ eigen-pairs of $m$

- *Reference:* none (ad hoc)

- *Function:* `make_pd`

- *Location:* `R/legacy_make_pd.R:61–21`

- *Notes:* Not a correlation matrix afterwards; rounding re-breaks PD (F1)


**V1-pleio-effects-14. Residual variance (context)**


$\sigma^2_{e}=v_g/h^2-v_g$, $v_g=\widehat{\mathrm{Var}}(\tilde g)$; multi-trait $\Sigma_e=\mathrm{diag}(\sigma_e)\,\text{cor\_res}\,\mathrm{diag}(\sigma_e)$


- *Variables:* —

- *Reference:* paper: trait variance = genetic + error (page unverified)

- *Function:* `phenotypes`

- *Location:* `R/legacy_Phenotypes.R:216,239,650`

- *Notes:* Post-transform $v_g$; seed $(s+j)\cdot\mathrm{round}(10h^2_1)$


**V1-pleio-effects-15. vQTL residual scale**


$\sigma_i=1+\sum_k v_k (x_{ik}+1)$; $k=\sqrt{(1/h^2-1)}/\mathrm{median}_i(\sigma_i)$; $e_{ij}\sim N(0,(k\sigma_i)^2)$; $y_{ij}=\mathrm{scale}(g)_i+e_{ij}+\mu$


- *Variables:* $v_k$ = `var_effect`, $x+1\in\{0,1,2\}$

- *Reference:* Murphy et al. 2022 Heredity 129:93–102 (form UNVERIFIABLE)

- *Function:* `vQTL`

- *Location:* `R/legacy_vQTL.R:89,62–71,76–89`

- *Notes:* $g$ standardised to unit variance; SD linear in dosage (not log-linear); $h^2$ exact at median $\sigma$; no guard on $\sigma\le0$


**V1-pleio-effects-16. vQTL sample heritability**


$\hat H^2=\mathrm{mean}_j\,1/\widehat{\mathrm{Var}}(y_{\cdot j})$


- *Variables:* —

- *Reference:* —

- *Function:* `vQTL`

- *Location:* `R/legacy_vQTL.R:251–193`

- *Notes:* Valid because $\mathrm{Var}(g)=1$; wrong for `nrow(h2)>1` (F10)


**V1-pleio-effects-17. Numeric normalisation**


if no $-1$ and some $2$: $x\leftarrow x-1$; NA $\leftarrow$ {Middle 0, Minor −1, Major 1, else 0}


- *Variables:* —

- *Reference:* —

- *Function:* `genotypes/numeric_df`

- *Location:* `R/legacy_Genotypes.R:79–62`

- *Notes:* 0/1-only data not shifted (F12); in-memory twin at `create_phenotypes.R:610–459` differs


**V1-pleio-effects-18. `maf_cutoff` filter**


keep $\min(p,1-p)\ge$ cutoff, $p=\mathrm{mean}((x+1)/2)$


- *Variables:* —

- *Reference:* —

- *Function:* `genotypes`

- *Location:* `R/legacy_Genotypes.R:163–111`

- *Notes:* Inclusive; `constraints` are strict


### Supplementary entries contributed by the Codex audit (verified by the reconciler)


**V1-pleio-effects-C1. Fully pleiotropic QTN-set construction**


$S_1=\cdots=S_T=S$, $S\sim$ sample w/o replacement from `index`; epistasis draws $w\,m_E$ physical loci and groups every $w$ consecutive loci


- *Variables:* $S$ shared causal set, $T$ traits, $w$ = `epi_interaction`, $m_E$ = `epi_QTN_num`

- *Reference:* Fernandes & Lipka 2020 *BMC Bioinformatics* 21:491 (architecture description; page unverified)

- *Function:* `qtn_pleiotropic`

- *Location:* `R/legacy_QTN_pleiotropic.R:169-201` (add), `313-345` (dom), `403-423` (var), `480-500` (epi)

- *Notes:* Effect-type sets drawn independently (not disjoint, F8); A/D share only under `same_add_dom_QTN`, A/V under `same_mv_QTN`; het loop requires $\geq$1 heterozygote anywhere in the set


**V1-pleio-effects-C2. Partial-pleiotropy locus construction**


$S_t=S_P\cup S_{T,t}$, intended $S_{T,t}\cap S_{T,u}=\varnothing$ ($t\neq u$); distinct-locus demand $\vert{}S_P\vert{}+\sum_t\vert{}S_{T,t}\vert{}$


- *Variables:* $S_P$ shared, $S_{T,t}$ trait-specific (drawn from `setdiff` pools)

- *Reference:* Fernandes & Lipka 2020 (partial pleiotropy; page unverified)

- *Function:* `qtn_partially_pleiotropic`

- *Location:* `R/legacy_QTN_partially_pleiotropic.R:165-233` (add+dom-loop), `294-360` (add only), `441-533` (dom), `616-684` (epi)

- *Notes:* Disjointness violated by the het re-sampling loop (F6/O8); preflight demands $\sum_t(\vert{}S_P\vert{}+\vert{}S_{T,t}\vert{})$ instead (O9); a single zero $\vert{}S_{T,t}\vert{}$ crashes (O4)


**V1-pleio-effects-C3. Epistatic physical-marker demand**


$N=w\,m_E$; partial: $N=w\,(\vert{}S_{P,E}\vert{}+\sum_t\vert{}S_{T,E,t}\vert{})$


- *Variables:* $w$ = `epi_interaction`

- *Reference:* none (implied by implementation)

- *Function:* `qtn_pleiotropic`, `qtn_partially_pleiotropic`

- *Location:* `R/legacy_QTN_pleiotropic.R:524-532`; `R/legacy_QTN_partially_pleiotropic.R:586-611` (also `:654`)

- *Notes:* Sampling uses the multiplier; the constraint preflight (`pleiotropic.R:105-109`, `partially.R:137-137`) does not (O9)


# Appendix A (at the audit commit, historical) — automated verification of cited locations

Every `file:line` reference above was checked against the worktree at HEAD `c511c6f`. Status meanings: **ok** = file exists, range inside file, and the named function occurs within the range or the 60 lines above it; **fn-in-file-only** = the function is defined in that file but not near the cited lines (usually a helper called from the cited lines); **fn-not-found** = the named function string does not occur in the cited file; **missing-file / beyond-eof** = the reference is wrong.


**272 of 395** cited ranges verified *ok*; the exceptions are listed below (they are informational: a helper cited by its caller, or a reference into a non-R file).

| Group | Topic | Function | Cited location | Status |
|---|---|---|---|---|
| v2-grammar/codex | Additive MAF effect scaling | `.pleio_draw` | `R/effects_pleioarch.R:143-159` | fn-in-file-only |
| v2-effects-arch | Trait-specific variance | `.pleio_draw`, `.pleio_unit_effects` | `effects_pleioarch.R:139-141` | fn-in-file-only |
| v2-effects-arch | Trait-specific variance | `.pleio_draw`, `.pleio_unit_effects` | `effects_pleioarch.R:536-536` | fn-in-file-only |
| v2-effects-arch | Allele -> genotype scaling | `.pleio_draw`; `.marker_maf_ref` | `effects_pleioarch.R:143-159` | fn-in-file-only |
| v2-effects-arch | Non-additive normaliser | `.pleio_unit_effects`; `.epi_unit_column` | `grammar_simulate_phenotype.R:452-455` | fn-not-found |
| v2-effects-arch | Non-additive normaliser | `.pleio_unit_effects`; `.epi_unit_column` | `grammar_realize.R:370-370` | fn-in-file-only |
| v2-effects-arch | Non-additive normaliser | `.pleio_unit_effects`; `.epi_unit_column` | `grammar_realize.R:408-408` | fn-in-file-only |
| v2-effects-arch | Effective per-component target | `.pleio_unit_effects` | `grammar_realize.R:510-515` | fn-in-file-only |
| v2-effects-arch | Total-correlation target | `.pleio_total_cor_check` | `grammar_realize.R:605-610` | fn-not-found |
| v2-effects-arch | Direct pair choice | `.draw_qtn_ld` | `arch_ld.R:112-122` | fn-in-file-only |
| v2-effects-arch | Direct pair choice | `.draw_qtn_ld` | `arch_ld.R:114-114` | fn-in-file-only |
| v2-effects-arch | Indirect flanks | `.draw_qtn_ld` | `arch_ld.R:123-153` | fn-in-file-only |
| v2-crossing | Heterozygote phasing | `as_population` | `R/cross_population.R:118-119` | fn-in-file-only |
| v2-crossing | Dosage from strands | `dosages` | `genome.rs:195-200` | fn-not-found |
| v2-crossing | Cross / self / DH | `mate_haplotypes` | `R/cross_mating.R:102-102` | fn-in-file-only |
| v2-crossing | Synthetic map rate | `synthetic_map` | `R/cross_map.R:157-157` | fn-in-file-only |
| v2-crossing | Synthetic map integration | `synthetic_map` | `R/cross_map.R:161-166` | fn-in-file-only |
| v2-crossing | Retention fractions (doc only) | roxygen | `R/cross_breed.R:68-72` | fn-not-found |
| v2-crossing | Pedigree key | `.mating_pedigree`, `.stable_key` | `hash.rs:10-21` | fn-not-found |
| v2-rust-core | Crossover count (count-location) | `.draw_meiosis` | `Genetics.cpp:44-44` | fn-not-found |
| v2-rust-core | Chiasma locations | `.draw_meiosis` | `Genetics.cpp:53-53` | fn-not-found |
| v2-rust-core | Chiasma locations | `.draw_meiosis` | `Genetics.cpp:86-97` | fn-not-found |
| v2-rust-core | Strand raffle | `.draw_meiosis` / `chromosome_mask` | `Genetics.cpp:55-55` | fn-not-found |
| v2-rust-core | Breakpoint rank | `breaks_at` | `Genetics.cpp:67-67` | fn-not-found |
| v2-rust-core | Gamete assembly | `recombine` | `Genetics.cpp:356-356` | fn-not-found |
| v2-rust-core | Cross | `mate_haplotypes` | `cross_mating.R:102-102` | fn-in-file-only |
| v2-rust-core | Genotype projection | `write_genotype`, `dosages` | `cross_population.R:154-158` | fn-in-file-only |
| v2-rust-core | Founder phasing | `as_population` | `cross_population.R:118-119` | fn-in-file-only |
| v2-rust-core/codex | Genetic-model post-transform (split out of Fable's single "codes" row) | `numericalize_core()` | `src/rust/src/numeric.rs:95-95` | fn-in-file-only |
| v2-rust-core/codex | Genetic-model post-transform (split out of Fable's single "codes" row) | `numericalize_core()` | `src/rust/src/numeric.rs:96-96` | fn-in-file-only |
| v2-rust-core/codex | Genetic-model post-transform (split out of Fable's single "codes" row) | `numericalize_core()` | `src/rust/src/numeric.rs:103-103` | fn-in-file-only |
| v2-rust-core/codex | Genetic-model post-transform (split out of Fable's single "codes" row) | `numericalize_core()` | `src/rust/src/numeric.rs:110-110` | fn-in-file-only |
| v2-rust-core/codex | Genetic-model post-transform (split out of Fable's single "codes" row) | `numericalize_core()` | `src/rust/src/numeric.rs:117-117` | fn-in-file-only |
| v2-selection | Selection differential | `select_ind` | `R/select_ind.R:281-281` | fn-in-file-only |
| v2-selection | Realized intensity | `select_ind`, `.select_culling` | `R/select_ind.R:282-283` | fn-in-file-only |
| v2-selection | Independent culling | `.select_culling` | `R/select_ind.R:712-726` | fn-in-file-only |
| v2-selection | Index : culling : tandem | roxygen; test-culling.R:71–89 | `R/select_ind.R:97-104` | fn-not-found |
| v2-selection/codex | Selfed / DH family BV correlation | `select_ind()` docs (`family_relationship`) | `R/select_ind.R:69-73` | fn-in-file-only |
| v2-selection/codex | Sequential culling | `.select_culling` | `R/select_ind.R:712-719` | fn-in-file-only |
| v2-selection/codex | Scheme RNG threading | all four wrappers, `.self_each`, `.intermate` | `R/select_schemes.R:89-89` | fn-in-file-only |
| v2-ocs-usefulness-marker | OCS objective | `optimum_contribution`, `.frank_wolfe` | `R/select_ocs.R:100-103` | fn-in-file-only |
| v2-ocs-usefulness-marker/codex | Ridge blend of G | `g_matrix` | `R/select_ocs.R:25-26` | fn-in-file-only |
| v2-prediction/codex | Known-variance mixed model (the model BLUP solves) | `predict_ebv` | `R/select_blup.R:90-97` | fn-in-file-only |
| v2-transcriptome | Expression model | `simulate_transcriptome` | `R/transcriptome_simulate.R:490-492` | fn-in-file-only |
| v2-transcriptome | cis score | `simulate_transcriptome` | `R/transcriptome_simulate.R:381-391` | fn-in-file-only |
| v2-transcriptome | trans factor | `simulate_transcriptome` | `R/transcriptome_simulate.R:392-392` | fn-in-file-only |
| v2-transcriptome | Additive genetic score | `simulate_transcriptome` | `R/transcriptome_simulate.R:428-457` | fn-in-file-only |
| v2-transcriptome | Epistatic blend | `simulate_transcriptome` | `R/transcriptome_simulate.R:402-451` | fn-in-file-only |
| v2-transcriptome | Residual | `simulate_transcriptome` | `R/transcriptome_simulate.R:464-474` | fn-in-file-only |
| v2-transcriptome | Realized heritability | `simulate_transcriptome`, `predict` | `R/transcriptome_simulate.R:770-773` | fn-in-file-only |
| v2-transcriptome | Effective coefficients | `simulate_transcriptome` | `R/transcriptome_simulate.R:503-513` | fn-in-file-only |
| v2-transcriptome | Effective coefficients | `simulate_transcriptome` | `R/transcriptome_simulate.R:532-532` | fn-in-file-only |
| v2-transcriptome | Effective coefficients | `simulate_transcriptome` | `R/transcriptome_simulate.R:539-539` | fn-in-file-only |
| v2-transcriptome | Effective coefficients | `simulate_transcriptome` | `R/transcriptome_simulate.R:550-550` | fn-in-file-only |
| v2-transcriptome | Genetic budget | `simulate_transcriptome` | `R/transcriptome_simulate.R:518-525` | fn-in-file-only |
| v2-transcriptome | Realized cis fraction | `simulate_transcriptome` | `R/transcriptome_simulate.R:526-528` | fn-in-file-only |
| v2-transcriptome | Mimic affine | `simulate_transcriptome` | `R/transcriptome_simulate.R:481-490` | fn-in-file-only |
| v2-io-formats | Allele frequency / MAF | `filter_geno` | `R/qc_filter_geno.R:203-210` | fn-in-file-only |
| v2-io-formats | Monomorphic / MAF filters | `filter_geno` | `R/qc_filter_geno.R:221-221` | fn-in-file-only |
| v2-io-formats | Monomorphic / MAF filters | `filter_geno` | `R/qc_filter_geno.R:224-224` | fn-in-file-only |
| v2-io-formats | Monomorphic / MAF filters | `filter_geno` | `R/qc_filter_geno.R:228-228` | fn-in-file-only |
| v2-io-formats | Heterozygosity | `filter_geno` | `R/qc_filter_geno.R:211-211` | fn-in-file-only |
| v2-io-formats | Heterozygosity | `filter_geno` | `R/qc_filter_geno.R:232-232` | fn-in-file-only |
| v2-io-formats | Heterozygosity | `filter_geno` | `R/qc_filter_geno.R:234-234` | fn-in-file-only |
| v2-io-formats | Pruning decision | `.ld_sweep` | `R/qc_filter_geno.R:498-498` | fn-in-file-only |
| v2-io-formats | Pruning decision | `.ld_sweep` | `R/qc_filter_geno.R:500-503` | fn-in-file-only |
| v2-io-formats | Window / step | `.ld_prune`, `.ld_sweep` | `R/qc_filter_geno.R:440-452` | fn-in-file-only |
| v2-io-formats | Window / step | `.ld_prune`, `.ld_sweep` | `R/qc_filter_geno.R:572-617` | fn-in-file-only |
| v2-io-formats | Numeric coding | `numericalize_core` | `src/rust/src/numeric.rs:95-117` | fn-in-file-only |
| v1-core-linkage | Epistatic value | `genetic_effect` | `legacy_genetic_effect.R:77-89` | fn-in-file-only |
| v1-core-linkage | Centering | `genetic_effect` | `legacy_genetic_effect.R:93-93` | fn-in-file-only |
| v1-core-linkage | Degree of dominance | `check_in` | `legacy_check_in.R:252-252` | fn-in-file-only |
| v1-core-linkage | Geometric series | `check_in` | `legacy_check_in.R:786-788` | fn-in-file-only |
| v1-core-linkage | Residual variance | `phenotypes` | `legacy_Phenotypes.R:216-216` | fn-in-file-only |
| v1-core-linkage | Single-trait phenotype | `phenotypes` | `legacy_Phenotypes.R:669-674` | fn-in-file-only |
| v1-core-linkage | Multi-trait residual | `phenotypes` | `legacy_Phenotypes.R:238-245` | fn-in-file-only |
| v1-core-linkage | Null trait (h2=0) | `phenotypes` | `legacy_Phenotypes.R:533-534` | fn-in-file-only |
| v1-core-linkage | Realized heritability | `phenotypes` | `legacy_Phenotypes.R:248-251` | fn-in-file-only |
| v1-core-linkage | Per-QTN PVE | `phenotypes` | `legacy_Phenotypes.R:260-268` | fn-in-file-only |
| v1-core-linkage | Residual seed | `phenotypes` | `legacy_Phenotypes.R:234-234` | fn-in-file-only |
| v1-core-linkage | LD-QTN seed | `qtn_linkage` | `legacy_QTN_linkage.R:110-110` | fn-in-file-only |
| v1-core-linkage | LD measure | `qtn_linkage` | `legacy_QTN_linkage.R:150-151` | fn-in-file-only |
| v1-core-linkage | Indirect selection | `qtn_linkage` | `legacy_QTN_linkage.R:125-212` | fn-in-file-only |
| v1-core-linkage | Direct selection | `qtn_linkage` | `legacy_QTN_linkage.R:1229-1308` | fn-in-file-only |
| v1-core-linkage | Direct selection | `qtn_linkage` | `legacy_QTN_linkage.R:978-1062` | fn-in-file-only |
| v1-core-linkage | MAF | `constraint`, `qtn_linkage`, `qtn_from_user` | `QTN_linkage.R:311-314` | fn-in-file-only |
| v1-core-linkage | MAF | `constraint`, `qtn_linkage`, `qtn_from_user` | `qtn_from_user.R:237-240` | fn-in-file-only |
| v1-core-linkage/codex | Reported residual "sample correlation" | `phenotypes()` | `legacy_Phenotypes.R:483-491` | fn-in-file-only |
| v1-core-linkage/codex | Reported residual "sample correlation" | `phenotypes()` | `legacy_Phenotypes.R:484-484` | fn-in-file-only |
| v1-core-linkage/codex | Reported residual "sample correlation" | `phenotypes()` | `legacy_Phenotypes.R:1336-1336` | fn-in-file-only |
| v1-core-linkage/codex | Reported residual "sample correlation" | `phenotypes()` | `legacy_Phenotypes.R:488-488` | fn-in-file-only |
| v1-core-linkage/codex | Reported residual "sample correlation" | `phenotypes()` | `legacy_Phenotypes.R:1340-1340` | fn-in-file-only |
| v1-core-linkage/codex | QTN anchor sampling without replacement | `qtn_linkage()` | `legacy_QTN_linkage.R:109-113` | fn-in-file-only |
| v1-core-linkage/codex | QTN anchor sampling without replacement | `qtn_linkage()` | `legacy_QTN_linkage.R:109-109` | fn-in-file-only |
| v1-core-linkage/codex | QTN anchor sampling without replacement | `qtn_linkage()` | `legacy_QTN_linkage.R:113-113` | fn-in-file-only |
| v1-core-linkage/codex | User marker resolution and order | `qtn_from_user()` | `legacy_qtn_from_user.R:1358-1371` | fn-in-file-only |
| v1-core-linkage/codex | User marker resolution and order | `qtn_from_user()` | `legacy_qtn_from_user.R:1447-1460` | fn-in-file-only |
| v1-core-linkage/codex | User marker resolution and order | `qtn_from_user()` | `legacy_qtn_from_user.R:1517-1530` | fn-in-file-only |
| v1-core-linkage/codex | User marker resolution and order | `qtn_from_user()` | `legacy_qtn_from_user.R:1589-1602` | fn-in-file-only |
| v1-core-linkage/codex | User marker resolution and order | `qtn_from_user()` | `legacy_qtn_from_user.R:1361-1361` | fn-in-file-only |
| v1-core-linkage/codex | Epistatic interaction diagnostic expansion | `qtn_from_user()` | `legacy_qtn_from_user.R:1531-1554` | fn-in-file-only |
| v1-core-linkage/codex | Epistatic interaction diagnostic expansion | `qtn_from_user()` | `legacy_qtn_from_user.R:1532-1532` | fn-in-file-only |
| v1-core-linkage/codex | Epistatic interaction diagnostic expansion | `qtn_from_user()` | `legacy_qtn_from_user.R:1541-1541` | fn-in-file-only |
| v1-core-linkage/codex | Epistatic interaction diagnostic expansion | `qtn_from_user()` | `legacy_qtn_from_user.R:1544-1544` | fn-in-file-only |
| v1-core-linkage/codex | *(reconciler-added)* Direct-LD acceptance sentinel | `qtn_linkage()` | `legacy_QTN_linkage.R:1231-1237` | fn-in-file-only |
| v1-core-linkage/codex | *(reconciler-added)* Direct-LD acceptance sentinel | `qtn_linkage()` | `legacy_QTN_linkage.R:123-125` | fn-in-file-only |
| v1-core-linkage/codex | *(reconciler-added)* Direct-LD acceptance sentinel | `qtn_linkage()` | `legacy_QTN_linkage.R:976-983` | fn-in-file-only |
| v1-core-linkage/codex | *(reconciler-added)* Direct-LD acceptance sentinel | `qtn_linkage()` | `legacy_QTN_linkage.R:1467-1473` | fn-in-file-only |
| v1-pleio-effects | Geometric effect series | `check_in` | `R/legacy_check_in.R:776-837` | fn-in-file-only |
| v1-pleio-effects | Dominance genetic value | `genetic_effect` | `R/legacy_check_in.R:252-252` | fn-not-found |
| v1-pleio-effects | Epistatic genetic value | `genetic_effect` | `R/legacy_genetic_effect.R:78-91` | fn-in-file-only |
| v1-pleio-effects | Total genetic value | `genetic_effect` | `R/legacy_genetic_effect.R:92-96` | fn-in-file-only |
| v1-pleio-effects | MAF of a QTN | `qtn_pleiotropic`, `qtn_partially_pleiotropic`, `constraint` | `R/legacy_QTN_pleiotropic.R:161-165` | fn-in-file-only |
| v1-pleio-effects | QTN sampling seeds (pleiotropic) | `qtn_pleiotropic` | `R/legacy_QTN_pleiotropic.R:130-130` | fn-in-file-only |
| v1-pleio-effects | Whitening | `base_line_multi_traits` | `R/legacy_Base_line_multi_traits.R:103-109` | fn-in-file-only |
| v1-pleio-effects | Colouring & rescale | `base_line_multi_traits` | `multi_traits.R:110-120` | fn-in-file-only |
| v1-pleio-effects | Residual variance (context) | `phenotypes` | `R/legacy_Phenotypes.R:216-216` | fn-in-file-only |
| v1-pleio-effects | vQTL sample heritability | `vQTL` | `R/legacy_vQTL.R:184-193` | fn-in-file-only |
| v1-pleio-effects/codex | Fully pleiotropic QTN-set construction | `qtn_pleiotropic` | `R/legacy_QTN_pleiotropic.R:125-157` | fn-in-file-only |
| v1-pleio-effects/codex | Partial-pleiotropy locus construction | `qtn_partially_pleiotropic` | `R/legacy_QTN_partially_pleiotropic.R:120-190` | fn-in-file-only |
| v1-pleio-effects/codex | Epistatic physical-marker demand | `qtn_pleiotropic`, `qtn_partially_pleiotropic` | `R/legacy_QTN_pleiotropic.R:480-488` | fn-in-file-only |
| v1-pleio-effects/codex | Epistatic physical-marker demand | `qtn_pleiotropic`, `qtn_partially_pleiotropic` | `R/legacy_QTN_partially_pleiotropic.R:616-641` | fn-in-file-only |
| v1-pleio-effects/codex | Epistatic physical-marker demand | `qtn_pleiotropic`, `qtn_partially_pleiotropic` | `R/legacy_QTN_partially_pleiotropic.R:654-654` | fn-in-file-only |

# Appendix B — quick index (topic, function, location)

| Ver | Group | Topic | Function | Location |
|---|---|---|---|---|
| V2 | v2-grammar | Layer scaling to prop | `.genetic_matrix` | R/grammar_realize.R:139-161 (146) |
| V2 | v2-grammar | Additive component | `.component_raw` | R/grammar_realize.R:498-509 (500), centred 432 |
| V2 | v2-grammar | Orthogonal genotypic value | `.component_raw` | R/grammar_realize.R:501-507 (506) |
| V2 | v2-grammar | Dominance component | `.component_raw` | R/grammar_realize.R:510 |
| V2 | v2-grammar | Epistasis unit column | `.epi_unit_column` | R/grammar_realize.R:579-589 (586, 588) |
| V2 | v2-grammar | Epistasis component | `.component_raw` | R/grammar_realize.R:511-529 (524) |
| V2 | v2-grammar | Geometric effect series | `.effect_series` | R/effects_series.R:23-72 (64) |
| V2 | v2-grammar | Repulsion phase | `.apply_phase` | R/grammar_layers.R:1482-1490 (1487) |
| V2 | v2-grammar | Residual budget | `.realize_phenotype` | R/grammar_realize.R:55-59 (59) |
| V2 | v2-grammar | Residual draw | `.draw_residual` | R/effects_series.R:85-102 (99) |
| V2 | v2-grammar | Phenotype | `.realize_phenotype` | R/grammar_realize.R:96 |
| V2 | v2-grammar | Requested budget identity | `.resolve_prop`, `.add_layer`, `.check_h2_complete` | R/grammar_layers.R:1085-1098, 1102-1115; R/grammar_realize.R:1153-1166 |
| V2 | v2-grammar | One-call split | `.build_one_call` | R/grammar_simulate_phenotype.R:767-780 (769) |
| V2 | v2-grammar | Realized H² | `.realized_h2` | R/grammar_realize.R:1078-1106 (716-720) |
| V2 | v2-grammar | Average effect of substitution | `.avg_effect` | R/grammar_realize.R:182-184 (183) |
| V2 | v2-grammar | Breeding value | `.breeding_value_matrix` | R/grammar_realize.R:282-303 (227-235) |
| V2 | v2-grammar | Orthogonal split | `.orthogonal_var_split`, `.variance_budget` | R/grammar_realize.R:778-794 (646-650), 536-551 |
| V2 | v2-grammar | HWE theory (check only) | — (audit check) | — |
| V2 | v2-grammar | vQTL loading | `.apply_vqtl` | R/grammar_realize.R:605-618 (616) |
| V2 | v2-grammar | vQTL residual | `.apply_vqtl` | R/grammar_realize.R:627-641 |
| V2 | v2-grammar | Transcriptome score | `.tx_raw`, `.transcriptome_matrix` | R/grammar_realize.R:454-476 (470, 473, 475), 308-338 |
| V2 | v2-grammar | Mediation split | `.mediation_budget` | R/grammar_realize.R:738-750 (604-606) |
| V2 | v2-grammar | Complex combination | `complex_phenotypes` | R/grammar_complex.R:122-134 (123, 127), 96-107 (99) |
| V2 | v2-grammar | MAF | `.marker_maf_ref` | R/grammar_simulate_phenotype.R:1022-1052 (1040, 1044) |
| V2 | v2-grammar | Sub-seed | `.layer_seed` | R/grammar_simulate_phenotype.R:1079-1094 (1083, 1092) |
| V2 | v2-grammar | Per-QTN variance | `.qtn_var` | R/io_write.R:1288-1312 (217-218) |
| V2 | v2-grammar | Fixed-scale phenotype | `phenotype_value` | R/cross_population.R:1004-1104 (1061) |
| V2 | v2-grammar | Genetic correlation (plot) | `.plot_cor` | R/grammar_plot.R:144-162 (150) |
| V2 | v2-grammar (codex) | PleioArch shared covariance | `.pleio_draw`, `.pleio_nonadditive_draw` | R/effects_pleioarch.R:65-93 (sigma at 78-80), 345-358 (356-357) |
| V2 | v2-grammar (codex) | PleioArch effect allocation (non-additive) | `.pleio_unit_effects`, `.draw_mvnorm` | R/effects_pleioarch.R:583-690 (535-537), 780-800 |
| V2 | v2-grammar (codex) | PleioArch attainability | `.check_pleio_feasible` | R/effects_pleioarch.R:870-917 |
| V2 | v2-grammar (codex) | Total pleiotropic correlation target | `.pleio_total_cor_check` | R/effects_pleioarch.R:724-788 (roxygen 543-574) |
| V2 | v2-grammar (codex) | Additive MAF effect scaling | `.pleio_draw` | R/effects_pleioarch.R:201-217 (144-149) |
| V2 | v2-grammar (codex) | LD window (two-trait linked distinct loci) | `.draw_qtn_ld` | R/arch_ld.R:62-225 (window 69-95, sampling 105-175) |
| V2 | v2-effects-arch | Pleiotropic covariance | `.pleio_draw`, `.pleio_nonadditive_draw` | `effects_pleioarch.R:90-91`, `:504-505` |
| V2 | v2-effects-arch | Trait-specific variance | `.pleio_draw`, `.pleio_unit_effects` | `:197-199`, `:685` |
| V2 | v2-effects-arch | Major/minor split | `.pleio_draw`, `.draw_mvnorm` | `:195-196`, `:934` |
| V2 | v2-effects-arch | MVN draw | `.draw_mvnorm` | `:929-939` |
| V2 | v2-effects-arch | Univariate draw | `.draw_univariate` | `:944-949` |
| V2 | v2-effects-arch | Shared-unit count | `.pleio_partition` | `:291-292, 320` |
| V2 | v2-effects-arch | Attainability (2 traits) | `.check_pleio_feasible` | `:876-890` |
| V2 | v2-effects-arch | Attainability (n traits) | `.check_pleio_feasible` | `:898-915` |
| V2 | v2-effects-arch | Allele -> genotype scaling | `.pleio_draw`; `.marker_maf_ref` | `:201-217`; `grammar_simulate_phenotype.R:1040-1044` |
| V2 | v2-effects-arch | Non-additive normaliser | `.pleio_unit_effects`; `.epi_unit_column` | `:601-604, 686`; `grammar_realize.R:579-589`; dominance design `:518` vs `grammar_realize.R:510` |
| V2 | v2-effects-arch | Effective per-component target | `.pleio_unit_effects` | `:659-664` |
| V2 | v2-effects-arch | Total-correlation target | `.pleio_total_cor_check` | `:754-759` |
| V2 | v2-effects-arch | Layer rescale (context) | `.genetic_matrix` | `grammar_realize.R:144-146` |
| V2 | v2-effects-arch | Few-unit attenuation | doc claim | `:26-28`; `grammar_simulate_phenotype.R:318-320` |
| V2 | v2-effects-arch | Complete-LD ensemble mean | doc claim | `:28-30`; DECISIONS.md 023 Scope |
| V2 | v2-effects-arch | Geometric effect series | `.effect_series` | `effects_series.R:55, 64` |
| V2 | v2-effects-arch | Repulsion phase | `.apply_phase` | `grammar_layers.R:1482-1490` |
| V2 | v2-effects-arch | Residual draw | `.draw_residual` | `effects_series.R:85-102` |
| V2 | v2-effects-arch | LD measure | `window_partners`, `r2_pair` | `arch_ld.R:111, 124-126` |
| V2 | v2-effects-arch | Direct pair choice | `.draw_qtn_ld` | `arch_ld.R:149-160` (`:151-152`); `partner` `:73-74` |
| V2 | v2-effects-arch | Indirect flanks | `.draw_qtn_ld` | `arch_ld.R:161-197` |
| V2 | v2-effects-arch | Layer sub-seed | `.layer_seed` | `grammar_simulate_phenotype.R:1083-1093` |
| V2 | v2-effects-arch | Gabriel blocks tag | `.gabriel_blocks` | `qc_ld_methods.R:24-29, 38-43` |
| V2 | v2-effects-arch (codex) | Zero-variance correlation guard | `.pleio_check_zero_var()` | `R/effects_pleioarch.R:246-275` |
| V2 | v2-effects-arch (codex) | Shared/specific unit count (clamped) | `.pleio_partition()` | `R/effects_pleioarch.R:288-329` (clamp `:291-292`) |
| V2 | v2-effects-arch (codex) | One shared unit | `.pleio_single_unit_consequence()` | `R/effects_pleioarch.R:412-441` |
| V2 | v2-effects-arch (codex) | Correlation input expansion | `.pleio_cor_matrix()` | `R/effects_pleioarch.R:957-986` |
| V2 | v2-effects-arch (codex) | Pleiotropic-share argument mapping | `.pleio_pi_vector()` | `R/effects_pleioarch.R:829-860` |
| V2 | v2-effects-arch (codex) | Architecture-specific QTN sampling | `.draw_qtn()` | `R/arch_independent.R:21-45` |
| V2 | v2-effects-arch (codex) | Distinct-chromosome allocation | `.draw_qtn_distinct_chr()` | `R/arch_independent.R:56-75` |
| V2 | v2-effects-arch (codex) | Epistatic-set sampling | `.draw_qtn_pairs()` | `R/arch_independent.R:81-103` |
| V2 | v2-effects-arch (codex) | Candidate-marker rule | `.candidate_markers()` | `R/arch_independent.R:114-134` |
| V2 | v2-effects-arch (codex) | Layer occurrence index | `.type_occurrence()` | `R/arch_independent.R:139-141` |
| V2 | v2-crossing | Heterozygote phasing | `as_population` | `R/cross_population.R:176-177` |
| V2 | v2-crossing | Dosage from strands | `dosages` | `R/cross_population.R:623`; `genome.rs:273-278` |
| V2 | v2-crossing | Crossover count | `.draw_meiosis` | `R/cross_mating.R:53-56` (dispatch), Poisson branch `:77-78`; gamma model `:215-275`; default `.INTERFERENCE_DEFAULT` `:117` (DECISION-047); cM$\rightarrow$M `:339,341` |
| V2 | v2-crossing | Chiasma positions | `.draw_meiosis` | `R/cross_mating.R:79-87` |
| V2 | v2-crossing | Strand choice | `.draw_meiosis` | `R/cross_mating.R:90` |
| V2 | v2-crossing | Ancestry mask | `chromosome_mask`, `breaks_at` | `src/rust/src/meiosis.rs:49-51, 63-72` |
| V2 | v2-crossing | Gamete | `recombine` | `meiosis.rs:224-237` |
| V2 | v2-crossing | Cross / self / DH | `mate_haplotypes` | `meiosis.rs:350-371`; events per progeny `R/cross_mating.R:345` |
| V2 | v2-crossing | Recombination fraction (implied, verified) | — (property of the Poisson process) | verified e4.R/e5.R; `tests/testthat/test-cross.R:179-202` (file runs under `simplePHENOTYPES.interference = "poisson"`, `:13`) |
| V2 | v2-crossing | Heterozygosity under selfing | `selfcross` (doc) | `R/cross_mating.R:627-629`; verified e1.R E3 |
| V2 | v2-crossing | Synthetic map rate | `synthetic_map` | `R/cross_map.R:157` |
| V2 | v2-crossing | Synthetic map integration | `synthetic_map` | `R/cross_map.R:168-173`, `:133` |
| V2 | v2-crossing | Additive value (fixed scale) | `additive_value` | `R/cross_population.R:746` |
| V2 | v2-crossing | Genotypic value | `genotypic_value` | `R/cross_population.R:821` |
| V2 | v2-crossing | Residual from h² | `phenotype_value` | `R/cross_population.R:1061` |
| V2 | v2-crossing | Breed composition | `breed_composition` | `R/cross_breed.R:37-43` |
| V2 | v2-crossing | Realized heterosis | `heterosis` | `R/cross_breed.R:135-139, 152` |
| V2 | v2-crossing | Expected F1 heterosis | `heterosis` + `.expected_cross_means` | `R/cross_breed.R:140-150`; `R/select_combining.R:244-255` |
| V2 | v2-crossing | Per-locus cross mean | `.expected_cross_means` | `R/select_combining.R:249-254` |
| V2 | v2-crossing | Retention fractions (doc only) | roxygen | `R/cross_breed.R:68-72`; `docs/SPEC-block3b.md` §7 |
| V2 | v2-crossing | Rotation sire sequence | `crossbreed` | `R/cross_breed.R:302-305` |
| V2 | v2-crossing | Pedigree key | `.mating_pedigree`, `.stable_key` | `R/cross_pedigree.R:25-37` (`.stable_key`), `:42-59` (`.key_part`), `:205-207` (mating key); `hash.rs:12-23` |
| V2 | v2-crossing | Generation | `.mating_pedigree` | `R/cross_pedigree.R:208-209` |
| V2 | v2-crossing | A-matrix (other group, consistency only) | `a_matrix` | `R/select_blup.R:69-90` |
| V2 | v2-crossing (codex) | Average-effect vs genotypic-value distinction (doc) | roxygen of `genotypic_value` (implemented in `select_ind(on = "bv")`, other group) | `R/cross_population.R:774-780` (verified: α text at 776-777) |
| V2 | v2-crossing (codex) | Full-sib family identity | `families` | `R/cross_pedigree.R:307-336` (verified: `families <-` at 307, closes at 336) |
| V2 | v2-crossing (codex) | Mating designs (counts) | `mating_design` | `R/cross_mate.R:360-486` (verified: `mating_design <-` at 229; file is 350 lines) |
| V2 | v2-crossing (codex) | Within-family relationship validation (test) | test helper `.icc` | `tests/testthat/test-family-relationship.R:1-58` (verified: 58 lines, `.icc` at 26, targets at 52-57) |
| V2 | v2-rust-core | Crossover count (count-location) | `.draw_meiosis` | R `R/cross_mating.R:77-78`; isqg `Genetics.cpp:44,84`; Rust: none (input `counts`) |
| V2 | v2-rust-core | Chiasma locations | `.draw_meiosis` | R `:79-87`; isqg `Genetics.cpp:86-97` |
| V2 | v2-rust-core | Strand raffle | `.draw_meiosis` / `chromosome_mask` | R `:90`; Rust `meiosis.rs:68-70`; isqg `:74-75` |
| V2 | v2-rust-core | Breakpoint rank | `breaks_at` | Rust `meiosis.rs:49-51`; isqg `Genetics.cpp:67` |
| V2 | v2-rust-core | Ancestry mask | `chromosome_mask`, `Bits::toggle_from`, `flip_all` | Rust `meiosis.rs:63-72`, `genome.rs:69-87`; isqg `:58-70` (XOR loop), `:74-75` (flip) |
| V2 | v2-rust-core | Gamete assembly | `recombine` | Rust `meiosis.rs:224-237`; isqg `Genetics.cpp:356` |
| V2 | v2-rust-core | Cross | `mate_haplotypes` | Rust `meiosis.rs:359-369`; R `.mate` `cross_mating.R:345,453,460-467` |
| V2 | v2-rust-core | Self | `selfcross` $\rightarrow$ `.mate(parent, parent)` | R `cross_mating.R:647-652`; Rust same path |
| V2 | v2-rust-core | Doubled haploid | `mate_haplotypes` (Dh) | Rust `meiosis.rs:263-269,351-358`; R `cross_mating.R:686-691` |
| V2 | v2-rust-core | Genotype projection | `write_genotype`, `dosages` | Rust `genome.rs:273-284`; R `cross_population.R:623` (`cis + trans - 1`) |
| V2 | v2-rust-core | Founder phasing | `as_population` | R `cross_population.R:176-177` |
| V2 | v2-rust-core | Map units | `.mate` | R `cross_mating.R:339,341` |
| V2 | v2-rust-core | Haldane map function (consequence, tested only) | — | test `tests/testthat/test-cross.R:179-202`; auditor `T2–T4` |
| V2 | v2-rust-core | Numericalization — orientation | `compute_flip` | R `io_detect_format.R:362-364`; reference mode `:354-359` |
| V2 | v2-rust-core | Numericalization — codes | `numericalize_core` | Rust `numeric.rs:76-88,135-180` |
| V2 | v2-rust-core | Pedigree key hash | `fnv1a_128`, `.stable_key` | Rust `hash.rs:12-23`; R `cross_pedigree.R:22-59` |
| V2 | v2-rust-core (codex) | Public ingestion dispatch | `as_numeric()` | `R/io_as_numeric.R:192-207` (verified: `if (is.character(x) && is.null(dim(x)))` at :198, `format_conversion(...)` at :206) |
| V2 | v2-rust-core (codex) | Genetic-model post-transform (split out of Fable's single "codes" row) | `numericalize_core()` | `src/rust/src/numeric.rs:149-180` (verified: `match model` at :156, `"Dom"` :157, `"Left"` :164, `"Right"` :171, default :178) |
| V2 | v2-selection | Truncation count from prop | `.resolve_keep` | `R/select_ind.R:574` |
| V2 | v2-selection | Selection intensity (large N) | `.resolve_keep` | `R/select_ind.R:585–587` |
| V2 | v2-selection | Selection differential | `select_ind` | `R/select_ind.R:460, 482` |
| V2 | v2-selection | Realized intensity | `select_ind`, `.select_culling` | `R/select_ind.R:461–478, 1089–1090` |
| V2 | v2-selection | Response (documented, tested here) | roxygen only | `R/select_ind.R:16–21` |
| V2 | v2-selection | Average effect | `.avg_effect` | `R/grammar_realize.R:182–184` |
| V2 | v2-selection | Breeding value | `.breeding_value_matrix` | `R/grammar_realize.R:236–305` (loop 282–303) |
| V2 | v2-selection | Smith–Hazel index | `.index_score`, `.index_weights` | `R/select_ind.R:769–782; 799–842` |
| V2 | v2-selection | QGSI | `.quadratic_index_score` | `R/select_ind.R:737, 744–747` |
| V2 | v2-selection | Lush combined index | `.combined_score` | `R/select_ind.R:867, 872–877, 889–902` |
| V2 | v2-selection | Within-family allocation | `.sel_within_family` | `R/select_ind.R:939–949` |
| V2 | v2-selection | Among-family | `.sel_among_family` | `R/select_ind.R:968–978` |
| V2 | v2-selection | Independent culling | `.select_culling` | `R/select_ind.R:1067–1081` |
| V2 | v2-selection | Index : culling : tandem | roxygen; test-culling.R:71–89 | `R/select_ind.R:157–164` |
| V2 | v2-selection | Bulk pool | `bulk` | `R/select_schemes.R:156–169`; `.bulk_counts` `:470–472` |
| V2 | v2-selection | Pedigree family size | `pedigree` | `R/select_schemes.R:304–307` |
| V2 | v2-selection | Tandem schedule | `pedigree`, `recurrent_selection` | `R/select_schemes.R:291, 398` |
| V2 | v2-selection | Intermating | `.intermate` | `R/select_schemes.R:542–545` |
| V2 | v2-selection (codex) | General index response (context, not implemented) | roxygen only | `R/select_ind.R:16-21` |
| V2 | v2-selection (codex) | Singular Smith–Hazel extension | `.index_weights` | `R/select_ind.R:799-830` |
| V2 | v2-selection (codex) | Selfed / DH family BV correlation | `select_ind()` docs (`family_relationship`) | `R/select_ind.R:94-98` |
| V2 | v2-selection (codex) | Sequential culling | `.select_culling` | `R/select_ind.R:1067-1074` |
| V2 | v2-selection (codex) | Population pooling | `c.Population` | `R/select_schemes.R:25-68` |
| V2 | v2-selection (codex) | Single seed descent | `single_seed_descent` | `R/select_schemes.R:102-118` |
| V2 | v2-selection (codex) | Recurrent-selection cycle | `recurrent_selection`, `.intermate` | `R/select_schemes.R:387-410; 354-370` |
| V2 | v2-selection (codex) | Scheme RNG threading | all four wrappers, `.self_each`, `.intermate` | `R/select_schemes.R:111, 155, 283, 384; 521; 544` |
| V2 | v2-ocs-usefulness-marker | Genomic relationship | `g_matrix` | `R/select_ocs.R:71–92` |
| V2 | v2-ocs-usefulness-marker | Genomic inbreeding | `g_matrix` (doc) | `R/select_ocs.R:16, 42–43` |
| V2 | v2-ocs-usefulness-marker | OCS objective | `optimum_contribution`, `.frank_wolfe` | `R/select_ocs.R:105–108, 601–624` |
| V2 | v2-ocs-usefulness-marker | Group coancestry | `optimum_contribution`, `.tune_lambda` | `R/select_ocs.R:319, 666` |
| V2 | v2-ocs-usefulness-marker | FW gradient / duality gap | `.frank_wolfe` | `R/select_ocs.R:601–607` |
| V2 | v2-ocs-usefulness-marker | Away-step and line search | `.frank_wolfe` | `R/select_ocs.R:609–625` |
| V2 | v2-ocs-usefulness-marker | Penalty tuning | `.tune_lambda` | `R/select_ocs.R:652–718` |
| V2 | v2-ocs-usefulness-marker | Parent sampling | `sample_parents` | `R/select_ocs.R:429–441`; `.allocate_slots` `:473–494` |
| V2 | v2-ocs-usefulness-marker | Usefulness | `cross_usefulness` | `R/select_usefulness.R:142–144` |
| V2 | v2-ocs-usefulness-marker | Selection intensity | `.intensity_from_p` | `R/select_usefulness.R:161–163` |
| V2 | v2-ocs-usefulness-marker | Fixed additive score | `.additive_model`, `.additive_gv` | `R/select_usefulness.R:204–209, 228` |
| V2 | v2-ocs-usefulness-marker | Average effect | `.avg_effect` | `R/grammar_realize.R:182–184` |
| V2 | v2-ocs-usefulness-marker | MAS feasibility | `marker_select` | `R/select_marker.R:116–119` |
| V2 | v2-ocs-usefulness-marker | MARS/oracle index | `additive_value` (documented in `marker_select`) | `R/cross_population.R:735–747`; doc `R/select_marker.R:20–33` |
| V2 | v2-ocs-usefulness-marker | Recurrent-allele score | `.mabc_founders` | `R/select_mabc.R:311` |
| V2 | v2-ocs-usefulness-marker | Recovery | `.mabc_recovery` | `R/select_mabc.R:464` |
| V2 | v2-ocs-usefulness-marker | Expected recovery | doc only | `R/select_mabc.R:25–31, 227–229` |
| V2 | v2-ocs-usefulness-marker | Interval weights | `.mabc_weights` | `R/select_mabc.R:520–529` |
| V2 | v2-ocs-usefulness-marker | Foreground / recombinant / background staging | `mabc_select` | `R/select_mabc.R:170–177, 194` |
| V2 | v2-ocs-usefulness-marker (codex) | Ridge blend of G | `g_matrix` | `R/select_ocs.R:25-26, 44-50, 93-95` |
| V2 | v2-ocs-usefulness-marker (codex) | MABC foreground feasibility | `mabc_select` | `R/select_mabc.R:168-174` |
| V2 | v2-prediction | Tabular relationship (off-diagonal) | `a_matrix` | `R/select_blup.R:73-79` |
| V2 | v2-prediction | Inbreeding / diagonal | `a_matrix` | `R/select_blup.R:80-88` |
| V2 | v2-prediction | Variance ratio | `predict_ebv`, `.blup_variances` | `R/select_blup.R:365`, `:676-681` |
| V2 | v2-prediction | GLS mean | `predict_ebv` | `R/select_blup.R:376-386` |
| V2 | v2-prediction | BLUP of u | `predict_ebv` | `R/select_blup.R:387-388` |
| V2 | v2-prediction | PEV / reliability | `predict_ebv` | `R/select_blup.R:395-398` |
| V2 | v2-prediction | GBLUP marker back-solve | `.gblup_marker_effects` | `R/select_blup.R:695-703` |
| V2 | v2-prediction | Accuracy / bias | `prediction_accuracy` | `R/select_blup.R:752-761` |
| V2 | v2-prediction | Expected cross mean (one locus) | `.expected_cross_means` | `R/select_combining.R:244-255` |
| V2 | v2-prediction | Tester-referenced average effect | `combining_ability` (doc) | `R/select_combining.R:30-37` |
| V2 | v2-prediction | Topcross / factorial GCA, SCA | `.decompose_ca` | `R/select_combining.R:321-324` |
| V2 | v2-prediction | Diallel GCA (method 4) | `.decompose_ca` | `R/select_combining.R:312-316` |
| V2 | v2-prediction | Diallel SCA (method 4) | `.decompose_ca` | `R/select_combining.R:317-318` |
| V2 | v2-prediction | Broad-sense residual for simulated crosses | `phenotype_value` via `.simulate_cross_means` | `R/cross_population.R:1061`; `R/select_combining.R:290-293` |
| V2 | v2-prediction | Realized-scale template effects | `.layer_scaled_effects` $\leftarrow$ `template_effects` | `R/grammar_realize.R:364`; `R/select_combining.R:400-406` |
| V2 | v2-prediction | Progeny-test expected mean | `progeny_test` (doc) | `R/select_progeny.R:17-35` |
| V2 | v2-prediction | Progeny-test accuracy | `progeny_test` (doc) | `R/select_progeny.R:37-53` |
| V2 | v2-prediction | Mate sampling | `progeny_test` | `R/select_progeny.R:121-135` |
| V2 | v2-prediction | Family mean | `progeny_test` | `R/select_progeny.R:149-150` |
| V2 | v2-prediction (codex) | Known-variance mixed model (the model BLUP solves) | `predict_ebv` | `R/select_blup.R:174-181` (roxygen), `:300-321` (inputs) |
| V2 | v2-prediction (codex) | Supplied-covariance validation (correlation-scale PSD test) | `.check_relationship` | `R/select_blup.R:587-650` (scaling `:630-632`, eigen `:640`, tolerance `:644-645`) |
| V2 | v2-prediction (codex) | VanRaden method-1 genomic relationship (consumed by GBLUP) | `predict_ebv` $\rightarrow$ `g_matrix` | `R/select_blup.R:337-343` (Codex wrote 168-173); `R/select_ocs.R:44-98` (Codex wrote 44-92; core at `:83, 90-92`) |
| V2 | v2-prediction (codex) | Progeny-test RNG draw order (reproducibility contract) | `progeny_test` | `R/select_progeny.R:113-114` (seed), `:132` (mates), `:138` (meioses), `:144-145` (residual) |
| V2 | v2-transcriptome | Expression model | `simulate_transcriptome` | `R/transcriptome_simulate.R:616-618` |
| V2 | v2-transcriptome | Reference-centered dosage | `simulate_transcriptome` | `:308-310`; predict `:845-847` |
| V2 | v2-transcriptome | cis score | `simulate_transcriptome` | `:501-511` |
| V2 | v2-transcriptome | trans factor | `simulate_transcriptome` | `:451-471`, `:512` |
| V2 | v2-transcriptome | Standardization | `z1` | `:492-496` |
| V2 | v2-transcriptome | Additive genetic score | `simulate_transcriptome` | `:548-578` |
| V2 | v2-transcriptome | Epistatic blend | `simulate_transcriptome` | `:522-572` |
| V2 | v2-transcriptome | Residual | `simulate_transcriptome` | `:585-595` |
| V2 | v2-transcriptome | Realized heritability | `simulate_transcriptome`, `predict` | `:619-628`; `:943-948` |
| V2 | v2-transcriptome | Effective coefficients | `simulate_transcriptome` | `:636-646`, `:666`, `:673`, `:684` |
| V2 | v2-transcriptome | Genetic budget | `simulate_transcriptome` | `:651-658` |
| V2 | v2-transcriptome | Realized cis fraction | `simulate_transcriptome` | `:659-661`; predict `:950` |
| V2 | v2-transcriptome | Cross-population genetic value | `predict.transcriptome_sim` | `:829-848`, `:860-886` |
| V2 | v2-transcriptome | Mimic affine | `simulate_transcriptome` | `:602-616` |
| V2 | v2-transcriptome | GRM | `.tx_grm` | `R/transcriptome_mimic.R:24-30` |
| V2 | v2-transcriptome | REML objective | `.greml_h2` | `R/transcriptome_mimic.R:76-84` |
| V2 | v2-transcriptome | GREML heritability | `.greml_h2` | `R/transcriptome_mimic.R:94-97` |
| V2 | v2-transcriptome | Factor count | `.tx_estimate_factors` | `R/transcriptome_mimic.R:111-119` |
| V2 | v2-transcriptome | kappa proxy | `.tx_estimate_kappa` | `R/transcriptome_mimic.R:149-172` (roxygen `:120-148`); called at `R/transcriptome_simulate.R:395` |
| V2 | v2-transcriptome | NB observation | `observe_counts` | `R/transcriptome_counts.R:75-80`, `:98-102` |
| V2 | v2-transcriptome | Layer score | `.tx_raw`, `.transcriptome_matrix` | `R/grammar_realize.R:454-476`, `:414-440` |
| V2 | v2-transcriptome | Mediation split | `.tx_raw(which="genetic")`, `.mediation_budget`, `mediation_split` | `R/grammar_realize.R:458-476`, `:730-761`; `R/io_write.R:1267-1270` |
| V2 | v2-transcriptome | Realized H2 | `.genetic_value_matrix`, `.realized_h2` | `R/grammar_realize.R:394-400`, `:1078-1106` |
| V2 | v2-transcriptome | Per-gene share | `.tx_qtn_var` | `R/io_write.R:1328-1345` |
| V2 | v2-transcriptome | Layer sub-seed | `.layer_seed`, `transcriptome` | `R/grammar_simulate_phenotype.R:1079-1092`; `R/transcriptome_layer.R:257-276` |
| V2 | v2-transcriptome (codex) | Single-component LMM (model statement) | `.greml_h2` | `R/transcriptome_mimic.R:32-53` (roxygen), fit `:58-97` |
| V2 | v2-transcriptome (codex) | Calibration metrics | benchmark 01 | `benchmarks/01_h2_calibration.R:61-62` (h2), `:93-94` (cis) |
| V2 | v2-transcriptome (codex) | Marginal eQTL scan | `scan_one` | `benchmarks/02_eqtl_recovery.R:61-80` ($r$ at `:75`, best rank `:80`) |
| V2 | v2-transcriptome (codex) | Co-expression rank test | benchmark 03 | `benchmarks/03_coexpression_fp_control.R:55-70` |
| V2 | v2-transcriptome (codex) | Benchmark-04 mediation target | benchmark 04 | `benchmarks/04_mediation_recovery.R:9-12, 89-90, 114-116` |
| V2 | v2-transcriptome (codex) | TWAS Pearson test | `cor_p` | `benchmarks/05_twas_power.R:103-112` ($t$ at `:110`, $p$ at `:111`) |
| V2 | v2-transcriptome (codex) | Bonferroni threshold | benchmark 05 | `benchmarks/05_twas_power.R:114, 146-147` |
| V2 | v2-io-formats | Allele frequency / MAF | `filter_geno` | `R/qc_filter_geno.R:290-297` |
| V2 | v2-io-formats | Monomorphic / MAF filters | `filter_geno` | `:311`, `:314`, `:318` |
| V2 | v2-io-formats | Heterozygosity | `filter_geno` | `:301`, `:322`, `:324` |
| V2 | v2-io-formats | Composite (genotypic) $r^2$ | `.ld_sweep` | `:556-565` |
| V2 | v2-io-formats | Pruning decision | `.ld_sweep` | `:532`, `:663`, `:665-668` |
| V2 | v2-io-formats | Window / step | `.ld_prune`, `.ld_sweep` | `:492-506`, `:605-617` (`win_at`), `:737-782` (slide) |
| V2 | v2-io-formats | Haplotypic $r^2$ | `.plink_hap_rsq`, `.plink_em_hethet`, `.plink_cubic_roots` | `:796-809` (r² at `:808`), `:822-879`, `:902-943` |
| V2 | v2-io-formats | VIF | `.ld_sweep` (vif branch) | `:683-736` (`solve` `:697`,`:721`; test `:725`) |
| V2 | v2-io-formats | Gabriel block MAF floor / tag | `.gabriel_blocks` | `R/qc_ld_methods.R:25`, `:40` |
| V2 | v2-io-formats | Major-allele flip (frequency) | `compute_flip` | `R/io_detect_format.R:362-364` |
| V2 | v2-io-formats | Reference flip | `compute_flip` | `:359` |
| V2 | v2-io-formats | Numeric coding | `numericalize_core` | `src/rust/src/numeric.rs:76-76`, `:138-145`, `:156-178` |
| V2 | v2-io-formats | HWE | — | — |
| V2 | v2-io-formats | Missing-rate filter | — | — |
| V2 | v2-io-formats (codex) | Raw genotype parsing and coding dispatch | `parse_hapmap_chars_to_raw()`, `.apply_coding()` | `R/io_detect_format.R:223-326`; `R/io_read_formats.R:18-194` |
| V2 | v2-io-formats (codex) | Exact two-locus ML cubic | `.plink_em_hethet()` | `R/qc_filter_geno.R:822-879` |
| V2 | v2-io-formats (codex) | Two-locus haplotype log-likelihood | `.plink_calc_lnlike()` | `R/qc_filter_geno.R:884-894` |
| V2 | v2-io-formats (codex) | Cubic real roots | `.plink_cubic_roots()` | `R/qc_filter_geno.R:902-943` |
| V2 | v2-io-formats (codex) | D' likelihood surface and class codes | `.plink_calc_lnlike_quantile()`, `.plink_blocks_classify()` | `R/qc_filter_geno.R:952-966`, `:978-1053` |
| V2 | v2-io-formats (codex) | Gabriel/Haploview block acceptance | `.plink_blocks_chrom()`, `.gabriel_blocks()` | `R/qc_filter_geno.R:1066-1180`; `R/qc_ld_methods.R:23-46` |
| V2 | v2-io-formats (codex) | Per-QTN marginal realized variance share | `.qtn_var()`, `qtn_table()` | `R/io_write.R:1288-1312`; `:1400-1403` |
| V2 | v2-io-formats (codex) | Per-gene transcriptome marginal share | `.tx_qtn_var()` | `R/io_write.R:1328-1345` |
| V2 | v2-io-formats (codex) | PLINK fixture allele serialization | fixture generator | `data-raw/plink_parity_fixtures.R:26-34` |
| V1 | v1-core-linkage | Additive genetic value | `genetic_effect` | `legacy_genetic_effect.R:109-112` |
| V1 | v1-core-linkage | Dominance value | `genetic_effect` | `:61-75` |
| V1 | v1-core-linkage | Epistatic value | `genetic_effect` | `:77-89` |
| V1 | v1-core-linkage | Centering | `genetic_effect` | `:93` |
| V1 | v1-core-linkage | Degree of dominance | `check_in` | `legacy_check_in.R:274` |
| V1 | v1-core-linkage | Geometric series | `check_in` | `:786-788, 796-798, 810-812, 823-825` |
| V1 | v1-core-linkage | Residual variance | `phenotypes` | `legacy_Phenotypes.R:216, 650, 1065, 1511` |
| V1 | v1-core-linkage | Single-trait phenotype | `phenotypes` | `:669-674, 1517-1522` |
| V1 | v1-core-linkage | Multi-trait residual | `phenotypes` | `:238-245, 1080-1086` |
| V1 | v1-core-linkage | Null trait (h2=0) | `phenotypes` | `:533-534, 86-92` |
| V1 | v1-core-linkage | Realized heritability | `phenotypes` | `:248-251, 336, 675-676, 876` |
| V1 | v1-core-linkage | Per-QTN PVE | `phenotypes` | `:260-268, 685-693` |
| V1 | v1-core-linkage | Residual seed | `phenotypes` | `:234, 665, 1076, 1513; 82, 529, 924, 1374` |
| V1 | v1-core-linkage | LD-QTN seed | `qtn_linkage` | `legacy_QTN_linkage.R:110, 383, 645, 971, 1220, 1243, 1455, 1479` |
| V1 | v1-core-linkage | LD measure | `qtn_linkage` | `:150-151, 176-177, 1025-1028, 1270-1273` |
| V1 | v1-core-linkage | Indirect selection | `qtn_linkage` | `:125-212` |
| V1 | v1-core-linkage | Direct selection | `qtn_linkage` | `:1229-1308` (fixed), `:978-1062, 1465-1546` (unfixed) |
| V1 | v1-core-linkage | MAF | `constraint`, `qtn_linkage`, `qtn_from_user` | `legacy_constraint.R:47-50`; `QTN_linkage.R:314-317`; `qtn_from_user.R:237-240` |
| V1 | v1-core-linkage | Genetic correlation (given `cor`) | `base_line_multi_traits` (not in group) | `legacy_Base_line_multi_traits.R:167-183` |
| V1 | v1-core-linkage (codex) | Heterozygote eligibility indicator | `constraint()` (correct); `create_phenotypes()` guard (wrong code) | `legacy_constraint.R:28-42`; `legacy_create_phenotypes.R:100-102, 848-935` |
| V1 | v1-core-linkage (codex) | Reported residual "sample correlation" | `phenotypes()` | `legacy_Phenotypes.R:483-491, 1335-1343` (verified: `matrix(0,...)` at `:484/:1336`, `cov2cor` at `:488/:1340`) |
| V1 | v1-core-linkage (codex) | QTN anchor sampling without replacement | `qtn_linkage()` | `legacy_QTN_linkage.R:109-113, 382-386, 644-648, 970-974, 1219-1223, 1454-1458` (verified: `while (s <= 10 & border)` at `:109,382,644,970,1219,1454`; `sample(index, ...)` at `:113,386,648,974,1223,1458`) |
| V1 | v1-core-linkage (codex) | User marker resolution and order | `qtn_from_user()` | `legacy_qtn_from_user.R:1358-1388` (add), `:1447-1460` (dom), `:1517-1530` (epi), `:1589-1602` (var) (verified: `selected_snps <- genotypes$snp %in% i` at `:1361, 1450, 1520, 1592`) |
| V1 | v1-core-linkage (codex) | Epistatic interaction diagnostic expansion | `qtn_from_user()` | `legacy_qtn_from_user.R:1582-1605` (verified: `e_len <- lengths(QTN_list$epi)/epi_interaction` at `:1532`, `each = epi_interaction` at `:1541`, `each = e_len[i]` at `:1544`) |
| V1 | v1-core-linkage (codex) | *(reconciler-added)* Direct-LD acceptance sentinel | `qtn_linkage()` | `legacy_QTN_linkage.R:1269-1269` (fixed add branch); same pattern `:123-125` (indirect add), `:976-983`, `:1467-1473` |
| V1 | v1-pleio-effects | Additive genetic value | `genetic_effect` | `R/legacy_genetic_effect.R:111–57` |
| V1 | v1-pleio-effects | Geometric effect series | `check_in` | `R/legacy_check_in.R:905–837` |
| V1 | v1-pleio-effects | Dominance genetic value | `genetic_effect` | `R/legacy_genetic_effect.R:120–77`; `R/legacy_check_in.R:274` |
| V1 | v1-pleio-effects | Epistatic genetic value | `genetic_effect` | `R/legacy_genetic_effect.R:151–91` |
| V1 | v1-pleio-effects | Total genetic value | `genetic_effect` | `R/legacy_genetic_effect.R:186–96` |
| V1 | v1-pleio-effects | Component variances | `genetic_effect` | `R/legacy_genetic_effect.R:115,59,69,76,87,90` |
| V1 | v1-pleio-effects | MAF of a QTN | `qtn_pleiotropic`, `qtn_partially_pleiotropic`, `constraint` | `R/legacy_QTN_pleiotropic.R:205–165`; `…partially….R:165–169`; `R/legacy_constraint.R:47–43` |
| V1 | v1-pleio-effects | QTN sampling seeds (pleiotropic) | `qtn_pleiotropic` | `R/legacy_QTN_pleiotropic.R:174,318,408,485` |
| V1 | v1-pleio-effects | QTN sampling seeds (partially) | `qtn_partially_pleiotropic` | `…partially….R:125,153,446,475,621,637` |
| V1 | v1-pleio-effects | Het re-sample loop | both QTN samplers | `R/legacy_QTN_pleiotropic.R:180–145,324–333`; `…partially….R:131–140,160–170,452–461,481–493` |
| V1 | v1-pleio-effects | Whitening | `base_line_multi_traits` | `R/legacy_Base_line_multi_traits.R:151–109,178–184` |
| V1 | v1-pleio-effects | Colouring & rescale | `base_line_multi_traits` | `…multi_traits.R:160–120,185–194` |
| V1 | v1-pleio-effects | PD repair | `make_pd` | `R/legacy_make_pd.R:61–21` |
| V1 | v1-pleio-effects | Residual variance (context) | `phenotypes` | `R/legacy_Phenotypes.R:216,239,650` |
| V1 | v1-pleio-effects | vQTL residual scale | `vQTL` | `R/legacy_vQTL.R:89,62–71,76–89` |
| V1 | v1-pleio-effects | vQTL sample heritability | `vQTL` | `R/legacy_vQTL.R:251–193` |
| V1 | v1-pleio-effects | Numeric normalisation | `genotypes/numeric_df` | `R/legacy_Genotypes.R:79–62` |
| V1 | v1-pleio-effects | `maf_cutoff` filter | `genotypes` | `R/legacy_Genotypes.R:163–111` |
| V1 | v1-pleio-effects (codex) | Fully pleiotropic QTN-set construction | `qtn_pleiotropic` | `R/legacy_QTN_pleiotropic.R:169-201` (add), `313-345` (dom), `403-423` (var), `480-500` (epi) |
| V1 | v1-pleio-effects (codex) | Partial-pleiotropy locus construction | `qtn_partially_pleiotropic` | `R/legacy_QTN_partially_pleiotropic.R:165-233` (add+dom-loop), `294-360` (add only), `441-533` (dom), `616-684` (epi) |
| V1 | v1-pleio-effects (codex) | Epistatic physical-marker demand | `qtn_pleiotropic`, `qtn_partially_pleiotropic` | `R/legacy_QTN_pleiotropic.R:524-532`; `R/legacy_QTN_partially_pleiotropic.R:586-611` (also `:654`) |

# Appendix C — line remapping from `c511c6f` to `e666a2c`

Every reference in the entries and in Appendix B was moved to `e666a2c` with `git diff -U0 -M c511c6f` line mapping (Appendix A is kept as it was at the audit commit). 35 references were unchanged, 396 shifted with unchanged content, and 101 point into lines that were edited after the audit (listed below: the equation there may have changed and must be re-checked against the code before the port relies on it). 33 references name a file outside this repository or an ambiguous shorthand and were left as written.

| Audit reference | Current reference |
|---|---|
| `R/grammar_realize.R:447-455` | `R/grammar_realize.R:578-588` |
| `R/grammar_realize.R:409-425 (422)` | `R/grammar_realize.R:510-528 (523)` |
| `R/effects_series.R:20-53 (52)` | `R/effects_series.R:23-72 (64)` |
| `R/grammar_realize.R:36-40 (40)` | `R/grammar_realize.R:54-58 (58)` |
| `R/grammar_realize.R:704-723` | `R/grammar_realize.R:1064-1092` |
| `R/grammar_simulate_phenotype.R:649-670` | `R/grammar_simulate_phenotype.R:1000-1030` |
| `R/effects_pleioarch.R:53-81` | `R/effects_pleioarch.R:65-93` |
| `R/effects_pleioarch.R:143-159` | `R/effects_pleioarch.R:201-217` |
| `grammar_simulate_phenotype.R:666-669` | `grammar_simulate_phenotype.R:1018-1022` |
| `grammar_realize.R:447-455` | `grammar_realize.R:578-588` |
| `grammar_realize.R:91-93` | `grammar_realize.R:143-145` |
| `grammar_simulate_phenotype.R:684-691` | `grammar_simulate_phenotype.R:1061-1071` |
| `R/cross_mating.R:51-52` | `R/cross_mating.R:77-78` |
| `R/cross_mating.R:53` | `R/cross_mating.R:79` |
| `R/cross_mating.R:55` | `R/cross_mating.R:81` |
| `R/cross_mating.R:102` | `R/cross_mating.R:345` |
| `R/cross_map.R:157` | `R/cross_map.R:157` |
| `R/cross_map.R:161-166` | `R/cross_map.R:168-173` |
| `R/cross_breed.R:68-72` | `R/cross_breed.R:68-72` |
| `meiosis.rs:232-234` | `meiosis.rs:359-369` |
| `cross_mating.R:102` | `cross_mating.R:345` |
| `cross_mating.R:217-222` | `cross_mating.R:647-652` |
| `cross_mating.R:98` | `cross_mating.R:337` |
| `io_as_numeric.R:8-12` | `io_as_numeric.R:8-12` |
| `numeric.rs:41-52` | `numeric.rs:76-88` |
| `io_detect_format.R:259-272` | `io_detect_format.R:286-297` |
| `R/select_ind.R:281` | `R/select_ind.R:460` |
| `R/select_ind.R:282` | `R/select_ind.R:461` |
| `R/select_ind.R:474` | `R/select_ind.R:769` |
| `R/select_schemes.R:130` | `R/select_schemes.R:165` |
| `R/select_ind.R:16-21` | `R/select_ind.R:16-21` |
| `R/select_ind.R:498-511` | `R/select_ind.R:799-830` |
| `R/select_schemes.R:22-61` | `R/select_schemes.R:25-68` |
| `R/select_schemes.R:86-94` | `R/select_schemes.R:102-118` |
| `R/select_schemes.R:89` | `R/select_schemes.R:106` |
| `R/select_ocs.R:315` | `R/select_ocs.R:424` |
| `R/select_blup.R:376-378` | `R/select_blup.R:752-761` |
| `R/select_combining.R:30-37` | `R/select_combining.R:30-37` |
| `R/select_blup.R:253-298` | `R/select_blup.R:427-650` |
| `R/select_progeny.R:89-90` | `R/select_progeny.R:113-114` |
| `R/transcriptome_simulate.R:490-492` | `R/transcriptome_simulate.R:616-618` |
| `R/transcriptome_mimic.R:101-110` | `R/transcriptome_mimic.R:149-172` |
| `benchmarks/03_coexpression_fp_control.R:52-67` | `benchmarks/03_coexpression_fp_control.R:55-70` |
| `benchmarks/04_mediation_recovery.R:9-12` | `benchmarks/04_mediation_recovery.R:9-12` |
| `R/qc_filter_geno.R:203-210` | `R/qc_filter_geno.R:290-297` |
| `src/rust/src/numeric.rs:41-45` | `src/rust/src/numeric.rs:76-76` |
| `legacy_QTN_linkage.R:1231-1237` | `legacy_QTN_linkage.R:1269-1269` |
| `R/legacy_vQTL.R:57` | `R/legacy_vQTL.R:89` |
| `R/legacy_vQTL.R:184` | `R/legacy_vQTL.R:251` |
| `R/legacy_QTN_partially_pleiotropic.R:120-190` | `R/legacy_QTN_partially_pleiotropic.R:165-233` |
| `partially.R:104-109` | `partially.R:137-137` |

