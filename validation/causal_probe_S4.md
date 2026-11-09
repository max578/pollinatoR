# Causal probes for Stage 4 (Step 0)

2026-11-10. dagitty 0.3.4, glmmTMB 1.1.14, pollinatoR 0.4.0 (branch `stage4-causal` at
2a5748b). Script: `validation/causal_probe_S4.R` (synthetic data only). Run from the package
root with `Rscript --vanilla validation/causal_probe_S4.R`.

## (a) dagitty on the field-study graph

| Graph | Effect | Adjustment set returned |
|---|---|---|
| Minimal: `trt -> veg -> yield`, `trt -> yield`, `rain -> veg`, `rain -> yield` | total | `{}` |
| Minimal | direct | `{ rain, veg }` |
| Field, branch H (veg at harvest; latent `fertility -> veg`, `fertility -> yield`; `forest -> fertility`, `forest -> bees -> yield`) | total | `{}` |
| Field, branch H | direct | **none** (empty list: not identifiable) |
| Field, branch H, exposure `bees` | total | `{ forest }` |
| Field, branch P (no `trt -> veg` edge) | total | `{}`; descendants of `trt` are `yield` only |

**Finding 1 (changes the plan).** Once a latent site variable such as soil fertility drives both
vegetative weight and yield, the direct effect of exclusion net of vegetative weight has no
adjustment set. Conditioning on `veg` opens `trt -> veg <- fertility -> yield`, and no measured
variable blocks it. The adjusted 31.5% share then estimates neither the direct effect nor the
total effect; it is biased by mediator-outcome confounding, the standard objection to "adjusting for
a mediator". The minimal graph returns `{ rain, veg }` only because it has no such confounder. Plan
v0.1 requires the adjust set to contain a dagitty direct-effect set when a graph is given; the probe
shows that check must also handle the empty case: `pollin_estimand(effect = "direct", dag = ...)`
refuses when `adjustmentSets(effect = "direct")` is empty and says why (an unblocked
mediator-outcome path, naming the latent node). Without a graph, the assumption row for a direct
effect must state "no unmeasured confounder of mediator and outcome" in addition to "no
treatment-by-mediator interaction".

**Finding 2.** A DAG cannot draw effect modification. The who-delivers question (bees modify the
exclusion contrast) is drawn as `bees -> yield`, and the set returned for it (`{ forest }`) is the
set for the bees-to-yield effect. That is the right confounding set for the moderation term under
the plan's assumption (no unmeasured site variable that changes the contrast and correlates with
the exposure), but the package must say the set is borrowed from the main-effect graph.

**Finding 3.** `localTests(type = "cis")` runs on the site-level subgraph (44 simulated sites, data
generated from that graph): three implied independencies, all with 95% intervals covering zero
(estimates -0.14, -0.23, -0.25; p 0.36, 0.13, 0.10). With 44 sites the intervals are about +/-0.29
wide, so the check detects only strong violations, as the plan assumed. Latent nodes must be
declared `[latent]` or dagitty tests independencies that involve unmeasured columns.

## (b) `pollin_moderation()` on a site-level-only exposure

The 0.4.0 behaviour is not a crash. On an exposure that is constant within site:

- glmmTMB drops `bees_within` and `trtOpen:bees_within` as rank-deficient, with a message;
- the result table still carries a `trtOpen:bees_within` row, labelled `within (robust)`, with NA
  estimate, standard error and p-value;
- the between term is estimated (0.012, SE 0.048; planted 0.053 per SD of the site means).

The non-positive-definite Hessian seen in the earlier probe is **not** caused by the zero column.
Rerunning with the site-intercept SD set to 0, 0 / 0.3, 0.3 (slope SD 0 or 0.2) gives pdHess
FALSE, FALSE, TRUE, TRUE: it comes from simulating no site intercept variance, a boundary fit of
the `(trt | Site)` term. The fix sketch (fit without the within column) reproduces the between
estimate exactly (0.01235) and the same boundary warning on that data. The shipped roxygen example
(within variation present) fits with pdHess TRUE.

**Corrected defect statement.** On a site-level-only exposure, `pollin_moderation()` returns an NA
row labelled robust, plus glmmTMB's rank-deficiency message, instead of saying the exposure has no
within-unit variation. Fix (Step 5): before fitting, drop any within component whose variance is
zero, message "`<ex>` varies between units only; no within-unit (robust) term", and label the
between row as the only evidence. The between term on 44 sites is then the whole of the paper_2
driver estimate, which is why the negative control and the adjustment set matter there.

## (c) Oracle grounding for C7

Grounded 2026-11-10 by web fetch (journal pages, Europe PMC, Crossref, OSF API). The shell
sandbox has no outbound network, so nothing below was re-read by a script here; science.org and
Wiley returned 403, so the Dainese supplementary tables and the Bishop and Nakagawa full text were
not read. Numbers below are as printed on the pages that were read.

| Candidate | Data, public and scriptable | Licence | Target number | Fit to S4 |
|---|---|---|---|---|
| Zou et al. 2017, *BMC Ecology* 17: 6, doi:10.1186/s12898-017-0116-1 | Additional file 1 (`12898_2017_116_MOESM1_ESM.xlsx`), Springer static URL, no login; 134 plants (72 open-top cage, 62 closed cage), 18 fields; seed set, 1000-seed weight, straw biomass; no yield column | Not stated on the record (BMC default CC BY 4.0, not confirmed) | Table 1, log 1000-seed weight, field random: t -3.25, df 111, P 0.002; straw biomass (sqrt) t -3.72 | Open-top versus closed cage: exclusion with the cage held roughly constant, i.e. the `reference = "sham"` route; no uncaged arm, so no cage-effect estimate |
| Bishop and Nakagawa 2021, *J. Appl. Ecol.* 58: 1030-1042, doi:10.1111/1365-2664.13830 | OSF B3E9U, `vf_metaanalysis_data_main.csv` (373 rows, raw means/SD/n, enclosure-control column `bagcontrol`), no login | None (`node_license: null`) | Enclosure moderator estimate: not read (full text blocked) | Meta-analytic moderator, not a package estimand |
| Dainese et al. 2019, *Sci. Adv.* 5: eaax0121, doi:10.1126/sciadv.aax0121 | Database S1 (`aax0121_Database_S1.xlsx`) in the supplement only; scripted download not confirmed | CC BY-NC 4.0 (article) | "Roughly a third" of the landscape effect on pollination mediated by richness (table S5, not read) | Path analysis; non-commercial clause sits badly with an MIT package reading it |
| Madin and Nelson 2023, *PLOS ONE* 18: e0289799 | Reanalysis of Dainese's database, not a new data set | CC BY | | Second reading of Dainese only |

**D9 ruling (Max, 2026-11-10): Zou et al. 2017 plus simulation.** Recommendation as made: Table 1 log 1000-seed weight (t -3.25, df 111) by
`pollin_estimand(reference = "sham")` on the open-top arm, plus the seed-set row as a second target.
It is the only candidate with public plant-level data a run-time read can fetch, and the model is
fully specified by the paper. It does not oracle the `cage` family itself: none of the candidates has
an uncaged arm. The adjustment route stays oracled by planted truth only (C3-C5).

Two things the record adds to the programme beyond the oracle:

- Closed-cage plants carried 39% (+/-11%) more straw biomass than open-cage plants (Results; the
  abstract says 18%, an internal inconsistency, so the Table 1 t value is the target, not the
  percentages). Exclusion raising vegetative biomass is direct evidence that vegetative weight is
  post-treatment in canola, i.e. branch H, the branch where the adjusted 31.5% is not identified
  (Finding 1).
- University of Reading Research Data Archive record 141 (Garratt and Degani, oilseed rape
  pollinator exclusion, CC BY 4.0, three CSVs) was found but not assessed. If it has an uncaged arm
  it is the oracle for the `cage` family; Mike Garratt is on the programme team.
