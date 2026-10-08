# 頭花の防御は「侵入口」を塞ぐのか、それとも別の入口へ誘導するのか？

2026-10-08. Supplemental mechanism/falsification note for PR #44. **Head only**: orientation, real phyllaries/spines, floret portal, stage, guild. No leaf trait work, new field authorization, fitted Cirsium effects, or altered frozen manuscript.

## Gap left by the first access ledger

One approach episode has a world approach direction and a contacted head interface, but a final observed entry route alone does **not** tell us whether a visitor (i) initially chose that portal, (ii) reached it after being stopped at another, or (iii) found the preceding route simply unavailable due to stage. Selection on successfully entered insects hides failed approaches. Furthermore, a photographed image-vertical axis is **not** a gravity-calibrated head orientation, nor is outline projection a real spine.

The smallest additional observation is the **time-ordered contact attempt sequence within the same independently identified approach episode**. It records first contact, blockage (including uncertain), subsequent contacts, access, and complete coverage. A failed attempt without complete continuous coverage cannot be classified as abandonment. A different last portal from first portal does not prove causal redirection without comparison against a control. Maintain head and plant as biological replication units.

## Minimal two-portal process null

Let `p` = fraction of all identifiable arrivals initially contacting outer involucre (rather than floret portal); `a` = probability of reaching a reproductive zone from that initial outer contact; `b` = the corresponding probability from the accessible floret portal; `s` = probability of attempting the floret portal **after failing** at the outer portal, all conditional on the same taxon, head stage and observation protocol.

`S = p [a + (1-a) s b] + (1-p)b`.

Under **fixed arrivals, p, b, s**, lowering outer access `a` cannot increase successful entry: `dS/da = p(1-sb) >= 0`. A blocked outer portal can have essentially no apparent fitness effect when nearly everyone reroutes successfully (`s*b ≈ 1`), but *post-failure rerouting of a fixed frequency by itself cannot reverse the entry effect*. This monotonicity is an algebraic restricted-model result, not Cirsium field evidence. Stage-specific rewiring changes `b` and invalidates a common effect assumption.

Only a change in arrival pressure, initial portal probabilities `p`, rerouting propensity `s`, alternative-portal access `b`, density/time constraints, or stage-specific guild composition permits a reversal under this model. Do not claim that any particular change actually occurs on thistles.

### Distinct predictions

1. **True blocking**: same-guild, same-stage initial outer contacts have more physically failed attempts under naturally stronger armature; entry per arriving candidate falls even after incorporating subsequent attempts. A randomized study needs separate approvals.
2. **Within-episode rescue / bypass**: first-contact failures rise, subsequently reached floret portals rise, but overall entry per arrival may be unchanged. During **bud stage**, the floret portal may be absent/closed, so the same outside architecture could function more strongly than at anthesis. This is not evidence for evolutionary reversal.
3. **Pre-entry reallocation**: a change in head cues/orientation could alter the first contacted portal `p`. In that case the sign of total entry may oppose the direct physical outside barrier effect. Because `p` is potentially post-treatment, fit/plot overall entry per randomized head alongside conditional processes rather than treating conditional coefficients as causal mediation.
4. **False no-benefit**: a consumer still enters after a blocked first attempt yet later fails to oviposit, so final success must be measured independently (egg evidence and seed fate); conversely pollinators may switch contact routes but not deposit pollen.

### Illustrative synthetic boundary (not field measurements)

With all first contacts at the outer involucre (`p=1`), opening access `a=0.7` versus barrier `a=0.3`, and rerouting after failure `s=0.8`:
- closed floret route (bud: `b=0`): total entry `0.70 -> 0.30`;
- exposed floret route (anthesis: `b=0.95`): `0.928 -> 0.832`; the barrier effect is attenuated but not reversed.

A hypothetical change from `p=0.8` to `p=0.2` with `a=0.5 -> 0.1`, `b=0.95`, `s=0` changes entry `0.59 -> 0.78`; that *requires* a changed initial choice, not simply rerouting with otherwise fixed behavioral parameters.

The supplied program `analysis/check_capitulum_portal_switching_bounds_v1.py` recomputes all identities and audits an **empty** add-on contact-attempt CSV. Synthetic examples are never treated as field observations.

## Observation integration and go/no-go

Append-only add-on `data/intake/capitulum_contact_attempt_sequence_v1.csv`: one row for each filmed attempted contact **within the same approach_episode_id** that is already present in the existing event ledger; `attempt_index` is 1,2,... ordered by first observed contact time. It is not a dataset of inferred intent or of independent insects. Keep missing/uncertain obstruction as `NA`. Heads with no arrivals remain in the original valid continuous-video effort ledger, not this attempt table. For new real rows, cross-check the existing event/effort join keys and footages before model fitting (this script checks sequence-only internal invariants; cross-ledger joins are a separately required gate).

If full approach space or between-portal motion is missing, report unknown route switching; an unobserved return is not a confirmed abandoned insect. For the hypothesis about real spines, require botanical measurements of length, rigidity, direction and phyllary interstice with gravity-referenced orientation and stage. Camera insect-order classifications are not adequate for a seed-feeder label. Keep pre-entry guild identification independent of the observed endpoint. Stage, host race, sex of head, terminal-versus-secondary head rank, weather, and alternative *Cirsium* flowering resources are explicit rival causes. Plants/heads (not attempts) define independent replication. Final filled/viable achenes plus independently verified pollen receipt/egg success/early parasitoid killing remain the fitness ceiling.

## Prior work and originality limit

- Thomas 2003, DOI 10.22543/0090-0222.2085: neutralizing *C. discolor* sticky traps did **not** raise seed predation or lower seed production; the null result directly cautions against inferring function from sticky surfaces.
- Thomas 2007, DOI 10.22543/0090-0222.2188: pollinators avoided sticky traps, some seed predators also bypassed traps.
- Leonard et al. 2013, DOI 10.1371/journal.pone.0055914: experimentally changed legitimate versus robbing *bumblebee* handling on artificial flowers with guides; flexible portal choice is **prior art outside Cirsium**, not a new discovered behavioral principle.
- Gijsman et al. 2020, DOI 10.1016/j.gecco.2020.e00945: *C. pitcheri* head rank and `Larinus` infestation matter; not a demonstrated spine-route switch.
- Russell & Louda 2005, DOI 10.1007/s00442-005-0204-3: alternative congeneric flowers can predict egg distribution independent of the focal head's geometry.

**Current outcome:** FORMAL_MECHANISTIC_PREDICTIONS_AND_EMPTY_INTAKE_ONLY. No real Cirsium switching, stage effect, genotype–phenotype covariance or adaptive trait reassembly identified.
