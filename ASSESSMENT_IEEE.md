# Honest assessment: chances at IEEE/ACM *Transactions on Computational Biology
and Bioinformatics* (TCBB)

**Status:** written after the synthetic-benchmark repair (full-dimensional re-embedding,
480/480 conditions, valid embedding arm), the BiB + IEEE manuscript updates, and a clean
audit. Every number quoted here is measured, not projected. Read together with
`ASSESSMENT.md` (the Briefings in Bioinformatics assessment); this document re-asks the
same questions for a different venue with a different bar.

**Why TCBB and not another IEEE journal.** TCBB is the flagship IEEE venue for exactly
this kind of work: computational methods for molecular biology, evaluated on public
benchmarks, with released code. *TBME* and *JBHI* are clinical/applied journals — a
genomic-island detector with no clinical cohort would be desk-rejected for scope.
*TPAMI/TKDE* would reject on field. TCBB is the right home and the only top-tier IEEE
journal where this paper belongs.

---

## 1. Where the paper stands after the synthetic repair

The repair changed the paper's evidentiary structure, so the venue calculus moved with
it. What the re-run established, measured over 2 hosts x 24 donors x 5 insert lengths
x 2 seeds (480 conditions, zero failures):

| Claim | Evidence | Strength |
|---|---|---|
| Embedding channel validates with exact ground truth | DNABERT-S-only recovery \synS{} overall vs \synComp{} composition; ahead at 20--40 kb | Strong — this arm was withdrawn a week ago |
| Boundary conclusions rest on exact coordinates too | MABE \synSMABE{} (embedding) vs \synCompMABE{} (composition) | Strong; the curated-data-only caveat is gone |
| Species-awareness transfers under exact ground truth | DNABERT-2-only \synDtwo{} overall, trails both channels | Strong; mirrors the IslandPick attribution |
| Small inserts challenge every channel | \synSShort{} / \synCompShort{} at 5--10 kb, near the 10 kb window scale | Honest limit, stated in the paper |
| Comp-only reproduced bitwise (69.1/37.9) | determinism check passed | Strong; pipeline hygiene referees notice |

Everything in `ASSESSMENT.md` §1 still holds (F1 **60.3 vs 55.2**, +5.2, p = 8.6x10^-5;
DNABERT-S - DNABERT-2 = **+0.105 AUC**; NT-v2 at 0.547; TreasureIsland 82.9/0.0 split;
strict-precision analysis; three reported negative results). The repair adds the one
thing that assessment listed as missing: boundary evidence that does not depend on
curated labels.

### What the repair did and did not buy

It bought the removal of the paper's most visible open wound. A TCBB reviewer reading
the pre-repair manuscript would have found a boundary-resolution claim whose "primary
evidence" section ended in a withdrawal — the single easiest major-revision demand to
write ("repeat with full-dimensional embeddings and resubmit"). That demand is now
un-writable; the experiment exists and supports the claim.

It did not buy a bigger headline. The embedding ties composition on synthetic overall
(68.9 vs 69.1) rather than beating it, and both channels fail most 5 kb inserts. The
paper reports this plainly, which is correct, but nobody should expect the repair to
move the acceptance needle by more than a few points on its own. Its value is
defensive: it closes the attack surface rather than opening new territory.

---

## 2. Realistic probability of acceptance at TCBB

TCBB is a step up in selectivity from BiB's protocol track: broader scope, no
venue-fit advantage (SSG-LUGIA was a BiB paper, not a TCBB one), a higher novelty
bar, and reviewers who routinely demand a second benchmark. Calibrated estimates:

| Outcome | Probability |
|---|---|
| Desk reject (no review) | ~10% |
| Reject after review | ~35% |
| Major revision → eventual accept | ~40% |
| Minor revision → accept | ~15% |
| **Eventual acceptance** | **≈50-55%** |

The realistic path is the same as at BiB — **major revision followed by acceptance**
— but the number is lower, for three concrete reasons. First, TCBB reviewers apply a
novelty test the BiB protocol track does not: the pipeline composes existing pieces
(tile pooling, PCA projection, MCD envelopes, CUSUM), and a hostile referee will call
the method "engineering around a pretrained model" rather than a new algorithm. The
controlled three-model comparison is the defence, and it is a good one, but it has to
carry that weight explicitly. Second, the 2008 benchmark is older than most TCBB
reviewers will tolerate without the 2018 revised set alongside it — at BiB this is a
"would be nice"; at TCBB it is closer to an expectation. Third, absolute performance
(F1 ≈ 0.6) reads worse in a transactions journal, where reviewers compare against the
whole ML literature rather than the GI-detection literature.

**Factors that help at TCBB specifically:**

- Reproducibility well above the transactions norm: public code, data, 22 tests,
  deterministic pipeline (Comp-only bitwise identical across re-runs), every number
  generated from results rather than typed.
- The three-model attribution (two independent general-purpose controls below
  hand-crafted statistics, species-aware above) is the kind of clean experimental
  design transactions referees reward, and it is now supported *twice*: detector-free
  AUC on curated data and detection/MABE on exact synthetic ground truth.
- Three reported negative results plus the withdrawn-then-repaired synthetic arm read
  as scientific seriousness, and the repair itself demonstrates exactly the
  revise-and-resubmit behaviour editors want to see.
- The manuscript is already in IEEE two-column format at 11 pages, inside TCBB's
  typical length envelope.

**Factors that hurt at TCBB specifically:**

- +5.2 F1 is solid but not dramatic, and the synthetic tie (68.9 vs 69.1) lets a
  referee frame the whole contribution as "matches hand-crafted features at far
  higher compute cost" — the paper must keep winning that framing fight with the
  boundary-error and no-training-data arguments.
- No wet-lab or literature-validated islands: the 2018 benchmark's 80 validated
  islands remain the strongest un-done experiment, and at TCBB it matters more.
- Single-group study with no independent replication; the supervisor co-authored the
  beaten baseline (a strength for credibility, but TCBB referees cannot be drawn
  from that circle either).
- Runtime: embedding extraction needs a GPU session; CPU-only users face ~70 h.
  Stated honestly in the paper, but a referee may still ask for a lighter operating
  point.

---

## 3. The reviews I expect at TCBB

Written as three plausible transactions referees. Deliberately unsparing.

### Reviewer 1 — computational biology methods (likely: *favourable, major revision*)

> This manuscript applies frozen genome foundation-model embeddings to unsupervised
> genomic island detection and evaluates them carefully: a reproduced baseline (within
> 0.09 F1 of published), paired non-parametric statistics with CIs, grouped
> cross-validation, a supervised-competitor analysis that actually re-runs the
> competitor, and a synthetic benchmark with exact ground truth that now covers the
> embedding channel too. The DNABERT-S/DNABERT-2/NT-v2 comparison is well designed
> and the conclusion — that species-aware pre-training, not genomic pre-training in
> general, transfers — is supported two independent ways. This is above the bar for
> methodological care.
>
> **Major points.**
> 1. The 2008 IslandPick benchmark should be joined by the 2018 revised set
>    (Bertelli & Brinkman, 104 genomes, 1,845 islands, 80 literature-validated).
>    The authors already cite it as future work; for a transactions paper it should
>    be present work, at least on the validated subset.
> 2. The detector itself (two-stage MCD) is taken as given. Please ablate it: does
>    a simple robust distance on the same host-projected features reach the same
>    F1? If the embedding is doing the work, the paper should show the detector is
>    not.
> 3. Report embedding-dimensionality and tile-length sensitivity. d=24 is stated to
>    be "not delicate" — show it.
>
> **Minor.** Runtime table is appreciated; add peak memory. State the random seeds
> governing the synthetic placements (they are recorded — say where).

### Reviewer 2 — microbial genomics domain expert (likely: *sceptical, major revision or reject*)

> The engineering is competent and the writing is unusually honest, which I want to
> acknowledge first: the authors report the thinner enrichment of their added
> regions (2.96x vs 5.85x agreed), the virulence lag, the failed detrending, and a
> withdrawn-then-repaired experiment, all in the main text. That honesty is why I am
> reviewing the claims rather than rejecting the paper.
>
> My concern is biological, not statistical. F1 ≈ 0.6 with ~18 kb boundary error
> means these predictions cannot be used without curation — the authors say so, and
> I agree with them, which leaves the question of what the advance buys a working
> microbiologist. The mobility enrichment (5.95x vs 5.25x baseline) is real but
> modest, and on virulence markers the 20-year-old baseline still wins. The
> synthetic benchmark is a genuine improvement over curated-only evaluation, but
> synthetic inserts from distant donors are the easy case; the hard case is an
> island from a close relative, compositionally near the host, and the paper should
> say how performance degrades along that axis rather than averaging over it.
>
> I also want to see behaviour on (a) a genome with no islands and (b) at least one
> metagenome-assembled genome, since applicability to unannotated data is the
> advertised motivation. Without some evidence there, the MAG claim is aspirational.
>
> The 80 literature-validated islands of the 2018 set would answer my deepest
> objection — that every label here is itself a prediction — and I join Reviewer 1
> in asking for them.

### Reviewer 3 — machine learning for genomics (likely: *positive, minor-to-major*)

> The central finding — two independent general-purpose genomic foundation models
> fall below hand-crafted k-mer statistics while the species-aware one clears them,
> shown both as detector-free AUC and as end-to-end detection on exact synthetic
> ground truth — is a clean, quotable result about training objectives, and the
> genre ("foundation models help task X") rarely manages attribution this clean.
> Reporting the fusion ablation against the proposed architecture is handled
> correctly, including the parsimony-not-superiority framing.
>
> Two things I want. First, the no-fine-tuning stance is principled but
> unquantified: fine-tune DNABERT-S on a held-out subset and report what is left
> on the table. If fine-tuning gains little, the frozen claim is strengthened; if
> it gains a lot, readers deserve to know the price of the principle. Second, the
> mean-minus-median direction (AUC 0.84 vs 0.76 on six genomes) is either tested on
> the full benchmark now that full-dimensional embeddings exist for 31 genomes, or
> cut. Six genomes is not evidence, and the full-dim archives the authors just
> built are the obvious vehicle.
>
> The detrending negative result is well handled and should stay.

**Common thread:** nobody disputes the measurements, and the synthetic repair removed
the one objection that would previously have been a one-line major revision. The
remaining pressure converges on two experiments — the 2018 benchmark (Reviewers 1+2)
and a quantified fine-tuning headroom (Reviewer 3) — plus ablations of the detector
and hyperparameters that are cheap given the infrastructure. None is fatal; all cost
time. The 2018 set is the one I would do before submitting; the fine-tuning run can
ride the same GPU sessions.

**One risk the repair introduces.** The paper now reports a tie, not a win, on the
only exact ground truth it has (68.9 vs 69.1). A skimming editor or referee may read
"matches composition at 100x the compute" and stop. The abstract and introduction
must keep the framing the evidence supports: the win is on the curated benchmark
(+5.2 F1, both metrics together), the boundary error (6.5 vs 8.8 kb), and the
no-training-data setting — the synthetic result corroborates rather than headlines.

---

## 4. What to do before submitting to TCBB, in priority order

### 4.1 Evaluate on the 2018 revised benchmark — **highest priority, 2–4 days**

Bertelli & Brinkman 2018: 104 genomes, 1,845 islands, and 80
literature-validated islands — the only ground truth in this area that is not
itself an algorithm's output. Two of three expected referees will ask for it; doing
it first converts their strongest demand into a reported result. Data are public
(pathogenomics.sfu.ca/islandviewer/download/). Cost is a further GPU embedding pass
for genomes not already covered, plus no pipeline changes (everything downstream is
genome-agnostic).

### 4.2 Quantify the fine-tuning headroom — **medium-high priority, ~1 GPU session + 1 day**

Reviewer 3's request is legitimate and cheap: fine-tune DNABERT-S on a held-out
genome subset, report the delta. Either outcome helps the paper (strengthens the
frozen claim or prices the principle). Do it on the same GPU sessions as 4.1.

### 4.3 Ablate the detector and the hyperparameters — **medium priority, CPU-only, 1 day**

Reviewer 1's asks need no GPU: (a) simple robust distance vs two-stage MCD on the
same features; (b) d and tile-length sensitivity curves from existing embeddings
(d) — tile-length needs re-embedding at 500/2000 bp, so fold it into the 4.1 GPU
pass or report d-sensitivity alone.

### 4.4 MAG + island-free genome behaviour — **medium priority, ~half day**

Run the frozen pipeline on one island-free chromosome and one MAG; report that it
behaves sanely (few/no calls, no crashes). This is a smoke test, not a benchmark,
and it defuses Reviewer 2's applicability challenge at low cost.

### 4.5 Presentation for a transactions venue — **cheap, do all**

- Keep the abstract framed as §3's risk note demands: curated win + boundary +
  no-training-data first, synthetic tie as corroboration.
- Add peak memory to the runtime paragraph; state synthetic seeds location
  (Reviewer 1's minor — both are one sentence).
- TCBB requires author biographies with photos at acceptance: prepare them now,
  not during the revision scramble.
- The 11-page IEEE manuscript has room (TCBB regular papers run ~12–14 pages);
  spend it on 4.1–4.3, not prose.

---

## 5. Submission mechanics worth knowing

- **No word limit binds the IEEE version.** TCBB judges by pages (~14 max for
  regular papers); at 11 pages there is comfortable headroom for the 2018 results.
- **ORCID** is still a placeholder in the BiB version and will be needed for any
  IEEE submission too — register before submitting anywhere.
- **Preprint posture:** a bioRxiv/medRxiv-style preprint is standard in this field
  and TCBB permits prior preprints; posting one at submission time establishes
  priority while review runs its 3–6 months.
- **If TCBB declines,** the natural next homes are *Bioinformatics Advances*,
  *NAR Genomics and Bioinformatics*, or *BMC Bioinformatics* — the work as it
  stands clears all three comfortably, and the TCBB reviews would directly
  strengthen a resubmission. Knowing that should make it easier to aim at TCBB
  first rather than below it.

---

## 6. Bottom line

You have a real result with unusually good hygiene, and the synthetic repair closed
the one hole a reviewer could previously have driven a major revision through. The
paper is submittable to TCBB today at roughly a coin flip — **≈50–55% eventual
acceptance, most likely via major revision** — with the proviso that the 2018
benchmark is the single addition that would move the number most. Two to four days
on §4.1 (+ the fine-tuning run in §4.2 on the same GPU sessions) is the highest-value
work remaining; everything else is scope and polish.

The honest one-sentence summary for a supervisor: *species-aware pre-training, and
not genomic pre-training in general, is what makes a foundation model useful for
finding horizontally acquired DNA — shown on curated labels, on exact synthetic
ground truth, and against two independent general-purpose controls — and the
resulting frozen representation replaces two decades of feature engineering
outright, with better boundaries and no training data.*
