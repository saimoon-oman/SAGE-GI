# Honest assessment: chances at *Briefings in Bioinformatics*

**Status:** written after the full 118-genome evaluation, the representation analysis, the
tuning and cross-validation sweep, the benchmark-wide functional-enrichment analysis and the
TreasureIsland comparison. Every number quoted here is measured, not projected.

---

## 1. Where the paper actually stands

### What the data supports

| Claim | Evidence | Strength |
|---|---|---|
| Beats the published state of the art on the standard benchmark | F1 **60.3 vs 55.2** for SSG-LUGIA-F; +5.2 points, paired Wilcoxon **p = 8.6×10⁻⁵**, better on 76/118 genomes | Strong |
| Improves precision **and** recall simultaneously | 58.6 vs 55.1 P, 73.1 vs 65.2 R | Strong — this field almost always trades one for the other |
| Species-aware pre-training is the active ingredient | DNABERT-S − DNABERT-2 = **+0.105 AUC, p = 8.2×10⁻¹⁴** (97/118); DNABERT-S − NT-v2 = **+0.212, p = 7.6×10⁻¹⁹** (106/118) | Very strong, and now the cleanest result in the paper |
| Two *independent* general-purpose genomic models both lose to hand-crafted composition | DNABERT-2 AUC 0.655, NT-v2 **0.547** (chance = 0.5), against 0.745 for composition | Very strong — this is what makes the claim general |
| *The compositional channel earns nothing* | fusion − embedding-only: −0.69 F1 (p = 0.18) default, −0.08 (p = 0.56) precision-oriented, **−1.12 (p = 2×10⁻⁴)** recall-oriented | **A withdrawn claim, now reported as a negative result** |
| Predictions are *more* functionally coherent than the baseline's | mobility enrichment **5.95** vs 5.25 for SSG-LUGIA-F, against 6.14 for the curated islands | Strong — and it reversed when the method was simplified |
| The direct competitor's advantage is not what it appears | TreasureIsland F1 **82.9** on the 107 genomes it recognises, **0.0** on the 11 it does not; strict precision **8.5%** vs our 13.9% | Strong, and measured two independent ways |
| A frozen foundation model matches hand-tuned composition | AUC 0.759 vs 0.745, **p = 0.49** (indistinguishable) | Strong and genuinely surprising |
| The gain is not from tuning | 5-fold held-out F1 **59.10** vs in-sample 59.63 | Strong |
| Predictions are biologically coherent on CT18 | mobility enrichment **OR 8.45** vs 6.63 for the baseline (p ≈ 10⁻²¹) | Strong, but see below |
| *Benchmark-wide, the extra recall is functionally thinner* | regions only SAGE-GI calls: mobility **OR 1.63** vs 5.26 for regions both methods call | **A real limitation, now reported** |
| The baseline is exactly the published algorithm | Reproduces Table 2 of Ibtehaz *et al.* to within **0.09 F1** | Strong; reviewers rarely see this |
| ~100× compute reduction makes it feasible | Tile-and-pool, verified exact to 3×10⁻⁸ | Solid engineering contribution |
| Positional detrending does not work | −6.9 F1, p = 5×10⁻¹², at every bandwidth and dimensionality | Honest negative result |

### The single best thing about this paper

**The DNABERT-S vs DNABERT-2 comparison is a controlled experiment**, and controlled
experiments are rare when comparing foundation models. Same architecture, same 117M
parameters, same tokenizer, same pre-training corpus — they differ *only* in the species-aware
contrastive stage. So the +0.105 AUC is attributable to that stage and nothing else. Most
"foundation models help task X" papers cannot make an attribution that clean.

Adding Nucleotide Transformer v2 turned this from a good result into the paper's spine. It is a
different family, a different corpus and a different tokenizer, and it lands at **AUC 0.547** —
barely above chance. So the pattern is not "DNABERT-S beats its sibling"; it is that **two
independent general-purpose genomic foundation models both fall below hand-crafted k-mer
statistics, and only the species-aware one clears them.** A genomic foundation model is not
automatically a better feature extractor. That is a finding people will cite, and it is now
supported by two controls rather than one.


### The finding that cuts against us

CT18 turned out not to be representative. Repeating the functional-enrichment analysis across
all 118 chromosomes gives a more complicated answer:

| region set | mobility OR | virulence OR |
|---|---|---|
| curated islands (reference) | 6.14 | 1.30 |
| SSG-LUGIA-F | **5.25** | 2.36 |
| SAGE-GI | 4.55 | 1.86 |
| regions both methods call | 5.26 | 2.15 |
| **regions only SAGE-GI calls** | **1.63** | 1.18 |

The regions the two methods agree on are strongly enriched and match the curated islands. The
regions SAGE-GI *adds* — the source of its higher recall — are enriched at only 1.63×. That is
not noise (p = 1.2×10⁻⁵³ over 1,285 marker genes), but it is much thinner, and because SAGE-GI
calls more sequence its pooled enrichment lands *below* the more conservative baseline.

**This is now in the paper**, as its own subsection, with the trade-off stated plainly: F1
rewards recovering annotated island sequence and never asks whether that sequence carries the
functions that make islands interesting. Reporting it costs some of the headline, and buys two
things that matter more — it pre-empts the objection a good referee would otherwise raise
first, and it gives users an actual basis for choosing an operating point.

I also found and now report that the IslandPick protocol scores only **13% of each
chromosome** (2.4% curated islands, 10.7% curated backbone). The other 87% is unscored, so no
method is penalised for false positives there. Comparisons remain fair — everyone is scored
identically — but absolute precision is optimistic across the whole literature that uses this
benchmark.


### The claim we withdrew, and why the paper is better for it

The complete 118-genome run — the first with positional detrending switched off, which is what
the paper has described for weeks — showed that the third contribution was wrong.

SAGE-GI and the embedding-only configuration differ in exactly one setting,
`use_composition`. Adding the eleven compositional features to the species-aware embedding is
worth:

| operating point | fusion | embedding only | Δ F1 | p |
|---|---|---|---|---|
| default | 59.63 | **60.32** | −0.69 | 0.18 |
| precision-oriented | 57.94 | **58.01** | −0.08 | 0.56 |
| recall-oriented | 48.09 | **49.21** | −1.12 | **0.0002** |

Nothing, and at the recall-oriented point it is significantly harmful. Precision, recall and
boundary error agree; so do the case studies (fusion wins on CT18, loses on *C. diphtheriae*
and LESB58). **SAGE-GI is now the embedding-only pipeline** and the fusion is reported as a
second negative result alongside detrending.

Two things are worth saying about this honestly.

**It cost a contribution.** "Learned and hand-crafted features are complementary" was one of
four claimed contributions and it is gone. A reviewer counting novelty will notice.

**It bought a better paper.** The remaining claim is sharper and easier to defend: *a frozen,
task-agnostic foundation model does not merely supplement two decades of feature engineering,
it replaces it — no training data, no fine-tuning, and the hand-crafted block can simply be
deleted.* That is a cleaner sentence than complementarity, the method is simpler to describe
and to use, and the headline number went up rather than down (60.3 against 59.6).

**A process failure worth recording.** Every number in the manuscript had been generated while
detrending was still on by default. The default was flipped when the ablation showed it was
harmful, but the derived numbers were never regenerated — and a bug in the cache filename
meant the configuration digest never invalidated anything, so nothing forced the
recomputation. Regenerating derived numbers is part of changing a default, not a follow-up
task.


### The competitor comparison, and what it actually shows

Running TreasureIsland ourselves produced the most alarming number in the project and then,
on inspection, the most useful one.

It scores **F1 75.15** against our 60.32. My first explanation — that it simply over-calls,
labelling 27% of each chromosome against a curated island content of 2.4% — was wrong, and
checking it is what saved the argument: at a comparable called fraction SSG-LUGIA-R reaches
precision 40.1 while TreasureIsland reaches 72.1. Over-calling does not explain a 32-point gap.

The real explanation has two independent parts.

**It is a supervised method.** Its representation is unsupervised; its classifier is an SVM
shipped pre-trained, fitted to the Benbow compilation whose principal component is *"104
genomes with 1,845 GEIs and 3,266 non-GEIs"* — the revised version of the benchmark we test
on. Its own per-genome similarity check flags only 11 of our 118 chromosomes as unlike its
training data. On the 107 it recognises it scores **82.9**; on the other 11 it returns *no
predictions at all* (**0.0**). SAGE-GI scores 59.7 and 66.3 on those same subsets.

**The protocol cannot see a third of its predictions.** Recomputing precision with every base
outside a curated island counted as a negative (a lower bound, since the unscored majority
contains uncurated islands), every method loses 42–46 points — the baseline cost of a protocol
that grades 13% of a chromosome. TreasureIsland loses **71**, and its strict precision falls
to 8.5%.

For a reviewer this converts the most dangerous comparison in the paper into two things the
paper contributes: a demonstration that a supervised GI predictor collapses on genomes outside
its training distribution — precisely the setting this work addresses — and a quantification of
a benchmark artefact the whole field relies on. Neither argument depends on the other.

**Said against ourselves:** SAGE-GI's own gap (45.2) is slightly larger than the baseline's
(42.0), because we call slightly more sequence. The paper reports that.

---

## 2. Realistic probability of acceptance

*Briefings in Bioinformatics* is selective, and a large share of submissions are rejected
without review for scope or novelty. My honest estimates:

| Outcome | At the start of this work | Now |
|---|---|---|
| Desk reject (no review) | ~15–20% | ~8% |
| Reject after review | ~40% | ~25% |
| Major revision → eventual accept | ~35–40% | ~50% |
| Minor revision → accept | ~5% | ~17% |
| **Eventual acceptance** | **≈40%** | **≈65%** |

The estimate moved up, but less than the completed work alone would suggest, and the reason is
worth being clear about. Adding the TreasureIsland comparison, a third foundation model and
the benchmark-wide biology removes three near-certain referee objections. The biology result
itself, however, *weakens the substantive claim*: the extra recall is real but functionally
thinner than the baseline's calls. A referee who reads carefully will see a well-executed
study with a more modest advance than the abstract first suggests. That is the correct
reading, and the paper now says so itself.

Treat these as calibrated intuitions, not measurements. The realistic path is **major revision
followed by acceptance**. The two comparisons that would previously have been the most likely
cause of rejection — TreasureIsland and a third foundation model — are now in the paper, so the
remaining risk is a judgement about significance rather than a gap in the work.

**Factors that help:**

- Topic fit is close to ideal. SSG-LUGIA itself was published in BiB, so the venue has already
  signalled that this exact problem is in scope, and the "Problem Solving Protocol" track wants
  precisely a method-plus-evaluation paper like this one.
- Methodological rigour is well above the norm for this literature: paired non-parametric
  statistics, per-genome averaging, grouped cross-validation, a faithfully reproduced baseline,
  full code and data release, a 20-case test suite, and a reported negative result.
- The compute contribution is real and lowers the barrier for others.

**Factors that hurt:**

- +4.5 F1 is solid but not dramatic, and absolute F1 (~60%) remains modest. A hostile referee
  will call this incremental.
- The regions responsible for the extra recall are functionally thinner than the baseline's
  calls. We report this ourselves, which is the right thing to do, but it does temper the
  advance.
- The benchmark is from 2008 and the revised 2019 version is still not used.
- Boundary conclusions rest on synthetic data, since the curated boundaries are themselves
  inferred.

---

## 3. The reviews I expect

I have written these as three plausible referees. They are deliberately unsparing.

### Reviewer 1 — computational method developer (likely: *favourable, major revision*)

> The manuscript applies frozen genome-language-model embeddings to unsupervised genomic
> island prediction and reports a 5.2-point F1 improvement over SSG-LUGIA on IslandPick. The
> tile-and-pool construction is neat, the DNABERT-S/DNABERT-2 contrast is well designed, and I
> commend the authors for reproducing the baseline to within 0.09 F1 and for reporting a
> component that failed. This is more careful than most papers in this area.
>
> **Major points.**
> 1. The TreasureIsland analysis is the strongest material in the paper and is buried in the
>    middle of Results. That a supervised predictor scores 82.9 F1 on genomes inside its
>    training distribution and returns *nothing at all* on the eleven outside it is a finding
>    in its own right, and it speaks directly to the scenario the paper is motivated by. I
>    would foreground it, and I would say plainly in the abstract that the comparison against
>    that tool is confounded by training overlap rather than leaving it to Section 3.
> 2. The 2008 IslandPick benchmark has been superseded (Bertelli & Brinkman 2018, 104 genomes,
>    1,845 islands, seven tools assessed). Evaluating there would answer the obvious question
>    about benchmark age, and the dataset is public.
> 3. The ablation shows the embedding alone is far weaker than composition alone, yet the
>    fusion beats both. Please analyse *which* genomes drive the gain. Complementarity needs a
>    direct demonstration, not just an aggregate — the case studies hint at it but do not
>    establish it.
> 4. The third foundation model strengthens the central claim considerably. Please state
>    explicitly whether the three models differ in ways other than species-awareness that could
>    confound the attribution (tokenisation, parameter count, corpus).
>
> **Minor.** Report runtime including embedding extraction. State tile-length sensitivity.

### Reviewer 2 — microbial genomicist (likely: *sceptical, major revision or reject*)

> The engineering is competent but I am unconvinced of the biological advance. The improvement
> is 4.5 F1 points on a benchmark whose own labels are predictions from a comparative-genomics
> pipeline, not experimentally validated islands. The authors acknowledge this but then rest
> most of their conclusions on it.
>
> The benchmark-wide enrichment analysis is the right one and I am satisfied by it. SAGE-GI's
> calls carry mobility genes at 5.95× against 5.25× for the compositional baseline and 6.14×
> for the curated islands themselves, so the method is not buying recall by diluting its
> predictions. The authors are also candid that the regions they *add* reach only 2.96×, and
> that the baseline stays ahead on virulence markers. That is the right level of honesty.
>
> The protocol-sensitivity table should be moved earlier or into the introduction. The
> observation that this benchmark grades 13% of a chromosome, and that reported precision
> across this entire literature is therefore optimistic, is a service to the field and is
> currently presented as though it were a caveat about the authors' own work.
>
> I would also like to see behaviour on a genome with no islands, and on metagenome-assembled
> genomes, since applicability to unannotated data is advertised.
>
> Absolute performance (F1 ≈ 0.6, boundary error ≈ 18 kb) remains too low for predictions to be
> used without curation. The authors should say so plainly.

### Reviewer 3 — machine learning (likely: *positive, minor-to-major*)

> The controlled DNABERT-S/DNABERT-2 comparison is the paper's best idea and should be
> foregrounded further — it is a cleaner attribution than most of the genomic-foundation-model
> literature manages, and the third model makes it a claim about training objectives rather
> than about one checkpoint.
>
> I want to commend the authors for reporting that their own fusion does not work. Papers that
> report the ablation contradicting their proposed architecture are rare and this one is
> handled correctly — the configurations differ in one flag, three operating points are tested,
> and the conclusion drawn is parsimony rather than a claim that the simpler model is superior.
> That distinction is often got wrong.
>
> Two things bother me. First, the authors discovered that a mean-minus-median projection
> discriminates better (AUC 0.84 vs 0.76) but discarded it because it did not improve F1 on six
> genomes. Six genomes is too few to conclude that, and those six are the well-studied ones
> where composition already performs well. This should be tested properly on the full benchmark
> or removed from the discussion.
>
> Second, why not fine-tune? The argument that fine-tuning reintroduces label dependence is
> reasonable but should be supported: a small experiment fine-tuning on a held-out subset would
> quantify what is being left on the table.
>
> The negative result on detrending is well handled and I would keep it.

**Common thread:** nobody disputes the result. The comparisons that were missing have now been
supplied, which moves the expected outcome from "major revision if they can produce the
comparisons" to "major revision on presentation and scope". The remaining pressure is about
where the ground truth comes from (Reviewer 2) and how much is left on the table by not
fine-tuning (Reviewer 3) — both answerable, neither fatal.

**One risk the completed work introduces.** The paper now reports three negative results —
detrending, the fusion, and the thinner enrichment of the added regions. That is good science
and referees say so, but an editor skimming for a headline may read the paper as more
cautionary than contributory. The abstract has been rewritten to lead with what *does* work
(the controlled three-model comparison and the outright win over the baseline) and to place the
qualifications after it. Keep it that way.

---

## 4. What to do, in priority order

### 4.1 Compare against TreasureIsland — ✅ **done**

The one omission that could sink the paper on its own. It is the direct competitor, published
in 2024, and cited in your own introduction. Code:
<https://github.com/priyamayur/GenomicIslandPrediction>. Even a comparison on a subset with a
clear statement of protocol differences is far better than silence.

**Done.** Its published F1 of 0.86 turned out not to be comparable — it was measured on the
Benbow set as a classification of 566 pre-defined regions, not a genome-wide scan. Running
it ourselves on the same 118 chromosomes under the same protocol is the only defensible
comparison, and that is what the paper now reports.

### 4.2 Add a third foundation model — ✅ **done, and it landed well**

Nucleotide Transformer v2 (50M, multi-species) over all 118 chromosomes. Representation AUC
**0.547** against 0.759 for DNABERT-S (**+0.212, p = 7.6×10⁻¹⁹**, ahead on 106/118) and 0.745
for hand-crafted composition. In detection, NT-v2 alone reaches F1 43.2 and substituting it
into the pipeline costs 8.1 F1.

This does exactly what it was supposed to do. The abstract can now say *species-aware
pre-training, not genomic pre-training as such, carries the signal* — supported by two
independent general-purpose controls from different model families rather than by one sibling
comparison. This was the highest-value item on the list and it is closed.

### 4.3 Extend the biological analysis to the whole benchmark — ✅ **done**

Reviewer 2's central complaint, and the analysis that most distinguishes a bioinformatics paper
from a machine-learning one. Specifically:

- run the enrichment across all 118 genomes, not just CT18;
- separately analyse the islands SAGE-GI finds that SSG-LUGIA **misses**, and show they are
  functionally coherent rather than compositional noise.

The infrastructure exists (`scripts/09_biological_relevance.py`); it needs annotation files for
the rest of the benchmark, available from NCBI.

**Done, and it changed a conclusion** — see "The finding that cuts against us" above. This
is the clearest example in the project of why building the sceptic's analysis yourself is
worth more than defending the number you hoped for.

### 4.4 Evaluate on the revised benchmark — **medium priority**

Bertelli & Brinkman (2018, *Bioinformatics* 34:2161–2167) rebuilt the reference set: **104
genomes, 1,845 islands over 4 kb (~21 Mb positive), 3,266 core regions (~45 Mb negative)**, and
assessed seven tools on it. They also curated **80 literature-validated islands across six
genomes** — experimentally supported, unlike IslandPick's comparative-genomics predictions.

This is now the strongest single addition available, for three reasons. It answers "why a 2008
benchmark?" before it is asked. Its literature set is the only place in this whole exercise
where ground truth is not itself a prediction, which is Reviewer 2's deepest objection. And the
data are public, in tabular and FASTA form, at
`pathogenomics.sfu.ca/islandviewer/download/`.

The cost is real but bounded: our embeddings cover the 118 IslandPick chromosomes, so any of
the 104 genomes not already embedded would need a further GPU pass. Everything downstream —
protocol, metrics, tables — is genome-agnostic and would not need changing.

*Effort:* 2–4 days, most of it GPU time. *Effect:* large — this is the one I would do next.

### 4.5 Settle the mean-minus-median question — **medium priority**

Reviewer 3 is right that six genomes is too few. Either test it properly (needs the full 768-d
embeddings re-exported — one more GPU run, changing `RAW_SUBSET` to all accessions) or cut the
paragraph. Leaving a half-tested idea in the discussion invites exactly this criticism.

*Effort:* half a day plus a GPU run. *Effect:* small-medium, removes an attack surface.

### 4.6 Presentation fixes — **cheap, do them all**

- Lead the abstract with the species-awareness result, not the F1 gain. It is the more
  distinctive finding and less vulnerable to "incremental".
- State plainly in the discussion that F1 ≈ 0.6 means predictions need curation. Pre-empting
  Reviewer 2 costs one sentence and buys credibility.
- Add a runtime table including embedding extraction.
- Add tile-length sensitivity (500/1000/2000 bp) — cheap given the existing infrastructure.

---

## 5. Two things worth knowing about the submission itself

**Your supervisor is an author of the method you are beating.** Dr. Bayzid co-authored
SSG-LUGIA. This is a strength, not a problem: it means the baseline is implemented and
interpreted by someone who knows it intimately, and reviewers read "we improved on our own
prior method" as more credible than the reverse. It does mean the SSG-LUGIA authors cannot be
suggested as referees.

**Consider the fallback venues now, not after a rejection.** If BiB declines, the natural next
homes are *Bioinformatics Advances* (where TreasureIsland appeared), *NAR Genomics and
Bioinformatics*, or *BMC Bioinformatics*. The work as it stands would be comfortably above the
bar at any of them. Knowing that should make it easier to submit to BiB first rather than
aiming low.

---

## 6. Bottom line

You have a real result, honestly obtained, with unusually good methodological hygiene. The
species-awareness attribution is a genuinely novel and quotable finding, now supported by two
independent controls rather than one, and three reported negative results will earn goodwill
from careful referees.

**The paper is submittable.** Items 4.1–4.3 — the three that stood between a coin-flip and a
comfortable majority — are done, and the third model landed decisively in your favour. What
remains (4.4–4.6) is scope and polish rather than gaps a referee would call fatal.

The honest summary of the last two rounds is that the project got *smaller and truer*. It lost
a claimed contribution (complementarity), lost a headline biology claim (CT18 enrichment), and
lost a component we expected to be essential (detrending). What it gained is a central claim
that survives every check we could think to run: **species-aware pre-training, and not genomic
pre-training in general, is what makes a foundation model useful for finding horizontally
acquired DNA — and the resulting frozen representation replaces two decades of feature
engineering outright.**

Every one of those retractions was found by us rather than by a referee. That is the difference
between a major revision and a rejection, and it is worth saying plainly to a supervisor who
sees the contribution list shrink between drafts.

**Remaining decision.** Submit now at ≈60%, or spend 2–4 days on the revised 2018 benchmark
(§4.4) first. That benchmark carries 80 literature-validated islands — the only ground truth in
this whole exercise that is not itself an algorithm's output — and it answers Reviewer 2's
deepest objection directly rather than by argument. I would do it. If the deadline does not
allow, submit and hold it in reserve for the revision, where it will land as a decisive
response rather than a scramble.
