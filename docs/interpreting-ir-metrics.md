# Interpreting the two IR metrics: junction vs depth

IRfinder-mdl reports two intron-retention ratios per intron, `ir_ratio_junction`
and `ir_ratio_depth`. They answer different questions, and the **depth** ratio is
deliberately the more *conservative* estimate of bona-fide retention. This note
explains the mechanics of why, so you can read the two numbers — and especially
their disagreement — correctly.

## Definitions

For an annotated intron `[s, e)` (donor boundary `s`, acceptor boundary `e`):

```
IR_junction = (R_L + R_R) / (R_L + R_R + S_L + S_R)

IR_depth    = D / (D + max(S_L, S_R)),   D = intron_coverage_bp / intron_length
```

| symbol | TSV column | meaning |
| --- | --- | --- |
| `S_L`, `S_R` | `splice_left`, `splice_right` | reads that **splice** the donor / acceptor — a CIGAR `N` op anchored at that boundary (within `--jitter`, with `--anchor` matched bp on the exonic side) |
| `R_L`, `R_R` | `retain_left`, `retain_right` | reads that cross the donor / acceptor as **continuous matched sequence** (≥ `--anchor` bp on both the exonic and intronic side, no `N` in the anchor window) |
| `D` | `intron_coverage_bp / intron_length` | mean per-base unspliced **depth** across the intron body |

The crucial detail is how `intron_coverage_bp` is accumulated
(`classify_read_vs_intron` in `irfinder_mdl/quantify.py`): a read contributes its
matched bp inside `[s, e)` **only if its CIGAR has no `N` op intersecting
`[s, e)`**. Any read that splices any part of this intron contributes **zero** to
the depth signal.

```python
n_intersects_intron = any(ns < e and ne > s for ns, ne in n_skips)
intronic_bp = 0 if n_intersects_intron else matched_bp_in(s, e, matched_intervals)
```

Both ratios are reported as `.` (undefined) when there is no splice evidence to
divide against (`splice == 0`): a relative retention measure at a locus with no
observed splicing is uninformative.

## Why depth is more conservative

A read earns "retention" credit very differently under the two metrics. Two
mechanisms make depth the stricter view.

### Mechanism 1 — depth excludes partially-spliced reads; junction rewards them

A junction "retain" read only has to cross **one** boundary as matched sequence;
it says nothing about the rest of the intron. That lets an **alternative
splice-site** read look like retention.

Consider a read that uses a cryptic acceptor 200 bp inside a 2,766 bp intron —
it splices `s → (e-200)`, then runs as matched sequence across the annotated
acceptor `e` into the downstream exon:

```
intron   s|========================================|e        2766 bp (spliced out in the canonical isoform)
read      [======  N: spliced s -> (e-200)  ======][MMMMM|MMMMM]
                                                    └200bp┘└ downstream exon
                                              matched across e
```

Scored against the annotated intron `[s, e)`:

| | junction | depth |
| --- | --- | --- |
| 5' boundary `s` | `N` starts at `s` → **`S_L`** (splice) | — |
| 3' boundary `e` | matched across `e`, no `N` in window → **`R_R`** (retain!) | — |
| intron body | not measured | `N` op intersects `[s, e)` → **0** coverage bp |

The *same read* lands in the junction **numerator** (`R_R`, pushing IR up) but is
**excluded** from the depth numerator and instead loads the depth **denominator**
via `S_L` (pushing IR down). It physically kept only 200 / 2,766 = 7 % of the
intron — it is alternative splicing, not retention. Junction calls it retention;
depth refuses to.

### Mechanism 2 — depth gives fractional credit by length; junction gives full credit per boundary

Even for a *genuine* unspliced read, junction awards a whole `R` count for merely
crossing a boundary with the minimum anchor (e.g. 8 bp into the intron). Depth
credits only the bp actually deposited: that read adds `8 / 2766 ≈ 0.003` of one
depth unit. To earn a full unit of `D` — the equivalent of one splice read in the
denominator — a read must blanket the **entire** intron.

So depth only climbs when reads physically cover the whole intron body, which is
the literal definition of intron retention.

## When the two agree

They **converge** when retention reads fully span the intron. If 100 reads each
span `[s, e)` unspliced and 100 reads splice it:

- `D = 100`, `max(S_L, S_R) = 100` → `IR_depth = 100 / 200 = 0.5`
- `R_L = R_R = 100`, `S_L = S_R = 100` → `IR_junction = 200 / 400 = 0.5`

Divergence appears **only** when the reads scored as boundary-retention do not
blanket the intron — because they are alternative-splice-site reads, or because
they are short/partial reads crossing just a boundary. In both cases
`IR_depth < IR_junction`.

## A worked example (SBX pilot)

The RNF185 intron `chr22:31192702-31195468` (2,766 bp) observed on SBX cDNA data:

| metric | value |
| --- | --- |
| `ir_ratio_junction` | **0.66** |
| `ir_ratio_depth` | **0.04** |
| `retain_right` | 138 |
| mean depth `D` | ~1.5 |

138 reads cross the 3' boundary as matched sequence (large `R_R` → junction says
66 % retained), but they are alt-acceptor reads that spliced out most of the
intron, so the body is nearly empty (`D ≈ 1.5` over 2.7 kb → depth says 4 %).
Depth reports the physical truth: the intron is mostly *gone*, not retained.

## How to read them together

- **Both high** → genuine intron retention (the body is unspliced *and* the
  boundaries behave like retention).
- **Junction high, depth low** → alternative splice-site usage near a boundary,
  or incomplete reads — not bona-fide retention. The *gap* between the two
  ratios is itself the diagnostic.
- **Both low** → constitutively spliced.

Lead with `ir_ratio_junction` for sensitivity; use `ir_ratio_depth` to confirm
that a junction signal reflects real retention rather than alternative splicing.

## Caveat: depth is coverage-hungry on long introns

Depth's length normalization makes it **biased against long introns when reads
are short relative to the intron**. A fully-retained 320 bp intron is easily
"filled" by single SBX reads (~400 bp median), but an 11 kb intron needs ~28
overlapping unspliced reads to reach `D = 1`. So at modest coverage depth can
*under*-call retention in long introns.

This is the price of its specificity:

- **junction** — sensitive, but inflated by alternative splicing;
- **depth** — specific, but coverage-hungry on long introns.

Reporting both, and watching where they disagree, is what lets you tell true
retention from alternative splicing.
