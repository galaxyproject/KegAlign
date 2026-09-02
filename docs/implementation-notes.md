# KegAlign Implementation Notes

A map of the parts of KegAlign where the reasoning is not obvious from the code, written for whoever has to change them next.

**This is not a changelog.** It is organised by mechanism, because someone arriving cold needs to know how scoring works, not what a particular pull request did.

Two conventions, both borrowed from LASTZ's own source. **Every claim names the file and line that backs it, or the command that checks it** — a claim you cannot check is a claim that will rot. And **known problems are written down where they live**, not omitted; LASTZ marks its own with `$$$`.

The authoritative explanation of any mechanism is the comment beside the code. Where this document and a source comment disagree, **the source comment wins and this document is stale.**

**Which tree this describes.** `main` as of the v0.3.1 release. On any earlier commit the symbols named here — `MIN_SEED_WEIGHT`, `MAX_SEED_WEIGHT`, `MAX_SEED_SPAN`, `drop_duplicate_hsps()`, `die_after()`, `sorted_commands()`, `lastz_command_sort_key()` — will not all exist, and line numbers will not match.

1. [The compressed alphabet](#1-the-compressed-alphabet)

2. [Seeds and the seed pattern](#2-seeds-and-the-seed-pattern)

3. [The scoring matrix](#3-the-scoring-matrix)

4. [Masking](#4-masking)

5. [Entropy](#5-entropy)

6. [HSPs and segments](#6-hsps-and-segments)

7. [Checking parity against LASTZ](#7-checking-parity-against-lastz)

8. [The orchestration layer](#8-the-orchestration-layer)

9. [Known problems](#9-known-problems)

## 1. The compressed alphabet

### Two representations, live at the same time

A base exists in two forms for the whole run, and almost every confusion in this codebase comes from not knowing which one you are looking at.

|  | raw | compressed |
|---|---|---|
| what it is | the FASTA byte, as read | one of eight codes, 0–7 |
| where it lives | host DRAM, `ref_DRAM` / `query_DRAM` | GPU global memory, `d_ref_seq` / `d_query_seq` |
| who reads it | seeding, and the seed-position table | extension, and only extension |
| case | preserved | collapsed into `L_NT` |

The eight codes and what they score are in §3. What matters here is that **neither form is derived from the other**. The raw bytes are uploaded to the GPU verbatim and mapped in place by a compression kernel, so the two forms are two mappings of the same bytes made by different code for different consumers — not a conversion you could follow from one to the other.

Be careful about which invariant that implies. The host and GPU mappings _are designed to disagree_: lower case becomes `N_NT` for seeding (no seed) and `L_NT` for extension (`bad_score`), and that disagreement _is_ the two-stage masking mechanism in §4. Reconciling them would destroy masking. The invariant that does have to hold is narrower, and it is the subject of the trap below: **the two GPU kernels must agree with each other**, and nothing checks that they do.

### Three character mappings, which must agree

Three functions turn a byte into something else — one `switch` and two `if`/`else if` chains, so grepping for `switch` finds only the first. A change to the alphabet touches all three:

| function | where | feeds | default case |
|---|---|---|---|
| `NtChar2Int()` | host, `ntcoding.cpp:11` | seeding | `N_NT` |
| `compress_string()` | GPU, `seed_filter_interface.cu:19` | reference extension | `X_NT` |
| `compress_string_rev_comp()` | GPU, `seed_filter.cu:114` | query extension, both strands | `X_NT` |

The host one is six cases wide and lumps everything else together:

```c
inline uint32_t NtChar2Int (char nt) {
    switch(nt) {
        case 'A': return A_NT;
        case 'C': return C_NT;
        case 'G': return G_NT;
        case 'T': return T_NT;
        case 'N': return N_NT;
        default:  return N_NT;      // lower case, IUPAC, '&', junk — all the same here
    }
}
```

That coarseness is deliberate and correct: `GetKmerIndexAtPos()` rejects the whole span the moment it sees `N_NT`, so lower case, an ambiguity code and a stray byte all produce the same outcome — _no seed_ — and there is nothing to gain by telling them apart. The GPU kernels are seven cases wide and _do_ tell them apart, because there the difference is a score.

> **Trap**
>
> **The two GPU kernels are near-duplicates, and adding a code means editing both.** They are in different files, one under `common/` and one under `src/`, and neither mentions the other. A ninth code — the fix §8 records for the literal `X` — means four edits: `parameters.h`, both kernels, and `build_substitution_matrix()`. Miss one kernel and the reference and query disagree about what a byte means, which shows up as HSPs that are wrong rather than absent.

### Why the query gets its own kernel

Not duplication for its own sake. The query kernel writes _two_ outputs from one read:

```c
dst_seq[i]            = dst;       // forward
dst_seq_rc[len-1-i]   = dst_rc;    // reverse complement, same pass
```

LASTZ's convention is to search both strands of the query against one strand of the reference, so the reference is never reverse-complemented and needs only one output. Fusing the complement into the compression pass means `src_seq[i]` is read once and the complement is a second table lookup rather than a second kernel launch over the whole sequence. The reference kernel is the same code with that half deleted.

Note how the complement is taken: `A_NT` ↔ `T_NT`, `C_NT` ↔ `G_NT`, and **every other code maps to itself** — `L_NT`, `N_NT`, `X_NT` and `E_NT` are all their own complement. That is what makes the fused write correct, and it is also why collapsing the four lower-case bases into one code costs nothing _here_ even though it costs the `[unmask]` option in §4.

### There are two reverse complements, for the same reason there are two representations

The GPU one above serves extension. Seeding needs its own, over raw bytes, and gets it from a 128-entry ASCII table (`ntcoding.cpp:82`):

```c
const char *rev_comp =
    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
    "NNNNNN&NNNNNNNNNNNNNNNNNNNNNNNNN"
    "NTVGHNNCDNNMNKNNNNYSAABWNRNNNNNN"
    "NtvghnncdnnmnknnnnysaabwnrnNNNNN";
```

Three things are worth reading out of that table. It complements the **IUPAC ambiguity codes properly** (`V`→`B`, `H`→`D`, `R`→`Y`), which the GPU version does not bother with because they all land in `X_NT` anyway — though not _every_ code: a literal `X` maps to `N` rather than to itself. That is harmless, because this table feeds only seeding and seeding rejects `X` and `N` alike — but know it before trusting the table wholesale. It **preserves case** for every letter, which is what makes the minus strand mask identically to the plus strand. And at index 38 it maps `&` to itself, so block separators survive the reversal in place.

### The separator is excluded twice, by unrelated mechanisms

Sequence blocks are concatenated with a single `&` between them (`main.cpp:399`, `main.cpp:519`). An HSP must never span a boundary, and two independent things guarantee it:

| stage | mechanism |
|---|---|
| seeding | `NtChar2Int('&')` falls to `default` → `N_NT` → `INVALID_KMER`. No seed's span may contain one. |
| extension | `E_NT` scores `-10 * xdrop` — a wall, not a penalty. See §3. |

These are not redundant, and it is worth being precise about why. The seeding arm guarantees only that no _seed span_ contains a separator; extension starts at a seed and walks outward through arbitrary sequence, so **the only thing stopping an HSP from crossing a boundary is the wall**. Do not treat seeding as a backstop for it — weaken the wall and nothing catches the case.

Note also that the seeding arm is accidental: `&` reaches `N_NT` through the same `default` arm that catches junk bytes, not through a case naming it. It works, and it is not written down anywhere in the source. Now it is.

### Span and weight are bounded separately, and both bounds are load-bearing

`GetKmerIndexAtPos()` reads the whole span into a fixed stack array before consulting `shape_pos[]`:

```c
uint32_t nt[MAX_SEED_SPAN];                   // ntcoding.cpp:64
for(int i = 0; i < seed_size; i++){           // seed_size is the SPAN, not the weight
    nt[i] = NtChar2Int(sequence[pos+i]);
```

So the span needs its own bound. The weight bound in §2 does not imply it — the two are independent, and `--seed` with twelve `T`s scattered over a hundred positions has a legal weight of 12 and a span of 100. Until September 2026 neither was checked before the arrays were filled; `MAX_SEED_SPAN` and `MAX_SEED_WEIGHT` now name the two limits, and both are checked before `GenerateShapePos()` writes anything.

## 2. Seeds and the seed pattern

### Weight is not span

A seed pattern has two sizes, and confusing them is the fastest way to get this wrong.

| pattern | shape | span | weight |
|---|---|---|---|
| `12of19` | `TTT0T00TT00T0T0TTTT` | 19 | **12** |
| `14of22` | `TTT0T0TT00TT00T0T0TTTT` | 22 | **14** |

**Span** is the length of the shape string and lives in `cfg.seed.size`. **Weight** is the count of match positions — the `T`s — and lives in `cfg.seed.kmer_size`, returned by `GenerateShapePos()` (`common/ntcoding.cpp:27`). The `0`s are don't-care positions: they sit inside the span but contribute no bits to the k-mer.

In a shape, `T` marks a match position that _permits a transition_ and `1` one that does not. `main.cpp` normalises a user's custom pattern so every match position becomes `T`, which is why weight and transition count are equal in practice — but the code keeps them separate, and so should you.

### How many seeds a position produces

`src/seeder.cpp` walks query positions. For each position whose k-mer is valid it pushes one seed for the position itself, then one more for every match position at which a transition is allowed:

```
seeds per valid position  =  1 + count(IsTransitionAtPos(t)), t < kmer_size
                          =  13 for 12of19
                          =  15 for 14of22
```

A position is valid only if _every_ base across the **whole span** — don't-care positions included — is upper-case `ACGT`. `GetKmerIndexAtPos()` loops `i < seed_size` and rejects on the first `N_NT`, and `NtChar2Int()` sends everything outside `ACGTN` — lower case included — to `N_NT`. LASTZ does the same, by restarting its word collection on any illegal base, so both tools exclude the same positions.

> **Trap**
>
> `SeedAndFilter()` copies a whole chunk's seeds to the GPU in one `cudaMemcpy`, into a buffer sized `MAX_SEEDS`. That size used to be a hardcoded `13 * wga_chunk_size`.
>
> 13 is `1 + 12` — the default `12of19` seed's weight plus one. It fits that one pattern _exactly_, with zero margin, and nothing else. `14of22` emits 15 per position and overran it:
>
> ```
> Error: cudaMemcpy of 26650680 bytes for seed_offsets failed
>        with error " invalid argument "
> ```
>
> 26650680 / 8 = **3,331,335** seeds into a 3,250,000-seed buffer. That is **13.325 seeds per position**, which sits strictly between 13 and 15 — arithmetically impossible under `12of19`, which cannot exceed 13, and expected under `14of22`. The ratio identifies the cause, not merely the symptom.
>
> The bound is now `MaxSeedsPerChunk()` in `common/seed_capacity.h`. It has two arms: under `--notransition` the seeder pushes one seed per valid position and the bound is just `wga_chunk_size`; otherwise it is `(1 + num_transitions) * wga_chunk_size`. Either way it is _exact_, not conservative — a chunk covers at most `wga_chunk_size` positions — and it evaluates to the old constant for `12of19`, so the default path is unchanged.
>
> It is driven by the **transition count**, not the weight. Those are equal for both built-in patterns and for any normalised custom one, but a shape containing a literal `1` has weight without the matching transition, and the capacity must follow what the seeder actually pushes. Note that the CLI cannot produce such a shape — `main.cpp:203-208` rewrites every `1` to `T` before `GenerateShapePos()` sees it — so the only caller that exercises the distinction is the unit test, which drives the function directly. The distinction is still worth preserving: it is what makes the bound follow the seeder rather than the weight.

### `--step` applies to the reference only, and that matches LASTZ

`cfg.step` reaches `GenerateSeedPosTable()` (`main.cpp:612`), which builds the reference seed-position table. The query loop in `seeder.cpp` advances one position at a time regardless. That asymmetry looks like a bug and is not one.

LASTZ does the same: `currParams->step` is passed only to `build_seed_position_table()` for the _target_ (`lastz.c:1203`, `:1214`, `:1229`, `:1305`, `:1311`). In `seed_search.c` the value appears only as `pt->step`, read back out of the table to reconstruct target coordinates (`pos1 = adjStart + step*pos`) — never to advance the query scan.

> **How to check**
>
> **Measured, and the parameter is live.** A parity run at `--step=2` on 2026-09-02 (laila, 4× A100) agreed exactly, and the counts moved:
>
> ```
> step=1   unmasked 1225 = 1225      masked 469 = 469
> step=2   unmasked 1207 = 1207      masked 461 = 461
> ```
>
> The drop matters as much as the agreement. Had KegAlign been ignoring the flag on a side LASTZ honours, its counts would have stayed at 1225/469 while LASTZ fell to 1207/461. Both moved, together.

### Why the guard is not an assert

There was an `assert(num_seeds <= MAX_SEEDS)` in exactly the right place. conda-forge builds with `-DCMAKE_BUILD_TYPE=Release`, which defines `NDEBUG`, which deletes it. Every released build had no guard at all, so the overrun surfaced as a bare `" invalid argument "` from `cudaMemcpy` hundreds of lines from its cause. It is now a real runtime check that survives `NDEBUG`.

### Why weight is capped at 15, and floored at 4

`MAX_SEED_WEIGHT` is 15, and the tighter-than-obvious bound is the load-bearing one. `GetKmerIndexAtPos()` packs two bits per match position into a `uint32_t`, and `INVALID_KMER` is `1<<31`. At weight 16 a legitimate k-mer can equal that sentinel and be silently discarded as invalid; above 16 it wraps outright. Staying at or below 15 also keeps `shape_pos[32]` and `transition_pos[32]` in bounds, and those matter more than they look. Each array sits immediately before a scalar it overruns into: `shape_pos[32]` aliases `shape_size` — the weight itself — and `transition_pos[32]` aliases `num_transitions`, the count `MaxSeedsPerChunk()` sizes the seed buffer from. An over-weight pattern corrupts the very numbers used to reject it and to size the buffer. Measured, not reasoned: on the pre-fix code `GenerateShapePos()` given forty `T`s returns **39**, because writing `shape_pos[32]` overwrote the loop counter mid-loop.

> **How to check**
>
> **How to check any of this.** Links the real `ntcoding.cpp`, so it exercises the actual shape parser rather than a copy of it, and pins the capacity against a loop mirroring the seeder's push sequence. Needs no GPU.
>
> ```bash
> bash tests/test_seed_capacity.bash
> ```

### The floor is a deleted assert

`MIN_SEED_WEIGHT` is 4 because `GenerateSeedPosTable()` asserts `kmer_size > 3`. That assert is deleted by `NDEBUG` in a Release build — the same way the seed-buffer assert was, which is how the `14of22` overflow reached a user as a bare `cudaMemcpy` failure. A weight-3 pattern was observed dying on a GPU with `thrust::system::system_error: device free failed: cudaErrorIllegalAddress`.

Be careful what you attribute that crash to. It is **not** the index table: at weights 1, 2 and 3 the table holds 5, 17 and 65 entries and every k-mer the seeder can produce indexes inside it. Enumerated against the real `ntcoding.cpp`. The bound is the upstream precondition, enforced; the mechanism of the crash is not established here.

### The index table quadruples with the weight

`index_table_size = (1 << 2*kmer_size) + 1` `uint32`s, allocated per GPU and per reference block:

| weight | index table |
| --- | --- |
| 12 (default `12of19`) | 64 MiB |
| 14 (`14of22`) | 1 GiB |
| 15 (the legal maximum) | **4 GiB** |

`MAX_SEED_WEIGHT` is a correctness bound, not a resource one. Nothing warns that the top of the legal range costs 4 GiB of device memory — plus the same again on the host, since the table is `calloc`'d before it is copied.

### Known problem — one position is seeded twice at interval boundaries

The chunk loop uses `e = min(i + wga_chunk_size, q_inter_end + 1)`, so it seeds position `q_inter_end` — which is the _first_ position of the next interval. That `+ 1` is not a mistake to delete: intervals are built half-open up to `seq_block_len - seed.size`, and the last valid seed position is that endpoint itself, so without the `+ 1` every block would lose its final seed. The overlap is real at _interior_ boundaries only, and produces duplicate HSPs in two different segments files. A correct fix widens the interval range and drops the `+ 1` together; doing only the latter is a regression.

## 3. The scoring matrix

### What it is

Extension scores a reference character against a query character. Both have already been compressed to one of eight codes (`common/parameters.h:25-32`):

| code | value | what reaches it |
|---|---|---|
| `A_NT C_NT G_NT T_NT` | 0–3 | upper-case `A C G T` |
| `L_NT` | 4 | lower-case `a c g t` — **all four collapsed into one** |
| `N_NT` | 5 | `n` and `N` |
| `X_NT` | 6 | everything else — every IUPAC ambiguity code, and any junk byte |
| `E_NT` | 7 | `&`, the separator between sequence blocks |

So the matrix is 8×8 = 64 ints (`NUC2`), living in `cfg.sub_mat` and copied to shared memory by `find_hsps`.

A LASTZ scoring file supplies only the 4×4 ACGT block. Everything else is derived. `build_substitution_matrix()` in `common/scoring.c` does that derivation, and it is the whole of the policy:

```
ACGT x ACGT              the 4x4 from the file, or HOXD70 by default
L_NT  vs ACGT, L_NT      bad_score            (-1000)
N_NT  vs ACGT, L_NT      bad_score  -- or ambiguous_penalty under --ambiguous
N_NT  vs N_NT            bad_score  -- or ambiguous_reward  under --ambiguous
X_NT  vs ACGT            fill_score           (-100)
X_NT  vs L_NT, N_NT      bad_score
X_NT  vs X_NT            fill_score
E_NT  vs anything        -10 * xdrop          (a wall, see below)
```

Read the rows as ranges, not as “anything”: they are written by loops that overwrite each other in order, so the last one to touch a cell wins. `L_NT` against `E_NT` is `-10 * xdrop`, not `bad_score`, because the separator loop runs last.

And `--ambiguous` takes _two_ values, which the table above only half-covers. `--ambiguous=n` rewrites the `N_NT` row as shown. `--ambiguous=iupac` additionally rewrites the whole `X_NT` row — every ambiguity code against a base becomes `ambiguous_penalty`, and `X_NT` against `X_NT` becomes `ambiguous_reward`. That second flag matters well beyond this section: it is what turns an aligned pair of ambiguity codes from something merely absorbed inside an HSP into something that actively extends one, which is the strongest trigger for the counter bug in §5.

### Why the separator is −10 × xdrop and not merely bad

Blocks are concatenated with `&` between them. An HSP must never span two blocks, so the separator's score has to be large enough that x-drop terminates on it _unconditionally_, whatever the peak score was. `bad_score` is a fixed −1000 and could in principle be climbed past by a long enough high-scoring run; `-10 * xdrop` is defined relative to the very threshold that decides termination, so it cannot be. This is a wall, not a penalty.

> **Trap**
>
> This is the one thing in this file most likely to mislead you, because the misleading evidence is real code sitting in this repository.
>
> `common/dna_utilities.c:1222-1234` — vendored, genuine LASTZ — copies the upper-case rows onto the lower-case characters:
>
> ```c
> // if the rows are DNA, make lower case rows equivalent to upper case
> memcpy (/*to*/ xss->ss.sub[r+'a'-'A'], /*from*/ xss->ss.sub[r], sizeof(scorerow));
> ```
>
> Read that alone and you will conclude that LASTZ scores a soft-masked base as its upper-case equivalent, that KegAlign's `L_NT = bad_score` diverges from it, and that the divergence should be closed. **All three conclusions are wrong.**
>
> Seven hundred lines earlier in the same file, `masked_score_set()` (`common/dna_utilities.c:518`, doc comment: _“Create a copy of a score set, with all lower case entries given ‘bad’ scores”_) builds a second matrix, stamping `badScore` across every non-upper-case row and column **and over `'n'` and `'X'`** — and over `'N'` too, but only when `N` is not itself a row character (`if (!nIsARow)`, `dna_utilities.c:561`). With the stock `ACGT` matrix it is not, so the parity argument below holds; a scoring file that gives `N` a row of its own changes that, and it is the same axis the `N_NT` claim runs along. LASTZ keeps both, and its own header says which is used where:
>
> ```c
> // lastz.h:201-205
> scoreset* maskedScoring;  // scoring set with lowercase penalized;  in
>                           //   general, we treat lower case as bad
>                           //   during the search for HSPs, then treat
>                           //   upper/lower case as equivalent during
>                           //   the gapped alignment stage
> ```
>
> The wiring, on the default `gfexXDrop` path:
>
> ```c
> lastz.c:2934        hpInfo->scoring = params->maskedScoring;   // gap-free extension
> seed_search.c:2539  scoreset* scoring = hp->scoring;           // xdrop_extend_seed_hit
> lastz.c:3421        gapped_extend(... params->scoring ...)     // the OTHER one
> ```
>
> Gap-free extension — the stage `find_hsps` reimplements — uses the **masked** matrix. KegAlign's `L_NT = bad_score` and `N_NT = bad_score` therefore match LASTZ's ungapped stage as read.
>
> Two honest limits on that. It is established _by reading LASTZ_, and the measurement in §7 only corroborates the `L_NT` half — the test data is pure upper-case `ACGT` with no `N` and no ambiguity codes, so `N_NT` and `X_NT` are never constructed in either parity case. And "always were" is a claim about history, not about the run: nothing here re-measures older revisions.
>
> `lastz.c`, `lastz.h` and `seed_search.c` are not vendored here; fetch them from `raw.githubusercontent.com/lastz/lastz/master/src/` if you need to re-check.
>
> > **This trap was not hypothetical.** In September 2026 a `--lastz_masking` flag was written, compiled, tested and very nearly shipped on the strength of reading the first block and not the second. It would have _removed_ parity behind a flag named for it. It was caught by an adversarial review and retracted; the patch is preserved at `bench/data/kegalign-stage4-retracted/stage4.patch`. The differential test in §7 would have refuted it in five minutes.

### What a scoring file actually supplies

`--scoring` takes a LASTZ scoring file and reads it with LASTZ's own parser (`common/scoring.c:load_scoring_matrix`). Three settings are read, and one that you might expect is not:

| | |
|---|---|
| **the 4×4 ACGT matrix** | read directly |
| `bad_score` | recovered from `sub[badRow][badCol]`, since LASTZ stores it by writing it across the bad character's row and column rather than keeping the scalar |
| `fill_score` | recovered from any pair the reader never wrote over, because LASTZ initialises all 256×256 entries to it before laying the matrix, the lower-case copies and the bad row/column on top |
| `gap_open_penalty` `gap_extend_penalty` | **not used here.** They belong to gapped extension, which is lastz's job, and the same file is handed to it via `--scores=` |

Before September 2026 only the 4×4 was used and the other two were silently discarded, so a file's `bad_score` and `fill_score` applied in the gapped stage but not the ungapped one — one file, two interpretations, inside one pipeline.

> **How to check**
>
> **How to check any of this.** Drives the real LASTZ reader through a scoring file whose every value differs from the built-in defaults, and asserts the values reach the matrix rather than only the reader. Needs no GPU.
>
> ```bash
> bash tests/test_scoring.bash
> ```

### Known problem — literal X

LASTZ gives the literal character `X` `badScore` in **both** its matrices; KegAlign folds `X` into the `X_NT` catch-all at `fill_score`. Since `X_NT` also holds the IUPAC ambiguity codes, which LASTZ _does_ give `fillScore`, the two cases are conflated and cannot both be right with the alphabet as it stands. Rare in genomic FASTA, and unaddressed.

## 4. Masking

### Two stages, two mechanisms, one outcome

A soft-masked base is handled twice, by unrelated code, and it is worth keeping them apart.

| stage | where | mechanism |
|---|---|---|
| **seeding** | host, over the raw sequence | `NtChar2Int()` sends everything outside `ACGTN` — lower case included — to `N_NT`, and `GetKmerIndexAtPos()` then returns `INVALID_KMER` if any base in the whole span is `N_NT`. No seed is produced. |
| **extension** | GPU, over the compressed sequence | lower case compresses to `L_NT`, which scores `bad_score` = −1000 against everything. At the default `xdrop` of 910, one masked base ends the extension. |

The two stages read _different representations_ — seeding walks the host DRAM buffer as raw ASCII, extension reads the GPU's compressed codes — and that separation matters: it is what would let the two policies be changed independently, if there were ever a reason to.

Case survives the whole path. The GPU stage reads the FASTA directly through `kseq`, not the `.2bit` built for the downstream lastz call, and nothing upper-cases on load. `RevComp()`'s table preserves case, so the minus strand is masked identically.

### Both match LASTZ

LASTZ excludes lower-case words from its seed table, and its gap-free extension uses `maskedScoring`, in which lower case, `N`, `n` and `X` are all `badScore`. The evidence, and the trap in reading it, are in §3. **KegAlign matches LASTZ on both stages**, established by reading LASTZ and corroborated in §7: 469 HSPs against 469, byte-identical, with half the bases masked.

How much that measurement carries is worth stating plainly, because it is easy to spend twice. The parity test masks target and query at _identical_ coordinates, so a seeding-only divergence — LASTZ seeding a lower-case word that KegAlign refuses — has nowhere to show itself: a seed inside an aligned masked pair scores `bad_score` and x-drops, producing no HSP either way. So 469 = 469 confirms the composite outcome of the two stages, not each of them independently. The per-stage claim rests on the source reading above.

### What masking costs, measured

It is not a small effect. On a real WindowMasker-softmasked pair (4× A100, `--seed 14of22 --strand both --num_threads 32 --debug`), the same run with the query upper-cased:

|  | masked | upper-cased | ratio |
|---|---|---|---|
| `#seeds` | 15,269,532,690 | 44,193,436,080 | 2.89× |
| `#seed hits` | 22,746,550,414 | 89,723,031,347 | 3.94× |
| `#HSPs` | 212,871,696 | 403,443,154 | 1.90× |
| pipeline | 56 sec | 117 sec | 2.09× |

WindowMasker was excluding **65.4%** of that query from seeding. Un-masking roughly triples the seeds, quadruples the hits, and doubles the HSPs — and the secondary ratios show why: hits per seed rise from 1.49 to 2.03 while HSPs per hit halve, which is what recovering repetitive sequence looks like.

> **Read that table as the cost of masking, not the cost of a bug.** An earlier reading of it treated the 1.90× as the price of a KegAlign divergence from LASTZ. There is no divergence — LASTZ pays the same price on the same input, because it excludes and penalises the same bases. The counters are real; the interpretation was wrong. `--debug` is required to print them at all.

### The connection to the seed buffer

Note that both the table above and the arithmetic below are `14of22` figures — 15 seeds per valid position, not the default `12of19`'s 13. Re-derive them if you change the pattern.

Masking is also why the `14of22` overflow in §2 took so long to appear. A chunk can only overrun the seed buffer if more than **86.7%** of its 250,000 positions yield a valid k-mer, and masked positions yield none. On a repeat-masked plant genome most chunks fall well under that; the one that finally tripped it was 88.8% clean. A crash that needs unmasked sequence, on input that is mostly masked, is a crash that arrives thousands of chunks in.

### Known problem — there is no way to turn it off

LASTZ exposes masking as a choice: the `[unmask]` sequence action converts lower case to upper case and puts those bases back into seeding. KegAlign has no equivalent and no CLI option. The only lever is to upper-case the FASTA before it goes in — which is a divergence _from_ LASTZ, not toward it, and worth being deliberate about.

Note also that `L_NT` collapses all four lower-case bases into one code, so once a base is compressed its identity is gone. Any future option that wanted masked bases scored as their upper-case equivalent would have to fold case in the compression kernels, not in the substitution matrix — the matrix cannot represent it.

## 5. Entropy

### What it is for

A run of low-complexity sequence can clear `hspthresh` on sheer length while carrying almost no information. The entropy adjustment demotes those. It is applied in one narrow band (`src/seed_filter.cu`):

```
if(total_score >= hspthresh && total_score <= 3*hspthresh && !noentropy)
```

At the defaults that is `[3000, 9000]`. Outside it nothing happens — an HSP scoring 28941 is never entropy-adjusted, which is worth remembering when attributing a surprising HSP to this machinery.

Inside the band, per-base match counts feed Shannon entropy over the four bases, normalised so a uniform composition gives 1.0:

```
entropy = -SUM( p[i] * ln p[i] ) / ln 4,   p[i] = count[i] / (extent+1)

d_hsp[hid].score = (int)(total_score * entropy)
```

Two consequences follow. The adjusted score is what gets _written into the segments file_ and handed to lastz, not merely a filter input. And the HSP is dropped when `total_score * entropy` falls back below `hspthresh` — so a wrong entropy changes both which HSPs survive and what score the survivors carry.

There is a floor: if fewer than 20 matched bases were counted, entropy is 1.0 and the adjustment is a no-op.

### The LASTZ ancestry

This is a port of `compute_entropy()` in `common/dna_utilities.c` — same formula, same `< 20` gate. LASTZ's version counts blindly and reads selectively:

```c
int count[256];
count['A'] = count['C'] = count['G'] = count['T'] = 0;
for (ix=0; ix<len; ix++)
    { if (s[ix] == t[ix]) count[s[ix]]++; }        // any byte, harmlessly
cA = count['A']; cC = count['C']; cG = count['G']; cT = count['T'];
```

A match of two `N`s increments `count['N']`, which is never read. The 256 entries exist precisely so that costs nothing.

> **Trap**
>
> KegAlign kept the blind increment and shrank the array:
>
> ```c
> short count[4];
> short count_del[4];
> ...
> if(r_chr == q_chr){          // r_chr is a compressed code in 0..7
>     if(...) count[r_chr] += 1;
>     else    count_del[r_chr] += 1;
> }
> ```
>
> `L_NT`=4, `N_NT`=5, `X_NT`=6 and `E_NT`=7 all index past the end. Indexing a local array out of range is undefined behaviour, so what those writes hit is a property of one build, not of the code. What can be said without disassembling anything is enough: the two arrays are adjacent locals in the same frame, `count_del` _is_ read back a few lines later by `count[i] = count[i] + count_del[i]`, and every landing site inside that frame is either a live counter or another live variable. So a non-nucleotide match either inflates a real base count or corrupts something else — and which of those it does can change with a compiler or toolkit bump.
>
> A specific layout — `count` at +0 and `count_del` at +8, making `count[4..7]` exactly `count_del[0..3]` — was reported for one build, and this document previously stated it as fact. It is not reproducible in the environment these notes were written in, so treat it as a plausible instance rather than the mechanism. The case for the fix does not need it.
>
> It is reachable with **stock flags and no masking**. Every IUPAC ambiguity code compresses to `X_NT`, so any two of them in one column compare equal, and `X_NT` against `X_NT` scores `fill_score` = −100 by default — absorbed inside an HSP rather than triggering the 910 x-drop.
>
> The guard is now `r_chr == q_chr && r_chr < L_NT`, which is memory-safe and follows LASTZ in counting only real bases — but _not_ provably identical to it. LASTZ discards the extra increments rather than excluding them, and its `entropy_lower_ok()` variant folds lower case back into `A/C/G/T`, which KegAlign structurally cannot do because `compress_string` has already collapsed all four into `L_NT`. That gap is real and unclosable without a wider alphabet. `r_chr` and `q_chr` are also initialised per tile to deliberately _unequal_ values, because they were assigned only inside the in-range guard yet read outside it — so a lane past the end of either sequence read them uninitialised on the first tile and stale afterwards.

> **Soft-masking cannot exercise this, and a masked test that comes back clean proves nothing about it.** Masked columns score `bad_score` against an xdrop of 910, so they x-drop in the same tile they are scanned, always take the `count_del` arm, and land in the half whose corruption falls outside the frame and is never read back. A differential run with 50% of bases masked was byte-identical to LASTZ — exactly what you would see whether or not the bug mattered. IUPAC codes, or `--ambiguous=n`, are what reach it.

### Where parity stops being provable

LASTZ has two entropy functions. Its HSP filter calls `entropy()`, which ignores lower case; `entropy_lower_ok()` folds `count['a'..'t']` back into `A/C/G/T`. KegAlign cannot offer the second, because `compress_string` collapses all four lower-case bases into the single code `L_NT` and the identity of the base is gone by then. For the filter LASTZ actually uses, the behaviours agree.

### Known problem — the counters are `short`

The warp scan sums 32 lanes' per-base counts into a `short`. A single-base count above 32767 wraps negative, `log()` of a negative yields NaN, and `(int)(total_score * NaN)` drops the HSP. It needs an HSP longer than ~32k bases that still scores at or below `3*hspthresh` — a very long, near-break-even alignment. No convincing genomic instance has been constructed, so this is recorded rather than fixed. Widen to `int` if these arrays are touched again.

## 6. HSPs and segments

### What actually crosses the boundary

Everything the GPU produces is a list of these (`src/graph.h:23-28`):

```c
typedef struct segmentPair {
    uint32_t ref_start;
    uint32_t query_start;
    uint32_t len;
    int      score;
} segmentPair;
```

Sixteen bytes, and _one_ length for both sequences — which is the structural signature of a gap-free stage. A gapped alignment needs two extents and an edit script; an HSP is a diagonal run, so `ref_start − query_start` is constant and one number describes both sides. Several things downstream depend on that: the diagonal is the sort key in `hspComp`, and it is what makes containment testable in a single comparison.

> **Trap**
>
> **`len` is an extent, not a count of bases.** The segment covers `len + 1` bases. Both places that consume it agree, and both look like off-by-one errors until you know that:
>
> ```c
> // src/seed_filter.cu — entropy denominator (see §5)
> p = count[i] / (extent + 1);
>
> // src/segment_printer.cpp:90 — printed end coordinate
> std::to_string(seg_r_start + e.len + 1 - r_chr_start[r_index])
> ```
>
> That is the same `+1` twice, and it is correct: LASTZ's `general` format uses origin-one inclusive coordinates, so `end − start + 1` must be the base count. Confirmed by measurement, not by reading — the parity test compares the printed columns against LASTZ's own with no shift applied and gets 1225 of 1225 byte-identical.
>
> If you ever change how `extent` accumulates, both sites move together or neither does.

### Three coordinate systems

A position is expressed three different ways on its path out, and the printer converts between all three in a single expression:

| system | origin | who uses it |
|---|---|---|
| block-local | 0 | the GPU — `segmentPair` holds these |
| concatenated-genome | 0 | the host, after `+ r_block_start` / `+ q_block_start` |
| chromosome-local | 1 | the `.segments` file, after an `upper_bound` lookup and `+1 − chr_start` |

The chromosome lookup is a binary search over `r_chr_start` / `q_chr_start`. For the query the printer caches the current chromosome and only re-searches when a segment leaves it (`segment_printer.cpp:81-88`) — HSPs arrive sorted by query position, so that cache almost always hits. The reference gets no such cache, because the reference-side order is not monotonic.

### The minus strand is walked backwards, on purpose

```c
for(int r = rc_hsps.size()-1; r >= 0; r--){    // segment_printer.cpp:132
```

This is not a mistake and does not need fixing. The reverse complement has its own DRAM buffer and its own chromosome table, reflected within the block:

```c
rc_q_chr_start = 2*seq_block_start + seq_block_len
                   - q_chr_start[c] - q_chr_len[c];   // main.cpp:359
```

Coordinates on that buffer run opposite to the forward query, so the ascending `hspCompLastz` order the GPU produced is _descending_ in the order LASTZ wants to read. Iterating backwards puts it right. The forward loop and this one are otherwise the same code.

### Ordering is a contract, not a convenience

Three sorts run inside `SeedAndFilter()`, and each has a distinct job:

| # | comparator | key | why |
|---|---|---|---|
| 1 | `hspComp` | diagonal, then `ref_start`, `len`, `score`↓ | puts same-diagonal HSPs adjacent, which is the precondition for step 2 |
| 2 | `hspEqual` via `unique_copy` | same diagonal **and** one contained in the other | collapses redundant HSPs — note this is _containment_, not equality |
| 3 | `hspCompLastz` | `query_start`, then `ref_start`, `len`, `score`↓ | the order LASTZ expects in a `--segments=` file |

The third sort exists only to satisfy the downstream reader; its name says so, and that is why `drop_duplicate_hsps()` sorts a **permutation** and compacts in place rather than sorting the HSPs — re-sorting would silently change what LASTZ receives.

Be careful how far you take that, though. Because the sort runs per batch and the batches are concatenated (see the trap below), **the emitted file is a sequence of locally-sorted runs, not a globally sorted file**. The order is a convention the code takes trouble to preserve rather than a guarantee it establishes — and nothing tests it: the parity test pipes both sides through `sort` before comparing, so a change that scrambled the emitted order would leave it green.

> **Trap**
>
> **Deduplication is per `MAX_HITS`-bounded iteration, not per chunk and not per file.** This is finer-grained than it looks, and getting it wrong understates the problem. All three sorts run _inside_ a loop:
>
> ```c
> int num_iter = num_hits/MAX_HITS + 1;        // seed_filter.cu:764
> for(int i = 0; i < num_iter; i++){
>     ... stable_sort(hspComp) ... unique_copy(hspEqual) ... stable_sort(hspCompLastz)
> }
> // then, seed_filter.cu:845 -- the per-iteration results are concatenated
> ```
>
> So duplicates already survive _within a single chunk_ whenever `num_iter > 1`, and the emitted order is a concatenation of locally-sorted runs rather than one sorted file. `seeder.cpp` then concatenates chunks on top of that with a plain `insert`, and `segment_printer` writes what it is handed:
>
> ```c
> fw_hsps.insert(fw_hsps.end(), anchors.begin()+1, anchors.end());
> ```
>
> An HSP reachable from seeds in two different chunks passes both `unique_copy` calls — they never see each other — and is written twice. LASTZ then repeats the same gapped extension and reports the alignment twice.
>
> Found by the parity test, not by reading: 1226 segments against LASTZ's 1225, with all 1225 matching. The arithmetic fits — the test query is ~376 kb, two chunks at the default 250,000, and exactly one HSP landed on the seam.
>
> `drop_duplicate_hsps()` now runs once the chunks are joined. It matches on the exact four-field key and **deliberately does not reproduce `hspEqual()`**: applying containment at this level could drop an HSP that LASTZ, working chunk-unaware, would have reported. It also decrements `num_hsps`, which is incremented per chunk, so the `--debug` count keeps matching what the files hold.

> **Trap**
>
> **The partition depends on the GPU's memory size, and so does the output.**
>
> ```c
> cudaGetDeviceProperties(&deviceProp, 0);
> float global_mem_gb = deviceProp.totalGlobalMem / 1073741824.0f;
> MAX_HITS = MAX_HITS_PER_GB * global_mem_gb;      // seed_filter.cu:863-868
> ```
>
> `MAX_HITS` sets `num_iter`, which sets which HSPs are ever compared with each other. And step 2 is `hspEqual`, a _containment_ test that deletes a real HSP rather than a copy of one. Put those together: **two cards with different memory can emit different HSP sets**, not merely different numbers of duplicates.
>
> Nothing in the codebase records this and no test covers it. It also qualifies every measured claim in this document — the parity numbers in §7 were taken on 4× A100 and are not guaranteed to reproduce on a smaller card.

> **How to check**
>
> The parity test counts duplicates **per file**, and the comment there explains why that scope and not the concatenation:
>
> ```bash
> for f in "$dir"/*.segments; do
>     # total - distinct, not `uniq -d | wc -l`: that counts duplicated GROUPS,
>     # so a segment emitted three times would report as one.
>     dup=$(( dup + $(wc -l < "$f") - $(sort -u "$f" | wc -l) ))
> done
> ```
>
> Each interval writes its own `.segments` consumed by its own LASTZ call, so counting across files would flag cross-interval duplicates the fix does not claim to remove.

### Known problem — the seam between intervals

Intervals are separate `seeder_body` invocations writing separate files, so the deduplication above cannot reach across them. In practice this is the smaller exposure: at the defaults an interval holds 40 chunks, against a single boundary with its neighbour.

There is a one-position overlap at that boundary, and it is **load-bearing**:

```c
uint32_t e = std::min(i + cfg.wga_chunk_size, q_inter_end + 1);   // seeder.cpp, both strand loops
```

Read alone, that `+1` is an off-by-one that makes adjacent intervals both seed their shared endpoint. But intervals are laid out to `end_pos = seq_block_len - cfg.seed.size` (`main.cpp:372`), and a span of `seed.size` fits exactly at that position — so it is the last _valid_ seed start in the block, and without the `+1` every block loses its final seed. A code review proposed removing it as the root cause of cross-interval duplicates; it was checked against the interval arithmetic and rejected. Documented here rather than fixed.

The right fix, if it is ever worth making, is to make the interval bounds half-open everywhere and extend only the last one — not to delete the `+1`.

## 7. Checking parity against LASTZ

### Why this test is different from the others

Every other test in `tests/` checks a piece of KegAlign against its own arithmetic: does the seed buffer's capacity match what the seeder pushes, does a scoring file's values reach the matrix, does the CUDA architecture list come out right. Those catch a component contradicting itself.

None of them can catch KegAlign faithfully implementing the wrong thing. `tests/test_lastz_parity.bash` is the one that can: it runs KegAlign and LASTZ over the same pair with the same parameters, asks each for its HSPs in the same eight columns, and diffs them.

### How the two are made comparable

`--nogapped` is the lever. It stops LASTZ after gap-free extension — precisely the stage `find_hsps` replaces on the GPU — so both tools are answering the same question rather than adjacent ones.

The output formats already agree. KegAlign writes LASTZ's own `--segments` _input_ format, and LASTZ can be asked for exactly those fields:

```
lastz "target.fa[multiple]" query.fa \
  --nogapped --strand=both \
  --seed=12of19 --xdrop=910 --hspthresh=3000 --step=1 \
  --format=general:name1,start1,end1,name2,start2,end2,strand2,score

kegalign target.fa query.fa work/ --seed 12of19 --strand both    # -> *.segments
```

Two field-name details that are not guessable. The column is `strand2`, not `strand` — LASTZ rejects the latter outright. And `start1`/`start2` are 1-based closed, matching what KegAlign writes, while `zstart1`/`zstart2` are 0-based half-open. Run `lastz … --format=general` with no field list and it prints its own header, which is authoritative.

### What it found

| case | LASTZ | KegAlign | identical | duplicates |
|---|---|---|---|---|
| unmasked | 1225 | 1226 | 1225 | **1** |
| soft-masked 50% | 469 | 469 | 469 | 0 |

Zero missed and zero spurious HSPs in either case, coordinates and scores byte-identical. **KegAlign's gap-free extension agrees with LASTZ exactly.** The single surplus segment was a duplicate, not an extra HSP — which is how the cross-chunk deduplication gap in §6 was found.

It also killed two plausible hypotheses. There is no systematic coordinate off-by-one: the columns matched as printed, with no shift applied, so `start + len` is right and `len` is already an offset rather than a base count. (The script does compute a −1-shifted comparison, but only prints it when the shift matches _more_ lines than the unshifted one — so a clean run bounds that count without reporting it. The unshifted match is the evidence; the shift is a diagnostic for when it fails.) And the duplicate was not an entropy artefact: at score 28941 it sits far outside the `[3000, 9000]` band where the adjustment applies at all.

> **How to check**
>
> **Running it.** Two cases, the second soft-masking both inputs at identical forward coordinates. Needs a GPU and both tools on `PATH`; exits 77 (skip) when it cannot run.
>
> ```bash
> bash tests/test_lastz_parity.bash [target.fa query.fa]
> ```

### What it cannot exercise

A green run is not a blanket guarantee, and three limits are worth knowing before citing one as evidence.

**Soft-masking does not reach the entropy counters.** Masked columns score `bad_score` against the default xdrop, so they terminate extension in the tile they are scanned and never land inside a counted extent. The masked case is meaningful for seeding and for scoring, and proves nothing about §5. IUPAC codes, or `--ambiguous=n`, are what reach that path — and the test does not yet carry them.

**Multi-block inputs are untested, but should _not_ disagree.** An earlier draft of this section claimed KegAlign's `&` separators would make genome-scale runs differ from LASTZ at every seam. That is wrong, and worth correcting explicitly because it would have licensed dismissing a real defect as expected noise. Both the block splits and the separators are only ever placed _between whole FASTA records_ (`main.cpp:350-402`: the size test and the `memset` both run after a complete `kseq` record is appended), and LASTZ does not extend an alignment across a record boundary either. The default pair is ~372 kb and ~377 kb against a 500 Mb `seq_block_size`, so it is a single block and the case has simply never been exercised. Treat a genome-scale disagreement as a finding, not as expected.

**It cannot run in CI.** Everything else in `tests/` is deliberately host-side and GPU-free so the workflow can run it. This one needs a GPU, so it is a manual pre-release check — which means it only protects you if someone remembers to run it.

> **Why this test exists.** In September 2026 a masking feature was designed, written, compiled and tested on the strength of reading one block of a scoring-file reader — and was wrong, because a second matrix seven hundred lines earlier overrides it for exactly the stage in question (§3). Every host-side test passed throughout; none of them could have known. This test would have refuted the premise in five minutes. _Reading half of a scoring-file reader is not evidence about behaviour._ When a parity question comes up, measure it here first.

## 8. The orchestration layer

The GPU binary is only the middle of the pipeline. `scripts/` drives the rest, and every defect found in the September 2026 review after the release itself was in here rather than in the CUDA. It is worth a section because none of it is obvious from reading `main.cpp`.

### The shape of a run

```
run_kegalign  ->  kegalign (the CUDA binary)     writes tmp*.segments, prints lastz commands
              ->  diagonal_partition.py          splits large .segments files
              ->  lastz --segments=...           gapped extension, one process per command
              ->  package_output.py              tars the whole thing for Galaxy
```

`runner.py` is the conductor. Three `multiprocessing.Manager()` queues carry work between the stages, and each stage ends when it reads as many `SENTINEL_VALUE`s as there are workers — so the sentinel count and the worker count must agree, and both come from `--num_cpu`.

The queues are **manager proxies**, not `queue.Queue` objects, which is why they can cross a `ProcessPoolExecutor` boundary at all. Anything that changes how workers are spawned has to keep that true.

### Ordering is inherited from SegAlign, and was silently lost

The shell script this replaced ran lastz in parallel and then concatenated deterministically:

```bash
for i in tmp*.plus.*;  do echo $i; done | sort -V
for i in tmp*.minus.*; do echo $i; done | sort -V
```

Plus before minus, version-ordered within each. The Python port dropped it: both output paths — the MAF concatenation and the lastz command file that becomes `galaxy/commands.json` — used the order work happened to finish in.

That went unnoticed for a long time because a second bug hid it. `executor.submit(worker(...))` **called** the worker instead of passing it, so every partitioner ran serially in the parent process and completion order was accidentally input order. Fixing the parallelism is what made the missing sort visible.

`KegAlignSegment.__lt__` had encoded the intended order the whole time — `["strand", "tmp", "block", "r", "split"]`, `plus=0`, `minus=1` — and was never called. Its only would-be caller, `KegAlignSegments.__iter__`, returned `self` alongside a `__next__` that was itself a generator function, so iterating it looped forever yielding generator objects. That reads as a half-finished port of `sort -V`: comparison written, iteration broken, call site never wired.

Both paths now go through it, via `sorted_commands()` and `lastz_command_sort_key()`. Nothing is buffered that was not already being read in full; only the final write is ordered.

> **Trap**
>
> There are **two** output paths and they are easy to confuse. `sorted_commands()` orders the MAF concatenation; `lastz_command_sort_key()` orders the command file. Fixing one and declaring victory is exactly what happened between v0.3.1's first and second attempts — the MAF was ordered while `commands.json`, which is what Galaxy hands downstream, was not. If you change the ordering rule, change it in both.

### Chunk size is estimated twice, by two different copies

`runner.py:estimate_chunk_size()` and `diagonal_partition.py` each estimate independently, with near-identical code. Both feed `statistics.quantiles()`, which needs **two** data points on Python 3.10–3.12 and one on 3.13+; the conda-forge `kegalign` recipe pins 3.12.

Neither guarded it adequately. `runner.py` checked `< 7` and `diagonal_partition.py` checked `len(files) < 2` — but the latter keys by *prefix*, so the two split files `DELETE_AFTER_CHUNKING` produces are `len(files) == 2` with `len(fdict) == 1`. Either reaches the quantile with too few points and dies with an unhandled `StatisticsError`, after all the GPU work is finished.

> **Trap**
>
> `chunk_size == 0` is **meaningful**, not an error. `diagonal_partition.py` documents *"set `<max-segments>` = 0 to skip partitioning"* and implements it by printing the command unchanged. Do not clamp the estimate's lower bound — flooring it at 1 emits one output file and one LASTZ command per segment line, for every file in the run, because the estimate is computed once and applied to all of them.

### Errors from workers are easy to lose

A worker that fails calls `sys.exit()`. While the partitioners were accidentally running in the parent, that ended the run. In a child process it does not: `SystemExit` is marshalled back into the future, and if nothing collects the futures it is simply discarded — the pool shuts down and the function reports success.

Both pools now collect with `as_completed(...).result()`. Note that `SystemExit` is a `BaseException`, so a surrounding `except Exception` will still not catch it; the process exits non-zero carrying the worker's own message, which is the behaviour we want but not via the path the code appears to take.

### How to check any of this

```
python3 tests/test_chunk_size.py     # imports the real estimator, drives it on temp dirs
python3 tests/test_output_order.py   # shuffles arrival order, requires a stable result
```

Both run in CI, need no GPU, and pin Python 3.12 — the behaviour under test is interpreter-dependent, and on 3.13+ the single-data-point case passes even with the fix reverted.

> **Known problem — `scripts/mps-mig/` does not run**
>
> `NamedPopen.__init__` forwards `name=` to `subprocess.Popen`, which raises `TypeError` on the first process constructed. `GPU_queue.__len__` referred to a name that exists only as a local in `main()`. Both are fixed, but the harness has never been exercised end to end and is not installed by the conda recipe, so treat it as unverified rather than working.

## 9. Known problems

Everything below is _open_ as of v0.3.1. The defects fixed along the way — the seed-buffer overflow, the discarded scoring-file values, the duplicate segments within an interval, the out-of-range entropy counters, the unbounded seed pattern arrays, the `exit()` calls from TBB workers, and the orchestration bugs in §8 — are described in their sections as mechanism, not repeated here as complaints. This list is what a reader should still not be surprised by.

Ordered roughly by exposure: the first four can change what comes out, the next two are behaviour gaps worth knowing before you are surprised by them, and the last two bound how much the tests can tell you.

### 1. Duplicate segments across intervals

_§6 · src/seeder.cpp · changes output · not fixed_

`drop_duplicate_hsps()` runs at the end of one `seeder_body` call, which covers one interval. Intervals are separate invocations writing separate `.segments` files consumed by separate LASTZ calls, so an HSP found from both sides of an interval boundary is still emitted twice and LASTZ still repeats that gapped extension.

**Exposure:** any query longer than `lastz_interval_size` (10 Mb) — so every real genome, and none of the test data. One boundary per interval against 40 chunks inside it, which is why this is the smaller half of the problem the fix addressed.

**Measured once, and it did not appear.** Setting `lastz_interval_size=100000` forces 15 intervals and 14 boundaries onto the small test pair, without needing a 10 Mb query. On 4× A100 that produced 804 segments across 15 distinct files with **zero** duplicates, within files or across them. That is one input at one interval size, so it does not prove the case cannot arise — but it is the only evidence either way, and it argues the exposure is narrower than the mechanism suggests.

**What it would take:** deduplication after the intervals are joined, which means it can no longer live in the seeder. The natural home is the printer, but the printer receives one interval at a time too.

### 2. The emitted HSP set depends on the GPU's memory size

_§6 · src/seed_filter.cu:764, :822-826 · changes output · not fixed_

`MAX_HITS = MAX_HITS_PER_GB * global_mem_gb`, read from `cudaGetDeviceProperties()`. It sets `num_iter`, which partitions a chunk's hits into the batches that `hspEqual` is applied within. Since `hspEqual` is a _containment_ test — it deletes a real HSP, not a copy — which HSPs survive is a function of how the hits were partitioned, and therefore of how much memory the card has.

**Exposure:** every run on hardware unlike the one a result was produced on. This is a reproducibility property, not an occasional glitch, and it silently qualifies every measured claim in this document.

**What it would take:** first a measurement — run the parity test on two cards of different capacity and diff the segment files; the answer decides whether this is theoretical or routine. If it is real, the fix is to apply containment once over the joined result rather than per batch, which is the same restructuring problem as item 1.

### 3. The entropy counters are `short`

_§5 · src/seed_filter.cu:266-267 · changes output · not fixed_

`short count[4]` holds a per-base match count over an HSP's whole extent. Above 32767 it wraps negative; `log()` of a negative is NaN; `(int)(total_score * NaN)` then fails the threshold test and the HSP is dropped — silently, with no diagnostic. LASTZ uses `int count[256]`.

**Exposure:** needs a single base to match more than 32767 times inside one HSP, so it wants a very long low-complexity extent. Rare, but exactly the input the entropy adjustment exists to handle in the first place.

**What it would take:** widening the counters to `int`, which costs 16 more bytes of frame per thread — cheap, but it changes register pressure in the hottest kernel, so it wants measuring rather than assuming.

### 4. Literal `X` is scored as an ambiguity code

_§3 · common/parameters.h, common/scoring.c · changes output · not fixable as-is_

LASTZ gives the character `X` `badScore` in _both_ its score sets. KegAlign has no code for it: `X_NT` is the catch-all holding every IUPAC ambiguity code as well, and those LASTZ scores with `fillScore`. One code, two required answers.

**Exposure:** low. `X` is rare in genomic FASTA and mostly appears in protein-derived or hand-edited sequence.

**What it would take:** a ninth code in the compressed alphabet, and therefore a wider substitution matrix — not a scoring change. Worth doing only alongside some other alphabet work.

### 5. One position is seeded twice at every interval boundary

_§2, §6 · src/seeder.cpp · load-bearing · will not be fixed as stated_

The chunk loop's `q_inter_end + 1` makes adjacent intervals both seed their shared endpoint. It reads as an off-by-one and is not one: intervals are laid out to `seq_block_len - cfg.seed.size`, and a span of `seed.size` fits exactly there, so removing the `+1` costs every block its final seed.

**Exposure:** one position in 10,000,000, and it feeds problem 1 rather than causing anything on its own.

**What it would take:** half-open interval bounds everywhere with the last interval extended by `seed.size` — a change to the layout in `main.cpp`, not to the `+1`. A code review proposed deleting the `+1` alone; it was checked against the interval arithmetic and rejected. Recorded here so the same fix is not proposed a third time.

### 6. Masking cannot be turned off

_§4 · no CLI surface · behaviour gap · not fixed_

LASTZ has the `[unmask]` sequence action; KegAlign has no equivalent and no option. The only lever is to upper-case the FASTA beforehand, which diverges _from_ LASTZ rather than toward it, and silently — nothing in the output records that it was done.

**What it would take:** case-folding in the compression kernels, not in the substitution matrix. `L_NT` collapses all four lower-case bases into one code, so by the time scoring sees a base its identity is gone and no matrix entry can recover it.

### 7. The strongest test in the repository cannot run in CI

_§7 · tests/test_lastz_parity.bash · confidence limit_

Everything else under `tests/` is deliberately host-side so the workflow can run it on every push. The parity test needs a GPU and both binaries on `PATH`, so it exits 77 and skips. It is a manual pre-release step, which means it protects the project exactly as often as someone remembers to run it.

**What it would take:** a GPU runner, or a scheduled job on a machine that has one. Until then, treat "the host-side suite is green" as saying nothing at all about extension behaviour — the retracted masking feature passed every host-side test it had.

### 8. Parity is established at one point in parameter space

_§7 · tests/test_lastz_parity.bash · confidence limit_

The measured agreement — 1225/1225 and 469/469 — is real, and it is a single configuration: seed `12of19`, `xdrop` 910, `hspthresh` 3000, `step` 1, one ~375 kb pair, single block. Three things it therefore does not cover:

· **the entropy path**, because soft-masked columns terminate extension before they can be counted — IUPAC codes or `--ambiguous=n` are what reach it, and the test carries neither;  · **the alphabet** — the test data is pure upper-case `ACGT`, with no `N`, no ambiguity code and no native lower case, so `N_NT` and `X_NT` are never constructed and the scoring policy in §3 is measured only on its ACGT rows;  · **multi-block inputs**, which have simply never been run (§7 — an earlier draft claimed they would differ legitimately, which was wrong);  · **emitted ordering**, which §6 calls a convention worth preserving and which the test destroys by sorting before it compares;  · **any other seed pattern**, including the `14of22` that produced the crash in §2. (`--step=2` _has_ now been run — see §2 — and agreed exactly.)

**What it would take:** more cases in the same script — it already loops over them, and each is three lines. This is the cheapest item on the list and probably the most valuable.

> **How to check**
>
> **Adding to this list.** An entry earns its place by naming the file, saying whether it changes output or only bounds confidence, and stating what a fix would actually require. An entry that says only "this is ugly" belongs in a comment beside the code, not here — and a problem that gets fixed leaves this list and becomes mechanism in its own section, because this document explains the code as it stands rather than recording what it survived.

---

Line references are against `main` at v0.3.0. LASTZ references are against `lastz/lastz` at master; the vendored files are those under `common/`.
