# Course proposal (LaTeX) and bibliography

Date: 2026-10-02
Model: claude-opus-5-5 (Claude Science)
Settings: default
Plan: none (direct instruction; text only, no analysis code)
Commit(s): none yet (add after the user commits)
Transcript: not saved

## Prompt actually used

> can you complete this research proposal
> also give me a .bib file which I can use as reference for the mentioned papers

(The user pasted the course LaTeX template with a draft Background and Aim.)

## What it produced

- `doc/proposal/proposal.tex`: Background, Aim, Approach, Context and a
  two-week Timeline, written into the course template (1-page A4 limit for
  the first four sections).
- `doc/proposal/references.bib`: 21 entries; 20 built from Crossref records
  of their DOIs on 2026-10-02, plus `hagberg2008networkx` (no DOI).

Facts in the text were checked against the Ionescu PDF in `data/Ionescu/`
(the six TFs, the HMGCR/SASP claim) and `doc/datasets.md` (accessions,
line counts).

## Verification

- Every `\cite` key exists in the .bib; `love2014deseq2` is in the .bib but
  not cited.
- **Not compiled**: no working LaTeX in the environment that wrote it. A
  font-metric estimate put page 1 at about 99% full; compile and check that
  the Timeline starts on page 2.

## Revision (same day)

User prompt: "in the aim, as I had written before I mostly want to identify the
gaps of the current multi-omics integration at GRN level methods and see if they
are doing well or if there are somethings they miss. and for that I am
primarily using CORNETO or i can try other methods if there's time"

Changed: title, Background (adds the methods' known weak points), Aim
(method evaluation; the PMS iNSC data become the test case), Methodology
(iv) gap tests and (v) comparison with GENIE3, Software (GENIE3), Context,
and Timeline days 10-12. Added `garciaalonso2019dorothea` to the .bib.
Page-1 estimate unchanged at about 99% full; still not compiled.

## Revision 2: supervisor feedback (same day)

Input: the user's edited version with the supervisor's inline comments
(stronger argument for the evaluation; say why the dataset suits the project
and what its quality is; aims as a bullet list, ordered with the
route-recovery question last; per-modality robustness, late integration,
modality-specific differences linked to MOFA/SOFA, early/intermediate
integration with CORNETO/COSMOS; data as sub-bullets).

Changed: Background (why method comparison matters; dataset scope, size and
quality), Aim (methodological aim + A1-A5 list), Data as sub-bullets (adds
lipidomics ST003328 and DoRothEA), Methodology mapped to A1-A5, Timeline
mapped to A1-A5. Kept the user's own edits (title wording, name, ID, Context).
Added `capraz2026sofa` (Nat Commun 17:8725, doi 10.1038/s41467-026-74694-6,
from Crossref) to the .bib. Page-1 estimate: about 99% full; not compiled.

## Revision 3: second round of supervisor comments (same day)

Input: four comments from the supervisor (I. Mohorianu), forwarded by the user:

1. Separate statistical approaches from machine-learning ones, and be more precise about
   directionality. GENIE3 is a suitable reference here; MOFA is not.
2. Add an intermediate step to the argument: current approaches infer GRNs one modality at a time,
   which assumes the modalities are statistically independent (biologically unlikely); next, MOFA and
   SOFA reduce dimensions across modalities but do not reach a GRN; only then CARNIVAL and COSMOS.
3. Describe CORNETO as "consolidating signal across samples and modalities".
4. Organise the methodology by objective, with sub-lists for the to-dos.

Changed:
- **Background**, in the order the supervisor gave:
  - statistical approaches (correlation, undirected) are separated from ML approaches (GENIE3,
    directed only because the regulators are fixed in advance, and unsigned);
  - the independence assumption is stated;
  - MOFA/SOFA are described as dimension reduction without a GRN;
  - CARNIVAL/COSMOS on a signed, directed PKN follow;
  - CORNETO "consolidates signal across samples and modalities";
  - the MOFA citation was removed from the data-driven sentence.
- **Methodology** is now a sub-list: Inputs, then A1–A5, one line each.

The Data sub-list, Aim, Software, Context and Timeline are unchanged. No new .bib entries were
needed; every `\cite` key exists in the .bib.

Not compiled. A character-count estimate puts the edits at roughly 5–6 lines longer than
revision 2. Revision 2 was already about 99% full, so page 1 probably overflows and needs a
matching cut.

## Revision 4: fit page 1 after revision 3 (same day)

Revision 3 already addressed the four comments (statistical vs ML and
directionality; GENIE3 kept, MOFA citation moved; independence -> MOFA/SOFA ->
CARNIVAL/COSMOS order; "consolidates signal across samples and modalities";
methodology as sub-list by aim). It was estimated ~6 lines over page 1.
Cuts, with the supervisor's structure kept: Background tightened (same
sequence of points; `garridorodriguez2022` no longer cited), Aim lead-in,
Software (citations grouped) and Context shortened to one line each.
Font-metric estimate is back to revision-2 level (~99% of page 1). Not compiled.

## What was changed by hand afterwards

The user's final version (pasted 2026-10-02 19:09) replaced
`doc/proposal/proposal.tex` verbatim. Hand changes relative to revision 4 include:
margins 0.75in and a compact title/section format; Background rewritten in the
user's words; data-quality sentence and "All analyses are illustrated on the case
study" sentence removed; Aim labels "Aim 1-5"; OmniPath citation dropped; Software
line edited; Timeline reverted to the revision-1 wording; name and ID filled in.

## Outputs stamped

None.
