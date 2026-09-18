# Recheck checklist for the prospective gold standard

Second-reviewer follow-up, 2026-09-18. Work through the batches in order.
Batches 1, 4 and 6 need no external lookup. Batches 2, 3 and 5 do.

Record each decision as a new row in a file `prospective_recheck.csv` with columns
`drugbank_id, mondo_id, batch, question, source_checked, answer, action, checked_by, checked_on`.
Actions are one of: `keep`, `keep_flagged`, `remap`, `merge`, `drop`, `reinstate`.

Where to look, once:
- **FDA**: Drugs@FDA, `accessdata.fda.gov/scripts/cder/daf/`, search the NDA/BLA number, then
  "Approval History, Letters, Reviews and Related Documents". The approval letter gives the
  indication text and the date.
- **EMA**: the product's EPAR page on ema.europa.eu, "Procedural steps taken after authorisation"
  (the annex table). This gives the variation number, the CHMP opinion date and the
  Commission decision date, which are different dates.
- **PMDA**: the "New Drug Approvals" list PDF for the review period, plus the review report
  (shinsa houkokusho) for the same item. The public-knowledge flag appears in the list notes.
- **TGA**: ARTG entry by number, plus the AusPAR if one exists.
- **Health Canada**: the NOC database entry by NOC number.

---

## Batch 1. Policy call on the Japanese public-knowledge pathway (no lookup)

These seven pairs all came through PMDA's public knowledge-based application route, which
relicenses use that is already standard of care on published evidence rather than new sponsor
trials. They are the best-documented associations in the set and four of them already leak in
BioKG. Decide the policy once, then apply it to all seven.

- [ ] Decide the policy. Three defensible options:
  - **Exclude the class.** Add a rule clause: drop rows whose regulator record states the
    application was public knowledge-based. Mechanical, no judgement per pair.
  - **Keep and report both ways.** Add a `public_knowledge` column, report the headline metric
    on all 115 and a sensitivity metric on 108.
  - **Keep silently.** Not recommended. A reviewer who reads the PMDA notes will find them.
- [ ] Apply to these seven pairs:
  - carboplatin / endometrial cancer, DB00958, MONDO:0011962, `PMDA:2024-06-24:81`
  - cyclophosphamide / graft versus host disease, DB00531, MONDO:0013730, `PMDA:2024-02-09:14`
  - methylprednisolone / Kawasaki disease, DB00959, MONDO:0012727, `PMDA:2023-09-25:26`
  - mycophenolate mofetil / nephrotic syndrome, DB00688, MONDO:0005377, `PMDA:2025-09-19:22`
  - mycophenolate mofetil / pulmonary systemic sclerosis, DB00688, MONDO:0002047, `PMDA:2024-06-24:59`
  - rituximab / lupus nephritis, DB00073, MONDO:0005556, `PMDA:2023-08-23:4`
  - rituximab / autoimmune hemolytic anemia, DB00073, MONDO:0020108, `PMDA:2026-02-19:11`
- [ ] Check whether the class is larger than seven. Grep `candidates_reviewed.csv` for
  `Public knowledge-based` across all 1,288 rows and confirm how many reached S3. The four
  the write-up names as BioKG leaks are a subset of the seven above.
- [ ] If you exclude the class, note that **methylprednisolone / Kawasaki disease is one of only
  two Hetionet-evaluable pairs**. Excluding it takes Hetionet to n=1 and the axis becomes
  unreportable. Decide what you will say about the Hetionet axis before you run the filter,
  not after.

---

## Batch 2. Prior-indication lookups (external sources needed)

For each, answer one yes/no question. Yes means the ingredient already held the disease or a
parent of it before 2023-01-01, so the pair fails the rule and should be dropped.

### 2a. Palmoplantar pustulosis, three pairs

- [ ] apremilast / PPP, DB05676, MONDO:0015597, `PMDA:2025-03-27:67`
- [ ] brodalumab / PPP, DB11776, MONDO:0015597, `PMDA:2023-08-23:61`
- [ ] risankizumab / PPP, DB14762, MONDO:0015597, `PMDA:2023-05-25:58`
  - **Question**: does the pre-2023 Japanese label for each product already carry a psoriasis
    indication whose wording covers pustular psoriasis, and does the Japanese label treat
    palmoplantar pustulosis as a form of psoriasis or as a separate disease?
  - **Where**: the PMDA review report for the 2023 or 2025 approval will state the existing
    approved indications in its opening section. Otherwise the current Japanese package insert
    for Otezla, Lumicef and Skyrizi, indications section.
  - **Decision**: if the prior label says "pustular psoriasis" or the review report frames PPP as a
    psoriasis variant, drop all three. If the Japanese label treats PPP as a distinct disease and
    the prior indication is plaque psoriasis only, keep all three but add a `related_to_prior`
    note, because MONDO puts PPP two `is_a` steps from psoriasis and BioKG filters two of them.
  - **Do not decide these three separately.** Whatever you conclude applies to all three.

### 2b. Trametinib and thyroid cancer

- [ ] trametinib / differentiated thyroid carcinoma, DB08911, MONDO:0015447, `EMA:23-26-march-2026:Mekinist`
  - **Question**: what exactly did trametinib hold for thyroid cancer before 2023, and in which
    jurisdictions?
  - **Where**: Drugs@FDA NDA 204114 (Mekinist), the May 2018 approval letter for anaplastic
    thyroid cancer. Then the Mekinist EPAR procedural steps, to see whether the EU ever held an
    anaplastic thyroid indication before the March 2026 DTC variation.
  - **Decision**: if the EU held no thyroid indication before 2023, the pair is defensible on the
    EU record but still fails the global "any of the five" rule because of the 2018 FDA anaplastic
    approval. DTC and anaplastic thyroid carcinoma are sibling subtypes under thyroid carcinoma,
    so under the rule's "not a subtype of a disease the ingredient already treated" clause this is
    a drop. Drop unless you can argue anaplastic and differentiated thyroid carcinoma are not
    subtypes of a common parent, which the ontology does not support.
  - **Same check applies to** encorafenib / thyroid cancer (DB11718, `PMDA:2024-05-17:71`) and
    binimetinib / thyroid cancer (DB11967, `PMDA:2024-05-17:72`), but for the opposite reason:
    confirm neither ingredient held any thyroid indication before 2023 (I believe neither did),
    then fix the term, see Batch 3.

### 2c. Juvenile arthritis, three pairs

- [ ] sarilumab / polyarticular JIA, DB11767, MONDO:0018456,
      `FDA:BLA761037-S15; EMA:14-17-october-2024:Kevzara; EMA:11-14-november-2024:Kevzara`
- [ ] baricitinib / polyarticular JIA, DB11817, MONDO:0018456, `PMDA:2024-03-26:66`
- [ ] ixekizumab / enthesitis-related JIA, DB11569, MONDO:0019437, `EMA:21-24-july-2025:Taltz`
  - **Question**: this is not really a lookup, it is a definitional decision. The ingredients held
    the adult counterparts (RA for sarilumab and baricitinib, axSpA and PsA for ixekizumab) well
    before 2023. Under the rule's "population extension" clause these should have been excluded.
  - **What to check anyway**: confirm for each that the pre-2023 adult indication existed in at
    least one of the five. Sarilumab RA, Kevzara FDA BLA 761037, 2017. Baricitinib RA, Olumiant
    FDA NDA 207924, 2018, and EMA 2017. Ixekizumab PsA and axSpA, Taltz FDA BLA 125521, 2017 and
    2020.
  - **Decision**: pick one of two positions and write it into section 1.
    - **Keep**, on the argument that pJIA and ERA are separate disease nodes with separate trials
      and separate labels, so the model is predicting an edge to a node it has never seen
      connected to that drug. This is the stronger argument for a link-prediction paper.
    - **Drop**, on the argument that the rule says population extensions are out.
  - **Either way**: delete the section 1 sentence "no pair is a setting or population change of an
    existing indication", which these three falsify.
  - **Also check** the two EMA Kevzara cards in the evidence list are the same procedure. The
    November 2024 card looks like a duplicate of the October 2024 one, and the row rationale says
    so. If they are the same procedure, the earlier date is the right `approval_date`.

### 2d. Tumour-agnostic overlap, three pairs

- [ ] trametinib / low grade glioma, DB08911, MONDO:0021637
- [ ] dabrafenib / low grade glioma, DB08912, MONDO:0021637
  - **Question**: does the FDA tumour-agnostic BRAF V600E solid-tumour indication granted to
    dabrafenib plus trametinib in June 2022 cover paediatric low-grade glioma?
  - **Where**: Drugs@FDA NDA 202806 (Tafinlar) and NDA 204114 (Mekinist), the 2022-06-22 approval
    letter, then the current label's Indications section, section 1, for the age range and the
    "no satisfactory alternative treatment options" wording.
  - **Decision**: if the 2022 indication covers patients aged 6 and over with any solid tumour
    including LGG, then the 2023 LGG approval is an organ-specific carve-out of an indication the
    ingredient already held, and the pair fails. This matters beyond these two rows, because it is
    the same logic by which the set drops six tumour-agnostic rows at S2 and then admits their
    organ-specific descendants.
- [ ] pembrolizumab / ovarian cancer, DB09037, MONDO:0008170, `FDA:BLA125514-S186; FDA:BLA761467-S2; HC:37409`
  - **Question**: do the 2017 MSI-H/dMMR and 2020 TMB-H tumour-agnostic indications already cover
    biomarker-selected ovarian carcinoma?
  - **Where**: Drugs@FDA BLA 125514, the 2017-05-23 and 2020-06-16 letters.
  - **Decision**: if yes, decide whether biomarker-restricted tumour-agnostic coverage counts as
    "an approved indication for that disease". I would argue it does, and that the same answer
    should then be applied to the LGG pairs above. Whatever you decide, apply it consistently and
    state it in section 1, because at the moment the set treats tumour-agnostic rows as
    out-of-scope at S2 without saying what that implies for their descendants.

### 2e. Four single-pair lookups

- [ ] methylprednisolone / Kawasaki disease, DB00959, MONDO:0012727 (also in Batch 1)
  - **Question**: do the pre-2023 systemic methylprednisolone labels in any of the five already
    carry an indication broad enough to cover Kawasaki disease, for example severe acute
    inflammatory disease or vasculitis?
  - **Where**: the US Solu-Medrol label, NDA 011856, Indications section, and the Japanese
    Solu-Medrol package insert as it stood before September 2023.
  - **Decision**: if the prior label names vasculitides or "severe acute inflammatory processes"
    generically, this is a parent-of-prior-indication case and the pair should go.
- [ ] nivolumab / skin carcinoma, DB09035, MONDO:0002656, `PMDA:2024-02-09:92`
  - **Question**: what does "unresectable advanced or recurrent epithelial skin malignancies" cover
    in the Japanese label, and does nivolumab hold any pre-2023 non-melanoma skin cancer
    indication in any of the five?
  - **Where**: PMDA review report for the February 2024 approval.
  - **Decision**: I do not believe nivolumab held cSCC or BCC anywhere before 2023, so the pair
    probably survives. The real problem is the term, see Batch 3.
- [ ] bexarotene / adult T-cell leukemia/lymphoma, DB00307, MONDO:0019471, `PMDA:2024-06-24:76`
  - **Question**: the indication is "ATL with skin lesions". Does the pre-2023 Japanese Targretin
    label's cutaneous T-cell lymphoma indication already cover the cutaneous manifestation of ATL?
  - **Where**: the PMDA review report, which will state the existing CTCL indication and say
    whether ATL skin lesions were previously in or out of scope.
  - **Decision**: if the review report frames this as extending an existing cutaneous T-cell
    lymphoma indication to a new histology, keep. If it frames it as covering skin lesions already
    within the CTCL indication, drop.
- [ ] melatonin / jet lag, DB01065, MONDO:0024381, `TGA:immela-link-medical-products-pty-ltd-ta-link-pharmaceuticals`
  - **Question**: did melatonin hold any approved sleep-disorder indication in Australia or the EU
    before 2023, and does the set handle the two halves of the IMMELA indication consistently?
  - **Where**: the ARTG entry for Circadin in Australia, and the Circadin and Slenyto EPARs.
  - **Decision**: the IMMELA text contains two indications, jet lag in adults and sleep disorders
    in children with ASD or ADHD. The second is the Slenyto indication, EU 2018, so it was rightly
    not carried. Confirm the jet lag half was genuinely unapproved in all five before 2023, then
    keep. If Circadin or any Australian melatonin product carried jet lag, drop.

---

## Batch 3. Term and scope fixes (read the label text only, no approval history)

These do not need a history lookup. Read the regulator text already in `candidates_reviewed.csv`
and pick the MONDO term that matches it. Record the old and new term.

- [ ] tezepelumab, DB15090: change **MONDO:0006031 chronic rhinosinusitis** to the CRSwNP term.
  All three regulator texts say "chronic rhinosinusitis with nasal polyps". Check whether MONDO
  has a nasal-polyps-specific child; if it does not, say so in the methods rather than silently
  using the parent.
- [ ] somapacitan, DB15093: the FDA letter says "growth failure associated with Noonan syndrome".
  **MONDO:0018997 Noonan syndrome is wrong.** Either remap to a short-stature or growth-failure
  term, or drop the pair as not a treatment indication for the mapped disease. Note the same
  letter's SGA and ISS components were dropped, so dropping this one is the consistent choice.
- [ ] faricimab, DB15303: the PMDA text says "choroidal neovascularization associated with angioid
  streaks". **MONDO:0011782 angioid streaks is wrong**; the treated lesion is CNV, which faricimab
  already treats in nAMD. Remap to a CNV term and then re-run the prior-indication test, which it
  will probably fail, or drop.
- [ ] miglustat, DB00419: the TGA text says miglustat is "an enzyme stabiliser of cipaglucosidase
  alfa". Decide whether an enabler of a co-administered enzyme counts as a treatment indication for
  Pompe disease. My recommendation is to drop and to widen the rule's toxicity carve-out to cover
  adjuncts and enablers generally.
- [ ] pegcetacoplan / AMD, DB16694, MONDO:0005150: the TGA text is geographic atrophy secondary to
  AMD. Remap to a geographic-atrophy term if MONDO has one, otherwise keep the parent and flag it,
  because MONDO:0005150 also covers neovascular AMD.
- [ ] encorafenib and binimetinib / thyroid cancer, MONDO:0002108: the PMDA text grants two things,
  BRAF-mutant thyroid cancer after chemotherapy and anaplastic thyroid cancer. Decide whether to
  keep the parent term or split. Keeping the parent is defensible here; just say so.
- [ ] brentuximab vedotin, DB08870, MONDO:0018905: the letter says large B-cell lymphoma including
  DLBCL NOS and high-grade B-cell lymphoma. DLBCL is narrower than what was approved. Remap to an
  LBCL term or note the narrowing.
- [ ] pembrolizumab / ovarian cancer, DB09037: one term covers epithelial ovarian, fallopian tube
  and primary peritoneal carcinoma. Note it or split it.
- [ ] capivasertib, DB12218, MONDO:0008315: generic prostate cancer parent kept for PTEN-deficient
  mHSPC. Note it.
- [ ] oxybutynin, DB01062, MONDO:0007754: both regulators say primary **palmar** hyperhidrosis;
  the term says palmaris et plantaris. Remap or note.
- [ ] human immunoglobulin g, DB00028, MONDO:1010185: the indication is antibody-mediated rejection
  in six named organ transplants. The term is a generic parent and DB00028 is a class-level node.
  Decide whether a class-level immunoglobulin node belongs in a link-prediction gold standard at
  all; it is unusually high-degree in most graphs.
- [ ] cabozantinib, DB08875: the same decision grants extra-pancreatic NETs, which was not carried.
  Add the second pair or record why only pNET was kept.

---

## Batch 4. One decision, two pairs (no lookup, just a merge decision)

- [ ] pegcetacoplan, DB16694: **complement 3 glomerulopathy (MONDO:0018013)** and
  **IC-MPGN (MONDO:0014005)** come from one sNDA, `FDA:NDA215014-S11`, one EMA variation and one
  trial, APL2-C3G-310. The two terms are two `is_a` steps apart, inside your own
  `prospective_max_hops = 2` relatedness threshold. Merge to one pair, or justify why one
  regulatory decision yields two independent positives.
- [ ] sirolimus, DB00877: **lymphangioma (MONDO:0002013)** and **vascular malformation
  (MONDO:0024291)** come from one PMDA decision, `PMDA:2024-01-18:11`, which lists ten entities
  (lymphangioma, lymphangiomatosis, Gorham disease, lymphangiectasia, hemangioendothelioma, tufted
  angioma, venous malformation, blue rubber bleb nevus syndrome, mixed vascular malformation,
  Klippel-Trenaunay-Weber syndrome). Two were kept, and MONDO:0024291 is effectively the parent
  class of the other. Either keep one term, or keep all ten and say so. Keeping an arbitrary two is
  the one option a reviewer will not accept.
- [ ] Fixed combination partners: dabrafenib and trametinib contribute LGG and hairy cell leukemia;
  encorafenib and binimetinib contribute NSCLC and thyroid cancer. That is **eight pairs from four
  decisions**, for drugs that are never given alone. Clustering the bootstrap by drug does not fix
  this because they are different drugs. Add a `decision_id` column and cluster the bootstrap by
  it, or state the effective independent-event count alongside n=115.

---

## Batch 5. Exclusions to reverse (external sources needed)

- [ ] **zanubrutinib / CLL-SLL**, currently dropped as `ingredient_prior_indication`, rows
  `FDA:NDA213217-S7` (2023-01-19) and `TGA:brukinsa-beigene-aus-pty-ltd-2` (2023-03-01).
  - **Question**: on what date did the European Commission grant Brukinsa the CLL indication? The
    rationale cites "EMA 2022-11 (EC decision)", but November 2022 looks like the CHMP positive
    opinion date, with the Commission decision in January 2023.
  - **Where**: the Brukinsa EPAR, "Procedural steps taken after authorisation", find the CLL
    variation and read both the opinion date and the Commission decision date.
  - **Decision**: if the Commission decision is on or after 2023-01-01, zanubrutinib held no CLL
    indication in the five at the cutoff, and **this pair should be reinstated** with
    `approval_date` 2023-01-19 (FDA) or the CHMP opinion date, whichever your `date_type`
    convention requires. Note that your convention elsewhere records CHMP opinion dates, so pick
    one convention and apply it to the exclusion test as well as to `approval_date`.
- [ ] **tremelimumab / non-small cell lung carcinoma**, currently dropped as `line_extension`, row
  `EMA:19-22-june-2023:Imjudo`.
  - **Question**: did tremelimumab hold any indication other than unresectable HCC with durvalumab
    before June 2023?
  - **Where**: Drugs@FDA BLA 761270 (Imjudo), approval history, and the Imjudo EPAR procedural
    steps.
  - **Decision**: if HCC was the only prior indication, NSCLC is a different disease and not a
    line, stage, regimen or population of it, so **reinstate the pair**. The row's own rationale
    already concedes this: "new organ but combination-defined; borderline new_indication". Note
    that reinstating adds a pair whose approval is combination-only, which interacts with the
    `combination_type` problem in Batch 6.
- [ ] **valganciclovir / symptomatic congenital CMV**, currently dropped as `line_extension`, row
  `PMDA:2023-03-27:64`.
  - **Question**: is congenital CMV a distinct disease entity from the CMV retinitis and CMV
    prophylaxis indications valganciclovir already held in Japan?
  - **Where**: the PMDA review report for the March 2023 approval, which will state the existing
    indications, plus the MONDO entry for congenital cytomegalovirus infection.
  - **Decision**: if MONDO treats congenital CMV as a separate entity and the trial endpoints are
    neurodevelopmental and audiological rather than antiviral suppression in the same population,
    reinstate or move to `unclear`. Calling it a population of CMV retinitis is a stretch.
- [ ] **emapalumab / HLH-MAS in Still disease**, currently dropped as `line_extension`, row
  `FDA:BLA761107-S18`. Lower priority.
  - **Question**: does MONDO separate primary HLH (the 2018 indication) from secondary
    HLH/macrophage activation syndrome?
  - **Decision**: if it does, this is borderline and belongs in `unclear` rather than
    `line_extension`.
- [ ] **melatonin Slenyto / neurogenetic disorders**, row `EMA:22-25-july-2024:Slenyto`. Evidence
  only.
  - The shipped EPAR string is truncated mid-word at "treatment of neurogenetic disorders (e."
    Retrieve the full text from the Slenyto EPAR before leaving the row classified either way.
- [ ] **sodium thiosulfate / cisplatin ototoxicity**, row `TGA:pedmarqsi`. Routing only.
  - Outcome is right but it was excluded as a prior indication when the rule's "prevention of
    another drug's toxicity" clause is what should have caught it. Reclassify so that the rule is
    demonstrably doing the work.
- [ ] **cisplatin / urothelial carcinoma** (`PMDA:2024-12-27:107`) versus **carboplatin /
  endometrial cancer** (kept). Consistency only, no lookup.
  - Both are decades-old Japanese public-knowledge catch-ups for a platinum agent. One is out
    because the US cisplatin label has carried bladder cancer since 1978; one is in because the US
    carboplatin label happens to be ovarian-only. Resolve this when you decide Batch 1.

---

## Batch 6. Data and prose fixes (no lookup)

- [ ] **Section 2, S3 row**: change `dropped` from **120 to 128**. `prospective_funnel.csv` says
  128, and 243 minus 115 is 128.
- [ ] **Section 1**: delete or qualify "no pair is a setting or population change of an existing
  indication". Ten pairs carry `setting_change_noted = True` and three juvenile-arthritis pairs are
  population counterparts.
- [ ] **`combination_type` column**: currently "single agent" for all 115. Thirty pairs have
  combination or add-on language in at least one regulator text, and at least ten are approved only
  as part of a fixed regimen (tucatinib with trastuzumab, adagrasib with cetuximab, sotorasib with
  panitumumab, miglustat with cipaglucosidase alfa, zanubrutinib with obinutuzumab, tafasitamab
  with lenalidomide and rituximab, brentuximab with lenalidomide and rituximab, capivasertib with
  abiraterone, encorafenib with binimetinib, dabrafenib with trametinib). Recompute the column or
  remove it. The earlier set's `Combo_Flag` exclusion has no successor here; say whether that is
  deliberate.
- [ ] **Section 2, BioKG leak sentence**: soften "a leak there is a known association, not a prior
  approval, which is exactly what check 1 is meant to remove". Four of the six leaks are
  decades-old standard-of-care uses and the other two have psoriasis as the nearest prior term.
  The leak pattern is evidence against those pairs, not a defence of them.
- [ ] **Section 1, "2,258 recorded decisions"**: either vendor `review/review_decisions.csv`
  alongside `candidates_reviewed.csv` or drop the number. It is not checkable from what ships in
  `data/gold_standards/`.
- [ ] **Section 2, "Hetionet holds 1,552 compounds ... and 137 diseases"**: not re-derivable from
  the shipped tables, unlike every other number in that section. Either ship the views or cite the
  release file.
- [ ] **Section 3, `min_new_pairs = 10`**: common support is currently **1 pair** (roflumilast to
  atopic dermatitis). State plainly that common-support scoring is not evaluable at present rather
  than leaving the reader to notice.
- [ ] **Section 5 citation claim**: "every `ingredient_prior_indication` and `line_extension`
  rationale names the prior approval (regulator and month)" holds for 84 of 96
  `ingredient_prior_indication` rows but only **9 of 304** `line_extension` rows, 244 of which name
  no year at all. Narrow the claim to the `ingredient_prior_indication` bucket.
- [ ] **Stray reviewer instructions in shipped rationales**: search `class_rationale` for
  `confirm`, `check`, `borderline`, `TBD`. At least `tislelizumab / nasopharyngeal carcinoma`
  ends with "confirm" and a pembrolizumab gastric row ends with "— check".
- [ ] **upadacitinib / alopecia areata**, DB15091, MONDO:0005340: the row's `class_rationale` reads
  "EPAR steps: non-segmental vitiligo", which is the other pair backed by the same EMA card. Fix
  the rationale text.
- [ ] **lisocabtagene maraleucel, DB16582**: nine vaccines were dropped at S2 as product-level
  identifiers that no graph models as compounds. A CAR-T cell product is the same kind of entity.
  Either drop it or state why cell therapies are treated differently from vaccines.
- [ ] **`line_extension` versus `expanded_population`**: the two buckets are not separable.
  Letermovir (HSCT to solid organ), concizumab (inhibitor to non-inhibitor), valganciclovir (adult
  to congenital) and melatonin (ASD to neurogenetic) are filed as line extensions while being
  population extensions on their face. Either merge the two counts in section 2 or define the
  boundary.
- [ ] **Per-graph cutoff**: section 1 fixes a global 2023-01-01 cutoff that is five and a half
  years after Hetionet's content date and two and a half after BioKG's. Everything approved between
  a graph's content date and 2023 is discarded although it is a valid prospective target for that
  graph. Either move to a per-graph cutoff, which would raise the Hetionet-evaluable count above 2,
  or state in section 1 that 2023 is a source-availability boundary and not a validity criterion.
- [ ] **Rule versus check 3**: section 1's rule and section 3's two-hop check are not the same
  criterion. Under the rule, ulcerative colitis is new for a Crohn drug; under check 3, psoriasis
  and palmoplantar pustulosis are two hops so the PPP pairs would be filtered. State the two-hop
  rule as the definition or explain the difference.
- [ ] **Tislelizumab and NMPA**: four pairs are in on EU dates for indications China approved
  earlier. The ingredient list where NMPA matters is short. Either check it or report the metric
  with those four removed.
- [ ] **FDA source bias**: section 1 treats the missing original NDAs and BLAs (Syfovre, Iwilfin,
  Opfolda, Zoryve foam, Zepbound) as a date-accuracy problem. Reframe it as a systematic bias
  against the target class, since a new product containing an old ingredient is the canonical
  repurposing event and is structurally invisible in a CDER efficacy-supplement feed.

---

## Order of work if time is short

1. Batch 1 policy decision, because it moves seven pairs and decides what you can say about Hetionet.
2. Batch 5, the two reinstatements, because a reviewer who samples the exclusions will find them.
3. Batch 6 items 1 to 4, because they are errors in the text rather than judgement calls.
4. Batch 2c, the juvenile arthritis position, because it needs a written justification either way.
5. Everything else.
