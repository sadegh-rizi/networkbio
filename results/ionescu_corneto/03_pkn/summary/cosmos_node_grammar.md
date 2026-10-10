# COSMOS meta-PKN node grammar

The selected NetworkCommons resource was inspected before filtering.
Reaction indices in `Gene<index>__<suffix>` are not Entrez IDs. PKN
outputs preserve reaction and numeric model-metabolite node names.

| Pattern | Unique nodes | Examples |
|---|---:|---|
| `complex: underscore-joined symbols` | 323 | AAAS_AHCTF1_GLE1_NDC1_NUP107_NUP133_NUP153_NUP155_NUP160_NUP188_NUP205_NUP210_NUP214_NUP35_NUP37_NUP42_NUP43_NUP50_NUP54_NUP58_NUP62_NUP85_NUP88_NUP93_NUP98_POM121_RAE1_RANBP2_SEC13_SEH1L_TPR, ABI1_BRK1_CYFIP1_NCKAP1_WASF1, ABI1_BRK1_CYFIP1_NCKAP1_WASF2, ACD_POT1_TERF1_TERF2_TERF2IP_TINF2, ACTL6A_ARID1A_ARID1B_ARID2_SMARCA4_SMARCB1_SMARCC1_SMARCC2_SMARCD1_SMARCE1 |
| `enzyme_or_reaction: reaction index plus suffix` | 19530 | Gene1000000002__SLC7A6, Gene1000000003__SLC7A6, Gene100000002__SLC22A5, Gene100000002__SLC6A14, Gene10001__SLC7A6_TRANSPORTER1 |
| `gene: plain symbol` | 5788 | A2M, AAK1, AANAT, AATF, AATK |
| `metabolite: HMDB` | 3491 | Metab__HMDB0000002_c, Metab__HMDB0000002_e, Metab__HMDB0000005_c, Metab__HMDB0000005_e, Metab__HMDB0000005_m |
| `metabolite: canonical or model identifier` | 797 | Metab__11_cis_retfa_c, Metab__11_cis_retfa_e, Metab__13_cis_oretn_c, Metab__13_cis_oretn_n, Metab__13_cis_retnglc_c |
| `metabolite: model identifier` | 2515 | Metab__10fthf5glu_c, Metab__10fthf5glu_e, Metab__10fthf5glu_l, Metab__10fthf5glu_m, Metab__10fthf6glu_c |
| `metabolite: numeric model identifier` | 3970 | Metab__10000638_c, Metab__100016_c, Metab__10003572_c, Metab__10004518_c, Metab__10007466_c |
| `other` | 1 | RET/PTC2 |

Unique nodes: 36415
Rows inspected: 83557
