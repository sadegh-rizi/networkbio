# COSMOS meta-PKN node grammar

The 2020 resource was inspected before filtering. Node identifiers are
kept in `raw_node` provenance; PKN outputs use canonical names where a
gene or numeric HMDB identifier can be resolved.

| Pattern | Unique nodes | Examples |
|---|---:|---|
| `complex: underscore-joined symbols` | 323 | AAAS_AHCTF1_GLE1_NDC1_NUP107_NUP133_NUP153_NUP155_NUP160_NUP188_NUP205_NUP210_NUP214_NUP35_NUP37_NUP42_NUP43_NUP50_NUP54_NUP58_NUP62_NUP85_NUP88_NUP93_NUP98_POM121_RAE1_RANBP2_SEC13_SEH1L_TPR, ABI1_BRK1_CYFIP1_NCKAP1_WASF1, ABI1_BRK1_CYFIP1_NCKAP1_WASF2, ACD_POT1_TERF1_TERF2_TERF2IP_TINF2, ACTL6A_ARID1A_ARID1B_ARID2_SMARCA4_SMARCB1_SMARCC1_SMARCC2_SMARCD1_SMARCE1 |
| `enzyme_or_reaction: XGene Entrez plus reaction` | 19530 | Gene1000000002__SLC7A6, Gene1000000003__SLC7A6, Gene100000002__SLC22A5, Gene100000002__SLC6A14, Gene10001__SLC7A6_TRANSPORTER1 |
| `gene: plain symbol` | 5788 | A2M, AAK1, AANAT, AATF, AATK |
| `metabolite: canonical HMDB or model identifier` | 10773 | Metab__10000638_c, Metab__100016_c, Metab__10003572_c, Metab__10004518_c, Metab__10007466_c |
| `other` | 1 | RET/PTC2 |

Unique nodes: 36415
Rows inspected: 83557
