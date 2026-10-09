# Isotopologue convention

The deposited ST003331 and ST003332 mwTab files were checked before stage-01
implementation. Their ST, SP, TR, AN and MS sections document 24-hour
[U-13C6]glucose exposure and the extraction volumes, but contain no explicit
statement that plain metabolite rows are total pools. The data use plain names
such as `D-Glucose` alongside labelled names such as `Glucose 13C6`.

R0 finding: no statement found that plain rows are total pools. D1 is therefore
registered as an interpretation, not a deposit fact: plain rows are treated as
M+0, labelled sums are called `isotopologue_sum_deposited`, and incomplete
labelled coverage is reported. For glucose-reachable metabolites, exchange
values describe changes in the unlabelled fraction rather than abundance.

## Whitelist coverage
The whitelist is fixed by the revised plan; no IDs were added during implementation.

### ST003331
- Found: C00135, C00123, C00047, C00073, C00079, C00188, C00078, C00183, C00082, C02989, C05335, C00318, C00072, C05422, C00864, C03492, C01595, C02700, C00108, C10164, C03722, C00322, C00463, C05658, C00331, C00637, C05635, C02796, C00788
- Absent: none

### ST003332
- Found: C00135, C00123, C00047, C00073, C00188, C00078, C00183, C00082, C02989, C00318, C00072, C05422, C00864, C03492, C01595, C02700, C00463, C05658, C00331, C00637, C05635
- Absent: C00079, C05335, C00108, C10164, C03722, C00322, C02796, C00788
