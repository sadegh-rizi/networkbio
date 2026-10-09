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
