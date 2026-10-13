base = StringReplace[DirectoryName[DirectoryName[$InputFileName]], "\\" -> "/"] <> "/";
$MaxExtraPrecision = 500;
f = 1 - 2/r;
hp = 1 + 4/r - 2/(r-2);
sigmad = -2/r^2;
weight = 4 sigma^2 (1-sigma^2);
fs = Get[base <> "investigation/zsource_rwparents_2.wl"];
raw = fs . {p^2, p (I w hp p + sigmad dp),
                    (I w hp p + sigmad dp)^2};
src = Together[-(raw /. r -> 2/sigma)/(2 weight)];
coeff = Together /@ {Coefficient[src,p,2],
                      Coefficient[Coefficient[src,p],dp],
                      Coefficient[src,dp,2]};
Print["Exporting rational coefficients; endpoint cancellation is tested on the full source."];
inp = Import[base <> "investigation/overlap_inputs.json", "RawJSON"];
cx[z_] := ToExpression[z[[1]]] + I ToExpression[z[[2]]];
freq = cx[inp["modes"]["220"]["omega"]];
enc[z_] := {Re[N[z,17]], Im[N[z,17]]};
rat[e_] := Module[{v=Together[e /. w -> freq]},
 <|"numerator" -> (enc /@ CoefficientList[Numerator[v],sigma]),
   "denominator" -> (enc /@ CoefficientList[Denominator[v],sigma])|>];
Export[base <> "theory_direction1/schwarzschild_source.json",
 <|"description" -> "Even 220 self-coupling to regularized Zerilli 44, hyperboloidal acceleration source, identical-pair factor included",
   "coefficients" -> (rat /@ coeff)|>, "RawJSON"];
Put[coeff,base <> "theory_direction1/schwarzschild_source.wl"];
Exit[];
