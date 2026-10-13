base = StringReplace[DirectoryName[DirectoryName[$InputFileName]], "\\" -> "/"] <> "/";
$MaxExtraPrecision=500;
inp=Import[base<>"investigation/overlap_inputs.json","RawJSON"];
cx[z_] := SetPrecision[ToExpression[z[[1]]]+I ToExpression[z[[2]]],60];
enc[z_] := {Re[N[z,17]],Im[N[z,17]]};
w=cx[inp["modes"]["220"]["omega"]];
aa=cx/@inp["modes"]["220"]["a"];
daa=Range[Length[aa]-1] Rest[aa];
poly[a_,x_]:=Fold[#1 x+#2&,Last[a],Reverse[Most[a]]];
coeff=Get[base<>"theory_direction1/schwarzschild_source.wl"];
rows={};
Do[
 grid=N[(1-Cos[Pi (Range[n]+(-1/2))/n])/2,60];
 p=(poly[aa,1-#]&/@grid)/Total[aa];
 dp=-(poly[daa,1-#]&/@grid)/Total[aa];
 b=MapThread[Total[(coeff/.sigma->#1){#2^2,#2 #3,#3^2}]&,{grid,p,dp}];
 rr=2/grid; ff=1-grid;
 hh=4+12*(rr-2)/(rr^2*(4*rr+6));
 zz=((hh+2 I w*(1-2 grid^2))*p-grid^2*(1-grid)*dp)/(4+2 I w);
 AppendTo[rows,<|"n"->n,"sigma"->N[grid,17],"source"->(enc/@b),
   "rw_parent"->(enc/@p),"z_parent"->(enc/@zz)|>];
 Print["N=",n," source first,last: ",N[b[[{1,-1}]],12]];
,{n,{32,48,64,80,96,128}}];
Export[base<>"theory_direction1/schwarzschild_samples.json",
 <|"omega220"->enc[w],"rows"->rows|>,"RawJSON"];
Exit[];
