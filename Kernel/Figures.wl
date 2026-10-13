plotStyle={Frame->True,Axes->False,FrameStyle->GrayLevel[.25],BaseStyle->{FontFamily->"Times",FontSize->12},
 ImageSize->390,AspectRatio->.58,GridLines->Automatic,GridLinesStyle->Directive[GrayLevel[.9],Thin],PlotRangePadding->Scaled[.025],ImagePadding->{{62,16},{46,32}}};
QCSFigures[new_:False]:=Module[{pair,rec,t,ref,pred,peak,deltaPeak,top,bottom,plots,fig2,fig4,kerr,dir,names,colors,ix,traces,norm,
 spins,w2,w4,w6,b4,b6,dp,dc,da,gap4,gap6,kappa,curves,titles,labels,styles,grid,targets,near,maxSpin},
 pair=If[new,Import[path["results/paired_preparations.wxf"],"WXF"],QCSData["source/reduced_input/paired_preparations_fine.mat"]];
 rec=If[new,Import[path["results/recovery_waveforms.wxf"],"WXF"],QCSData["source/reduced_input/completion/recovery_waveforms.mat"]];
 t=rec["t"];ref=rec["reference_plus"]-rec["reference_minus"];pred=rec["prediction_plus"]-rec["prediction_minus"];
 peak=MapThread[Max,{Max/@Transpose[Abs[rec["reference_plus"]]],Max/@Transpose[Abs[rec["reference_minus"]]]}];deltaPeak=Max/@Transpose[Abs[ref]];
 colors={RGBColor[.09,.50,.48],RGBColor[.51,.34,.62]};
 top=Table[ix=Flatten[Position[pair["t"],x_/;x<=If[j==1,70.,25.]]];
 traces=Table[pair["radiation_"<>sign][[ix,3-j]],{sign,{"plus","minus"}}];norm=Max[Abs[traces]];
 ListLinePlot[Transpose[{pair["t"][[ix]],Re[#]/norm}]&/@traces,Evaluate[Sequence@@plotStyle],PlotRange->All,
 PlotStyle->{Directive[RGBColor[.07,.36,.57],Thickness[.004]],Directive[RGBColor[.78,.34,.18],Dashed,Thickness[.004]]},
 PlotLabel->If[j==1,"(a) Infinity radiation","(b) Horizon radiation"],FrameLabel->{"Time tau/M","Normalized radiation"},
 PlotLegends->Placed[{"X+","X-"},Above]],{j,2}];
 bottom=Table[ListLinePlot[Transpose[{t,100 Abs[ref[[All,j]]]/peak[[j]]}],Evaluate[Sequence@@plotStyle],PlotRange->{{0,70},{0,12}},
 PlotStyle->colors[[j]],PlotLabel->"("<>FromCharacterCode[98+j]<>") Cubic difference ("<>ToString[5+j]<>","<>ToString[5+j]<>")",
 FrameLabel->{"Time tau/M","Difference / common peak (%)"}],{j,2}];fig2=GraphicsGrid[{top,bottom},Spacings->{0,.3},ImageSize->820];
 top=Table[Show[ListLinePlot[Transpose[{t,Re[ref[[All,j]]]/peak[[j]]}],Evaluate[Sequence@@plotStyle],PlotStyle->GrayLevel[.25],
 PlotRange->All,FrameLabel->{"Time tau/M","Re(delta h) / common peak"},PlotLabel->"("<>FromCharacterCode[96+j]<>") Cubic difference ("<>ToString[5+j]<>","<>ToString[5+j]<>")"],
 ListPlot[Transpose[{t[[41;;-1;;100]],Re[pred[[41;;-1;;100,j]]]/peak[[j]]}],PlotStyle->colors[[j]],PlotMarkers->{"\[EmptyCircle]",7}]],{j,2}];
 bottom=Table[ListLogPlot[Transpose[{Rest[t],Rest[Abs[pred[[All,j]]-ref[[All,j]]]]/deltaPeak[[j]]}],Joined->True,Evaluate[Sequence@@plotStyle],
 PlotRange->{{0,70},{10^-8,2 10^-4}},PlotStyle->colors[[j]],FrameLabel->{"Time tau/M","Error / peak difference"},
 PlotLabel->"("<>FromCharacterCode[98+j]<>") Relative L2 error: "<>ToString[ScientificForm[(l2[pred-ref,t]/l2[ref,t])[[j]],3],StandardForm]],{j,2}];
 fig4=Column[{GraphicsGrid[{top,bottom},Spacings->{0,.3},ImageSize->820],Style["Solid: forward evolution. Open circles: reconstruction from boundary radiation.",12,FontFamily->"Times"]},Alignment->Center];
 maxSpin=Min[Table[Last[modeTable[ell,n]][[1]],{ell,{2,4,6}},{n,0,7}]];
 spins=Select[modeTable[2,0][[All,1]],.97<=#<=maxSpin&];w2=First[mode[2,0,spins]];{w4,b4}=mode[4,0,spins];{w6,b6}=mode[6,0,spins];
 dp=2 w2-w4;dc=w4+w2-w6;da=3 w2-w6;kappa=Sqrt[1-spins^2]/(2(1+Sqrt[1-spins^2]));
 targets=kerrTargets[];curves={{Abs[dp],Abs[dc],Abs[da]},{Abs[b4],Abs[b6],Abs[b4 b6/(dp da)]}};
 titles={"(a) Near-extremal Kerr detunings","(b) Linear pole factors (fixed convention)"};
 labels={{"Delta p","Delta c","Delta A"},{"B440","B660","|K0|"}};
 plots=Table[ListLogLogPlot[Transpose[{kappa,#}]&/@curves[[j]],Joined->True,Evaluate[Sequence@@plotStyle],PlotRange->All,
 FrameLabel->{"Surface gravity kappa M","Magnitude"},FrameTicks->{{Automatic,None},{{.001,.01,.05},None}},PlotLabel->titles[[j]],PlotLegends->Placed[labels[[j]],Above]],{j,2}];
 AppendTo[plots,ListLogPlot[Transpose[{Lookup[targets,"a"],Lookup[targets,#]}]&/@{"Delta_p_over_440_spacing","Delta_A_over_660_spacing"},
 Joined->True,Evaluate[Sequence@@plotStyle],PlotRange->All,FrameLabel->{"Spin chi","Detuning / overtone spacing"},PlotLabel->"(c) Pole isolation",
 PlotLegends->Placed[{"Production","Total"},Above],GridLines->{Automatic,{1}}]];
 AppendTo[plots,ListLogPlot[Transpose[{Lookup[targets,"a"],Lookup[targets,#]}]&/@{"term441_over_440_abs","term661_over_660_abs"},
 Joined->True,Evaluate[Sequence@@plotStyle],PlotRange->All,FrameLabel->{"Spin chi","Equal-vertex pole-term ratio"},PlotLabel->"(d) Overtone stress test",
 PlotLegends->Placed[{"|t441/t440|","|t661/t660|"},Above],GridLines->{Automatic,{1}}]];
 kerr=GraphicsGrid[Partition[plots,2],Spacings->{0,.3},ImageSize->820];dir=path["results/figures"];If[!DirectoryQ[dir],CreateDirectory[dir]];
 names={"fig2_preparation_difference","fig4_waveform_reconstruction","kerr_exchange"};
 MapThread[(Export[FileNameJoin[{dir,#1<>".pdf"}],#2,"PDF"];Export[FileNameJoin[{dir,#1<>".png"}],#2,"PNG",ImageResolution->160])&,{names,{fig2,fig4,kerr}}];
 FileNameJoin[{dir,#<>".pdf"}]&/@names];
