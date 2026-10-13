QCSCheckData[]:=Module[{pair,rec,t,ref,pred,err,peak,sep,percent,ix,rp,rm,mid,coeff,base,a,energies,report},
 pair=QCSData["source/reduced_input/paired_preparations_fine.mat"];
 rec=QCSData["source/reduced_input/completion/recovery_waveforms.mat"];t=rec["t"];
 ref=rec["reference_plus"]-rec["reference_minus"];pred=rec["prediction_plus"]-rec["prediction_minus"];
 err=l2[pred-ref,t]/l2[ref,t];peak=MapThread[Max,{Max/@Transpose[Abs[rec["reference_plus"]]],Max/@Transpose[Abs[rec["reference_minus"]]]}];
 sep=2 l2[ref,t]/(l2[rec["reference_plus"],t]+l2[rec["reference_minus"],t]);percent=100(Max/@Transpose[Abs[ref]])/peak;
 ix=First[weights[pair["t"],0,400]];rp=pair["radiation_plus"][[ix,2]];rm=pair["radiation_minus"][[ix,2]];mid=(rp+rm)/2;
 coeff=hiddenCoefficients[];base=ConstantArray[0.,36];base[[1]]=.25;base/=Sqrt[energy[base]];
 a=Last[json["source/reduced_input/pair_results.json"]]["hidden_amplitude"];energies=energy[.5 base+# a coeff]&/@{1,-1};
 report=<|"scope"->"Diagnostics recomputed in Wolfram Language from the retained arrays.",
 "infinity_mismatch_to_midpoint_0_400M"->l2[rp-mid,pair["t"][[ix]]]/l2[mid,pair["t"][[ix]]],
 "relative_complex_L2_difference_errors"->err,"common_cubic_peaks"->peak,"pointwise_difference_peak_percent"->percent,
 "complex_L2_separation_percent"->100 sep,"hidden_direction_energy"->energy[coeff],"preparation_energies"->energies|>;
 save["data_checks",report];report];
