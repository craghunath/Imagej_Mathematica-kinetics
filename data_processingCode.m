
mxt=50;
st=0;
fr=0.03;frame rate 
no=1;number of frame jumps
mbot=7*10^-3;
n=0.00000013;
n2=0.000004;
radius=16*10^-3;
g=9.8066;
pow=4;
km2=1.4002637*10^-99.8085913*10^-10;
Needs["ErrorBarPlots`"]
v01[path_,nn_,nn2_,pixel_,fr1_,er1_]:=Block[{path1=path,n=nn,n2=nn2,pix=0.01/pixel,fr2=fr1,er0=er1,plt2,rplt,vplt,errfrc,mnplt,ind},SetDirectory[path1];
Needs["ErrorBarPlots`"];
flst=FileNames[];
fr=1/fr2;
indexp=Mean[Transpose[Take[#,5]&/@rdat1*pix][[1]]];
mbot=0.00716;
dat=Import[#,"Data"]&/@flst;
adat=Map[{#[[4]],#[[5]]}&,Drop[#,1]&/@dat,{2}];
bdat=Differences[#,1,(#//Length)/2]&/@adat;
rdat=Apply[Plus,bdat^2,{2}]//Sqrt;Map[Sqrt[#[[1]]^2+#[[2]]^2]&,bdat,{2}]; This also works fine
ind=er0+ #[[1]]*pix&/@rdat//Mean;
min=Min[Length[#]&/@rdat];
co=Take[#,min]*pix&/@rdat;
vf1={Differences[#]&/@co/(fr)};vf=Table[MapThread[List,{Take[co[[#]],Length[vf1[[i,#]]]],vf1[[i,#]]}]&/@Range[1,Length[vf1[[i]]],1],{i,1,Length[vf1]}];
sz=Tiny;
plt2=vf[[1]];
mnfrc=MovingAverage[TrimmedMean[Differences[#]/(fr)&/@co,.05],3];
mnplt=MapThread[List,{Take[Mean[co],Length[mnfrc]],mnfrc}];
sdfrc=Sqrt[TrimmedVariance[Differences[#]/(fr)&/@co,.1]];
sdfrc=StandardDeviation[Differences[#]/(fr)&/@co];
errfrc=Partition[Riffle[mnplt,Map[ErrorBar,(sdfrc)]],2];/Sqrt[Length[plt2]]
rplt=MapThread[List,{Drop[Flatten[Table[  Evaluate[r[t]/.(NDSolve[{mbot r''[t]+10^-6 r'[t]+10^-6 (r'[t])^2-km2/(r[t]^pow)==0,r[0]==ind,r'[0]==0},r[t],{t,st,mxt}])]  ,{t,0,5,fr}   ]   ],-1] ,Differences[Flatten[Table[  Evaluate[r[t]/.(NDSolve[{mbot r''[t]+10^-6 r'[t]+10^-6 (r'[t])^2-km2/(r[t]^pow)==0,r[0]==ind,r'[0]==0},r[t],{t,st,mxt}])]  ,{t,0,5,fr}  ]   ] ]/(fr)}];
vplt=MapThread[List,{Drop[Flatten[Table[  Evaluate[r[t]/.(NDSolve[{mbot r''[t]+n r'[t]+n2 (r'[t])^2-km2/(r[t]^pow)==0,r[0]==ind,r'[0]==0},r[t],{t,st,mxt}])]  ,{t,0,mxt,fr}   ]   ],-1] ,Differences[Flatten[Table[  Evaluate[r[t]/.(NDSolve[{mbot r''[t]+n r'[t]+n2 (r'[t])^2-km2/(r[t]^pow)==0,r[0]==ind,r'[0]==0},r[t],{t,st,mxt}])]  ,{t,0,mxt,fr}   ]   ]]/(fr)}];
a2=ListPlot[plt2,PlotRange->{{0,0.25},{0,0.05}},PlotStyle->PointSize[0.003]];
a1=ListLinePlot[{rplt,vplt},PlotRange->{{0,0.25},{0,0.05}},PlotStyle->{{Darker[Green,0.26],Thickness[0.0025]},{Blue,Thickness[0.002]}}];

a3=ErrorListPlot[errfrc,PlotStyle->{PointSize[Tiny],Darker[Red,0.15]},PlotRange->{{0,0.25},{0,0.05}}];
a4=ListLinePlot[mnplt,PlotRange->{{0,0.25},{0,0.05}},PlotStyle->{PointSize[0.002],Black}];

Show[a3,a2,a4,a1,PlotRange->{{0.02,.3},{0,0.05}},GridLines->{{0.035}~Join~Range[0,0.3,0.01], Automatic},Ticks->{Range[0,0.3,0.01], Automatic},Axes->True,AxesOrigin->{0.1,0},AxesStyle->Directive[Darker[Blue,0.3],Bold,12],AxesLabel->{"Distance (m)","Speed (\!\(\*SuperscriptBox[\(ms\), \(-1\)]\))"},PlotLabel->{ ToString[n]<>"(n)",ToString[ n2]<>"(\!\(\*SubscriptBox[\(n\), \(2\)]\))", ToString[N[Round[ind,0.0001]]]<>"(initial_distance)" },LabelStyle->Directive[Bold,Red,14]]
];
adr={"/path/to/datafile/folder1/","/path/to/datafile/folder2/","/path/to/datafile/folder3/"};
ppcm={53,47,47};ppcm is pixel per cm data from calibration
fhu=v01[adr[[#]],0.0006,.1,ppcm[[#]],60,0]&/@{1,2,3};
fhu//Show



