function extract_erp_measures(matrixDir, phenotypeCsv, outputDir)
% Extract the final study measures from archived subject-level bin averages.
% Does not rerun filtering, ICA, manual edits or artifact rejection.
% Input matrices: ERPtotal_{CON,PEP}_{MMN,CNV_PINV}_0p05Hz.mat.
% Output directory must not exist. No input file is ever modified.
assert(~isfolder(outputDir) && ~isfile(outputDir),'Choose a new output directory.');
channels=[3 4 5 9 10];sites=["Fz","FCz","Cz","FC3","FC4"];
weights=[150 25 150 25 25 25];stdBins=[1 3 5 6];devBins=[2 4];
tShort=[];tLongS4=[];tLongS1=[];
W=struct();ids=strings(0,1);groups=strings(0,1);
for group=["CON","PEP"]
 s=load(fullfile(matrixDir,"ERPtotal_"+group+"_MMN_0p05Hz.mat"));
 l=load(fullfile(matrixDir,"ERPtotal_"+group+"_CNV_PINV_0p05Hz.mat"));
 x=double(s.("ERPtotal_"+group+"_MMN"));y=double(l.("ERPtotal_"+group+"_CNV_PINV"));
 assert(size(x,1)==16 && size(x,3)==6);
 assert(size(y,1)==16 && size(y,3)==6);
 ts=erp_time_axis(s, "ERPtotal_"+group+"_MMN", 'mmn');
 tl=erp_time_axis(l, "ERPtotal_"+group+"_CNV_PINV", 'long');
 if isempty(tShort)
  tShort=ts;tLongS4=tl;tLongS1=tLongS4+900;
 else
  assert(numel(ts)==numel(tShort) && all(abs(ts-tShort)<1e-6), 'MMN axes differ across groups.');
  assert(numel(tl)==numel(tLongS4) && all(abs(tl-tLongS4)<1e-6), 'Long-epoch axes differ across groups.');
 end
 shortNames=string(s.subject_names(:));longNames=string(l.subject_names(:));
 [ok,order]=ismember(shortNames,longNames);assert(all(ok));y=y(:,:,:,order);
 assert(all(isfinite(x(:))) && all(isfinite(y(:))));
 num=str2double(extractBetween(shortNames,'Sujeto','_'));
 assert(all(isfinite(num)),'Adapt ID parsing to your deidentified matrix filenames.');
 prefix="C";label="CON";if group=="PEP",prefix="P";label="FEP";end
 ids=[ids;prefix+compose('%02d',num)];groups=[groups;repmat(label,numel(num),1)];
 % Rebaseline PINV per participant/channel/bin without changing accepted trials.
 yp=y-mean(y(:,tLongS4>=-100&tLongS4<=0,:,:),2);
 temp=struct();temp.mmn_std=combine(x,channels,stdBins,weights);
 temp.mmn_dev=combine(x,channels,devBins,weights);
 temp.mmn=temp.mmn_dev-temp.mmn_std;
 temp.all=combine(y,channels,1:6,weights);
 temp.pinv=combine(yp,channels,devBins,weights)-combine(yp,channels,stdBins,weights);
 for f=string(fieldnames(temp))'
  if ~isfield(W,f),W.(f)=temp.(f);else,W.(f)=cat(1,W.(f),temp.(f));end
 end
end
assert(numel(unique(ids))==numel(ids));
gm=squeeze(mean(W.mmn(:,2,:),1));gn=squeeze(mean(W.all(:,2,:),1));
idx=find(tShort>=100&tShort<=250);[~,k]=min(gm(idx));mmnPeak=tShort(idx(k));
idx=find(tLongS1>=50&tLongS1<=200);[~,k]=min(gn(idx));nPeak=tLongS1(idx(k));
idx=find(tLongS1>=150&tLongS1<=300);[~,k]=max(gn(idx));pPeak=tLongS1(idx(k));
% Common group-collapsed FCz windows, never selected on the group contrast.
mm=tShort>=mmnPeak-50&tShort<=mmnPeak+50;
nm=tLongS1>=nPeak-25&tLongS1<=nPeak+25;pm=tLongS1>=pPeak-25&tLongS1<=pPeak+25;
assert(nPeak+25<pPeak-25 && pPeak+25<300,'Review N1/P2 localizer.');
cn=tLongS1>=0&tLongS1<900;cl=tLongS1>=800&tLongS1<900;
pi=tLongS4>=300&tLongS4<=1000;pe=tLongS4>=300&tLongS4<=500;
T=table();
for i=1:numel(ids)
 for j=1:5
  m=squeeze(W.mmn(i,j,:));a=squeeze(W.all(i,j,:));p=squeeze(W.pinv(i,j,:));
  mt=tShort(mm);pt=tLongS4(pi);[mp,mk]=min(m(mm));[pp,pk]=min(p(pi));
  row=table(ids(i),groups(i),sites(j),mean(a(cn)),mean(m(mm)),mean(p(pi)),mean(a(nm)),mean(a(pm)),mean(a(cl)),mean(p(pe)),mp,mt(mk),pp,pt(pk),...
   'VariableNames',{'Subject','GROUP','Electrode','CNV','MMN','PINV','N1','P2','CNV_late','PINV_early','MMN_peak','MMN_latency','PINV_peak','PINV_latency'});
  T=[T;row];
 end
end
P=readtable(phenotypeCsv,'VariableNamingRule','preserve','TextType','string');
assert(numel(unique(P.Subject))==height(P),'Phenotype input must have one row per participant.');
[found,order]=ismember(T.Subject,P.Subject);assert(all(found),'Missing phenotype record.');
if ismember('GROUP',P.Properties.VariableNames),assert(all(P.GROUP(order)==T.GROUP));end
covs=setdiff(P.Properties.VariableNames,{'Subject','GROUP'},'stable');
T=[T P(order,covs)];mkdir(outputDir);
writetable(T,fullfile(outputDir,'analysis_input.csv'));
windows=table(["MMN";"N1";"P2"],[mmnPeak;nPeak;pPeak],[mmnPeak-50;nPeak-25;pPeak-25],[mmnPeak+50;nPeak+25;pPeak+25],...
 'VariableNames',{'Component','Peak_ms','Start_ms','End_ms'});
writetable(windows,fullfile(outputDir,'localizer_windows.csv'));
save(fullfile(outputDir,'waveforms.mat'),'W','ids','groups','sites','tShort','tLongS4','tLongS1','-v7');
disp(windows);
end

function a=combine(x,ch,bins,weights)
% Nominal bin weights reproduce the archived pipeline, not accepted-trial weights.
v=sum(x(ch,:,bins,:).*reshape(weights(bins),1,1,[],1),3)/sum(weights(bins));
a=permute(v,[4 1 2 3]);
end
