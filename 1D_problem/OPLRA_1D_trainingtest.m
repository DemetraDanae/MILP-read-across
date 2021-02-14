%%%%%Import files
data_file = 'FILENAME.csv'; 
data=readtable(data_file,'Delimiter',';','ReadRowNames',true);
totalsamples=size(data,1);
data_sc=scaling(data);
data_sc(:,1) = table2array(data(:,1)); 
sampleNames=data.Properties.RowNames;
variableNames=data.Properties.VariableNames;

%%%%%external validation
samples=round(0.75*totalsamples); 
test_samples=totalsamples-samples;
[training_samples,Test_samples]=kenstone(data_sc(:,2:end),samples);    
training_data=data_sc(training_samples,:);                  
test_data=data_sc(Test_samples,:);

exp_test=test_data(:,1);
Asf_test=test_data(:,2:end);            

Asf=training_data(:,2:end);
exp=training_data(:,1);
features=size(Asf,2);

regions=1;
U=10;
U2=abs(sum(exp));
l=0.01;
beta=0.05;
epsilon=0.05;
Z=0.5;

%%%%%simple linear regression for R=1 training 
simplelinear=fitlm(Asf,exp);
initialpredictions=predict(simplelinear,Asf);
MAE_in=sum(abs(initialpredictions-exp))/samples;
W=simplelinear.Coefficients.Estimate(2:end);   
REG_in=sum(abs(W));
z_in=MAE_in+l*REG_in;

ERRORcurrent=z_in;
ERRORold=Inf;
ERRORtmp=Inf;

%%%%%simple linear regression for R=1 test 
initialpredictions_test=predict(simplelinear,Asf_test);
MAE_in_test=sum(abs(initialpredictions_test-exp_test))/test_samples;
REG_in=sum(abs(W));
z_in_test=MAE_in_test+l*REG_in;

ERRORcurrent_test=z_in_test;
ERRORold_test=Inf;
ERRORtmp_test=Inf;

regions=regions+1;
fbest=0;
step=1'
ERRORcurtest_table(step,1)=ERRORcurrent_test;
step=step+1;

%%%%%R21 training
final1=[exp,initialpredictions];
R1=corrcoef(final1(:,2),final1(:,1));
R21=R1(1,2)^2;

%%%%%R21 test
final1_test=[exp_test,initialpredictions_test];
R21_test=1-sum((final1_test(:,1)-final1_test(:,2)).^2)/sum((final1_test(:,1)-mean(exp)).^2); %external%q2
 
variableNames_Asf=data(:,2:end).Properties.VariableNames;

for p_f=1:feature
Frs=binvar(regions,samples);
Br=sdpvar(regions,1);
Pred=sdpvar(regions,samples);
Xrf=sdpvar(regions,features);
Es=sdpvar(samples,1);
Ers=sdpvar(regions,samples);
Wrf=sdpvar(regions,features);
Wrfpos=sdpvar(regions,features);
MAE=sdpvar(1);
REG=sdpvar(1);

constraints=[];

constraints=[constraints,MAE==sum(Es)/samples];

for re=1:regions
    for fe=1:features       
        constraints=[constraints,Wrfpos(re,fe)>=Wrf(re,fe)];
        constraints=[constraints,Wrfpos(re,fe)>=-Wrf(re,fe)];
    end
end

constraints=[constraints, REG==sum(sum(Wrfpos,2))]; 

for re=2:regions
    constraints=[constraints,Xrf(re,p_f)>=Xrf(re-1,p_f)+epsilon];
end 

constraints=[constraints,Xrf(regions,p_f)==1];
constraints=[constraints,Xrf(1,p_f)>=epsilon];

for i=1:samples
    constraints=[constraints,sum(Frs(:,i))==1];
end 
  
for re=2:regions
    for i=1:samples
        constraints=[constraints,Xrf(re-1,p_f)+epsilon-U*(1-Frs(re,i))<= Asf(i,p_f)];
    end
end 

for re=1:regions-1 
    for i=1:samples
        constraints=[constraints,Asf(i,p_f)<=Xrf(re,p_f)-epsilon+U*(1-Frs(re,i))];
    end 
end 

for re=1:regions
   for i=1:samples
        constraints=[constraints,Pred(re,i)==sum(Asf(i,:).*Wrf(re,:))+Br(re,1)];
   end 
end 

for re=1:regions
    for i=1:samples 
      constraints=[constraints,Es(i,1)>=0];
      constraints=[constraints,Ers(re,i)<=U2*(Frs(re,i))];
      constraints=[constraints,Ers(re,i)>=exp(i,1)-Pred(re,i)-U2*(1-Frs(re,i))];
      constraints=[constraints,Ers(re,i)>=Pred(re,i)-exp(i,1)-U2*(1-Frs(re,i))];
      constraints=[constraints,Es(i,1)>=Ers(re,i)]; 
    end
end

z=MAE+l*REG;

options=sdpsettings('solver','mosek'); 

readAcross=optimize(constraints,z,options);

value(Xrf(:,p_f))
value(z);

if value(z)<ERRORtmp && value(z)~=0
    ERRORtmp=value(z);
    zopt=value(z);
    bestMAE=value(MAE);
    bestREG=value(REG);
    bestEs=value(Es);
    fbest=p_f;
    bestFrs=value(Frs);
    bestPred=value(Pred);
    bestXrf=value(Xrf);
    bestWrf=value(Wrf);
    bestBr=value(Br);
end

end

finalregions=regions;

z_test=OPLRA_1D_test(fbest,features,test_samples,exp_test,finalregions,Asf_test,l,bestXrf,bestWrf,bestBr);

ERRORold=ERRORcurrent;
ERRORcurrent=ERRORtmp;

ERRORold_test=ERRORcurrent_test;
ERRORcurrent_test=z_test;

ERRORcurtest_table(step,1)=ERRORcurrent_test;
step=step+1;

while ERRORcurrent_test<(1-beta)*ERRORold_test && ERRORcurrent_test~=0
    clear p_f
    clear Frs
    clear Xrf
    clear Br
    clear Pred
    clear Wrf
    clear Wrfpos
    clear Es
    clear Ers
    clear MAE
    clear REG 
    
    regions=regions+1;
    
Frs=binvar(regions,samples);
Br=sdpvar(regions,1);
Pred=sdpvar(regions,samples);
Xrf=sdpvar(regions,features);
Es=sdpvar(samples,1);
Ers=sdpvar(regions,samples);
Wrf=sdpvar(regions,features);
Wrfpos=sdpvar(regions,features);
MAE=sdpvar(1);
REG=sdpvar(1);

constraints=[];

constraints=[constraints,MAE==sum(Es)/samples];
for re=1:regions
    for fe=1:features       
        constraints=[constraints,Wrfpos(re,fe)>=Wrf(re,fe)];
        constraints=[constraints,Wrfpos(re,fe)>=-Wrf(re,fe)];
    end
end

constraints=[constraints, REG==sum(sum(Wrfpos,2))]; 

for re=2:regions
    constraints=[constraints,Xrf(re,fbest)>=Xrf(re-1,fbest)+epsilon],
end 

constraints=[constraints,Xrf(regions,fbest)==1];
constraints=[constraints,Xrf(1,fbest)>=epsilon];

for i=1:samples
    constraints=[constraints,sum(Frs(:,i))==1];
end 

for re=2:regions
    for i=1:samples
        constraints=[constraints,Xrf(re-1,fbest)+epsilon-U*(1-Frs(re,i))<= Asf(i,fbest)];
    end
end 

for re=1:regions-1 
    for i=1:samples
        constraints=[constraints,Asf(i,fbest)<=Xrf(re,fbest)-epsilon+U*(1-Frs(re,i))];
    end 
end 

for re=1:regions
   for i=1:samples
        constraints=[constraints,Pred(re,i)==sum(Asf(i,:).*Wrf(re,:))+Br(re,1)];
   end 
end 

for re=1:regions
    for i=1:samples        
        constraints=[constraints,Es(i,1)>=0];
        constraints=[constraints,Ers(re,i)<=U2*(Frs(re,i))];
        constraints=[constraints,Ers(re,i)>=exp(i,1)-Pred(re,i)-U2*(1-Frs(re,i))];
        constraints=[constraints,Ers(re,i)>=Pred(re,i)-exp(i,1)-U2*(1-Frs(re,i))];
        constraints=[constraints,Es(i,1)>=Ers(re,i)]; 
    end
end

z=MAE+l*REG;

options=sdpsettings('solver','mosek');

readAcross=optimize(constraints,z,options);

value(Xrf(:,fbest))

z_test=OPLRA_1D_test(fbest,features,test_samples,exp_test,regions,Asf_test,l,value(Xrf),value(Wrf),value(Br));

ERRORold=ERRORcurrent;
ERRORcurrent=value(z);

ERRORold_test=ERRORcurrent_test;
ERRORcurrent_test=z_test;

ERRORcurtest_table(step,1)=ERRORcurrent_test;
step=step+1;

if ERRORcurrent_test<(1-beta)*ERRORold_test && ERRORcurrent_test~=0
    bestFrs=value(Frs);
    bestPred=value(Pred);
    bestXrf=value(Xrf);
    bestWrf=value(Wrf);
    bestBr=value(Br);
    zopt=value(z);
    bestMAE=value(MAE);
    bestREG=value(REG);
    bestEs=value(Es);
    finalregions=regions;
end 

end 

Ps=value(bestFrs.*bestPred);
Ps=Ps';
Ps1=sum(Ps,2);

%%%%%R2 training
final=[exp,Ps1];
R=corrcoef(final(:,2),final(:,1));
R2=R(1,2)^2;

for re=1:finalregions
    if sum(bestFrs(re,:))<=1.5
        R2_re(re)=NaN;
    else 
    clear final2
    clear joint
    clear exp_reduced 
    exp_reduced=exp(:,1).*(bestFrs(re,:)');
    joint=[exp_reduced,Ps(:,re)];
    rowfinal2=1;
    for i=1:samples
    if bestFrs(re,i)==1 
       final2(rowfinal2,:)=joint(i,:);
       rowfinal2=rowfinal2+1;
    end 
    end 
    R_re=corrcoef(final2(:,2),final2(:,1));
    R2_re(re)=R_re(1,2)^2;
    end 
end  

%RMSE train
RMSE = sqrt(mean((final(:,1) - final(:,2)).^2));

Frs_test=zeros(finalregions,test_samples);

for i=1:test_samples
    for re=2:finalregions-1 
        if Asf_test(i,fbest)>=bestXrf(re-1,fbest) && Asf_test(i,fbest)<=bestXrf(re,fbest)
            Frs_test(re,i)=1;
        end
    end
end

for i=1:test_samples 
    re=1;
    if bestXrf(re,fbest )>Asf_test(i,fbest)
        Frs_test(re,i)=1;
    end
end

for i=1:test_samples
    re=finalregions;
    if bestXrf(re-1,fbest )<Asf_test(i,fbest) && Asf_test(i,fbest)<=bestXrf(re,fbest)
        Frs_test(re,i)=1;
    end
end

for re=1:finalregions
   for i=1:test_samples
        Pred_test(re,i)=sum(Asf_test(i,:).*bestWrf(re,:))+bestBr(re,1);
   end 
end 

Ps_test=value(Frs_test.*Pred_test);
Ps_test=Ps_test';
Ps1_test=sum(Ps_test,2);

for i=1:test_samples
  E(i,1)=abs(exp_test(i,1)-Ps1_test(i,1));      
end

MAE_test=sum(E)/test_samples;

for re=1:finalregions
    for fe=1:features       
        Wpos_t(re,fe)=abs(bestWrf(re,fe));
    end
end

REG_test=sum(sum(Wpos_t,2));
z_test=MAE_test+l*REG_test;

%%%%%q2 test 
final_test=[exp_test,Ps1_test];
q2 = 1-sum((final_test(:,1)-final_test(:,2)).^2)/sum((final_test(:,1)-mean(exp)).^2); %external%q2
q2_test = 1-sum((final_test(:,1)-final_test(:,2)).^2)/sum((final_test(:,1)-mean(final_test(:,1))).^2); %mean values of the test endpoint

%%%%%R2 test
R_test=corrcoef(final_test(:,2),final_test(:,1));
R2_test=R_test(1,2)^2;

for re=1:finalregions
    if sum(Frs_test(re,:))<=1.5
        R2_test_re(re)=NaN;
    else 
        clear final2_test
        clear joint_test
        clear exp_reduced_test
        exp_reduced_test=exp_test(:,1).*(Frs_test(re,:)');
        joint_test=[exp_reduced_test,Ps_test(:,re)];
        rowfinal2_test=1;
        for i=1:test_samples
            if Frs_test(re,i)==1 
                final2_test(rowfinal2_test,:)=joint_test(i,:);
                rowfinal2_test=rowfinal2_test+1;
            end 
        end 
        R_test_re=corrcoef(final2_test(:,2),final2_test(:,1));
        R2_test_re(re)=R_test_re(1,2)^2;
    end  
end

%RMSE test
RMSE_test = sqrt(mean((final_test(:,1) - final_test(:,2)).^2));
