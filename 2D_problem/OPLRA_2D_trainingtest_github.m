
%%%%%Import files
%data_file = 'filename.csv';
data=readtable(data_file,'Delimiter',';','ReadRowNames',true);
data=data(:,:);
totalsamples=size(data,1);
data_sc=scaling(data);
data_sc(:,1) = table2array(data(:,1)); %endpoints not normalized
sampleNames=data.Properties.RowNames;
variableNames=data.Properties.VariableNames;

%%%%%external validation
samples=round(0.66*totalsamples); 
test_samples=totalsamples-samples;
[training_samples,Test_samples]=kenstone(data_sc(:,2:end),samples);     %no endpoint
training_data=data_sc(training_samples,:);                  
test_data=data_sc(Test_samples,:);
trainingsampleNames=sampleNames(training_samples,1);

exp_test=test_data(:,1);
Asf1_test=test_data(:,2:40);          
Asf2_test=test_data(:,41:end);   
Asf_test=test_data(:,2:end);

Asf1=training_data(:,2:40);
Asf2=training_data(:,41:end);
Asf=training_data(:,2:end);
exp=training_data(:,1);
%exp=exp(randperm(length(exp)));        %y-scrambling
features1=size(Asf1,2);
features2=size(Asf2,2);
features=size(Asf,2);

regions1=2;
regions2=2;
U=10;
U2=abs(sum(exp));
l=0.005;
beta=0.05;
epsilon=0.05;
Z=0.5;

%%%%%simple linear regression for R=1 training
simplelinear=fitlm(Asf,exp);
initialpredictions=predict(simplelinear,Asf);
MAE_in=sum(abs(initialpredictions-exp))/samples;
W=simplelinear.Coefficients.Estimate(2:end);   %no intercept
REG_in=sum(abs(W));
z_in=MAE_in+l*REG_in;

ERRORcurrent=z_in;  
ERRORold=Inf;
ERRORtmp=Inf;

%%%%%simple linear regression for R=1 test
initialpredictions_test=predict(simplelinear,Asf_test);
MAE_in_test=sum(abs(initialpredictions_test-exp_test))/test_samples;
REG_in_test=sum(abs(W));
z_in_test=MAE_in_test+l*REG_in_test;
ERRORcurrent_test=z_in_test;  
ERRORold_test=Inf;
ERRORtmp_test=Inf;

fbest1=0;
fbest2=0;
step=1;
ERRORcur_table(step,1)=ERRORcurrent;
step=step+1;

%%%%%R2 training
final1=[exp,initialpredictions];
R1=corrcoef(final1(:,2),final1(:,1));
R21=R1(1,2)^2;

%%%%%R21 test 
final1_test=[exp_test,initialpredictions_test];
R21_test=1-sum((final1_test(:,1)-final1_test(:,2)).^2)/sum((final1_test(:,1)-mean(exp)).^2); %external%q2
 
variableNames_Asf=data(:,2:end).Properties.VariableNames;

reliability=domain(Asf,Asf_test,Z);

for p_f1=1:features1
    for p_f2=1:features2
   p_f1
   p_f2
Fr1s=binvar(regions1,samples);
Fr2s=binvar(regions2,samples);
Br1r2=sdpvar(regions1,regions2,'full');
Pred=sdpvar(regions1,regions2,samples,'full');
Xrf1=sdpvar(regions1,features1);
Xrf2=sdpvar(regions2,features2);
Ers1=sdpvar(regions1,samples);
Ers2=sdpvar(regions2,samples);
Es=sdpvar(samples,1);
Wrf=sdpvar(regions1,regions2,features,'full');
Wrfpos=sdpvar(regions1,regions2,features,'full');
MAE=sdpvar(1);
REG=sdpvar(1);

constraints=[];

constraints=[constraints,MAE==sum(Es)/samples];

for re1=1:regions1
    for re2=1:regions2
        for fe=1:features       
        constraints=[constraints,Wrfpos(re1,re2,fe)>=Wrf(re1,re2,fe)];
        constraints=[constraints,Wrfpos(re1,re2,fe)>=-Wrf(re1,re2,fe)];
        end
    end
end

constraints=[constraints, REG==sum(sum(sum(Wrfpos)))]; 

for re=2:regions1
    constraints=[constraints,Xrf1(re,p_f1)>=Xrf1(re-1,p_f1)+epsilon];
end 

constraints=[constraints,Xrf1(regions1,p_f1)==1];
constraints=[constraints,Xrf1(1,p_f1)>=epsilon];

for re=2:regions2
    constraints=[constraints,Xrf2(re,p_f2)>=Xrf2(re-1,p_f2)+epsilon];
end 

constraints=[constraints,Xrf2(regions2,p_f2)==1];
constraints=[constraints,Xrf2(1,p_f2)>=epsilon];

for i=1:samples
    constraints=[constraints,sum(Fr1s(:,i))==1];
end 

for i=1:samples
    constraints=[constraints,sum(Fr2s(:,i))==1];
end

for re=2:regions1
    for i=1:samples
        constraints=[constraints,Xrf1(re-1,p_f1)+epsilon-U*(1-Fr1s(re,i))<= Asf1(i,p_f1)];
    end
end 

for re=1:regions1-1 
    for i=1:samples
        constraints=[constraints,Asf1(i,p_f1)<=Xrf1(re,p_f1)-epsilon+U*(1-Fr1s(re,i))];
    end 
end 

for re=2:regions2
    for i=1:samples
        constraints=[constraints,Xrf2(re-1,p_f2)+epsilon-U*(1-Fr2s(re,i))<= Asf2(i,p_f2)];
    end
end 

for re=1:regions2-1 
    for i=1:samples
        constraints=[constraints,Asf2(i,p_f2)<=Xrf2(re,p_f2)-epsilon+U*(1-Fr2s(re,i))];
    end 
end 


for re1=1:regions1
    for re2=1:regions2
        for i=1:samples
         constraints=[constraints,Pred(re1,re2,i)==sum(Asf(i,:).*Wrf(re1,re2,:))+Br1r2(re1,re2)];
        end 
    end 
end 

for re1=1:regions1
   for re2=1:regions2
        for i=1:samples
            constraints=[constraints,Es(i,1)>=0];
            constraints=[constraints,Ers1(re1,i)<=U2*(Fr1s(re1,i))];
            constraints=[constraints,Ers2(re2,i)<=U2*(Fr2s(re2,i))];
            constraints=[constraints,Ers1(re1,i)>=exp(i,1)-Pred(re1,re2,i)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
            constraints=[constraints,Ers1(re1,i)>=Pred(re1,re2,i)-exp(i,1)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
            constraints=[constraints,Ers2(re2,i)>=exp(i,1)-Pred(re1,re2,i)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
            constraints=[constraints,Ers2(re2,i)>=Pred(re1,re2,i)-exp(i,1)-U2*(1-Fr1s(re1,i))-U2*(1-Fr2s(re2,i))];
            constraints=[constraints,Es(i,1)>=Ers1(re1,i)];
            constraints=[constraints,Es(i,1)>=Ers2(re2,i)];  
        end
    end
end

z=MAE+l*REG;

options=sdpsettings('solver','mosek'),

readAcross=optimize(constraints,z,options);

value(Xrf1(:,p_f1))
value(Xrf2(:,p_f2))
value(Fr1s)
value(Fr2s)

if value(z)<ERRORtmp && value(z)~=0
    ERRORtmp=value(z);
    zopt=value(z);
    bestMAE=value(MAE);
    bestREG=value(REG);
    bestEs=value(Es);
    fbest1=p_f1;
    fbest2=p_f2;
    bestFr1s=value(Fr1s);
    bestFr2s=value(Fr2s);
    bestPred=value(Pred);
    bestXrf1=value(Xrf1);
    bestXrf2=value(Xrf2);
    bestWrf=value(Wrf);
    bestBr1r2=value(Br1r2);
end

end
end

finalregions1=regions1;
finalregions2=regions2;

[z_test, E_test]=OPLRA_2D_test_final(fbest1,fbest2,test_samples,exp_test,finalregions1,finalregions2,Asf_test,Asf1_test,Asf2_test,l,bestXrf1,bestXrf2,bestWrf,bestBr1r2);

ERRORold=ERRORcurrent;
ERRORcurrent=ERRORtmp;

ERRORold_test=ERRORcurrent_test;
ERRORcurrent_test=z_test;

ERRORcurtest_table(step,1)=ERRORcurrent_test;
step=step+1;

bestF=zeros(regions1,regions2,samples);

for i=1:samples
    for re1=1:regions1
        for re2=1:regions2
            if bestFr1s(re1,i)>=0.5 && bestFr2s(re2,i)>=0.5
                bestF(re1,re2,i)=1;
            end 
        end
    end
end

while ERRORcurrent_test<(1-beta)*ERRORold_test && ERRORcurrent_test~=0

clear MAE
clear REG
clear Es
clear Ers1
clear Ers2
clear Pred
clear Wrf
clear Wrfpos
clear Fr1s
clear Fr2s
clear Xrf1
clear Xrf2

r_a1=regions1+1;
r_a2=regions2;
[z_a,MAE_a,REG_a,Es_a,Fr1s_a,Fr2s_a,Br1r2_a,Wrf_a,Pred_a,Xrf1_a,Xrf2_a]=OPLRA_2D_final(fbest1,fbest2,exp,r_a1,r_a2,Asf,Asf1,Asf2,epsilon,U,U2,l);

r_b1=regions1;
r_b2=regions2+1;
[z_b,MAE_b,REG_b,Es_b,Fr1s_b,Fr2s_b,Br1r2_b,Wrf_b,Pred_b,Xrf1_b,Xrf2_b]=OPLRA_2D_final(fbest1,fbest2,exp,r_b1,r_b2,Asf,Asf1,Asf2,epsilon,U,U2,l);

r_c1=regions1+1;
r_c2=regions2+1;
[z_c,MAE_c,REG_c,Es_c,Fr1s_c,Fr2s_c,Br1r2_c,Wrf_c,Pred_c,Xrf1_c,Xrf2_c]=OPLRA_2D_final(fbest1,fbest2,exp,r_c1,r_c2,Asf,Asf1,Asf2,epsilon,U,U2,l);

minz=Inf;
if z_a<minz && z_a~=0
    minz=z_a;
    regions1=r_a1;
    regions2=r_a2;
    z=z_a;
    MAE=MAE_a;
    REG=REG_a;
    Es=Es_a;
    Fr1s=Fr1s_a;
    Fr2s=Fr2s_a;
    Pred=Pred_a;
    Xrf1=Xrf1_a;
    Xrf2=Xrf2_a;
    Wrf=Wrf_a;
    Br1r2=Br1r2_a;      
end

if z_b<minz && z_b~=0
    minz=z_b;
    regions1=r_b1;
    regions2=r_b2;
    z=z_b;
    MAE=MAE_b;
    REG=REG_b;
    Es=Es_b;
    Fr1s=Fr1s_b;
    Fr2s=Fr2s_b;
    Pred=Pred_b;
    Xrf1=Xrf1_b;
    Xrf2=Xrf2_b;
    Wrf=Wrf_b;
    Br1r2=Br1r2_b;    
end

if z_c<minz && z_c~=0
    minz=z_c;
    regions1=r_c1;
    regions2=r_c2;
    z=z_c;
    MAE=MAE_c;
    REG=REG_c;
    Es=Es_c;
    Fr1s=Fr1s_c;
    Fr2s=Fr2s_c;
    Pred=Pred_c;
    Xrf1=Xrf1_c;
    Xrf2=Xrf2_c;
    Wrf=Wrf_c;
    Br1r2=Br1r2_c;    
end

ERRORold=ERRORcurrent;
ERRORcurrent=minz;

F=zeros(regions1,regions2,samples);

for i=1:samples
    for re1=1:regions1
        for re2=1:regions2
            if Fr1s(re1,i)>=0.5 && Fr2s(re2,i)>=0.5
                F(re1,re2,i)=1;
            end 
        end
    end
end

[z_test, E_test]=OPLRA_2D_test_final(fbest1,fbest2,test_samples,exp_test,regions1,regions2,Asf_test,Asf1_test,Asf2_test,l,Xrf1,Xrf2,Wrf,Br1r2);

ERRORold_test=ERRORcurrent_test;
ERRORcurrent_test=z_test;

if ERRORcurrent_test<(1-beta)*ERRORold_test
    clear bestF
    clear zopt
    bestF=F;
    bestPred=Pred;
    zopt=z;
    bestMAE=MAE;
    bestREG=REG;
    bestEs=Es;
    bestXrf1=Xrf1;
    bestXrf2=Xrf2;
    bestWrf=Wrf;
    bestBr1r2=Br1r2;
    finalregions1=regions1;
    finalregions2=regions2;
end 

end 

Ps=value(bestF.*bestPred);
A=sum(sum(Ps,1));
Ps1(:,1)=A(1,1,:);

%%%%%R2 training
final=[exp,Ps1];
R=corrcoef(final(:,2),final(:,1));
R2=R(1,2)^2;

samplestoregions=zeros(samples,1);
for re1=1:finalregions1
    for re2=1:finalregions2
    if sum(bestF(re1,re2,:),3)<=1
        R2_re(re1,re2)=NaN;
    else 
    clear final2
    clear joint
    clear exp_reduced 
    samplestoregions(:,1)=bestF(re1,re2,:);
    exp_reduced=exp(:,1).*samplestoregions(:,1);
    Ps2(:,1)=Ps(re1,re2,:);
    joint=[exp_reduced,Ps2(:,1)];
    rowfinal2=1;
    for i=1:samples
    if bestF(re1,re2,i)==1 
       final2(rowfinal2,:)=joint(i,:);
       rowfinal2=rowfinal2+1;
    end 
    end 
    R_re=corrcoef(final2(:,2),final2(:,1));
    R2_re(re1,re2)=R_re(1,2)^2;
    end 
    end 
end 

%RMSE train
RMSE = sqrt(mean((final(:,1) - final(:,2)).^2));

F_test=zeros(finalregions1,finalregions2,test_samples);
Fr1s_test=zeros(finalregions1,test_samples);
Fr2s_test=zeros(finalregions2,test_samples);

for i=1:test_samples
    for re1=2:finalregions1-1 
         if Asf1_test(i,fbest1)>=bestXrf1(re1-1,fbest1) && Asf1_test(i,fbest1)<=bestXrf1(re1,fbest1)
            Fr1s_test(re1,i)=1;
        end
    end
end

for i=1:test_samples
    for re2=2:finalregions2-1 
         if Asf2_test(i,fbest2)>=bestXrf2(re2-1,fbest2) && Asf2_test(i,fbest2)<=bestXrf2(re2,fbest2)
            Fr2s_test(re2,i)=1;
        end
    end
end

for i=1:test_samples 
    re=1 ;
    if bestXrf1(re,fbest1)>Asf1_test(i,fbest1)
        Fr1s_test(re,i)=1;
    end
end

for i=1:test_samples 
    re=1 ;
    if bestXrf2(re,fbest2)>Asf2_test(i,fbest2)
        Fr2s_test(re,i)=1;
    end
end

for i=1:test_samples
    re=finalregions1;
    if bestXrf1(re-1,fbest1)<Asf1_test(i,fbest1) && Asf1_test(i,fbest1)<=bestXrf1(re,fbest1)
        Fr1s_test(re,i)=1;
    end
end

for i=1:test_samples
    re=finalregions2;
    if bestXrf2(re-1,fbest2 )<Asf2_test(i,fbest2) && Asf2_test(i,fbest2)<=bestXrf2(re,fbest2)
        Fr2s_test(re,i)=1;
    end
end

for re1=1:finalregions1
    for re2=1:finalregions2
      bestW(1,:)=bestWrf(re1,re2,:);
     for i=1:test_samples
        Pred_test(re1,re2,i)=sum(Asf_test(i,:).*bestW(1,:))+bestBr1r2(re1,re2);
     end 
    end
end 

for i=1:test_samples
    for re1=1:finalregions1
        for re2=1:finalregions2
          if Fr1s_test(re1,i)>=0.5 && Fr2s_test(re2,i)>=0.5
               F_test(re1,re2,i)=1;
          end 
       end
    end
end

Ps_test=value(F_test.*Pred_test);
B=sum(sum(Ps_test,1));
Ps1_test(:,1)=B(1,1,:);

for i=1:test_samples
  E(i,1)=abs(exp_test(i,1)-Ps1_test(i,1));      
end

MAE_test=sum(E)/test_samples;

for re1=1:finalregions1
    for re2=1:finalregions2
        for fe=1:features      
            Wpos_t(re1,re2,fe)=abs(bestWrf(re1,re2,fe));
        end
    end
end

REG_test=sum(sum(sum(Wpos_t)));
z_test=MAE_test+l*REG_test;

%%%%%q2 test 
final_test=[exp_test,Ps1_test];
q2 = 1-sum((final_test(:,1)-final_test(:,2)).^2)/sum((final_test(:,1)-mean(exp)).^2); %external%q2
q2_test = 1-sum((final_test(:,1)-final_test(:,2)).^2)/sum((final_test(:,1)-mean(final_test(:,1))).^2); %mean values of the test endpoint

%%%%%q2 test per region
samplestoregions=zeros(samples,1);
samplestoregions_test=zeros(test_samples,1);
for re1=1:finalregions1
    for re2=1:finalregions2
        if sum(F_test(re1,re2,:),3)<=1.5
            q2_re(re1,re2)=NaN;
        else 
            clear final2_test
            clear joint_test
            clear exp_reduced
            clear exp_reduced_test
            samplestoregions(:,1)=bestF(re1,re2,:);
            samplestoregions_test(:,1)=F_test(re1,re2,:);
            exp_reduced = exp(:,1).*samplestoregions(:,1);
            exp_reduced_test=exp_test(:,1).*samplestoregions_test(:,1);
            n = sum(bestF(re1,re2,:),3);
            Ps2_test(:,1)=Ps_test(re1,re2,:);
            joint_test=[exp_reduced_test,Ps2_test(:,1)];
            rowfinal2_test=1;
            for i=1:test_samples
                if F_test(re1,re2,i)==1 
                    final2_test(rowfinal2_test,:)=joint_test(i,:);
                    rowfinal2_test=rowfinal2_test+1;
                end 
            end 
            q2_re(re1,re2)=1-sum((final2_test(:,1)-final2_test(:,2)).^2)/sum((final2_test(:,1)-sum(exp_reduced)/n).^2);
            q2_re_test(re1,re2)=1-sum((final2_test(:,1)-final2_test(:,2)).^2)/sum((final2_test(:,1)-mean(final2_test(:,1))).^2);
        end
    end 
end 

%%%%%R2 test
R_test=corrcoef(final_test(:,2),final_test(:,1));
R2_test=R_test(1,2)^2;

samplestoregions_test=zeros(test_samples,1);
for re1=1:finalregions1
for re2=1:finalregions2
    if sum(F_test(re1,re2,:),3)<=1.5
        R2_test_re(re1,re2)=NaN;
    else 
        clear final2_test
        clear joint_test
        clear exp_reduced_test
        samplestoregions_test(:,1)=F_test(re1,re2,:);
        Ps2_test(:,1)=Ps_test(re1,re2,:);
        exp_reduced_test=exp_test(:,1).*samplestoregions_test(:,1);
        joint_test=[exp_reduced_test,Ps2_test(:,1)];
        rowfinal2_test=1;
        for i=1:test_samples
            if F_test(re1,re2,i)==1 
                final2_test(rowfinal2_test,:)=joint_test(i,:);
                rowfinal2_test=rowfinal2_test+1;
            end 
        end 
        R_test_re=corrcoef(final2_test(:,2),final2_test(:,1));
        R2_test_re(re1,re2)=R_test_re(1,2)^2;
    end  
end
end

%RMSE test
RMSE_test = sqrt(mean((final_test(:,1) - final_test(:,2)).^2));