function z_test=OPLRA_1D_test(fbest,features,test_samples,exp_test,finalregions,Asf_test,l,bestXrf,bestWrf,bestBr)

Frs_test=zeros(finalregions,test_samples);

for i=1:test_samples
    for re=2:finalregions-1 
        if Asf_test(i,fbest)>=bestXrf(re-1,fbest) && Asf_test(i,fbest)<=bestXrf(re,fbest)
            Frs_test(re,i)=1;
        end
    end
end

for i=1:test_samples 
    re=1 ;
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
