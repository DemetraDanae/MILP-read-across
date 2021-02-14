function [z_test, E_test]=OPLRA_2D_test_final(fbest1,fbest2,test_samples,exp_test,finalregions1,finalregions2,Asf_test,Asf1_test,Asf2_test,l,bestXrf1,bestXrf2,bestWrf,bestBr1r2)

features=size(Asf_test,2);
    
F_test=zeros(finalregions1,finalregions2,test_samples);

for i=1:test_samples
    for re1=2:finalregions1-1 
        for re2=2:finalregions2-1
         if Asf1_test(i,fbest1)>=bestXrf1(re1-1,fbest1) && Asf1_test(i,fbest1)<=bestXrf1(re1,fbest1) && Asf2_test(i,fbest2)>=bestXrf2(re2-1,fbest2) && Asf2_test(i,fbest2)<=bestXrf2(re2,fbest2)
            F_test(re1,re2,i)=1;
        end
    end
    end
end 


for i=1:test_samples 
    re1=1;
    re2=1;
    if bestXrf1(re1,fbest1 )>Asf1_test(i,fbest1) && bestXrf2(re2,fbest2 )>Asf2_test(i,fbest2)
        F_test(re1,re2,i)=1;
    end
end


for i=1:test_samples
    re1=finalregions1;
    re2=finalregions2;
    if bestXrf1(re1-1,fbest1 )<Asf1_test(i,fbest1) && Asf1_test(i,fbest1)<=bestXrf1(re1,fbest1) && bestXrf2(re2-1,fbest2 )<Asf2_test(i,fbest2) && Asf2_test(i,fbest2)<=bestXrf2(re2,fbest2)
        F_test(re1,re2,i)=1;
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

Ps_test=value(F_test.*Pred_test);
B=sum(sum(Ps_test,1));
Ps1_test(:,1)=B(1,1,:);

for i=1:test_samples
  E_test(i,1)=abs(exp_test(i,1)-Ps1_test(i,1));      
end

MAE_test=sum(E_test)/test_samples;

for re1=1:finalregions1
    for re2=1:finalregions2
        for fe=1:features       
            Wpos_t(re1,re2,fe)=abs(bestWrf(re1,re2,fe));
        end
    end
end


REG_test=sum(sum(sum(Wpos_t)));
z_test=MAE_test+l*REG_test;
