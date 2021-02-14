function scaled=scaling(prop)
samples=size(prop,1);
features=size(prop,2);

prop = table2array(prop);
for j=1:features
   for i=1:samples
      prop_sc(i,j) = (prop(i,j)-min(prop(:,j)))/(max(prop(:,j))-min(prop(:,j)));
   end
end

scaled = prop_sc(:,all(~isnan(prop_sc)));
end 
