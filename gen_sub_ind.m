function ind = gen_sub_ind(sizes,ratio)
% Generate uniform bool sub-indexes of a matrix with given density ratio
ind=zeros(sizes,'logical');
row=zeros(sizes(1),1);
for i=0:(sizes(1)-1)/ratio
row(i*ratio+1)=true;
end
for i=0:(sizes(2)-1)/ratio
    ind(:,i*ratio+1)=row;
end
end

