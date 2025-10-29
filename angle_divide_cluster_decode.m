function angle_sequence = angle_divide_cluster_decode(angle_series,nbin)
for time_points = 1:size(angle_series,2)
angle_sequence(:,time_points)=angle_divide_cluster(angle_series(:,time_points),nbin);
end

function angle_sequence=angle_divide_cluster(angle_topo,nbin)

angle_range=0:(2*pi/nbin):2*pi;
angle_sequence=zeros(size(angle_topo));
for j = 1:length(angle_topo)
    for i = 1:nbin
        if angle_range(i)<=angle_topo(j) && angle_range(i+1)>=angle_topo(j)
        angle_sequence(j)=i;
        end
    end
end
end

end