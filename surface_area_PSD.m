function f = surface_area_PSD(S,ff,index_PSD)
% calculate the surface area corresponding to the PSD
    aux = [];
    for j = 1:length(index_PSD)
        aux_A = ff(ff(:,1)==index_PSD(j),:);
        aux_B = ff(ff(:,2)==index_PSD(j),:);
        aux_C = ff(ff(:,3)==index_PSD(j),:);
        aux = [aux;aux_A;aux_B;aux_C];
    end
    aux = unique(aux);

    x_i = S(ff(aux,1),:);
    x_jm = S(ff(aux,2),:);
    x_j = S(ff(aux,3),:);

    N = cross(x_jm-x_i,x_j-x_i);
    f = sum(sqrt(N(:,1).^2+N(:,2).^2+N(:,3).^2))/2;
end