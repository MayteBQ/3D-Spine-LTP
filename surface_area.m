function f = surface_area(S,ff)
    x_i = S(ff(:,1),:);
    x_jm = S(ff(:,2),:);
    x_j = S(ff(:,3),:);

    N = cross(x_jm-x_i,x_j-x_i);
    f = sum(sqrt(N(:,1).^2+N(:,2).^2+N(:,3).^2))/2;
end