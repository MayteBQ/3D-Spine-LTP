function [S,ff,P] = Omega_profile(j,aux_f_t,max_tension,S,ff,S_track,P)
% create the Omega-profile in the spines
%     obtin the index of the mesh vertex with the highest tension
    aux_ind_tracking = find(aux_f_t{j-1}>= max_tension);
    [~,aux_ind] = min(sqrt((S_track(aux_ind_tracking,1)-S(:,1)).^2 + (S_track(aux_ind_tracking,2)-S(:,2)).^2 + ...
       (S_track(aux_ind_tracking,3)-S(:,3)).^2));

%     from that vertex obtain the neighbouring vertices

    aux_A = ff(ff(:,1)==aux_ind,:);
    aux_B = ff(ff(:,2)==aux_ind,:);
    aux_C = ff(ff(:,3)==aux_ind,:);
    aux_all = unique([aux_A;aux_B;aux_C]);

    aux_all2=[];
    for k =1:length(aux_all)
        aux_A = ff(ff(:,1)==aux_all(k),:);
        aux_B = ff(ff(:,2)==aux_all(k),:);
        aux_C = ff(ff(:,3)==aux_all(k),:);
        aux_all2 = [aux_all2;aux_A;aux_B;aux_C];
    end
    aux_all2 = unique(aux_all2);

    aux_all3=[];
    for k =1:length(aux_all2)
        aux_A = ff(ff(:,1)==aux_all2(k),:);
        aux_B = ff(ff(:,2)==aux_all2(k),:);
        aux_C = ff(ff(:,3)==aux_all2(k),:);
        aux_all3 = [aux_all3;aux_A;aux_B;aux_C];
    end
    aux_all3 = unique(aux_all3);

    aux_all4=[];
    for k =1:length(aux_all3)
        aux_A = ff(ff(:,1)==aux_all3(k),:);
        aux_B = ff(ff(:,2)==aux_all3(k),:);
        aux_C = ff(ff(:,3)==aux_all3(k),:);
        aux_all4 = [aux_all4;aux_A;aux_B;aux_C];
    end
    aux_all4 = unique(aux_all4);
%    allocate different weights according to the distince to the maximum tension vertex 
    aux_loc = zeros(size(S,1),1);
    aux_loc(aux_ind) = 1;
    aux_loc(setdiff(aux_all,aux_ind)) = 0.7;
    aux_loc(setdiff(aux_all3,[aux_all;aux_ind])) = 0.5;
    aux_loc(setdiff(aux_all4,[aux_all3;aux_all;aux_ind])) = 0.4;
%     push inwards the membrane as a gaussian around the vertex with the
%     highest tension
    sig = .1;
    alp = .003;
    x_diff = S(:,1)-S(:,1)';
    y_diff = S(:,2)-S(:,2)';
    z_diff = S(:,3)-S(:,3)';
    d = sqrt(x_diff.^2 + y_diff.^2 + z_diff.^2);
    W = (alp/(sig*sqrt(2*pi)))*exp(-(d.^2)./(2*sig^2));
    B =W*aux_loc;
    aux_dist = sqrt((-S(:,1)).^2 + (-S(:,2)).^2 + (-S(:,3)).^2);%);%
    aux_dir = [-S(:,1) -S(:,2) -S(:,3)];
    f = B.*(aux_dir./aux_dist);
    S_new = S + f;

%     move half-way the connecting points to the rest of the membrane
    aux_all5=[];
    for k =1:length(aux_all4)
        aux_A = ff(ff(:,1)==aux_all4(k),:);
        aux_B = ff(ff(:,2)==aux_all4(k),:);
        aux_C = ff(ff(:,3)==aux_all4(k),:);
        aux_all5 = [aux_all5;aux_A;aux_B;aux_C];
    end
    aux_all5 = unique(aux_all5);
    aux_change = setdiff(1:size(S,1),[aux_all5;aux_all4;aux_all3;aux_all2;aux_all]);
    for k = 1:length(aux_change)
        S_new(aux_change(k),:) = S(aux_change(k),:);
    end
    aux_change = setdiff(aux_all5,[aux_all4;aux_all3;aux_all2;aux_all]);
    for k = 1:length(aux_change)
        S_new(aux_change(k),:) = 0.5*(S(aux_ind,:) +S_new(aux_change(k),:));
    end

%     remesh
    [ff_new,S_new] = remeshing(ff, S_new, int32(P.index2), P.delta_S, int32(3));
    [S_new,P.index,P.index2] =  fix_points(S_new,ff_new,P);



S = S_new;
ff = ff_new;
% j = 100;
% aux_S{j} = S;
% aux_ff{j} = ff;
% aux_S_tracking{j}= S_track2;
% aux_ff_tracking{j} = ff_track2;
% aux_S_track_RE{j} = S_track_RE;
% save_volume(j) =  volume_sphere(ff,S);
% save_area(j) = surface_area(S,ff);
% [F_mem, f_P, f_t, f_k] = membrane_force_3D_iso(aux_S_tracking{j},aux_ff_tracking{j},P);
% aux_f_P{j} = sqrt(f_P(:,1).^2 +f_P(:,2).^2 +f_P(:,3).^2);
% aux_f_t{j} = sqrt(f_t(:,1).^2 +f_t(:,2).^2 +f_t(:,3).^2);
% aux_f_k{j} = sqrt(f_k(:,1).^2 +f_k(:,2).^2 +f_k(:,3).^2);
% j2 = j;
% j = j+1;
% 
% % P.t_end = 20*60;
% % aux_t= P.t_initial:P.delta_t:P.t_end;
% P.a_points = [P.a_points;0 -.43 .2];
end