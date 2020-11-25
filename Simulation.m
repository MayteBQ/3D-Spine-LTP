close all;
clear all;

load('resting_shape.mat')
% load resting shape mesh, S corresponds to the vertices positions, ff to
% the triangulation and P to the paramenter structure
% Here P.tau = sigma : surface tension
P.s_scale = 0.99;

% % select tracking points for more than one actin polymerization focus
% P.new_mesh = 8;%for 22 foci
% [P.a_points,S_track,ff_track] =select_foci_tracking(S,ff,P);

% select tracking points for a single actin polymerization focus
[P.a_points,S_track,ff_track] =select_foci_tracking_dist(S,ff,P,1262);%for the focus near to the PSD

% model parameters
P.n_fil = 70;%number of actin filaments in the spine head
P.phi = P.n_fil/(size(S,1)*size(P.a_points,1));
P.alpha = P.phi*3.8;
P.K = size(S,1);
P.zeta = 0.004;
P.d_max_int = 5e-4;
P.t_initial = 0; %s
P.t_end =9*60;% s 

aux_t = P.t_initial:P.delta_t:P.t_end;
aux_S = cell(length(aux_t),1);
aux_S_tracking = cell(length(aux_t),1);
aux_ff = cell(length(aux_t),1);%triangulation
aux_ff_tracking = cell(length(aux_t),1);
save_volume = zeros(length(aux_t),1);%save volume
aux_f_P = cell(length(aux_t),1);
aux_f_t = cell(length(aux_t),1);
aux_f_k = cell(length(aux_t),1);
save_area = zeros(length(aux_t),1);
% save initial shape
aux_S{1} = S;
aux_S_tracking{1} = S_track;
aux_ff{1} = ff;
aux_ff_tracking{1} = ff_track;
save_volume(1) =  volume_sphere(ff,S);
save_area(1) = surface_area(S,ff);
[F_mem, f_P, f_t, f_k] = membrane_force_3D_iso(S_track,ff_track,P);
aux_f_P{1} = sqrt(f_P(:,1).^2 +f_P(:,2).^2 +f_P(:,3).^2);
aux_f_t{1} = sqrt(f_t(:,1).^2 +f_t(:,2).^2 +f_t(:,3).^2);
aux_f_k{1} = sqrt(f_k(:,1).^2 +f_k(:,2).^2 +f_k(:,3).^2);


max_tension = 2.5*max(aux_f_t{1});
max_volume = 2.5*save_volume(1);
j=2; 
while j<=length(aux_t) && max(aux_f_t{j-1})< max_tension && save_volume(j-1)<max_volume
    [S, ff, P] = solve_system_threshold_3D_rk_iso(S,ff,1,P);
    P.K = size(S,1);
    aux_S{j} = S;
    aux_ff{j} = ff;
    save_volume(j) =  volume_sphere(ff,S);
    save_area(j) = surface_area(S,ff);
    [aux_S_tracking{j},aux_ff_tracking{j}] = track_points(S,aux_S_tracking{j-1},ff,aux_ff_tracking{j-1});
    
    [F_mem, f_P, f_t, f_k] = membrane_force_3D_iso(aux_S_tracking{j},aux_ff_tracking{j},P);
    aux_f_P{j} = sqrt(f_P(:,1).^2 +f_P(:,2).^2 +f_P(:,3).^2);
    aux_f_t{j} = sqrt(f_t(:,1).^2 +f_t(:,2).^2 +f_t(:,3).^2);
    aux_f_k{j} = sqrt(f_k(:,1).^2 +f_k(:,2).^2 +f_k(:,3).^2); 
    
    j = j+1;
end


figure(1)
for jj = 1:1/P.delta_t:j-1
    figure(1)
    clf
    lll = aux_S{jj}; 
    ff_ll = aux_ff{jj};
    aux_color = aux_f_t{jj};
    lll2 = aux_S_tracking{jj}; 
    ff_ll2 = aux_ff_tracking{jj};
    trimesh( ff_ll, lll(:,1), lll(:,2), lll(:,3),'edgecolor','k')
    hold on
    scatter3(lll2(:,1),lll2(:,2),lll2(:,3),15,aux_color,'filled')
    title(['time = ' num2str(aux_t(jj),'%.1f') ' s  '])
    set(gca,'fontsize',20);
    xlabel('x (\mu m)')
    ylabel('y (\mu m)')
    zlabel('z (\mu m)')   
    colormap(hot)
    caxis([0 max_tension])
    axis([-.8 .8 -.8 .8 -.8 .8])
    colorbar
end

figure;
plot(aux_t(1:j-1),save_volume(1:j-1),'linewidth',2)
set(gca,'fontsize',20)
xlabel('t (sec)')
ylabel('Volume (\mu{m}^3)')


% make the omega profile
jj = j-1;
S_track= aux_S_tracking{jj};
[S,ff,P] = Omega_profile(j,aux_f_t,max_tension,S,ff,S_track,P);

figure
trimesh(ff,S(:,1),S(:,2),S(:,3),'edgecolor','b')
alpha 0.5

aux_S{j} = S;
aux_ff{j} = ff;
save_volume(j) =  volume_sphere(ff,S);
save_area(j) = surface_area(S,ff);
j2 = j;
j = j+1;

% evolution of the fusion of the Omega-profile 
while j<=j2+30/P.delta_t
    [S, ff, P] = solve_system_threshold_3D_rk_iso(S,ff,1,P);
    P.K = size(S,1);
    aux_S{j} = S;
    aux_ff{j} = ff;
    save_volume(j) =  volume_sphere(ff,S);
    save_area(j) = surface_area(S,ff);
    j = j+1;
end
figure;
plot(aux_t(j2:j-1)-aux_t(j2),save_volume(j2:j-1),'linewidth',2)
set(gca,'fontsize',20)
xlabel('t (sec)')
ylabel('Volume (\mu{m}^3)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%changing the size of the PSD without actin force
P.delta_r_PSD = P.delta_S/75;
j3 = j-1;
aux = find(S(P.index2,3)>0);
aux_area_PSD = zeros(size(aux_t));
aux_area_PSD(1:j-1)= surface_area_PSD(S,ff,P.index2(aux));
P.index_PSD = P.index2(aux);

while j<= j3+7*60/P.delta_t && aux_area_PSD(j-1) <=aux_area_PSD(j3)*2
    [S,P.index_PSD,P] = increase_PSD(S,ff,j-j3,P);
    [S, ff, P] = solve_system_threshold_3D_rk_iso_exo(S,ff,0,P);
    P.K = size(S,1);
    aux_S{j} = S;
    aux_area_PSD(j) = surface_area_PSD(S,ff,P.index_PSD);
    aux_ff{j} = ff;
    save_volume(j) =  volume_sphere(ff,S);
    save_area(j) = surface_area(S,ff);
    
    j = j+1;
end

figure;
plot(aux_t(j3:j-1)-aux_t(j3),aux_area_PSD(j3:j-1),'linewidth',2)
set(gca,'fontsize',20)
xlabel('t (sec)')
ylabel('Volume (\mu{m}^3)')

for jj = 1:10/P.delta_t:j-1
    figure(1)
    clf
    lll = aux_S{jj}; 
    ff_ll = aux_ff{jj};
    trimesh( ff_ll, lll(:,1), lll(:,2), lll(:,3),'edgecolor','k')
    title(['time = ' num2str(aux_t(jj),'%.1f') ' s  '])
    set(gca,'fontsize',20);
    xlabel('x (\mu m)')
    ylabel('y (\mu m)')
    zlabel('z (\mu m)')   
    axis([-.8 .8 -.8 .8 -.8 .8])
end