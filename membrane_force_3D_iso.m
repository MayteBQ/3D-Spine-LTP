function [f,f_P, f_t, f_k] = membrane_force_3D_iso(vv,ff,P)
    
    f_C=zeros(length(vv),3);
    f_A=zeros(length(vv),3);
    f_V=zeros(length(vv),3);
    
%     for the curvature calculation
    for j=1:length(vv)
        aux_A = ff(ff(:,1)==j,:);
        aux_B = ff(ff(:,2)==j,:);
        aux_C = ff(ff(:,3)==j,:);
        aux_all = [aux_A;aux_B;aux_C];
        
        if ~isempty(aux_all)
            platelet = unique(aux_all);
            aux_platelet = setdiff(platelet,j);
            [N_i,~] = size(aux_platelet);
            A = zeros(N_i,1);
            B = zeros(N_i,1);
            C = zeros(N_i,1);

            jj=1;
            if ~isempty(aux_A)
                m = size(aux_A,1)+jj-1;
                A(jj:m) = aux_A(:,1);
                B(jj:m) = aux_A(:,2);
                C(jj:m) = aux_A(:,3); 
                jj = m+1; 
            end
            if ~isempty(aux_B)
                m = size(aux_B,1)+jj-1;
                A(jj:m) = aux_B(:,2);
                B(jj:m) = aux_B(:,3);
                C(jj:m) = aux_B(:,1); 
                jj = m + 1; 
            end
            if ~isempty(aux_C)
                m = size(aux_C,1)+jj-1;
                A(jj:m) = aux_C(:,3);
                B(jj:m) = aux_C(:,1);
                C(jj:m) = aux_C(:,2); 
            end
            D = zeros(N_i,1);
            for jj = 1:N_i
                D(jj) = C(B==C(jj));
            end
        
        
           
            x_i = vv(A,:);
            x_jm = vv(B,:);
            x_j = vv(C,:);
            x_jp = vv(D,:);

            a1 = x_i - x_jm;
            b1 = x_j - x_jm;
            l_i_jm = sqrt(a1(:,1).^2+a1(:,2).^2+a1(:,3).^2);
            l_j_jm = sqrt(b1(:,1).^2+b1(:,2).^2+b1(:,3).^2);
            cos_theta1 = (a1(:,1).*b1(:,1) + a1(:,2).*b1(:,2) + a1(:,3).*b1(:,3))./...
                (l_i_jm.*l_j_jm);
            cot_theta1 = cos_theta1./sqrt(1-cos_theta1.^2);

            a2 = x_i - x_jp;
            b2 = x_j - x_jp;
            l_i_jp = sqrt(a2(:,1).^2+a2(:,2).^2+a2(:,3).^2);
            l_j_jp = sqrt(b2(:,1).^2+b2(:,2).^2+b2(:,3).^2);
            cos_theta2 = (a2(:,1).*b2(:,1) + a2(:,2).*b2(:,2) + a2(:,3).*b2(:,3))./...
                (l_i_jp.*l_j_jp);
            cot_theta2 = cos_theta2./sqrt(1-cos_theta2.^2);

            aux_deltas_xi = sum(1/2*(cot_theta1 + cot_theta2)'.*(x_i'-x_j'),2);

            l_i_j = sqrt((x_i(:,1)-x_j(:,1)).^2+(x_i(:,2)-x_j(:,2)).^2+(x_i(:,3)-x_j(:,3)).^2);
            A_i = sum((1/8)*((cot_theta1+cot_theta2).*(l_i_j.^2)));


            N = cross(x_jm(1,:)-x_i(1,:),x_j(1,:)-x_i(1,:));
            l_N = sqrt(N*N');
            n= N/l_N;
            deltas_xi = aux_deltas_xi/A_i;
            aux_n = (deltas_xi - (deltas_xi'*n')*n');
            H = (1/2)*deltas_xi'*n';

            dcot_theta1_di = (1./((1-cos_theta1.^2).^(3/2))).*...
                    ((x_j-x_jm)./(l_i_jm.*l_j_jm)-cos_theta1.*(x_i-x_jm)./(l_i_jm.^2));
            dcot_theta1_djm = (1./((1-cos_theta1.^2).^(3/2))).*...
                    ((2*x_jm-x_j-x_i)./(l_i_jm.*l_j_jm)+cos_theta1.*...
                    ((x_i-x_jm)./(l_i_jm.^2)+(x_j-x_jm)./(l_j_jm.^2)));
            dcot_theta1_dj = (1./((1-cos_theta1.^2).^(3/2))).*...
                    ((x_i-x_jm)./(l_i_jm.*l_j_jm)-cos_theta1.*(x_j-x_jm)./(l_j_jm.^2));

            dcot_theta2_di = (1./((1-cos_theta2.^2).^(3/2))).*...
                    ((x_j-x_jp)./(l_i_jp.*l_j_jp)-cos_theta2.*(x_i-x_jp)./(l_i_jp.^2));
            dcot_theta2_dj = (1./((1-cos_theta2.^2).^(3/2))).*...
                    ((x_i-x_jp)./(l_i_jp.*l_j_jp)-cos_theta2.*(x_j-x_jp)./(l_j_jp.^2));
            dcot_theta2_djp = (1./((1-cos_theta2.^2).^(3/2))).*...
                    ((2*x_jp-x_j-x_i)./(l_i_jp.*l_j_jp)+cos_theta2.*...
                    ((x_i-x_jp)./(l_i_jp.^2)+(x_j-x_jp)./(l_j_jp.^2)));


            dAi_di = (1/8)*((dcot_theta1_di+dcot_theta2_di).*(l_i_j.^2)...
                    + 2*(x_i-x_j).*(cot_theta1 + cot_theta2));
            dAi_djm = (1/8)*(dcot_theta1_djm.*(l_i_j.^2));
            dAi_dj = (1/8)*((dcot_theta1_dj+dcot_theta2_dj).*(l_i_j.^2)...
                - 2*(x_i-x_j).*(cot_theta1 + cot_theta2));
            dAi_djp = (1/8)*(dcot_theta2_djp.*(l_i_j.^2));


            daux_deltas_xi_diI = zeros(3,3);
            daux_deltas_xi_dII = zeros(3,3);
            daux_deltas_xi_dIII = zeros(3,1);


            for jj=2:N_i

                daux_deltas_xi_diI = daux_deltas_xi_diI + (dcot_theta1_di(jj,:)+dcot_theta2_di(jj,:)).*(x_i(jj,:)'-x_j(jj,:)');
                aux_dII = (cot_theta1(jj) + cot_theta2(jj))*eye(3);
                daux_deltas_xi_dII = daux_deltas_xi_dII + aux_dII;
                aux_dIII = (cot_theta1(jj) + cot_theta2(jj))*(x_i(jj,:)'-x_j(jj,:)');
                daux_deltas_xi_dIII = daux_deltas_xi_dIII + aux_dIII;

                daux_deltas_xi_djmI = dcot_theta1_djm(jj,:).*(x_i(jj,:)'-x_j(jj,:)');
                daux_deltas_xi_djI = (dcot_theta1_dj(jj,:)+dcot_theta2_dj(jj,:)).*(x_i(jj,:)'-x_j(jj,:)');
                daux_deltas_xi_djpI = dcot_theta2_djp(jj,:).*(x_i(jj,:)'-x_j(jj,:)');

                deltas_xi_djm = (1/2)*(daux_deltas_xi_djmI/A_i - dAi_djm(jj,:).*aux_dIII/(A_i^2));
                deltas_xi_dj = (1/2)*(daux_deltas_xi_djI/A_i - aux_dII/A_i ...
                    - dAi_dj(jj,:).*aux_dIII/(A_i^2));
                deltas_xi_djp = (1/2)*(daux_deltas_xi_djpI/A_i - dAi_djp(jj,:).*aux_dIII/(A_i^2));

                dH_xjm = (1/2)*(n*deltas_xi_djm);
                dH_xj = (1/2)*(n*deltas_xi_dj);
                dH_xjp = (1/2)*(n*deltas_xi_djp);

                f_C(B(jj),:) = f_C(B(jj),:) +(2*H*dH_xjm*A_i+ (H^2)*dAi_djm(jj,:));
                f_C(C(jj),:) = f_C(C(jj),:) + (2*H*dH_xj*A_i+ (H^2)*dAi_dj(jj,:));
                f_C(D(jj),:) = f_C(D(jj),:) + (2*H*dH_xjp*A_i+ (H^2)*dAi_djp(jj,:));

            end

            jj=1;
            dAi_di = sum(dAi_di);

            daux_deltas_xi_diI = daux_deltas_xi_diI + (dcot_theta1_di(jj,:)+dcot_theta2_di(jj,:)).*(x_i(jj,:)'-x_j(jj,:)');
            aux_dII = (cot_theta1(jj) + cot_theta2(jj))*eye(3);
            daux_deltas_xi_dII = daux_deltas_xi_dII + aux_dII;
            aux_dIII = (cot_theta1(jj) + cot_theta2(jj))*(x_i(jj,:)'-x_j(jj,:)');
            daux_deltas_xi_dIII = daux_deltas_xi_dIII + aux_dIII;

            daux_deltas_xi_djmI = dcot_theta1_djm(jj,:).*(x_i(jj,:)'-x_j(jj,:)');
            daux_deltas_xi_djI = (dcot_theta1_dj(jj,:)+dcot_theta2_dj(jj,:)).*(x_i(jj,:)'-x_j(jj,:)');
            daux_deltas_xi_djpI = dcot_theta2_djp(jj,:).*(x_i(jj,:)'-x_j(jj,:)');

            deltas_xi_di = (1/2)*(daux_deltas_xi_diI/A_i + daux_deltas_xi_dII/A_i ...
                - dAi_di.*daux_deltas_xi_dIII/(A_i^2));
            deltas_xi_djm = (1/2)*(daux_deltas_xi_djmI/A_i - dAi_djm(jj,:).*aux_dIII/(A_i^2));
            deltas_xi_dj = (1/2)*(daux_deltas_xi_djI/A_i - aux_dII/A_i ...
                - dAi_dj(jj,:).*aux_dIII/(A_i^2));
            deltas_xi_djp = (1/2)*(daux_deltas_xi_djpI/A_i - dAi_djp(jj,:).*aux_dIII/(A_i^2));

            dH_xi = (1/2)*(n*deltas_xi_di + cross(aux_n,x_j(jj,:)-x_jm(jj,:))/l_N);
            dH_xjm = (1/2)*(n*deltas_xi_djm + cross(aux_n,x_i(jj,:)-x_j(jj,:))/l_N);
            dH_xj = (1/2)*(n*deltas_xi_dj + cross(aux_n,x_jm(jj,:)-x_i(jj,:))/l_N);
            dH_xjp = (1/2)*(n*deltas_xi_djp);


            f_C(A(jj),:) = f_C(A(jj),:) + (2*H*A_i*dH_xi+ (H^2)*dAi_di);
            f_C(B(jj),:) = f_C(B(jj),:) +(2*H*dH_xjm*A_i+ (H^2)*dAi_djm(jj,:));
            f_C(C(jj),:) = f_C(C(jj),:) + (2*H*dH_xj*A_i+ (H^2)*dAi_dj(jj,:));
            f_C(D(jj),:) = f_C(D(jj),:) + (2*H*dH_xjp*A_i+ (H^2)*dAi_djp(jj,:));
        end
    end

    x_i = vv(ff(:,1),:);
    x_jm = vv(ff(:,2),:);
    x_j = vv(ff(:,3),:);

    N = cross(x_jm-x_i,x_j-x_i);
    l_N = sqrt(N(:,1).^2+N(:,2).^2+N(:,3).^2);
    n = N./l_N;%n = N;%
% for the area

    f_A(ff(:,1),:) = cross(n,x_j-x_jm);
    f_A(ff(:,2),:) = f_A(ff(:,2),:)+cross(n,x_i-x_j);
    f_A(ff(:,3),:) = f_A(ff(:,3),:)+cross(n,x_jm-x_i);
% for the volume
    f_V(ff(:,1),:) = cross(x_jm,x_j);
    f_V(ff(:,2),:) = f_V(ff(:,2),:)+cross(x_j,x_i);
    f_V(ff(:,3),:) = f_V(ff(:,3),:)+cross(x_i,x_jm);
        
      
    f =-(2*P.kappa)*f_C -(P.P/6)*f_V-(P.tau/2)*f_A;%
    f_P = -(P.P/6)*f_V;%force generated by membrane pressure
    f_t = -(P.tau/2)*f_A;%force generated by membrane tension
    f_k = -(2*P.kappa)*f_C;%force generated by membrane curvature
    
end