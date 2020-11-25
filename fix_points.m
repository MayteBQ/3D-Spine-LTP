function [S_new,index,index2] = fix_points(S,ff,P)
    S_new = S;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux = find(S(:,3)>=P.h_PSD & sqrt(S(:,1).^2 +S(:,2).^2)<=P.r_PSD);
    S_new(aux,3)=P.h_PSD;
    
    aux2 = find(S(:,3)<=P.h_neck &  sqrt(S(:,1).^2 +S(:,2).^2)<=P.r_neck);
    S_new(aux2,3)=P.h_neck;
    
%     find the closest neighboring vertices to those of neck
%     because the membrane continues to the neck instead of closing, these
%     points are also fixed, but their y-coordinate is unchanged. 
    aux3 = [];
    for j=1:length(aux2)
        aux_A = ff(ff(:,1)==aux2(j),:);
        aux_B = ff(ff(:,2)==aux2(j),:);
        aux_C = ff(ff(:,3)==aux2(j),:);
        aux3 = [aux3;aux_A;aux_B;aux_C];
    end
    aux3 = unique(aux3);
    index2 =  [aux;aux2;aux3];    
    index2 = unique(index2);
    index = setdiff(1:length(S),index2);

end