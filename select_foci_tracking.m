% Select tracking points and the nucleation position of the actin 
% polymerization foci 
function [a_points,S_track,ff_track] = select_foci_tracking(S,ff,P) 
%     crete a new mesh for the nucleation postion
    [~,S_aux] = remeshing(ff, S*P.s_scale, int32([]), P.delta_S*P.new_mesh, int32(3));
%     select the positions that do not correspond to the PSD or neck
    P.h_PSD = P.s_scale*P.h_PSD;
    P.h_neck = P.s_scale*P.h_neck;
    [S_aux,P.index,P.index2] = fix_points(S_aux,ff,P);
    aux_ind = setdiff(P.index,P.index2);
    a_points = S_aux(aux_ind,:);  
%     select the tracking points
    [ff_aux,S_aux] = remeshing(ff, S*P.s_scale, int32([]), P.delta_S*2, int32(3));
    [S_track,ff_track] = track_points(S_aux,S,ff_aux,ff);
end