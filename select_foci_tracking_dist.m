% tracking points for a single actin polymerization focus
function [a_points,S_track,ff_track] = select_foci_tracking_dist(S,ff,P,ind)
% ind: index of the vertex related to the nucleation position of the actin
% polymerization focus
% selecting the starting position of the focus
    a_points = P.s_scale*S(ind,:);
% select the tracking points
    [ff_aux,S_aux] = remeshing(ff, S*P.s_scale, int32([]), P.delta_S*2, int32(3));
    [S_track,ff_track] = track_points(S_aux,S,ff_aux,ff);
end