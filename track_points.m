function [S_track,ff_track] = track_points(S_new,S_old,ff_new,ff_old)
% update the position of the tracking points in terms of the new mesh    
    if size(S_new,1) >= size(S_old,1)
        S_track = zeros(size(S_old));
        for k=1:size(S_old,1)
            dist = sqrt((S_new(:,1)-S_old(k,1)).^2+...
                (S_new(:,2)-S_old(k,2)).^2+...
                (S_new(:,3)-S_old(k,3)).^2);
            ind = find(dist == min(dist),1);
            S_track(k,:) = S_new(ind,:);
            ff_track = ff_old;
        end
    else
        S_track = S_new;
        ff_track = ff_new;
    end
end
        
        
        