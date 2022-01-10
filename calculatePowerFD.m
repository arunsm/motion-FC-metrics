% Function to calculate framewise displacement (FD) using the Power
% formulation.
% Input: 6-parameter realignment parameters as (t x 6) matrix, where t is
% the number of time points; rotational parameters expected to be in
% radians
% Output: Framewise displacement (FD) at each time point

function FD = calculatePowerFD(RP)
    nTimePoints = size(RP, 1);
    FD = zeros(1, nTimePoints);
    for i = 1:nTimePoints
        if i == 1
            FD(i) = 0;
        else
            FD(i) = sum(abs(RP(i-1, 1:3) - RP(i, 1:3))) + sum(50*(abs(RP(i-1, 4:6) - RP(i, 4:6))));
        end
    end


%%
% nTimePoints = size(RP, 1);
% rmax = 80;
% affineParameters = zeros(nTimePoints, 12);
% 
% % create matrix containing affine transform parameters
% for i = 1:nTimePoints
%     R = calculateRotationalAffine(RP(i, 4), RP(i, 5), RP(i, 6));
%     affineParameters(i, 1:9) = reshape(R, 1, 9);
%     affineParameters(i, 10:12) = RP(i, 1:3);
% end
% 
% pm = zeros(nTimePoints, 16);
% pm(:, 1:12) = affineParameters;
% pm(:, 13:16) = repmat([0, 0, 0, 1], nTimePoints, 1);
% T_rb_prev = eye(4);
% relativeRMSmotion = zeros(1, nTimePoints);
% 
% for i = 1:nTimePoints
%     T_rb = reshape(pm(i, :), [4, 4]); % reshaping 12 motion parameters for current time point into 4x4 matrix
%     if i == 1
%         relativeRMSmotion(i) = 0;
%     else
%         M = dot(T_rb, inv(T_rb_prev)) - eye(4);
%         A = M(1:4, 1:4);
%         b = M(1:4, 4);
%         relativeRMSmotion(i) = sqrt((rmax*rmax/5) * trace(dot(A', A)) + dot(b', b));
%     end
%     T_rb_prev = T_rb;
% end
% meanRelativeRMSmotion = mean(relativeRMSmotion);

end