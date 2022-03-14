function fiberKinematics =...
        calcFixedWidthPennatedFiberKinematics(fiberLengthAlongTendon,...
                                        fiberVelocityAlongTendon,...
                                        optimalFiberLength,...
                                        pennationAngleAtOptimalFiberLength)
%%
%This function uses the fixed width pennation model to transforms the 
%length and velocity of the fiber along the tendon into a fiber length, 
%pennation angle, fiber velocity, and pennation angular velocity. 
%
%This model very simply assumes that the fibers remain parallel that their  
%length and pennation angle vary such that the width and height of the
%parallelogram remain constant. This is a simple approximation to the
%constant volume property of muscle: the area of this parallelogram will
%be constant, and so if you give this a constant depth the volume of the
%parallelpiped remains constant.
%
%                       
%                          /----------/=====  -
%                         /          /        | h
%                   =====/----------/         -
%                       |<--  w -- >|
% 
%  = : tendon
%  / : fiber
%
% Note that this the terms 'height' and 'thickness' are used interchangably
% with 'width' when referring to this model. 
%
% @param fiberLengthAlongTendon  (m)
% @param fiberVelocityAlongTendon (m/s)
% @param optimalFiberLength (m)
% @param pennationAngleAtOptimalFiberLength (radians)
%
% @returns fiberKinematics, a structure containing the fields
% 
%     .fiberLength    (m)
%     .fiberVelocity   (m/s)
%     .pennationAngle  (radians)
%     .pennationAngularVelocity  (radians/sec)
%
%%

lceAT     = fiberLengthAlongTendon;    % lce*cos(aPen)
dlceAT     = fiberVelocityAlongTendon;  % (dlce/dt)*cos(aPen) - lce*sin(aPen)*(daPen/dt)


lce    = [];
dlce   = [];
alpha  = [];
dalpha = [];


if(pennationAngleAtOptimalFiberLength > eps^0.5)
    lopt      = optimalFiberLength;         
    alphaOpt  = pennationAngleAtOptimalFiberLength;

    %%
    %Length information
    %%
    
    %%
    %1b. Insert the correct expression for the fiber pennation angle alpha
    %    and the fiber length lce
    %%    
    h     = lopt*sin(alphaOpt); %the height/thickness of the pennated fiber, which is constant
    alpha = atan2(h,lceAT); %alphaOpt
    lce   = sqrt( h*h + lceAT*lceAT ); % lceOpt

    %%
    %Velocity information: obtained by solving:
    % [sin(alpha) ,  lce*cos(alpha)] (dlce/dt  ) = 0        [1]
    % [cos(alpha) , -lce*sin(alpha)] (dalpha/dt) = vceAT    [2]
    %           matrix                  vector     vector        
    %
    % Eqn 1 is just dh/dt = d/dt (lce*sin(alpha)) = 0
    % Eqn 2 is just d/dt (lce*cos(alph)) = fiberVelocityAlongTendon
    %%
    %%
    %1c. Solve the system of equations described above for
    %    the fiber velocity (dlce) and the fiber velocity along the tendon
    %    dlceAT.
    %%
    A = [sin(alpha), lceAT;... % [1 0; 0 1];%
        cos(alpha), -h   ]; 

    b = [0; ...   %[0;0];
         dlceAT];

    if(isnan(dlceAT) == 0)
        if( max(max(isnan(A))) > 0 || max(max(isinf(A))) >0 )
           here=1; 
        end
        x = pinv(A)*b;

        dlce   = x(1);
        dalpha = x(2);
    else
        dlce = NaN;
        dalpha=NaN;
    end
    
    
    
else
    lce    = lceAT;
    dlce   = dlceAT;
    alpha  =0;
    dalpha = 0;    
end

fiberKinematics.fiberLength    = lce;
fiberKinematics.fiberVelocity  = dlce;
fiberKinematics.pennationAngle = alpha;
fiberKinematics.pennationAngularVelocity = dalpha;



