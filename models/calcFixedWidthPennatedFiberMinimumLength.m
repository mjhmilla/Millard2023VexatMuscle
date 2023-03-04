function minFiberKinematics =...
    calcFixedWidthPennatedFiberMinimumLength(...
                         minimumActiveFiberNormalizedLength,...
                         maximumPennationAngle,...
                         optimalFiberLength,...
                         pennationAngleAtOptimalFiberLength)
%%
% This function calculates the minimum length the fiber can contract to
% without exceeding a given maximum pennation angle, nor falling below the
% minimum physical limit.
%
% @param minimumActiveFiberNormalizedLength approximately 0.5 * lceOpt (m)
% @param maximumPennationAngle usually slightly less than pi/2 (radians)
% @param optimalFiberLength (m)
% @param pennationAngleAtOptimalFiberLength (radians)
%
% @return minFiberKinematics, a structure with the fields
%
%         .minimumFiberLength            
%         .minimumFiberLengthAlongTendon 
%         .pennationAngleAtMinimumFiberLength          

%
%%
lceOpt      = optimalFiberLength;         
lceActMin   = lceOpt*minimumActiveFiberNormalizedLength;


minFiberLength             = [];
minFiberLengthAlongTendon  = [];
pennationAngleAtLceMin     = [];

epsRoot = eps^0.5;
if(pennationAngleAtOptimalFiberLength > epsRoot) 
    alphaOpt    = pennationAngleAtOptimalFiberLength;
    h           = lceOpt*sin(alphaOpt); %the height/thickness of t
                                        %he pennated fiber, which is constant    
    lcePenMin                 = h/sin(maximumPennationAngle);

    minFiberLength            = max([lcePenMin,lceActMin]);

    assert( minFiberLength > epsRoot, ...
            ['Minimum fiber length is too close to 0!',...
             'A fiber length of 0 will cause singularities',...
             ' in the pennation model']);
   
    pennationAngleAtLceMin    = asin( h/minFiberLength );
    minFiberLengthAlongTendon = minFiberLength*cos(pennationAngleAtLceMin); 

else
    minFiberLength            = lceActMin;
    minFiberLengthAlongTendon = lceActMin;
    pennationAngleAtLceMin    = 0;
end



assert( pennationAngleAtLceMin < pi/2 - epsRoot, ...
    ['Maximum pennation angle too close to pi/2!',...
     'If reached this will cause singularities',...
     ' in the pennation model']);

assert( minFiberLength > epsRoot, ...
    ['Minimum fiber length is too close to 0!',...
     'A fiber length of 0 will cause singularities',...
     ' in the pennation model']);

 
 
minFiberKinematics.minimumFiberLength       = minFiberLength;
minFiberKinematics.minimumFiberLengthAlongTendon ...
                                            = minFiberLengthAlongTendon;
minFiberKinematics.pennationAngleAtMinimumFiberLength ...
                                            =  pennationAngleAtLceMin;

