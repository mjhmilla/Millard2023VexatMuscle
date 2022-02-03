function dy = stateDerivativeDwellDependency(t, y, pathFcn,...
    omega, zeta,timeConstant)

if(t > 0.1)
    here=1;
end

path = pathFcn(t);
l  = path(1);
dl = path(2);

x  = y(1);
dx = y(2);

f = (dl-dx)/timeConstant;

ddx = f - (2*zeta*omega)*dx -(omega*omega)*x;

dy = [dx;ddx];

%t0 = dl/velocityStick;
%s  = exp(-t0*t0);
%dy = dl*(1-s) - s*(lsrs/timeConstant);



