function f = calcWormLikeChainModelDer(z,L,Tkb_div_A,Tkb_div_ASq,derOrder)
%%
%
% @param z : end-to-end length
% @param A : persistence length 
% @param L : the contour length which is the end-to-end length of the chain
% when stretched with infinite force
% @param T : absolute temperature of the WLC in Kelvin
% @param kb: Boltzmann's constant.
% @param f : output
%
%% 

%assert(length(derOrder)==3);

%T  = 1;%273.2 + 38.0; %Body temperature
%kb = 1;%1.380649e-23;%Boltzman's constant in J/K

%                  z               A                L
derCase = derOrder(1)*1 + derOrder(2)*10 + derOrder(2)*100 ; 


f = 0;
zL = z/L;
switch(derCase)
  case 0    
    %F
    t0 =(1-zL);
    t1 = t0*t0;
    f  = (Tkb_div_A*(1/(4*t1)+zL-0.25));    
  case 1
    % dF/dz
    t0 =(1-zL);
    t1 = t0*t0; % (1-z/L)^2
    t2 = t1*t0; % (1-z/L)^3
    f  = (Tkb_div_A*(1/(2*L*t2)+1/L));
  case 10
    %dF/dA
    t0 =(1-zL);
    t1 = t0*t0;    
    f = -(Tkb_div_ASq*(1/(4*t1)+t0-0.25));
  case 100
   %dF/dL
    t0 =(1-zL);
    t1 = t0*t0; % (1-z/L)^2    
    t2 = t1*t0; % (1-z/L)^3
    s0 = L*L;
    f = (Tkb_div_A*(-z/(2*s0*t2)-z/s0));   
  otherwise 
    assert(0,'Higher derivatives not yet implemented');
end
