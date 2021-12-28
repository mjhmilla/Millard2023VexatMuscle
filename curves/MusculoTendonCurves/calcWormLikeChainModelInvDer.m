function output = calcWormLikeChainModelInvDer(fstar,zguess,L,Tkb_div_A,Tkb_div_ASq,derOrder)
%%
%
% @param z : end-to-end length
% @param A : persistence length 
% @param L : the contour length which is the end-to-end length of the chain
% when stretched with infinite force
% @param T : absolute temperature of the WLC in Kelvin
% @param kb: Boltzmann's constant.
% @param z : output length
%
%% 

%assert(length(derOrder)==1);
output=0;

zL = zguess/L;
flag_solveRootsOfCubicPolynomial=1;


if(flag_solveRootsOfCubicPolynomial ==0)
    iter=0;
    iterMax=100;
    tol=1e-10;
    err=1;
    
    while(abs(err) > tol && iter < iterMax )
      f = calcWormLikeChainModelDer(zL*L,L,Tkb_div_A,Tkb_div_ASq,[0,0,0]);
      err = f-fstar;  
      derr = calcWormLikeChainModelDer(zL*L,L,Tkb_div_A,Tkb_div_ASq,[1,0,0])*L;
      dzL = -err/derr;
      zL = zL+dzL;
      if(zL > 1)
        zL = 1-tol;
      end
      if(zL < 0)
        zL = 0;
      end
      
      iter=iter+1;
    end
    
    if(abs(err) >= tol && abs(dzL) >= tol)
      here=1;
    end
    
    assert(abs(err)<= tol || abs(dzL) <= tol);
    
    derCase = derOrder(1)*1 + derOrder(2)*10 + derOrder(2)*100 ; 
    
    if(derCase == 0)
      output=zL*L;
    elseif (derCase == 1)
      df = calcWormLikeChainModelDer(zL*L,L,Tkb_div_A,Tkb_div_ASq,derOrder);
      output = 1/df;
    else
      assert(0,'This derivative case has not been implemented');
    end

end

if(flag_solveRootsOfCubicPolynomial==1)

  m = Tkb_div_A;

  %The WLC is described by:
  %
  % f = (T*kb*(1/(4*(1-z/L)^2)+z/L-0.25))/A
  %
  % Setting
  %  m = kbT/A
  %  s = z/L
  %
  % yields the polynomial
  %
  % m*(0.25/(1-s)^2 + s - 0.25)-f = 0
  % m*(0.25 + s*(1-s)^2 - 0.25*(1-s)^2)-f*(1-s)^2 = 0
  % 
  % Since
  % (1-s)^2 = s^2 -2*s + 1
  %
  % we have
  % m*(0.25 + s*( s^2 -2*s + 1) - 0.25*( s^2 -2*s + 1) )-f*( s^2 -2*s + 1) = 0
  %
  % m s^3 + (-2m-0.25m-f)s^2 + (1.5m+2f)s + (0.25m-0.25-f)
  %
  %
  % For the polynomial
  % ax^3 + bx^2 +cx + d = 0

  a = m;
  b = -(2.25)*m-fstar;
  c = 1.5*m+2*fstar;
  d = -fstar;

  %%
  % I could use the formula to solve a cubic, but Matlab has a function
  % that will (probably) do exactly the same thing. The formulat can be found
  % here:
  %
  %https://math.vanderbilt.edu/schectex/courses/cubic/
  p = -b/(3*a);
  q = p*p*p + (b*c-3*a*d)/(6*a*a);
  r = c/(3*a);
  %
  % x =  (q + ( (q*q + (r-p^2)^3 )^(1/2)) )^(1/3)
  %    + (q - ( (q*q + (r-p^2)^3 )^(1/2)) )^(1/3)
  %    + p
  %
  % tmp0 = r-p*p;
  % tmp1 = sqrt( q*q + tmp0*tmp0*tmp0);
  % 
  % x =  (q + tmp1)^(1/3) ...
  %    + (q - tmp1)^(1/3)...
  %    + p;
  %
  % This gives one root. Then by synthetic division we could reduce
  % the cubic system to a quadratic and use the quadratic formula.
  % I suspect there are some scaling issues with this direct formula
  % so I'm going to stick with the Matlab function for now

  w = roots([a,b,c,d]);

  for i=1:1:length(w)
    if(abs(imag(w(i,1))) < eps)
      s = w(i,1);
    end
  end

  output = s*L;

end




