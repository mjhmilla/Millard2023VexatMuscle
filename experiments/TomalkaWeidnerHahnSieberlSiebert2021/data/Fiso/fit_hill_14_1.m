%---------------------------------------------------------------
% fit hill 
% fv=(1-v_norm)./(1+v_norm/curv);   curv=a/Fiso (wie in Lit)
%---------------------------------------------------------------

% initial parameters
x0 = [1 0.1]; % bei normierten vmax Startwert [1, curv] _ curv=a/Fiso (wie in Lit), abhängig von Experimentaltemp.& Muskel

xdata=[ velocity2 velocity7 velocity5 velocity3 velocity1 velocity4 velocity6 velocity8 -.9]*-1; %v [mm/s]_EDL_02.09.15_isotonic_force-step_fiber length - lopt [mm]:1.55
ydata=[ relF2 relF7 relF5 relF3 relF1 relF4 relF6 relF8 .12]; %F/Fmax;  
fiso=0.33;

% set optim options
options = optimset('maxfunevals', 100)
p = lsqcurvefit(@hill_funktion,x0,xdata,ydata); 
force = hill_funktion(p, xdata);

%Hill parameter a,b
%a=p(2)*fiso;  %a=curv*Fiso; !!!!!!!!!!! Kontrolle
a=p(2)*1;  %wenn F auf Fiso normiert
%a=p(2)*ydata(1);  %a=curv*Fiso
b=p(2)*p(1);  %b=curv*vmax
c=b*(fiso+a)
pmax=a*b+c-2*sqrt(a*b*c);

%v_norm=0:0.1:1;
v_norm=-0.2:0.1:6;
%fv=(p(1)-v_norm)./(p(1)+v_norm/p(2));
fv= hill_funktion(p, v_norm);

%Parameter
fiso
vmax=p(1)
curv=p(2)
a
b
c=b*(fiso+a)
pmax
%%
figure(103)
hold on
plot(xdata, force, 'ro-')   % Fit
plot(xdata, ydata, 'k.-') %Messdaten
plot(v_norm,fv,'g*') % v in mm/s

xlabel('velocity_a_b_s [Lopt/s]')
ylabel('Force [F/F_I_M]')

xx=[0 -a;10*p(1) -a]; % -a [N] = Asymptote zur x-Achse im Unendlichen 
plot(xx(:,1),xx(:,2),'b-')
yy=[-b 0;-b 20*fiso]; % -b [mm/s] = Asymptote zur y-Achse im Unendlichen 
plot(yy(:,1),yy(:,2),'m-')
% Koordinatenachsen
x=[0 0;10*p(1) 0];
plot(x(:,1),x(:,2),'k-')
y=[0 0;0 20*fiso];
plot(y(:,1),y(:,2),'k-')


legend('fitted data','experimental data','Hill-fit')

xlim ([0 3])
ylim ([-0.1 1.01])

grid on
hold off