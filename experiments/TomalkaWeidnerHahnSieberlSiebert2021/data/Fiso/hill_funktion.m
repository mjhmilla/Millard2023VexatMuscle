function force = fun(p,xdata)

% fv=(1-v_norm)./(1+v_norm/curv);  curv=a/Fiso (wie in Lit)
% p(1) = vmax;
% p(2) = curv ;  % v_norm=0:.002:1; % auf vmax normierte Geschwindigkeit
% p(3) = fiso

%Variante A) auf f normiert aber nicht auf vmax normiert 
force=(p(1)-xdata)./(p(1)+xdata/p(2)); % Verwenden wenn Kraft auf 1 normiert ist

% % %Variante B) nicht auf vmax normiert und nicht auf f normiert
% % fiso=evalin('base','fiso')
% % force=(p(1)-xdata)./(p(1)+xdata/p(2))*fiso;

%Variante C) auf vmax normiert
%force=(1-xdata)./(1+xdata/p(1));