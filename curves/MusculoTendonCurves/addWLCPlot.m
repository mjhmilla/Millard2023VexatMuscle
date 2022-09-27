function [figH,lengthV,forceV,forceWLCV] ...
        = addWLCPlot(figH, forceLengthTitinCurve,...
            lengthAtOneFpeHalf,lengthContour,normForceFailure,...
            curveColor,curveName)

x    = lengthAtOneFpeHalf;
y    = calcBezierYFcnXDerivative(x,forceLengthTitinCurve,0);
dydx = calcBezierYFcnXDerivative(x,forceLengthTitinCurve,1);

lo = x-(y/dydx);
lengthContourUpd=lengthContour-lo;
% Solve for WLC coefficient s.t. the two curves are equal at fiso
ca = 1;
fPaWLC=calcWormLikeChainModelDer(lengthAtOneFpeHalf-lo,...
                                  lengthContourUpd,...
                                  ca,ca*ca,[0,0]);
fPa = calcBezierYFcnXDerivative(lengthAtOneFpeHalf,...
                forceLengthTitinCurve,0);

ca = fPa/fPaWLC;
% Solve for WLC length at 5.5*fiso, the length

lengthFailure = calcWormLikeChainModelInvDer(normForceFailure,lengthContourUpd,... 
                                      ca,ca*ca,[0,0]);
lengthFailure = lengthFailure + lo;

%fNTestP=calcWormLikeChainModelDer(lengthFailure-lo,...
%    lengthContour,ca,ca*ca,[0,0]);

n=100;
n01 = [0:(1/(n-1)):1]';    
lengthV = n01.*(lengthFailure);
forceV = zeros(size(lengthV));
forceWLCV = zeros(size(lengthV));


for i=1:1:n
    forceV(i,1) = calcBezierYFcnXDerivative(lengthV(i,1),...
                forceLengthTitinCurve,0);

    if(lengthV(i,1)>lo)
        forceWLCV(i,1)= calcWormLikeChainModelDer(lengthV(i,1)-lo,...
                                      lengthContourUpd,...
                                      ca,ca*ca,[0,0]);
    else
        forceWLCV(i,1)=0;
    end
end

zLength = lengthV./lengthContourUpd;

xTxt = lengthV(end);
yTxt = forceV(end);
text(xTxt,yTxt,curveName,'Color',curveColor,...
    'VerticalAlignment','bottom','HorizontalAlignment','left');
hold on;
xTxt = lengthV(end);
yTxt = forceV(end);
text(xTxt,yTxt,'WLC','Color',[1,1,1].*0.5,...
    'VerticalAlignment','top','HorizontalAlignment','left');
hold on;

plot(lengthV,forceWLCV,'--','Color',[1,1,1].*0.5,'LineWidth',2);
hold on;
plot(lengthV,forceV,'-','Color',curveColor);
xlabel('Norm. Length ($$\ell/\ell_{o}^{M}$$)')
ylabel('Norm. Force');
%title(curveName);
