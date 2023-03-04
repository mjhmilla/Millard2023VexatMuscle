function fig = plotStructOfBezierSplines( curveStruct ,...
                                          ignoreStructsWithTheseKeyWords)


curveNames =fieldnames(curveStruct);

fig  = [];

for i=1:1:length(curveNames)

  flag_ignore=0;
  for j=1:1:length(ignoreStructsWithTheseKeyWords)
    idxKeyWord = strfind(curveNames{i},ignoreStructsWithTheseKeyWords{j});
    if(isempty(idxKeyWord)==0)
        flag_ignore=1;
    end
  end
  if(isempty(curveStruct.(curveNames{i})) == 0 ...
     && flag_ignore==0)
    
    fig.(curveNames{i}) = figure;

    curveSample = calcBezierYFcnXCurveSampleVector(...
                    curveStruct.(curveNames{i}), 200,[]);

    xmin = min(curveSample.x);
    xmax = max(curveSample.x);
    ymin = min(curveSample.y);
    ymax = max(curveSample.y);

    xV   = curveSample.x;
    yV   = curveSample.y;
    y1V  = curveSample.dydx;
    y2V  = curveSample.d2ydx2;

    subplot(2,2,1);
    plot(curveSample.x, curveSample.y,'k');
      hold on;        

      xlabel('x');
      ylabel('y');
      title(curveNames{i});
      grid on;
      axis square;
      box off;            
      xlim([xmin,xmax]);

    subplot(2,2,2);        
    plot(curveSample.x, curveSample.dydx,'r');
      hold on;            
      xlabel('x');
      ylabel('dy/dx');
      grid on;
      axis square;
      box off;            
      xlim([xmin,xmax]); 

    subplot(2,2,3);        
    plot(curveSample.x, curveSample.d2ydx2,'b');
      hold on;            
      xlabel('x');
      ylabel('d2y/dx2');
      grid on;
      axis square;
      box off;            
      xlim([xmin,xmax]); 
    if(isempty(curveStruct.(curveNames{i}).integral)==0)  
        intYdx = curveSample.intYdx;

        subplot(2,2,4)
            plot(curveSample.x, intYdx,'g');
            hold on;
            xlabel('x');
            ylabel(['int(y)']);
            axis square;
            box off;                
            grid on;
            hold on;
            xlim([xmin,xmax]);  

    end
  end
  
end


