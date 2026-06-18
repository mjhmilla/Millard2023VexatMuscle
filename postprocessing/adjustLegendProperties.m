function lgdH = adjustLegendProperties(lgdH)


lgdHeight = lgdH.Position(4) * 1.05;
lgdH.Position(2) = lgdH.Position(2) + 8*(lgdHeight-lgdH.Position(4));
lgdH.Position(4) = lgdHeight;
lgdH.Position(3) = lgdH.Position(3)*0.5;