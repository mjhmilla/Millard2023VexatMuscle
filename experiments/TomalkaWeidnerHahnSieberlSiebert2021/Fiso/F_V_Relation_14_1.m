clear;
close all;
clc

% Data Import

cd 'C:\Users\bwpc\Documents\Data\F-V-Relation\Fiso'
Data1 = dlmread(['Fiso_018mN_3_exp.dat'],'',4,0);

Data2 = dlmread(['Fiso_028mN_4_exp.dat'],'',4,0);

Data3 = dlmread(['Fiso_022mN_5_exp.dat'],'',4,0);

Data4 = dlmread(['Fiso_013mN_6_exp.dat'],'',4,0);

Data5 = dlmread(['Fiso_022mN_7_exp.dat'],'',4,0);

Data6 = dlmread(['Fiso_009mN_8_exp.dat'],'',4,0);

Data7 = dlmread(['Fiso_025mN_9_exp.dat'],'',4,0);

Data8 = dlmread(['Fiso_003mN_10_exp.dat'],'',4,0);




% Defining time columns to plot

time=Data1(:,1)-66.3;
time2=Data2(:,1)-66.3;
time3=Data3(:,1)-66.3;
% time4=Data7(:,1)-66.4;
% time5=Data9(:,1)-66.4;
% time6=Data11(:,1)-66.4;
% time7=Data13(:,1)-66.4;

% Defining length columns to plot

length1=Data1(:,2);
length2=Data2(:,2);
length3=Data3(:,2);
length4=Data4(:,2);
length5=Data5(:,2);
length6=Data6(:,2);
length7=Data7(:,2);
length8=Data8(:,2);




% Calculation of forces and mean maximal isometric force and mean plateau
% isotonic force

force1=Data1(:,3);
force2=Data2(:,3);
force3=Data3(:,3);
force4=Data4(:,3);
force5=Data5(:,3);
force6=Data6(:,3);
force7=Data7(:,3);
force8=Data8(:,3);


meanF1 = mean(force1(28000:28500));
meanF2 = mean(force2(28000:28500));
meanF3 = mean(force3(28000:28500));
meanF4 = mean(force4(28000:28500));
meanF5 = mean(force5(28000:28500));
meanF6 = mean(force6(28000:28500));
meanF7 = mean(force7(28000:28500));
meanF8 = mean(force8(28000:28500));


meanISOf1 = mean (force1(28820:28950));
meanISOf2 = mean (force2(28820:28950));
meanISOf3 = mean (force3(28820:28950));
meanISOf4 = mean (force4(28820:28950));
meanISOf5 = mean (force5(28820:28950));
meanISOf6 = mean (force6(28820:28950));
meanISOf7 = mean (force7(28820:28950));
meanISOf8 = mean (force8(28820:28950));



% Calculate percentage of max isometric force

relF1 = meanISOf1/meanF1;
relF2 = meanISOf2/meanF2;
relF3 = meanISOf3/meanF3;
relF4 = meanISOf4/meanF4;
relF5 = meanISOf5/meanF5;
relF6 = meanISOf6/meanF6;
relF7 = meanISOf7/meanF7;
relF8 = meanISOf8/meanF8;



% Calculate shortening velocity

velocity1= (length1(28950)-length1(28000))/(time(28950)-time(28000));

velocity2= (length2(28950)-length2(28000))/(time(28950)-time(28000));
velocity3= (length3(28950)-length3(28000))/(time(28950)-time(28000));
velocity4= (length4(28950)-length4(28000))/(time(28950)-time(28000));
velocity5= (length5(28950)-length5(28000))/(time(28950)-time(28000));
velocity6= (length6(28950)-length6(28000))/(time(28950)-time(28000));
velocity7= (length7(28950)-length7(28000))/(time(28950)-time(28000));
velocity8= (length8(28950)-length8(28000))/(time(28950)-time(28000));


% Off-setting force Data

offset1=0-force1(1);
force1=force1+offset1;

offset2=0-force2(1);
force2=force2+offset2;

offset3=0-force3(1);
force3=force3+offset3;

offset4=0-force4(1);
force4=force4+offset4;

offset5=0-force5(1);
force5=force5+offset5;

offset6=0-force6(1);
force6=force6+offset6;

offset7=0-force7(1);
force7=force7+offset7;

offset8=0-force8(1);
force8=force8+offset8;




% Putting all SSC Trials into one array

% SSCall=[SSC1,SSC2,SSC4,SSC5,SSC7,SSC8,SSC10,SSC11,SSC12,SSC13,SSC14,SSC15,SSC16,SSC17,SSC18,SSC19,SSC20,SSC21,SSC22,SSC23,SSC24,SSC25,SSC26,SSC27,SSC28,SSC29,SSC30,SSC31,SSC32,SSC33,SSC34,SSC35,SSC36,SSC37,SSC38,SSC39,SSC40,SSC44,SSC42,SSC43,SSC45,SSC46,SSC47];
% 
% SSCmean= mean(SSCall,2);
% SSCstd= std(SSCall,1,2);
% 
% 
% SSCmeanPstd(:,1)= SSCmean(:,1) + SSCstd(:,1);
% SSCmeanMstd(:,1)= SSCmean(:,1) - SSCstd(:,1);



% Glätten der Messdaten

iN = 10; % Länge des Filters
vg1 = filter(ones(1,iN)/iN, 1, force1);
iN = 10;
vg2 = filter(ones(1,iN)/iN, 1, force2);
iN = 10;
vg3 = filter(ones(1,iN)/iN, 1, force3);
iN = 10;
vg4 = filter(ones(1,iN)/iN, 1, force4);
iN = 10;
vg5 = filter(ones(1,iN)/iN, 1, force5);
iN = 10;
vg6 = filter(ones(1,iN)/iN, 1, force6);
iN = 10;
vg7 = filter(ones(1,iN)/iN, 1, force7);
iN = 10;
vg8 = filter(ones(1,iN)/iN, 1, force8);
iN = 10;


% Plotting the Data

figure (1)

% subplot(2,1,1)
hold on
grid on
plot(time,vg1,'-','linewidth',1);
plot(time2,vg2,'-','linewidth',1);
plot(time3,vg3,'-','linewidth',1);
plot(time,vg4,'-','linewidth',1);
plot(time,vg5,'-','linewidth',1);
plot(time,vg6,'-','linewidth',1);
plot(time,vg7,'-','linewidth',1);
plot(time,vg8,'.','linewidth',1);

% plot(time,SSC10,'.','linewidth',1);


hold off
title({'Plot of EDL-SSCs with same amount of stretch and diffrent velocities';});
set(gca,'fontsize',14,'linewidth',1)
xlabel('Time [s]');
ylabel('force [F_0]');
legend('40% V_m_a_x','1% V_m_a_x','2,5% V_m_a_x','5% V_m_a_x','1,75% V_m_a_x')

figure (2)

hold on;
grid on;

plot(time,length1);
plot(time,length2);
plot(time,length3);
plot(time,length4);
plot(time,length5);
plot(time,length6);
plot(time,length7);
plot(time,length8);


hold off;

title({'Plot of SSCs mean with std (n=43)';});
set(gca,'fontsize',14,'linewidth',1)
xlabel('Time [s]');
ylabel('force [F_0]');
legend('SSC 082 1 082 @ 40% V_m_a_x','mean minus std','mean plus std')


