% load Si.mat

PolySi = importdata('PolySi.mat');

wl = PolySi(:,1)/1000;
N= PolySi(:,2);
MSiO2 =[wl,N];

f = fit(wl,N,"poly2")
curveFitter

% x_PolySi = PolySi(:,1)/1000 %in um
% y_PolySi = PolySi(:,2);
% myfittype = fittype('a + b/(x^2)',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a','b'}) 
% 
% myfittype = fittype('sqrt(a+b*x.^2/(x.^2-c)+d*x.^2/(x.^2-e))',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a','b','c','d','e'}) 
% 
% mfp = sqrt(1+0.286044141+1.07044083*x.^2/(x.^2-0.0100585997)+1.10202242*x.^2/(x.^2-100))
% 
% myfit = fit(x,y,myfittype)
% % 
% plot(x,y) %, hold on
% plot(x,mfp), hold on
% legend('Original','Equation','wtf')
% 
% cauchycoeff = coeffvalues(myfit)
% cauchycoeff(2)*1000000000

% https://refractiveindex.info/?shelf=main&book=SiO2&page=Radhakrishnan-o
% https://refractiveindex.info/?shelf=main&book=SiO2&page=Ghosh-o
%fittype https://se.mathworks.com/help/curvefit/fittype.html

a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + 
              a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2) + 
              a5*exp(-((x-b5)/c5)^2) + a6*exp(-((x-b6)/c6)^2) + 
              a7*exp(-((x-b7)/c7)^2)