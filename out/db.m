function [pension] = db(AIME)
db(1) =    0.92791401924319428 ;%   0.89546759290755285;%2.5341578829431910 ;%0.51638944535881093;
db(2) =   -6.5225371224241505E-006;%-6.9447e-06;%-7.9562281381743094E-002 ;% -1.2816038704094608E-008;% -9.6197733540772514E-006;
bound = -db(1) /(2.0*db(2));
 
pension =  (db(2))*AIME.^2 + db(1)*AIME;
pension(AIME > bound) =(db(2))*bound.^2 + db(1)*bound;
% % delta(3) = 8.2270; %mean([9.083298    8.5256042   7.906901    7.392151]);
% % delta(2)=0.0767; %!mean([0.029333  0.061582208 0.094386    0.121426]);
% % delta(1)=-7.2503e-04; %mean([-0.00023 -0.0005901  -0.00091    -0.00117]);75
% % 
%  t=1:60;
% % %figure()plot(t,polyval(delta,t))
% % 
% %fc = table2array(readtable('1bias.txt'));
% fc =zeros(1,60)';
% %plot(t,exp(polyval(delta,t)-fc(1:60)'))
% %figure()
% %plot(t,exp(polyval(delta,t)))
% hold on
% delta(3) = 9.083298; %8.2270; %mean([9.083298    8.5256042   7.906901    7.392151]);
% delta(2)=0.029333; %0.0767; %!mean([0.029333  0.061582208 0.094386    0.121426]);
% delta(1)=-0.00033; %-7.2503e-04; %mean([-0.00023 -0.0005901  -0.00091    -0.00117]);75
% 
% 
% plot(21:80,exp(polyval(delta,t)-fc(1:60)'))
% [~,loc]=max(exp(polyval(delta,t)-fc(1:60)'))
% hold on
% fc = table2array(readtable('1bias.txt'));
% plot(21:80,exp(polyval(delta,t)-fc(1:60)'))
% [~,loc]=max(exp(polyval(delta,t)-fc(1:60)'))
% grid on
% legend('DATA','No fc','fc');
%     
% % 
% % plot(AIME,pension)
% % hold on
% % pension =  -(6.9447e-06)*AIME.^2 + 0.8334*AIME;
% % pension(AIME > 60000) = 25000;
% % plot(AIME,pension)
% % legend('Estimated','Calibrated')
% 
% 
% % AIME = 1.2525e+04; 
% % db(1) = 0.893279806195867;
% % db(2) = -8.223557117746168E-002;
% % bound = -db(1) /(2.0*db(2));
% % pension =  (db(2))*AIME.^2 + db(1)*AIME;
% % max((db(2))*bound.^2 + db(1)*bound,pension)
% % 
% %  AIME = 1:70000;
% %  pension =  (db(2))*AIME.^2 + db(1)*AIME;
% %  pension(AIME > bound) =(db(2))*bound.^2 + db(1)*bound;
% %  
% %  plot(AIME,pension)
% % 

end