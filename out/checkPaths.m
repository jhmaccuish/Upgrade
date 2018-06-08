clear;
close all;

lpath = table2array(readtable('lpath'));
figure(); 
plot(21:(20+85),lpath);
title('Labour Supply')

cpath = table2array(readtable('Cpath'));
figure(); 
plot(21:(20+85),cpath);
title('Consumption')

apath = table2array(readtable('Apath'));
figure(); 
plot(21:(20+85),apath);
title('Assets')

% Ypath = table2array(readtable('YempPath'));
% figure();
% plot(21:(20+60),Ypath(1:60));
% title('Income')
% 
% AIMEpath = table2array(readtable('AIMEPath'));
% figure();
% plot(21:(20+85),AIMEpath);
% title('AIME')
% 
% 
% Incpath = zeros(1,85);
% for ixt=1:85
%     Incpath(ixt) = Ypath(ixt);
%     if ixt >= 40
%         Incpath(ixt) = Incpath(ixt) +  107.45*52; %params%pension
%         if  ixt >= 45
%             Incpath(ixt) = Incpath(ixt) +  107.45*52; %params%pension
%             %if (ixL==0) then
%             %    Y =   Y + dbPension(params,AIME);
%             %end if
%         else
%             Incpath(ixt) = Incpath(ixt) +  6235.89;% params%spouseInc
%         end
%     else
%         Incpath(ixt) = Incpath(ixt) +  6235.89; %params%spouseInc
%     end
% end
% 
% benefit = [57.90*52*ones(5,1);73.10*52*ones(80,1)];
% 
% for ixt=1:85
%     if ixt >= 40
%         benefit(ixt) = benefit(ixt) +  107.45*52; %params%pension
%         if  ixt >= 45
%             benefit(ixt) = benefit(ixt) +  107.45*52; %params%pension
%             
%             benefit(ixt) =   benefit(ixt) + db(AIMEpath(45)); %dbPension(params,AIME);
%         else
%             benefit(ixt) = benefit(ixt) +  6235.89;% params%spouseInc
%         end
%     else
%         benefit(ixt) = benefit(ixt) +  6235.89; %params%spouseInc
%     end
% end
% 
% figure()
% plot(21:(20+85),Incpath);
% hold on
% plot(21:(20+85),benefit);
% legend('Employed','Unemployed')
% title('Income')



% Richpath = table2array(readtable('RichPath'));
% figure(); 
% %plot(21:(20+85),Richpath);
% plot(52:75,Richpath(32:32+23))
% title('Rich Labour Supply')
% 
% Poorpath = table2array(readtable('PoorPath'));
% figure(); 
% %plot(21:(20+85),Poorpath);
% plot(52:75,Poorpath(32:32+23))
% title('Poor Supply')



figure();
moments = xlsread('../../moments/Moments.xlsx');
%moments = mean(moments,1);
plot(52:75,moments,52:75,lpath(32:32+23))
title('Mean Partipcation Rate:Simulation vs Data Labour')
xlabel('Age')
ylabel('Mean Partipcation Rate')
legend('Data','Simulation')

% figure();
% Amoments = xlsread('../../moments/AssetMoments.xlsx');
% %moments = mean(moments,1);
% plot(1:24,Amoments,1:24,apath(32:32+23))
%title('Sim vs Data Assets')



%%%%%%%%%%%%%%%%%%
% figure();
% moments = xlsread('./inworl.xlsx');
% moments = mean(moments,1);
% plot(1:24,moments,1:24,lpath(32:32+23))
% title('Sim vs Data')
% 
% figure();
% plot(1:24,moments)
% title('Data')
% gmm = abs(lpath(32:32+23)'-moments)*abs(lpath(32:32+23)'-moments)';
% 
% %figure();
% paths = xlsread('./Lpath52to85.xlsx');
% %moments = mean(moments,1);
% %plot(1:24,moments,1:24,lpath(32:32+23))
% %title('Sim vs Data')
