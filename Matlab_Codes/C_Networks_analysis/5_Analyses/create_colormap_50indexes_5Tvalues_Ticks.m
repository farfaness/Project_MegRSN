% Script to create colormap by customizing the 64 index with triplet RGB

% %Index of colormap
% %from t=2.4450 to t=3 ; orange fonc?
% for k=1:10, mymap_5hotColor(k, :) = [0.9 0.5 0], end
% %from t=3 to t=3.5 ; orange clair
% for k=11:20, mymap_5hotColor(k, :) = [1 0.7 0], end
% %from t=3.5 to t=4 ; jaune p?le
% for k=21:30, mymap_5hotColor(k, :) = [1 0.9 0.2], end
% %from t=4 to t=4.5; jaune tr?s p?le
% for k=31:50, mymap_5hotColor(k, :) = [1 1 0.4], end
% %from t=4.5 to t=4.7805; blanc
% for k=41:50, mymap_5hotColor(k, :) = [1 1 1], end
% 
% Index_mymap_5hotColor = ['t=2.4450' ; 't=3.0   '; 't=3.5   '; 't=4.0   '; 't=4.7805']

%Index of colormap
%from t=2.4450 to t=3 ; orange fonc?
for k=1:12, mymap_5hotColor.color(k, :) = [0.9 0.5 0], end
%from t=3 to t=3.5 ; orange clair
for k=13:23, mymap_5hotColor.color(k, :) = [1 0.7 0], end
%from t=3.5 to t=4 ; jaune p?le
for k=24:33, mymap_5hotColor.color(k, :) = [1 0.9 0.2], end
%from t=4 to t=4.5; jaune tr?s p?le
for k=34:43, mymap_5hotColor.color(k, :) = [1 1 0.4], end
%from t=4.5 to t=4.7805; blanc
for k=44:50, mymap_5hotColor.color(k, :) = [1 1 1], end

mymap_5hotColor.Index = ['t=2.4450' ; 't=3.0   '; 't=3.5   '; 't=4.0   '; 't=4.5   '; 't=4.7805']
mymap_5hotColor.Ticks = [2.4450 ; 3.0; 3.5 ; 4.0 ; 4.5; 4.7805]

save('MyColormap','mymap_5hotColor', '-append')

colormapeditor %To see the data values associated to each index