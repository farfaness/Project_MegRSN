% Script to create colormap by customizing the 64 index with triplet RGB

%Index of colormap
%from t=2.10 to t=2.9 (p=0.05 to p= 0.01); green
for k=1:19, mymap_GBOY(k, :) = [0 1 0], end
%from t=2.9 to t=3.2 (p=0.01); blue
for k=20:26, mymap_GBOY(k, :) = [0 0 1], end
%from t=3.2 to t=3.9 (p=0.005); orange
for k=27:42, mymap_GBOY(k, :) = [1 0.5 0], end
%from t=3.9 to t=4.8; yellow
for k=43:64, mymap_GBOY(k, :) = [1 1 0], end

save('MyColormap','mymap_GBOY', '-append')


colormapeditor %To see the data values associated to each index