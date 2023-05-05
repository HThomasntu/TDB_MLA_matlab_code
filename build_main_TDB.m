% Copyright, M.Bencsik, H.Thomas, 2023

close all 
clear all 


[A S_R] = audioread('G:\Harriet Data\Varroa Vibration\petri_dish_video_analysis\varroa_brood_comb\sound_500.flac', [1 2]);
[A S_R] = audioread('G:\Harriet Data\Varroa Vibration\petri_dish_video_analysis\varroa_brood_comb\sound_500.flac', [1 round((40*60)*S_R)]);

[B S_R] = audioread('G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\500 hz bee file 27.07.21 2.wav', [1 2]);
[B S_R] = audioread('G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\500 hz bee file 27.07.21 2.wav', [1 round((((1*60)*60)+35*60)*S_R)]);

[C S_R] = audioread('G:\Harriet Data\Varroa Vibration\petri_dish_video_analysis\varroa_brood_comb\sound_500.flac', [1 2]);
[C S_R] = audioread('G:\Harriet Data\Varroa Vibration\petri_dish_video_analysis\varroa_brood_comb\sound_500.flac', [1 round(((40*60))*S_R)]);

[D S_R] = audioread('G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\500 hz bee file 27.07.21 2.wav', [1 2]);
[D S_R] = audioread('G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\500 hz bee file 27.07.21 2.wav', [1 round(((2*60)+5)*S_R)]);

[E S_R] = audioread('G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\varroa 4 20.09.21_500hz.wav', [1 2]);
[E S_R] = audioread('G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\varroa 4 20.09.21_500hz.wav', [1 round(((1*60))*S_R)]);

[F S_R] = audioread('E:\Additional brood videos\20.08.21 sample 1_500hz.wav', [1 2]);
[F S_R] = audioread('E:\Additional brood videos\20.08.21 sample 1_500hz.wav', [1 round(((28*60))*S_R)]);

A = A*0.02651;
B = B*0.02651;
C = C*0.02651;
E = E*0.02651;
F = F*0.02651;


% Signals of interest (in seconds):

mite_timings = [571 927.5 1184 1185 1326 1329 1329.5 1329.8 1330 1330.5 1330.8 1335 1336 1338.8 1382 1407.5 1436 1852 1931 1936 1938 1939 1945.5 1947 1948 1951.5 1996.5 2039 2040 2251 2252 2253 2327 2328 2329];
bee_timings = [3367 3368 3369 3370 3371 3376 3380 3388 3390.5 3390.8 3392 3393 3393.3 3394 3396 3397 3397.5 3397.8 3400 3400.8 3401 3402 3402.5 3402.8 3403 3403.3 3404 3406 3408 3408.8 3409 3410 3410.8 3411 3413 3414 3415 3417 3417.5 3421 3424 3429 3432 3438 3443 3454 3455 3459 3468 3480 3504 3810 3811 3812 3815 3817 3818 3819 3820 3821 3822 3823 3824 3839 3840 3853 3854 3855 3856 5415 5418 5420 5422 5426 5429 5430 5433 5434 5440 5441 5444 5445 5446 5450 5451 5452 5453];
mite_background = [5 6 7 39 40 41 45.5 53 54 55 56 57 58 65 66 67 71 79 80 83 86 87 88 89 92.5 93.5 94.5 96 97 98 99 100 101 104.5 105.5 106 109.5 110.5 111 112.5 113.5 115 117.5 118.5 1321.5 1322.5 1323.5 1324 1325.5 1341.5 1343.5 1344.5 1348 1349.5];
bee_background = [30 31 32 33 34 35 36 37 38 39 40 49 50 63 64 88 89 90 91 92 93 94 95 96 97 98 119 120];
V4_mite = [40.75 41 41.25 41.5 41.75 42 42.25 42.5 42.75 43 43.25 43.5 43.75 44 44.25 44.5 44.75 45 45.25 45.5 45.75 46 46.25 46.5 46.75 47 47.25 47.5 47.75 48 48.5];
V4_bg = [49.75 50 50.25 50.5 50.75];
bee_video = [1500.8 1501 1501.5 1501.8 1502 1502.3 1502.5 1502.8 1503 1503.3 1503.5 1503.8 1504.3 1504.5 1504.8 1505 1505.5 1505.8 1506 1506.3 1506.5 1506.8 1507 1507.3 1507.8 1508 1508.5 1508.8 1509 1509.8 1510 1510.3 1510.5 1510.8 1511 1511.3 1511.5 1511.8 1512 1512.3 1512.5 1512.8 1513 1513.3 1513.5 1515.8 1516 1516.8 1518 1519.5 1519.8 1520 1521.5 1522.3 1522.5 1522.8 1523.8 1524 1524.3 1524.5 1524.8 1525.5 1525.8 1526 1526.3 1526.5 1527 1527.3 1528 1528.3 1528.5 1528.8 1529.5 1529.8 1530.8 1531 1531.3 1531.5 1531.8 1532.3 1532.5 1532.8 1533 1533.3 1533.5 1534 1534.3 1534.8 1538.5 1539.5 1540.8 1541 1541.3 1541.5 1541.8 1542.8 1543 1543.3 1543.5 1543.8 1544 1544.3 1544.5 1544.8 1545 1545.3 1545.5 1545.8 1546 1546.3 1546.5 1546.8 1547 1547.3 1547.5 1547.8 1549 1550 1550.3 1550.5 1550.8 1551 1551.3 1551.5 1551.8 1552 1552.3 1552.5 1552.8 1553 1553.3 1553.5 1553.8 1554 1554.3 1554.5 1554.8 1555 1555.3 1555.5 1555.8 1556 1556.3 1556.5 1556.8 1557 1557.3 1557.5 1557.8 1558 1558.3 1558.5 1558.8 1559 1559.3];

% mite_timings = [571 927.5 1184 1185 1326 1326.3 1329.5 1329.8 1330 1330.8 1335 1336 1336.3 1338.8 1382 1407.5 1436 1852 1931 1936 1938 1939 1945.5 1947 1948 1951.5 1996.5 2039 2040 2251 2252 2253 2327 2328 2329];
% bee_timings = [3367 3368 3369 3370 3371 3376 3380 3388 3390.8 3392 3393.3 3394 3396 3397 3397.8 3400.5 3401 3202.5 3402.8 3403 3404 3406 3408.8 3409 3410.5 3410.8 3411 3411.3 3413 3414 3415 3417 3417.5 3421 3424 3429 3432 3438 3443 3454 3455 3459 3468 3480 3504 3810 3811 3812 3815 3817 3818 3819 3820 3821 3822 3823 3824 3839 3840 3853 3854 3855 3856 5415 5418 5420 5422 5426 5429 5430 5433 5434 5440 5441 5444 5445 5446 5450 5451 5452 5453];
% mite_background = [5 6 7 39 40 41 45.5 53 54 55 56 57 58 65 66 67 71 79 80 83 86 87 88 89 92 92.5 93.5 94 96 97 98 99 100 101 104 104.5 105.75 106 109 109.5 110.75 111 112.5 113 113.5 115 117.75 118.5 1321.8 1322 1322.3 1322.5 1322.8 1323 1323.3 1323.5 1324 1324.3 1325 1325.3 1325.5 1325.8 1341.8 1343.5 1343.8 1344 1344.3 1344.5 1344.8 1345 1348 1348.3 1349 1349.3];
% bee_background = [30 31 32 33 34 35 36 37 38 39 40 49 50 63 64 88 89 90 91 92 93 94 95 96 97 98 119 120];

% Convert signals into 2DFTs: 

window_length = 0.5;

MF = 4;
tr = 0.015;


for periodicity_section = 1:length(mite_timings)
    
LL = round((mite_timings(periodicity_section)-window_length)*S_R);
UL = round((mite_timings(periodicity_section)+window_length)*S_R);

mite_tdft(:,:,periodicity_section) = two_D_FT_Gaussian(A(LL:UL,1),MF,tr,S_R,size(A,1)/(2*S_R));

end 

for periodicity_section2 = 1:length(bee_timings)
    
LL2 = round((bee_timings(periodicity_section2)-window_length)*S_R);
UL2 = round((bee_timings(periodicity_section2)+window_length)*S_R);

bee_tdft(:,:,periodicity_section2) = two_D_FT_Gaussian(B(LL2:UL2,1),MF,tr,S_R,size(B,1)/(2*S_R));

end 

for periodicity_section3 = 1:length(mite_background)
    
LL3 = round((mite_background(periodicity_section3)-window_length)*S_R);
UL3 = round((mite_background(periodicity_section3)+window_length)*S_R);

mite_bg_tdft(:,:,periodicity_section3) = two_D_FT_Gaussian(C(LL3:UL3,1),MF,tr,S_R,size(C,1)/(2*S_R));

end 

for periodicity_section4 = 1:length(bee_background)
    
LL4 = round((bee_background(periodicity_section4)-window_length)*S_R);
UL4 = round((bee_background(periodicity_section4)+window_length)*S_R);

bee_bg_tdft(:,:,periodicity_section4) = two_D_FT_Gaussian(D(LL4:UL4,1),MF,tr,S_R,size(D,1)/(2*S_R));

end 

for periodicity_section5 = 1:length(V4_mite)
    
LL5 = round((V4_mite(periodicity_section5)-window_length)*S_R);
UL5 = round((V4_mite(periodicity_section5)+window_length)*S_R);

V4_mite_tdft(:,:,periodicity_section5) = two_D_FT_Gaussian(E(LL5:UL5,1),MF,tr,S_R,size(E,1)/(2*S_R));

end 

for periodicity_section6 = 1:length(V4_bg)
    
LL6 = round((V4_bg(periodicity_section6)-window_length)*S_R);
UL6 = round((V4_bg(periodicity_section6)+window_length)*S_R);

V4_bg_tdft(:,:,periodicity_section6) = two_D_FT_Gaussian(E(LL6:UL6,1),MF,tr,S_R,size(E,1)/(2*S_R));

end 

for periodicity_section7 = 1:length(bee_video)
    
LL7 = round((bee_video(periodicity_section7)-window_length)*S_R);
UL7 = round((bee_video(periodicity_section7)+window_length)*S_R);

new_bee_tdft(:,:,periodicity_section7) = two_D_FT_Gaussian(F(LL7:UL7,1),MF,tr,S_R,size(F,1)/(2*S_R));

end 

% interpolate the 2DFTs:

for uu = 1:size(mite_tdft,3)
    for vv = 1:size(mite_tdft,1)
        mite_tdft_int(vv,:,uu) = interp1(1:size(mite_tdft,2),mite_tdft(vv,:,uu),1:5:size(mite_tdft,2));
    end
%     subplot(2,2,1)
%     imagesc(log10(mite_tdft(:,:,uu)))
%     subplot(2,2,2)
%     imagesc(log10(mite_tdft_int(:,:,uu)))
%     pause(1)
end 

for uu2 = 1:size(bee_tdft,3)
    for vv2 = 1:size(bee_tdft,1)
        bee_tdft_int(vv2,:,uu2) = interp1(1:size(bee_tdft,2),bee_tdft(vv2,:,uu2),1:5:size(bee_tdft,2));
    end
%     subplot(2,1,1)
%     imagesc(log10(bee_tdft(:,:,uu2)))
%     subplot(2,1,2)
%     imagesc(log10(bee_tdft_int(:,:,uu2)))
%     pause(1)
end 

for uu3 = 1:size(mite_bg_tdft,3)
    for vv3 = 1:size(mite_bg_tdft,1)
        mite_bg_tdft_int(vv3,:,uu3) = interp1(1:size(mite_bg_tdft,2),mite_bg_tdft(vv3,:,uu3),1:5:size(mite_bg_tdft,2));
    end
%     subplot(2,1,1)
%     imagesc(log10(mite_bg_tdft(:,:,uu3)))
%     subplot(2,1,2)
%     imagesc(log10(mite_bg_tdft_int(:,:,uu3)))
%     pause(1)
end 

for uu4 = 1:size(bee_bg_tdft,3)
    for vv4 = 1:size(bee_bg_tdft,1)
        bee_bg_tdft_int(vv4,:,uu4) = interp1(1:size(bee_bg_tdft,2),bee_bg_tdft(vv4,:,uu4),1:5:size(bee_bg_tdft,2));
    end
%     subplot(2,1,1)
%     imagesc(log10(bee_bg_tdft(:,:,uu4)))
%     subplot(2,1,2)
%     imagesc(log10(bee_bg_tdft_int(:,:,uu4)))
%     pause(1)
end 

for uu5 = 1:size(V4_mite_tdft,3)
    for vv5 = 1:size(V4_mite_tdft,1)
        V4_mite_tdft_int(vv5,:,uu5) = interp1(1:size(V4_mite_tdft,2),V4_mite_tdft(vv5,:,uu5),1:5:size(V4_mite_tdft,2));
    end
%     subplot(2,1,1)
%     imagesc(log10(bee_bg_tdft(:,:,uu4)))
%     subplot(2,1,2)
%     imagesc(log10(bee_bg_tdft_int(:,:,uu4)))
%     pause(1)
end 

for uu6 = 1:size(V4_bg_tdft,3)
    for vv6 = 1:size(V4_bg_tdft,1)
        V4_bg_tdft_int(vv6,:,uu6) = interp1(1:size(V4_bg_tdft,2),V4_bg_tdft(vv6,:,uu6),1:5:size(V4_bg_tdft,2));
    end
%     subplot(2,1,1)
%     imagesc(log10(bee_bg_tdft(:,:,uu4)))
%     subplot(2,1,2)
%     imagesc(log10(bee_bg_tdft_int(:,:,uu4)))
%     pause(1)
end 

for uu7 = 1:size(new_bee_tdft,3)
    for vv7 = 1:size(new_bee_tdft,1)
        new_bee_tdft_int(vv7,:,uu7) = interp1(1:size(new_bee_tdft,2),new_bee_tdft(vv7,:,uu7),1:5:size(new_bee_tdft,2));
    end
%     subplot(2,1,1)
%     imagesc(log10(bee_bg_tdft(:,:,uu4)))
%     subplot(2,1,2)
%     imagesc(log10(bee_bg_tdft_int(:,:,uu4)))
%     pause(1)
end 

% Crop the 2DFTs appropriately:

cr_bee_tdft = bee_tdft_int(3:60,:,:);
cr_mite_tdft = mite_tdft_int(3:60,:,:);
cr_mite_bg_tdft = mite_bg_tdft_int(3:60,:,:);
cr_bee_bg_tdft = bee_bg_tdft_int(3:60,:,:);
cr_V4_mite_tdft = V4_mite_tdft_int(3:60,:,:);
cr_V4_bg_tdft = V4_bg_tdft_int(3:60,:,:);
cr_new_bee_tdft = new_bee_tdft_int(3:60,:,:);


% Scale the 2DFTs by the mean: 

% for counter = 1:size(cr_bee_tdft,3)   
%     scaled_bee(:,:,counter) = (1/mean(cr_bee_tdft(:,1,counter)))*cr_bee_tdft(:,:,counter); 
% end 
% 
% for counter = 1:size(cr_mite_tdft,3)
%     scaled_mite(:,:,counter) = (1/mean(cr_mite_tdft(:,1,counter)))*cr_mite_tdft(:,:,counter); 
% end 
% 
% for counter = 1:size(cr_mite_bg_tdft,3)
%     scaled_mite_bg(:,:,counter) = (1/mean(cr_mite_bg_tdft(:,1,counter)))*cr_mite_bg_tdft(:,:,counter); 
% end 
% 
% for counter = 1:size(cr_bee_bg_tdft,3)
%     scaled_bee_bg(:,:,counter) = (1/mean(cr_bee_bg_tdft(:,1,counter)))*cr_bee_bg_tdft(:,:,counter); 
% end 
% 
% scaled_bg = cat(3, scaled_bee_bg, scaled_mite_bg);



% Scale the 2DFTs by the max: 

for counter = 1:size(cr_bee_tdft,3)
    bee_mag = max(cr_bee_tdft(:,1,counter)); 
    
    new_matrix(:,:,counter) = bee_mag; 
end 

scaled_bee = cr_bee_tdft./new_matrix;


for counter = 1:size(cr_mite_tdft,3)
    mite_mag = max(cr_mite_tdft(:,1,counter)); 
    
    new_matrix2(:,:,counter) = mite_mag; 
end 

scaled_mite = cr_mite_tdft./new_matrix2;


for counter = 1:size(cr_mite_bg_tdft,3)
    mite_bg_mag = max(cr_mite_bg_tdft(:,1,counter)); 
    
    new_matrix3(:,:,counter) = mite_bg_mag; 
end 

scaled_mite_bg = cr_mite_bg_tdft./new_matrix3;


for counter = 1:size(cr_bee_bg_tdft,3)
    bee_bg_mag = max(cr_bee_bg_tdft(:,1,counter)); 
    
    new_matrix4(:,:,counter) = bee_bg_mag; 
end 

scaled_bee_bg = cr_bee_bg_tdft./new_matrix4;


for counter = 1:size(cr_V4_mite_tdft,3)
    V4_mite_mag = max(cr_V4_mite_tdft(:,1,counter)); 
    
    new_matrix5(:,:,counter) = V4_mite_mag; 
end 

scaled_V4_mite = cr_V4_mite_tdft./new_matrix5;


for counter = 1:size(cr_V4_bg_tdft,3)
    V4_bg_mag = max(cr_V4_bg_tdft(:,1,counter)); 
    
    new_matrix6(:,:,counter) = V4_bg_mag; 
end 

scaled_V4_bg = cr_V4_bg_tdft./new_matrix6;


for counter = 1:size(cr_new_bee_tdft,3)
    new_bee_mag = max(cr_new_bee_tdft(:,1,counter)); 
    
    new_matrix7(:,:,counter) = new_bee_mag; 
end 

scaled_new_bee = cr_new_bee_tdft./new_matrix7;


%scaled_bg = cat(3, scaled_bee_bg, scaled_mite_bg);
scaled_bg = cat(3, scaled_bee_bg, scaled_mite_bg, scaled_V4_bg);
scaled_mite_final = cat(3, scaled_mite, scaled_V4_mite);
scaled_bee_final = cat(3, scaled_bee, scaled_new_bee);

save IMPROVED_tdft_TDB_interpolated_biggerTDB.mat
