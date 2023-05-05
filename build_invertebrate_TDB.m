% Copyright, M.Bencsik, H.Thomas, 2023

clear all

% % load the accelerometer tracks: 
% load invertebrate_recordings.mat
% load woodlouse_recordings.mat
% load beetle_recordings.mat
% load varroa_recordings.mat
% 
% % load the background removal information:
% load background_analysis.mat
% load woodlouse_background_analysis.mat
% load beetle_background_analysis.mat
% load varroa_background_analysis.mat

[A S_R] = audioread('G:\Harriet Data\Varroa Vibration\Varroa Videos\Sound Booth\New input on Sony with recorder\Individual Varroa.wav');
[B S_R] = audioread('G:\Harriet Data\Varroa Vibration\Varroa Videos\Sound Booth\New input on Sony with recorder\2 Varroa.wav');
[C S_R] = audioread('G:\Harriet Data\Varroa Vibration\Varroa Videos\Sound Booth\New input on Sony with recorder\5 Varroa.wav');
[D S_R] = audioread('G:\Harriet Data\Varroa Vibration\Varroa Videos\Sound Booth\New input on Sony with recorder\6 Varroa.wav');
[E S_R] = audioread('G:\Harriet Data\Varroa Vibration\Varroa Videos\Sound Booth\New input on Sony with recorder\7 Varroa.wav');
[F S_R] = audioread('G:\Harriet Data\Varroa Vibration\Varroa Videos\Sound Booth\New input on Sony with recorder\8 Varroa.wav');
[G S_R] = audioread('G:\Harriet Data\Varroa Vibration\Varroa Videos\Sound Booth\New input on Sony with recorder\19 Varroa.wav');
[H S_R] = audioread('G:\Harriet Data\Varroa Vibration\Varroa Videos\Sound Booth\New input on Sony with recorder\07.10.19.wav');

[I S_R] = audioread('F:\Insect discrimination videos\Woodlouse A.wav');
[J S_R] = audioread('F:\Insect discrimination videos\Woodlouse B.wav');
[K S_R] = audioread('F:\Insect discrimination videos\Woodlouse C.wav');
[L S_R] = audioread('F:\Insect discrimination videos\Woodlouse D 1.wav');

[M S_R] = audioread('F:\Insect discrimination videos\Beetle A 1.wav');
[N S_R] = audioread('F:\Insect discrimination videos\Beetle A 2.wav');
[O S_R] = audioread('F:\Insect discrimination videos\Beetle B 1.wav');
[P S_R] = audioread('F:\Insect discrimination videos\Beetle B 2.wav');
[Q S_R] = audioread('F:\Insect discrimination videos\Beetle D 1.wav');
[R S_R] = audioread('F:\Insect discrimination videos\Beetle D 2.wav');
[S S_R] = audioread('F:\Insect discrimination videos\Beetle E 1.wav');
[T S_R] = audioread('F:\Insect discrimination videos\Beetle E 2.wav');

load all_average_background.mat

varroa_times = [62:67 69:71 73:85 87 89:91 308:310 312:314 320:334 336 354:357 372 373 562 563 568:570 790:792 820:822 848 850];
woodlouse_times = [26:30 39:42 45:49 55 70 71 73:84 105:110 143:147 165:168 218:223];
beetle_times = [35 36 51:55 58 59 64 65 74 76:78 80 83 84 90:92 103 106:108 110 111 113 114 127 134 136 140:142 168:171 179:183 185:187 193 194 199 208 215 220 228:230 234 247 285 300:302 307:309];
gnat_times = [26:40 53:57 60:76 78:115 120:146 229:240 250:270];

woodlouse2 = [22:26 28:33 35 47:49 58 60:70 78:98 102 108:130 132:142 144 146 147 149:164 166:170];
woodlouse3 = [47:63 66:69 80:95 110:126 128:140 160:172 174 175];

beetle2 = [28:37 41:45 48 49 58 59 83:85 87 125:127 145 146 153 162 164 165];
beetle3 = [48:50 53:55 62 64 65 71 72 100:102 127 167];

varroa2 = [134 423:429 451:453 457 359 460 473:476 532 534:538 811:815 820:822 824:828];
varroa3 = [53:67 79 86 87 100 103 104];


window_length = 0.5;

MF = 4;
TR = 0.015;

for varroa_section = 1:length(varroa_times)
    
    LL = round((varroa_times(varroa_section)-window_length)*S_R);
    UL = round((varroa_times(varroa_section)+window_length)*S_R);
    
    varroa_tdft(:,:,varroa_section) = abs(two_D_FT_Gaussian_calibrated(A(LL:UL,1),MF,TR,S_R,size(A(LL:UL,1),1)/(2*S_R)) - background_varroaind); 
    
end

for woodlouse_section = 1:length(woodlouse_times)
    
    LL2 = round((woodlouse_times(woodlouse_section)-window_length)*S_R);
    UL2 = round((woodlouse_times(woodlouse_section)+window_length)*S_R);
    
    woodlouse_tdft(:,:,woodlouse_section) = abs(two_D_FT_Gaussian_calibrated(I(LL2:UL2,1),MF,TR,S_R,size(I(LL2:UL2,1),1)/(2*S_R))-background_woodlouseA); 
    
end

for beetle_section = 1:length(beetle_times)
    
    LL3 = round((beetle_times(beetle_section)-window_length)*S_R);
    UL3 = round((beetle_times(beetle_section)+window_length)*S_R);
    
    beetle_tdft(:,:,beetle_section) = abs(two_D_FT_Gaussian_calibrated(S(LL3:UL3,1),MF,TR,S_R,size(S(LL3:UL3,1),1)/(2*S_R))-background_beetleE1); 
    
end

for woodlouse_section2 = 1:length(woodlouse2)
    
    LL4 = round((woodlouse2(woodlouse_section2)-window_length)*S_R);
    UL4 = round((woodlouse2(woodlouse_section2)+window_length)*S_R);
    
    woodlouse2_tdft(:,:,woodlouse_section2) = abs(two_D_FT_Gaussian_calibrated(K(LL4:UL4,1),MF,TR,S_R,size(K(LL4:UL4,1),1)/(2*S_R))-background_woodlouseC); 
    
end

for woodlouse_section3 = 1:length(woodlouse3)
    
    LL5 = round((woodlouse3(woodlouse_section3)-window_length)*S_R);
    UL5 = round((woodlouse3(woodlouse_section3)+window_length)*S_R);
    
    woodlouse3_tdft(:,:,woodlouse_section3) = abs(two_D_FT_Gaussian_calibrated(L(LL5:UL5,1),MF,TR,S_R,size(L(LL5:UL5,1),1)/(2*S_R))-background_woodlouseD1); 
    
end

for beetle_section2 = 1:length(beetle2)
    
    LL6 = round((beetle2(beetle_section2)-window_length)*S_R);
    UL6 = round((beetle2(beetle_section2)+window_length)*S_R);
    
    beetle2_tdft(:,:,beetle_section2) = abs(two_D_FT_Gaussian_calibrated(O(LL6:UL6,1),MF,TR,S_R,size(O(LL6:UL6,1),1)/(2*S_R))-background_beetleB1); 
    
end

for beetle_section3 = 1:length(beetle3)
    
    LL7 = round((beetle3(beetle_section3)-window_length)*S_R);
    UL7 = round((beetle3(beetle_section3)+window_length)*S_R);
    
    beetle3_tdft(:,:,beetle_section3) = abs(two_D_FT_Gaussian_calibrated(R(LL7:UL7,1),MF,TR,S_R,size(R(LL7:UL7,1),1)/(2*S_R))-background_beetleD2); 
    
end

for varroa_section2 = 1:length(varroa2)
    
    LL8 = round((varroa2(varroa_section2)-window_length)*S_R);
    UL8 = round((varroa2(varroa_section2)+window_length)*S_R);
    
    varroa2_tdft(:,:,varroa_section2) = abs(two_D_FT_Gaussian_calibrated(D(LL8:UL8,1),MF,TR,S_R,size(D(LL8:UL8,1),1)/(2*S_R))-background_varroa6); 
    
end

for varroa_section3 = 1:length(varroa3)
    
    LL9 = round((varroa3(varroa_section3)-window_length)*S_R);
    UL9 = round((varroa3(varroa_section3)+window_length)*S_R);
    
    varroa3_tdft(:,:,varroa_section3) = abs(two_D_FT_Gaussian_calibrated(H(LL9:UL9,1),MF,TR,S_R,size(H(LL9:UL9,1),1)/(2*S_R))-background_varroa07); 
    
end

varroa_tdft = varroa_tdft(4:85,1:100,:);
woodlouse_tdft = woodlouse_tdft(4:85,1:100,:);
beetle_tdft = beetle_tdft(4:85,1:100,:);

woodlouse2_tdft = woodlouse2_tdft(4:85,1:100,:);
woodlouse3_tdft = woodlouse3_tdft(4:85,1:100,:);

beetle2_tdft = beetle2_tdft(4:85,1:100,:);
beetle3_tdft = beetle3_tdft(4:85,1:100,:);

varroa2_tdft = varroa2_tdft(4:85,1:100,:);
varroa3_tdft = varroa3_tdft(4:85,1:100,:);


%scale by the mean

for counter = 1:size(varroa_tdft,3)
    
    scaled_varroa(:,:,counter) = (1/mean(varroa_tdft(:,1,counter)))*varroa_tdft(:,:,counter);

end 

for counter = 1:size(woodlouse_tdft,3)
    
    scaled_woodlouse(:,:,counter) = (1/mean(woodlouse_tdft(:,1,counter)))*woodlouse_tdft(:,:,counter);

end 

for counter = 1:size(beetle_tdft,3)
    
    scaled_beetle(:,:,counter) = (1/mean(beetle_tdft(:,1,counter)))*beetle_tdft(:,:,counter);

end 


for counter = 1:size(woodlouse2_tdft,3)
    
    scaled_woodlouse2(:,:,counter) = (1/mean(woodlouse2_tdft(:,1,counter)))*woodlouse2_tdft(:,:,counter);

end 


for counter = 1:size(woodlouse3_tdft,3)
    
    scaled_woodlouse3(:,:,counter) = (1/mean(woodlouse3_tdft(:,1,counter)))*woodlouse3_tdft(:,:,counter);

end 


for counter = 1:size(beetle2_tdft,3)
    
    scaled_beetle2(:,:,counter) = (1/mean(beetle2_tdft(:,1,counter)))*beetle2_tdft(:,:,counter);

end 


for counter = 1:size(beetle3_tdft,3)
    
    scaled_beetle3(:,:,counter) = (1/mean(beetle3_tdft(:,1,counter)))*beetle3_tdft(:,:,counter);

end 


for counter = 1:size(varroa2_tdft,3)
    
    scaled_varroa2(:,:,counter) = (1/mean(varroa2_tdft(:,1,counter)))*varroa2_tdft(:,:,counter);

end 


for counter = 1:size(varroa3_tdft,3)
    
    scaled_varroa3(:,:,counter) = (1/mean(varroa3_tdft(:,1,counter)))*varroa3_tdft(:,:,counter);

end 

scaled_varroa = cat(3, scaled_varroa, scaled_varroa2, scaled_varroa3);
scaled_woodlouse = cat(3, scaled_woodlouse, scaled_woodlouse2, scaled_woodlouse3);
scaled_beetle = cat(3, scaled_beetle, scaled_beetle2, scaled_beetle3);

%build into one TDB

varroa_array = [];

for counter = 1:size(scaled_varroa,3)
temp = scaled_varroa(:,:,counter);
varroa_array(:,end+1) = temp(:);
end 

varroa_array = varroa_array';


woodlouse_array = [];

for counter = 1:size(scaled_woodlouse,3)
temp = scaled_woodlouse(:,:,counter);
woodlouse_array(:,end+1) = temp(:);
end 

woodlouse_array = woodlouse_array';


beetle_array = [];

for counter = 1:size(scaled_beetle,3)
temp = scaled_beetle(:,:,counter);
beetle_array(:,end+1) = temp(:);
end 

beetle_array = beetle_array';

TDB = [varroa_array; woodlouse_array; beetle_array];

figure(1)
imagesc((abs(TDB)), [max(((1/1000)*TDB(:))) max((TDB(:)))])
title('\fontsize{20} Full training database')
xlabel('\fontsize{20} Data points per signal')
ylabel('\fontsize{20} Signal number')
colormap jet
set(gca,'ColorScale','log')
colorbar
h = colorbar;
set(get(h,'label'),'string','\fontsize{20} Acceleration magnitude (m/s^2)');
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)


save('harriet_trainingdatabase_UPDATED','TDB','scaled_varroa','scaled_woodlouse','scaled_beetle')