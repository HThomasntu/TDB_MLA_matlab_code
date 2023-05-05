% Copyright, M.Bencsik, H.Thomas, 2023

clear all 
close all 

% load IMPROVED_2DFT_DFA_outcome.mat
% load IMPROVED_masking_data.mat 

load discrimination_spectra.mat

file_name = '500 hz bee file 27.07.21 2.wav';

% BEE

% [A S_R] = audioread('G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\500 hz bee file 27.07.21 2.wav', [1 2]);
% [A S_R] = audioread(['G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\',file_name], [round((2*60)*S_R) round((((1*60)*60)+45*60)*S_R)]);
% A = single(A(:,1));
 
% MITE 

% [A S_R] = audioread('G:\Harriet Data\Varroa Vibration\petri_dish_video_analysis\varroa_brood_comb\sound_500.flac', [1 2]);
% [A S_R] = audioread(['G:\Harriet Data\Varroa Vibration\petri_dish_video_analysis\varroa_brood_comb\',file_name], [round((2*60)*S_R) round((34*60)*S_R)]);
% A = single(A(:,1));
 
% TEST VARROA4 VIDEO: 

% [A S_R] = audioread('G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\varroa 4 20.09.21_500hz.wav', [1 2]);
% [A S_R] = audioread(['G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\',file_name], [1 round((2*60)*S_R)]);

% TEST BEE VIDEOS (NO VARROA): 

% [A S_R] = audioread('E:\Additional brood videos\20.08.21 sample 1_500hz.wav', [1 2]);
% [A S_R] = audioread(['E:\Additional brood videos\',file_name], [round((2*60)*S_R) round((((1*60)*60)+40*60)*S_R)]);
% A = single(A(:,1));

% BROOD SAMPLE RECORDINGS: 

 [A S_R] = audioread('G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\500 hz bee file 27.07.21 2.wav', [1 2]);
 [A S_R] = audioread(['G:\Harriet Data\Varroa Vibration\Brood measurements\TDB\',file_name], [round((56*60)*S_R) round((((1*60)*60)+16*60)*S_R)]);
 A = single(A(:,1));


A = A*0.02651;

increment_shift = round(S_R*0.25);

df_scores_bee = [];
index_array = [];
threshold = 1.4;

for start_time = 1:increment_shift:length(A)-(1*S_R);
    audio_section = A(start_time:start_time+(1*S_R)-1,1);
    tdft = two_D_FT_Gaussian(audio_section(:,1),mf,tr,S_R,size(audio_section,1)/(2*S_R));
    for vv = 1:size(tdft,1)
        tdft_int(vv,:) = interp1(1:size(tdft,2),tdft(vv,:),1:5:size(tdft,2));
    end
    tdft_int = tdft_int(3:60,:);
%         imagesc(tdft)
%         pause(1)
    mean_audio = max(tdft_int(:,1));
    scaled_audio = tdft_int./mean_audio;
    
    ccp = sum(new_dfa(:).*scaled_audio(:));
    ccp2 = sum(new_dfa2(:).*scaled_audio(:));
    
    df_scores_bee(:,end+1) = [ccp; ccp2];   
    index_array(end+1) = start_time;
    
    100*start_time/(length(A)-(1*S_R))
    
end 

df_x = df_scores_bee(1,:);
df_y = df_scores_bee(2,:);

% Save the DF scores that have been computed for this particular file name:

save([file_name(1:end-4),'.mat'],'df_x','df_y','index_array')