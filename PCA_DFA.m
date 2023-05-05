% Copyright, M.Bencsik, H.Thomas, 2023

clear all
close all
path(path,'C:\Users\PHY3HALLH\OneDrive - Nottingham Trent University\Desktop\Harriet_matlab_code\useful_matlab_function')

tr = 0.015;
mf = 4;


load IMPROVED_tdft_TDB_interpolated_biggerTDB.mat

% convert the scaled and cropped matrixes into arrays for the TDB: 

bee_array = [];

for counter = 1:size(scaled_bee_final,3)
temp = scaled_bee_final(:,:,counter);
bee_array(:,end+1) = temp(:);
end 

bee_array = bee_array';


mite_array = []; 

for counter = 1:size(scaled_mite_final,3)
temp = scaled_mite_final(:,:,counter);
mite_array(:,end+1) = temp(:);
end 

mite_array = mite_array';


bg_array = []; 

for counter = 1:size(scaled_bg,3)
temp = scaled_bg(:,:,counter);
bg_array(:,end+1) = temp(:);
end 

bg_array = bg_array';


% mite_bg_array = []; 
% 
% for counter = 1:size(scaled_mite_bg,3)
% temp = scaled_mite_bg(:,:,counter);
% mite_bg_array(:,end+1) = temp(:);
% end 
% 
% mite_bg_array = mite_bg_array';
% 
% 
% bee_bg_array = []; 
% 
% for counter = 1:size(scaled_bee_bg,3)
% temp = scaled_bee_bg(:,:,counter);
% bee_bg_array(:,end+1) = temp(:);
% end 
% 
% bee_bg_array = bee_bg_array';
% 
% 
% TDB = [mite_array; bee_array; mite_bg_array; bee_bg_array];

TDB = [mite_array; bee_array; bg_array];


figure(1)
imagesc(log10(abs(TDB)), [max(log10((1/300)*TDB(:))) max(log10(TDB(:)))])
title('\fontsize{20} Training database for bee, mite & background discrimination')
xlabel('\fontsize{20} Data points per pulse')
ylabel('\fontsize{20} Pulse number (Mite(1:28), Bee(29:57), Background(57:end))')
colormap jet
colorbar
h = colorbar;
set(get(h,'label'),'string','\fontsize{20} Acceleration magnitude (m/s^2)');
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)

mite_wavs = size(scaled_mite_final,3);
bee_wavs = size(scaled_bee_final,3);
bg_wavs = size(scaled_bg,3);
% mite_bg_wavs = size(scaled_mite_bg,3);
% bee_bg_wavs = size(scaled_bee_bg,3);


% undertake PCA of the collection of spectrograms:
 
%TDB = [matrix_A(:,indices) matrix_B]';
% Centre the data set:
temp2 = mean(TDB);
centred_data_set2 = (TDB - ones(size(TDB,1),1)*temp2)';

% Calculate the PCA scores and eigenvectors:
L = centred_data_set2*centred_data_set2'; % L is the covariance matrix C=A*A'.
[V D] = eig(L); % Diagonal elements of D are the eigenvalues for both L=A'*A and C=A*A'.

scores = flipud(V'*centred_data_set2);
Eigenspectra = fliplr(V);

figure(2)
plot(scores(1,1:mite_wavs),scores(2,1:mite_wavs),'ro')
hold on
plot(scores(1,(mite_wavs+1):(mite_wavs+bee_wavs)),scores(2,(mite_wavs+1):(mite_wavs+bee_wavs)),'ko')
hold on
plot(scores(1,(mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)), scores(2,(mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)),'bo')
% plot(scores(1,(mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)), scores(2,(mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)),'bo')
% hold on
% plot(scores(1,(mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)), scores(2,(mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)),'co')


% Calculate the DFA scores and DF scores:

limit = 6;

% [U,V, eigenval] = dfa([scores(1:limit,:)]',[zeros(1,mite_wavs) ones(1,bee_wavs) 2*ones(1,mite_bg_wavs) 3*ones(1,bee_bg_wavs)],2);

[U,V, eigenval] = dfa([scores(1:limit,:)]',[zeros(1,mite_wavs) ones(1,bee_wavs) 2*ones(1,bg_wavs)],2);

% Calculate the DF scores on the entire collection of swarming and non
% swarming spectra:
DFA_spectrum_01 = sum(Eigenspectra(:,1:limit)'.*(V(:,1)*ones(1,size(Eigenspectra,2))));
DFA_spectrum_02 = sum(Eigenspectra(:,1:limit)'.*(V(:,2)*ones(1,size(Eigenspectra,2))));
%DFA_spectrum_03 = sum(Eigenspectra(:,1:limit)'.*(V(:,3)*ones(1,size(Eigenspectra,2))));
 validation_data = TDB;
% plot(DFA_spectrum_01,'r')
% hold on
% plot(DFA_spectrum_02,'k')
% figure
A_x = sum((validation_data).*(ones(size(validation_data,1),1)*DFA_spectrum_01),2);
A_y = sum((validation_data).*(ones(size(validation_data,1),1)*DFA_spectrum_02),2);
%A_z = sum((validation_data).*(ones(size(validation_data,1),1)*DFA_spectrum_03),2);

figure(3)
subplot(1,2,1)
plot(A_x(1:mite_wavs),A_y(1:mite_wavs),'ro')
% Calculate the centroid of the varroa:
centroid_01 = [mean(A_x(1:mite_wavs)) mean(A_y(1:mite_wavs))];
hold on
plot(A_x((mite_wavs+1):(mite_wavs+bee_wavs)),A_y((mite_wavs+1):(mite_wavs+bee_wavs)),'ko')
centroid_02 = [mean(A_x((mite_wavs+1):(mite_wavs+bee_wavs))) mean(A_y((mite_wavs+1):(mite_wavs+bee_wavs)))];
hold on
plot(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)),A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)),'bo')
centroid_03 = [mean(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs))) mean(A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)))];
% plot(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)),A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)),'bo')
% centroid_03 = [mean(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs))) mean(A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)))];
% hold on
% plot(A_x((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)),A_y((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)),'co')
% centroid_04 = [mean(A_x((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs))) mean(A_y((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)))];
xlabel('\fontsize{20} DF score 1')
ylabel('\fontsize{20} DF score 2')
title('\fontsize{20} DFA  outcome (scatterplot)')
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)
legend('walking','bee','background', 'Location', 'northwest')

subplot(1,2,2)
new_dfa = reshape(DFA_spectrum_01, size(scaled_bee,1), size(scaled_bee,2));
new_dfa2 = reshape(DFA_spectrum_02, size(scaled_bee,1), size(scaled_bee,2));
imagesc([0 0.5*mf/tr], [0 0.5*S_R], (new_dfa))
colormap jet
colorbar
h = colorbar;
set(get(h,'label'),'string','Acceleration magnitude (m/s^2)');
title('\fontsize{20} DFA 1 outcome (spectrum)')
xlabel('\fontsize{20} Spectral Repetition (Hz)')
ylabel('\fontsize{20} Frequency (Hz)')
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)

figure(4)
subplot(1,2,1)
plot(A_x(1:mite_wavs),A_y(1:mite_wavs),'ro')
% Calculate the centroid of the varroa:
centroid_01 = [mean(A_x(1:mite_wavs)) mean(A_y(1:mite_wavs))];
hold on
plot(A_x((mite_wavs+1):(mite_wavs+bee_wavs)),A_y((mite_wavs+1):(mite_wavs+bee_wavs)),'ko')
centroid_02 = [mean(A_x((mite_wavs+1):(mite_wavs+bee_wavs))) mean(A_y((mite_wavs+1):(mite_wavs+bee_wavs)))];
hold on
plot(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)),A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)),'bo')
centroid_03 = [mean(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs))) mean(A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)))];
% plot(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)),A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)),'bo')
% centroid_03 = [mean(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs))) mean(A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)))];
% hold on
% plot(A_x((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)),A_y((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)),'co')
% centroid_04 = [mean(A_x((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs))) mean(A_y((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)))];
xlabel('\fontsize{20} DF score 1')
ylabel('\fontsize{20} DF score 2')
title('\fontsize{20} DFA outcome (scatterplot)')
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)
legend('walking','bee', 'background', 'Location', 'northwest')

subplot(1,2,2)
imagesc([0 0.5*mf/tr], [0 0.5*S_R], (new_dfa2))
colormap jet
colorbar
h = colorbar;
set(get(h,'label'),'string','Acceleration magnitude (m/s^2)');
title('\fontsize{20} DFA 2 outcome (spectrum)')
xlabel('\fontsize{20} Spectral Repetition (Hz)')
ylabel('\fontsize{20} Frequency (Hz)')
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)


%save IMPROVED_2DFT_DFA_outcome.mat
save('discrimination_spectra.mat','new_dfa','new_dfa2','mf','tr','centroid_01','centroid_02','centroid_03')


% FOR THESIS AND PAPER: 

%scatterplot: 

figure(5)
plot(A_x(1:mite_wavs),A_y(1:mite_wavs),'ko','MarkerFaceColor','k')
% Calculate the centroid of the varroa:
centroid_01 = [mean(A_x(1:mite_wavs)) mean(A_y(1:mite_wavs))];
hold on
plot(A_x((mite_wavs+1):(mite_wavs+bee_wavs)),A_y((mite_wavs+1):(mite_wavs+bee_wavs)),'bo','MarkerFaceColor','b')
centroid_02 = [mean(A_x((mite_wavs+1):(mite_wavs+bee_wavs))) mean(A_y((mite_wavs+1):(mite_wavs+bee_wavs)))];
hold on
plot(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)),A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)),'co','MarkerFaceColor','c')
centroid_03 = [mean(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs))) mean(A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)))];
% plot(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)),A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)),'bo')
% centroid_03 = [mean(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs))) mean(A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)))];
% hold on
% plot(A_x((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)),A_y((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)),'co')
% centroid_04 = [mean(A_x((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs))) mean(A_y((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)))];
xlabel('\fontsize{20} DF 1')
ylabel('\fontsize{20} DF 2')
title('\fontsize{20} DFA scatterplot')
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)
legend('mite','bee', 'background', 'Location', 'northwest')

% scree plot: 
figure(6)
subplot(2,1,1)
plot(mean(abs(scores)'))
xlabel('\fontsize{22} Principal component number')
ylabel('\fontsize{22} Average modulus PC score')
title('\fontsize{22} a) Scree plot for PCA of mite, bee and background vibrations')
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)
xline(7, 'r','linewidth',1)

subplot(2,1,2)
plot(mean(abs(scores)'))
axis([0 400 0 1])
xlabel('\fontsize{22} Principal component number')
ylabel('\fontsize{22} Average modulus PC score')
title('\fontsize{22} b) Scree plot for PCA of mite, bee and background vibrations (zoomed)')
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)
xline(7, 'r','linewidth',1)

decay_curve = mean(abs(scores)');
PCA_deviations(decay_curve,6)

%df spectrum plot:
figure(7)
subplot(1,2,1)
%imagesc([0 0.5*mf/tr], [0 0.5*S_R], (new_dfa), [-0.2 0.2])
imagesc([0 0.5*mf/tr], [0 4000], (new_dfa), [-0.28 0.28])
colormap jet
colorbar
h = colorbar;
set(get(h,'label'),'string','Acceleration magnitude (a.u.)');
title('\fontsize{20} a) DF spectrum 1')
xlabel('\fontsize{20} Spectral Repetition (Hz)')
ylabel('\fontsize{20} Frequency (Hz)')
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)

subplot(1,2,2)
%imagesc([0 0.5*mf/tr], [0 0.5*S_R], (new_dfa2), [-0.2 0.2])
imagesc([0 0.5*mf/tr], [0 4000], (new_dfa2), [-0.28 0.28])
colormap jet
colorbar
h = colorbar;
set(get(h,'label'),'string','Acceleration magnitude (a.u.)');
title('\fontsize{20} b) DF spectrum 2')
xlabel('\fontsize{20} Spectral Repetition (Hz)')
ylabel('\fontsize{20} Frequency (Hz)')
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)


figure(8)
figure('position', [1 1 1920 1080])
subplot(2,1,1)
plot(A_x(1:mite_wavs),A_y(1:mite_wavs),'ko','MarkerFaceColor','k')
% Calculate the centroid of the varroa:
centroid_01 = [mean(A_x(1:mite_wavs)) mean(A_y(1:mite_wavs))];
hold on
plot(A_x((mite_wavs+1):(mite_wavs+bee_wavs)),A_y((mite_wavs+1):(mite_wavs+bee_wavs)),'bo','MarkerFaceColor','b')
centroid_02 = [mean(A_x((mite_wavs+1):(mite_wavs+bee_wavs))) mean(A_y((mite_wavs+1):(mite_wavs+bee_wavs)))];
hold on
plot(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)),A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)),'co','MarkerFaceColor','c')
centroid_03 = [mean(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs))) mean(A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+bg_wavs)))];
hold on
plot((centroid_01(:,1)), (centroid_01(:,2)), 'ro', 'markersize', 20, 'Linewidth', 3)
hold on
plot((centroid_02(:,1)), (centroid_02(:,2)), 'go', 'markersize', 20, 'Linewidth', 3)
hold on
plot((centroid_03(:,1)), (centroid_03(:,2)), 'mo', 'markersize', 20, 'Linewidth', 3)
% plot(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)),A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)),'bo')
% centroid_03 = [mean(A_x((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs))) mean(A_y((mite_wavs+bee_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs)))];
% hold on
% plot(A_x((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)),A_y((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)),'co')
% centroid_04 = [mean(A_x((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs))) mean(A_y((mite_wavs+bee_wavs+mite_bg_wavs+1):(mite_wavs+bee_wavs+mite_bg_wavs+bee_bg_wavs)))];
xlabel('\fontsize{20} DF score 1')
ylabel('\fontsize{20} DF score 2')
title('\fontsize{20} a) DFA scatterplot - TDB')
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)
legend('mite','bee', 'background', 'Location', 'northeast')
grid on
axis equal

subplot(2,2,3)
%imagesc([0 0.5*mf/tr], [0 0.5*S_R], (new_dfa), [-0.2 0.2])
imagesc([0 0.5*mf/tr], [0 4000], (new_dfa), [-0.28 0.28])
colormap jet
colorbar
h = colorbar;
set(get(h,'label'),'string','Acceleration magnitude (a.u.)');
title('\fontsize{20} b) Discriminant 2DFT 1')
xlabel('\fontsize{20} Spectral Repetition (Hz)')
ylabel('\fontsize{20} Frequency (Hz)')
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)

subplot(2,2,4)
%imagesc([0 0.5*mf/tr], [0 0.5*S_R], (new_dfa2), [-0.2 0.2])
imagesc([0 0.5*mf/tr], [0 4000], (new_dfa2), [-0.28 0.28])
colormap jet
colorbar
h = colorbar;
set(get(h,'label'),'string','Acceleration magnitude (a.u.)');
title('\fontsize{20} c) Discriminant 2DFT 2')
xlabel('\fontsize{20} Spectral Repetition (Hz)')
ylabel('\fontsize{20} Frequency (Hz)')
set(gcf,'color','w');
a = get(gca,'TickLabel');  
set(gca,'TickLabel',a,'fontsize',20)

print('Figure 5.tif', '-dtiff', '-r400')