% Copyright, M.Bencsik, H.Thomas, 2023

clear all 
load IMPROVED_2DFT_DFA_outcome.mat

mite_X = [];
bee_X = [];
bg_X = [];

for pulse = 1:size(A_x,1)
    
    if pulse < 67
        
        mite_X(end+1,:) = A_x(pulse,:);
        
    elseif pulse > 308
        
        bg_X(end+1,:) = A_x(pulse,:);
        
    else bee_X(end+1,:) = A_x(pulse,:);
        
    end 
end 

mite_Y = [];
bee_Y = [];
bg_Y = [];

for pulse2 = 1:size(A_y,1)
    
    if pulse2 < 67
        
        mite_Y(end+1,:) = A_y(pulse2,:);
        
    elseif pulse2 > 308
        
        bg_Y(end+1,:) = A_y(pulse2,:);
        
    else bee_Y(end+1,:) = A_y(pulse2,:);
        
    end 
end 

figure(1)
plot(mite_X, mite_Y,'ro')
hold on 
plot(bee_X, bee_Y,'ko')
hold on
plot(bg_X, bg_Y,'bo')

save('polygonal_DF_areas.mat','mite_X','bee_X','bg_X','mite_Y','bee_Y','bg_Y')
%save IMPROVED_masking_data.mat