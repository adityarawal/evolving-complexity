clear all; clc; close all;
rng('shuffle');
len = 300
% X = 0 + (1-0).*randn(len,1);
% Y = X+rand(len,1);
% Y = Y-mean(Y);
X = randn(len,1); %Normal distribution
Y = X+rand(len,1);
%Y = rand(len,1);
mean_X = mean(X)
var_X = var(X)
mean_Y = mean(Y)
var_Y = var(Y)

%Make sure the input is normalized to 0-1 (also preferrably zero mean
%and unit variance
max_X = max(X)
min_X = min(X)
max_Y = max(Y)
min_Y = min(Y)
X= (X-min(X))/(max(X)-min(X));
Y= (Y-min(Y))/(max(Y)-min(Y));
max_X_norm = max(X)
min_X_norm = min(X)
max_Y_norm = max(Y)
min_Y_norm = min(Y)

mean_X_norm = mean(X)
var_X_norm = var(X)
mean_Y_norm = mean(Y)
var_Y_norm = var(Y)

if (max_X_norm ~= 1 || max_Y_norm ~=1 || min_X_norm ~=0 || min_Y_norm ~=0)
    display('ERRROR: Values of X or Y outside 0-1 range');
end

% dlmwrite('del_Y.txt', Y, '\n'); %Write to file to be read by python/C++
% dlmwrite('del_X.txt', X, '\n'); %Write to file to be read by python/C++

c = corr([X Y])
I_exact = -0.5*log(1-c(1,2)^2)

k = 3;
[ I1, I2, points_knn dist_knn nx1 ny1] = KraskovMI( X, Y, k);
I1
[ I1_fast, I2_fast, fast_points_knn fast_dist_knn fast_nx1 fast_ny1] = fastKraskovMI( X, Y, k);
I1_fast

if(I1~=I1_fast)
    display('Error: Discrepancy in I1 and I1_fast');
else
    display('Correct; I1 and I1_fast match');
end

[X_entropy] = entropy_nearest_neighbor(X);
%[X Y points_knn fast_points_knn]
%[X Y dist_knn fast_dist_knn]

% %To avoid stupid rounding errors
% tol = 1.e-6
% error = sum(sum(abs(points_knn-fast_points_knn)));
% if (error > tol)
%     display('Error: Discrepancy in nearest neighbors');
% else
%     display('Correct: Nearest neighbors match'); 
% end
% %[nx1 ny1 fast_nx1 fast_ny1]
% if(sum(sum([(nx1-fast_nx1) (ny1-fast_ny1)])) > 0)
%     display('Error: Discrepancy in nx1, ny1');
% else
%     display('Correct: Nearest points match');
% end
% t=0;
% 
% f1 = @() KraskovMI( X, Y, k);
% timeit(f1)
% 
% f2 = @() fastKraskovMI( X, Y, k);
% timeit(f2)
% 
t=0;
