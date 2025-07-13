clear all;
close all;
clc
%% -------------------------------------------------------------
%  File: run_SFOA_VMD.m
%  Version: 1.0
%  Description: Main script for SFOA-based parameter optimization and 
%               VMD-based seismic data denoising.
%
%  This code is part of the paper submitted to Computers & Geosciences:
%  "Seismic Denoising Using SFOA-Optimized Variational Mode Decomposition"
%
%  Author: Anonymous (for peer review)
%  License: MIT (see LICENSE file in repository)
%  Repository: https://github.com/anonymous-researcher/SFOA-VMD-denoising
%
%  Requirements:
%     - MATLAB R2018a or later
%     - Signal Processing Toolbox
%
%  Usage:
%     Simply run this file after placing your signal in example_data/
%
% -------------------------------------------------------------

%% 数据导入

data123=textread('F:\随钻罗江3\83077811\Z83077811_13_1s.tsegy', '' , 'headerlines', 2) ;

%% SFOA参数初始化

Npop=5;
Max_it=20;
lb=[1,1];
ub=[50,50000];

%% 利用SFOA选取最佳参数

[xposbest,fvalbest,Curve]=SFOA_VMD_old(data_cl,Npop,Max_it,lb,ub);

%% 利用最佳参数进行VMD分解

[imf_normal,residual,info] = vmd(data_cl,'NumIMF',round(xposbest(1)),'PenaltyFactor',round(xposbest(2))); 

%% 计算每个IMF相关系数

for i=1:round(xposbest(1))
    cor_data(i)=corr(x,imf_normal2(:,i));
end

%% 选取合适的IMF重构信号

result_SFOAVMD=sum(imf_normal(:,:),2);

%% 评价参数

% 模拟信号 信噪比 SNR

SNR_signal_moni=data;
SNR_denoise_moni=result_SFOAVMD;

SNR = 10 * log10(mean(SNR_signal_moni.^2) / mean(SNR_denoise_moni.^2));

% 实际信号 信噪比 SNR

SNR_signal_shiji=data;
SNR_denoise_shiji=result_SFOAVMD;

SNR_residual = SNR_signal_shiji - SNR_denoise_shiji;  % 估计残余噪声
snr_est = 10 * log10(mean(SNR_denoise_shiji.^2) / mean(SNR_residual.^2));

% 模拟信号 均方根误 RMSE

RMSE_signal_moni=data;
RMSE_denoise_moni=result_SFOAVMD;

RMSE_moni=sqrt(mean((RMSE_signal_moni - RMSE_denoise_moni).^2));

% 实际信号 均方根误差

RMSE_signal_shiji=data;
RMSE_denoise_shiji=result_SFOAVMD;

RMSE_residual=RMSE_signal_shiji-RMSE_denoise_shiji;
RMSE_moni=sqrt(mean((RMSE_signal_moni - RMSE_denoise_moni).^2));

% 模拟信号 标准RMSE RRMSE

RRMSE_signal_moni=data;
RRMSE_denoise_moni=result_SFOAVMD;

RRMSE_moni=norm(RRMSE_signal_moni-RRMSE_denoise_moni)/norm(RRMSE_signal_moni);

% 归一化互相关 NCC

NCC_signal_moni=data;
NCC_denoise_moni=result_SFOAVMD;

NCC_moni=sum((NCC_signal_moni-mean(NCC_signal_moni)).*(NCC_denoise_moni-mean(NCC_denoise_moni)))/(norm(NCC_signal_moni-mean(NCC_signal_moni))*norm(NCC_denoise_moni-mean(NCC_denoise_moni)));

