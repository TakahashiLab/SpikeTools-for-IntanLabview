function [ALL,results_stats] = magnetStatsCSV2025b(species)
% 単一CSV 'speed_heading_data2025.csv' に
%  - sample 行（時系列サンプル）
%  - window 行（時間窓×魚の要約）
%  - intensity 行（時間窓平均の強度集計）
% を縦持ちで集約。対応（ペア）解析用の PairKey を付与。
%
% 依存: magnetBehav2.m（本ファイルと同一フォルダ想定）

% ------- 入力一覧 -------
% mat_files = {
%     '2025-05-28 10-01-55_vmc_checked.mat', 'beringSea',   1;
%     '2025-05-28 10-26-38_vmc_checked.mat', 'okhotskSea',  1;
%     '2025-05-28 10-47-04_vmc_checked.mat', 'controlSea1', 1;
%     '2025-05-28 11-05-32_vmc_checked.mat', 'beringSea',   2;
%     '2025-05-28 11-30-09_vmc_checked.mat', 'okhotskSea',  2;
%     '2025-05-28 11-49-48_vmc_checked.mat', 'controlSea1', 2;
%     '2025-05-28 12-10-14_vmc_checked.mat', 'beringSea',   3;
%     '2025-05-28 12-35-16_vmc_checked.mat', 'okhotskSea',  3;
%     '2025-05-28 12-54-39_vmc_checked.mat', 'controlSea1', 3;
%     '2025-05-28 13-14-15_vmc_checked.mat', 'beringSea',   4;
%     '2025-05-28 13-40-20_vmc_checked.mat', 'okhotskSea',  4;
%     '2025-05-28 14-00-01_vmc_checked.mat', 'controlSea1', 4;
%     '2025-05-28 14-20-03_vmc_checked.mat', 'controlT1',   5;
%     '2025-05-28 14-45-15_vmc_checked.mat', 'microT100',   5;
%     '2025-05-28 15-05-29_vmc_checked.mat', 'microT10',    5;
%     '2025-05-28 15-25-08_vmc_checked.mat', 'rotation90',  5;
%     '2025-05-28 15-43-24_vmc_checked.mat', 'controlT1',   6;
%     '2025-05-28 16-08-08_vmc_checked.mat', 'microT100',   6;
%     '2025-05-28 16-27-27_vmc_checked.mat', 'microT10',    6;
%     '2025-05-28 16-47-22_vmc_checked.mat', 'rotation90',  6;
%     '2025-05-28 17-06-11_vmc_checked.mat', 'controlT1',   7;
%     '2025-05-28 17-25-21_vmc_checked.mat', 'microT100',   7;
%     '2025-05-28 17-44-35_vmc_checked.mat', 'microT10',    7;
%     '2025-05-28 18-03-16_vmc_checked.mat', 'rotation90',  7
% };

%mat_files = {
% '2026-03-03 14-48-32DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 1;
% '2026-03-03 15-18-09DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 1;
% '2026-03-03 15-45-23DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 1;
% '2026-03-03 16-03-43DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 2;
% '2026-03-03 16-28-15DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 2;
% '2026-03-03 16-54-00DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 2;
% '2026-03-03 17-11-29DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 3;
% '2026-03-03 17-37-29DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 3;
% '2026-03-03 18-02-54DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 3;
% '2026-03-04 08-42-29DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 4;
% '2026-03-04 09-08-03DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 4;
% '2026-03-04 09-33-12DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 4;
% '2026-03-04 09-51-54DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 5;
% '2026-03-04 10-17-01DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 5;
% '2026-03-04 10-42-13DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 5;
% '2026-03-04 11-02-12DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 6;
% '2026-03-04 11-28-01DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 6;
% '2026-03-04 11-52-20DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 6;
% '2026-03-04 12-10-18DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT100', 7;
% '2026-03-04 12-34-16DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 7;
% '2026-03-04 12-54-00DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 7;
% '2026-03-04 13-13-14DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT100', 8;
% '2026-03-04 13-32-19DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 8;
% '2026-03-04 13-50-02DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 8;
% '2026-03-04 14-10-18DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 9;
% '2026-03-04 14-29-25DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 9;
% '2026-03-04 14-49-31DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat',  'microT100', 9;
% '2026-03-04 15-07-52DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 10;
% '2026-03-04 15-26-13DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 10;
% '2026-03-04 15-44-35DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat',  'microT100', 10;
% '2026-03-04 16-02-57DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 11;
% '2026-03-04 16-21-10DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 11;
% '2026-03-04 16-41-17DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat',  'microT100', 11;
% '2026-03-04 16-59-13DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat',  'microT100', 12;
% '2026-03-04 17-19-28DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 12;
% '2026-03-04 17-39-06DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 12;
%};

if nargin < 1 || isempty(species)
    species = "salmon";
end
species = lower(string(species));

%%trout
mat_files_trout = {
'2025-01-22 15-00-37DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'okhotskSea',  1,'5mo';
'2025-01-22 15-25-15DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'beringSea',   1,'5mo';
'2025-01-22 15-50-06DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlSea1', 1,'5mo';
'2025-01-22 16-15-14DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'okhotskSea',  2,'5mo';
'2025-01-22 16-40-20DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'beringSea',   2,'5mo';
'2025-01-22 17-05-12DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlSea1', 2,'5mo';
'2025-01-23 08-50-13DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT10',    3,'5mo';
'2025-01-23 09-15-17DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT100',   3,'5mo';
'2025-01-23 09-40-02DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlT1',   3,'5mo';
'2025-01-23 10-05-10DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT10',    4,'5mo';
'2025-01-23 10-30-24DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT100',   4,'5mo';
'2025-01-23 10-55-03DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlT1',   4,'5mo';
'2025-01-23 11-20-15DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'okhotskSea',  5,'5mo';
'2025-01-23 11-45-37DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'beringSea',   5,'5mo';
'2025-01-23 12-10-08DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlSea1', 5,'5mo';
'2025-01-23 13-04-49DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'okhotskSea',  6,'5mo';
'2025-01-23 13-29-43DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'beringSea',   6,'5mo';
'2025-01-23 13-55-01DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlSea1', 6,'5mo';
'2025-01-23 14-20-11DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'rotation0',  6,'5mo';
'2025-01-23 14-45-41DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'rotation180',  6,'5mo';
% '2025-03-13 13-00-45DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat',
'2025-03-13 13-25-10DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'beringSea',   7,'5mo';
'2025-03-13 13-50-10DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'okhotskSea',  7,'5mo';
'2025-03-13 14-15-06DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlSea1', 7,'5mo';
'2025-03-13 14-33-19DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'okhotskSea',  8,'5mo';
'2025-03-13 15-00-10DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'beringSea',   8,'5mo';
'2025-03-13 15-25-03DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlSea1', 8,'5mo';
'2025-03-13 15-45-02DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlT1',   9,'5mo';
'2025-03-13 16-10-09DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT100',   9,'5mo';
'2025-03-13 16-35-23DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT10',   9,'5mo';
'2025-03-13 16-55-02DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'rotation0',   9,'5mo';
'2025-03-13 17-00-15DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'rotation180',  9,'5mo';
'2025-03-14 08-49-36DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlT1',   10,'5mo';
'2025-03-14 09-15-21DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT100',   10,'5mo';
'2025-03-14 09-40-11DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT10',   10,'5mo';
'2025-03-14 10-05-36DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'rotation180',  10,'5mo';
'2025-03-14 10-25-02DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlT1',   11,'5mo';
'2025-03-14 10-50-10DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT100',   11,'5mo';
'2025-03-14 11-15-07DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT10',   11,'5mo';
'2025-03-14 11-40-11DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'rotation180',  11,'5mo';
'2025-03-14 12-29-50DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'controlT1',   12,'5mo';
'2025-03-14 12-55-38DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT100',   12,'5mo';
'2025-03-14 13-20-16DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'microT10',   12,'5mo';
'2025-03-14 13-44-41DLC_mobnet_100_trout_2025Oct17shuffle1_500000_vmc_checked_c.mat', 'rotation180',  12,'5mo';
};

%%salmon
mat_files_salmon = {
'2026-03-03 14-48-32DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 1, '2mo';
'2026-03-03 15-18-09DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 1, '2mo';
'2026-03-03 15-45-23DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 1, '2mo';
'2026-03-03 16-03-43DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 2, '2mo';
'2026-03-03 16-28-15DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 2, '2mo';
'2026-03-03 16-54-00DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 2, '2mo';
'2026-03-03 17-11-29DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 3, '2mo';
'2026-03-03 17-37-29DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 3, '2mo';
'2026-03-03 18-02-54DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 3, '2mo';
'2026-03-04 08-42-29DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 4, '2mo'
'2026-03-04 09-08-03DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 4, '2mo';
'2026-03-04 09-33-12DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 4, '2mo';
'2026-03-04 09-51-54DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 5, '2mo';
'2026-03-04 10-17-01DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 5, '2mo';
'2026-03-04 10-42-13DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 5, '2mo';
'2026-03-04 11-02-12DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'okhotskSea', 6, '2mo';
'2026-03-04 11-28-01DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlSea1', 6, '2mo';
'2026-03-04 11-52-20DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'beringSea', 6, '2mo';
'2026-03-04 12-10-18DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT100', 7, '2mo';
'2026-03-04 12-34-16DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 7, '2mo';
'2026-03-04 12-54-00DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 7, '2mo'
'2026-03-04 13-13-14DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT100', 8, '2mo';
'2026-03-04 13-32-19DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 8, '2mo';
'2026-03-04 13-50-02DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 8, '2mo';
'2026-03-04 14-10-18DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 9, '2mo';
'2026-03-04 14-29-25DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 9, '2mo';
'2026-03-04 14-49-31DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat',  'microT100', 9, '2mo';
'2026-03-04 15-07-52DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 10, '2mo';
'2026-03-04 15-26-13DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 10, '2mo'
'2026-03-04 15-44-35DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat',  'microT100', 10, '2mo';
'2026-03-04 16-02-57DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 11, '2mo';
'2026-03-04 16-21-10DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 11, '2mo';
'2026-03-04 16-41-17DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat',  'microT100', 11, '2mo';
'2026-03-04 16-59-13DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat',  'microT100', 12, '2mo';
'2026-03-04 17-19-28DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'microT10', 12, '2mo';
'2026-03-04 17-39-06DLC_mobnet_100_salmon2025-2Mar17shuffle1_500000_vmc.mat', 'controlT1', 12, '2mo';
'2025-05-28 10-01-55_vmc_checked.mat', 'beringSea',   13, '5mo';
'2025-05-28 10-26-38_vmc_checked.mat', 'okhotskSea',  13, '5mo';
'2025-05-28 10-47-04_vmc_checked.mat', 'controlSea1', 13, '5mo';
'2025-05-28 11-05-32_vmc_checked.mat', 'beringSea',   14, '5mo';
'2025-05-28 11-30-09_vmc_checked.mat', 'okhotskSea',  14, '5mo';
'2025-05-28 11-49-48_vmc_checked.mat', 'controlSea1', 14, '5mo';
'2025-05-28 12-10-14_vmc_checked.mat', 'beringSea',   15, '5mo';
'2025-05-28 12-35-16_vmc_checked.mat', 'okhotskSea',  15, '5mo';
'2025-05-28 12-54-39_vmc_checked.mat', 'controlSea1', 15, '5mo';
'2025-05-28 13-14-15_vmc_checked.mat', 'beringSea',   16, '5mo';
'2025-05-28 13-40-20_vmc_checked.mat', 'okhotskSea',  16, '5mo';
'2025-05-28 14-00-01_vmc_checked.mat', 'controlSea1', 16, '5mo';
'2025-05-28 14-20-03_vmc_checked.mat', 'controlT1',   17, '5mo';
'2025-05-28 14-45-15_vmc_checked.mat', 'microT100',   17, '5mo';
'2025-05-28 15-05-29_vmc_checked.mat', 'microT10',    17, '5mo';
'2025-05-28 15-25-08_vmc_checked.mat', 'rotation90',  18, '5mo';
'2025-05-28 15-43-24_vmc_checked.mat', 'controlT1',   18, '5mo';
'2025-05-28 16-08-08_vmc_checked.mat', 'microT100',   18, '5mo';
'2025-05-28 16-27-27_vmc_checked.mat', 'microT10',    19, '5mo';
'2025-05-28 16-47-22_vmc_checked.mat', 'rotation90',  19, '5mo';
'2025-05-28 17-06-11_vmc_checked.mat', 'controlT1',   19, '5mo';
'2025-05-28 17-25-21_vmc_checked.mat', 'microT100',   20, '5mo';
'2025-05-28 17-44-35_vmc_checked.mat', 'microT10',    20, '5mo';
'2025-05-28 18-03-16_vmc_checked.mat', 'rotation90',  20, '5mo';
};

switch char(species)
    case 'trout'
        mat_files = mat_files_trout;
    case 'salmon'
        mat_files = mat_files_salmon;
    otherwise
        error('species must be either ''salmon'' or ''trout''.');
end
% ------- パラメータ -------
fps = 30;
num_tracks = 9;
heading_jump_threshold_deg = 90;

% 時間窓（対応解析は窓単位で揃える）
time_windows = [0 5*60];   % 9分以降 / 20分以降
win_len_sec = inf;

% 条件→µT（必要に応じて調整）
defaultEarth = 50;
cond2uT = containers.Map( ...
    {'microT10','microT20','microT30','microT50','microT70','microT100', ...
     'controlT1','controlSea1','beringSea','okhotskSea','rotation90'}, ...
    [10, 20, 30, 50, 70, 100, ...
     defaultEarth, defaultEarth, defaultEarth, defaultEarth, defaultEarth] ...
);

% 条件→磁場北オフセット（deg）
cond2MagNorthOffset = containers.Map({'rotation90'}, {90});
defaultMagNorth = 0;

% ------- 出力（単一テーブル） -------
ALL = table();

for i = 1:size(mat_files,1)
    filename  = mat_files{i,1};
    condition = mat_files{i,2};
    cohort_id = mat_files{i,3};
    age_group = string(mat_files{i,4});
    order_in_session = sum([mat_files{1:i,3}] == cohort_id);

    % 磁場北オフセット
    if isKey(cond2MagNorthOffset, condition)
        magOff = cond2MagNorthOffset(condition);
    else
        magOff = defaultMagNorth;
    end

    % 強度（µT）
    if isKey(cond2uT, condition)
        inten_uT = cond2uT(condition);
    else
        inten_uT = defaultEarth;
    end

    % --- 解析（magnetBehav2） ---
    behav = magnetBehav2(filename, struct( ...
        'fps', fps, ...
        'magNorthDeg', magOff, ...
        'window_starts_sec', num2cell(time_windows), ...
        'win_len_sec', win_len_sec, ...
        'heading_jump_threshold_deg', heading_jump_threshold_deg ...
    ));

    speed = behav.speed;                % [9 x T]
    heading = behav.heading;            % [9 x T]
    hmag = behav.headingMag;            % [9 x T]
    turnR = behav.turnRate;             % [9 x (T-1)]
    curva = behav.curvature;            % [9 x (T-1)]
    T = size(speed,2);
    tvec = behav.meta.time(:);          % [T x 1]
    fish_global_ids = (cohort_id-1)*num_tracks + (1:num_tracks);

    % ========== 1) sample 行 ==========
    cosMag = cosd(hmag); sinMag = sind(hmag);
    rows_sample = table();
    for f = 1:num_tracks
        n = T;
        tmp = table( ...
            repmat("sample", n,1), ...
            repmat(string(filename), n,1), repmat(string(condition), n,1), ...
            repmat(cohort_id, n,1), repmat(age_group, n,1), repmat(fish_global_ids(f), n,1), ...
            repmat(f, n,1), repmat(order_in_session, n,1), tvec(1:n), ...
            speed(f,1:n)', heading(f,1:n)', hmag(f,1:n)', ...
            cosMag(f,1:n)', sinMag(f,1:n)', ...
            repmat(inten_uT, n,1), ...
            'VariableNames', {'Level','Session','Condition','Cohort','AgeGroup','Fish','Track','OrderInSession','Time', ...
                              'Speed','Heading','HeadingMag','cosMag','sinMag','Intensity_uT'});

        % turnRate / curvature は長さ T-1 → 2..end に入れる（長さ一致）
        tr = nan(n,1); cv = nan(n,1);
        if ~isempty(turnR)
            tr(2:end) = turnR(f,:)';
        end
        if ~isempty(curva)
            cv(2:end) = curva(f,:)';
        end
        tmp.TurnRate  = tr;         % deg/s
        tmp.Curvature = cv;         % 1/unit

        % window 要約の列は sample では NaN
        tmp.WinStartSec   = nan(n,1);
        tmp.Tortuosity    = nan(n,1);
        tmp.Straightness  = nan(n,1);
        tmp.meanAngle_deg = nan(n,1);
        tmp.R             = nan(n,1);
        tmp.RayleighZ     = nan(n,1);
        tmp.RayleighP     = nan(n,1);

        % 対応キー（n×1 で代入）
        keyFish = strcat("C", string(cohort_id), "_F", string(fish_global_ids(f)));
        tmp.PairKey_Fish    = repmat(keyFish, n, 1);
        tmp.PairKey_FishWin = strings(n,1); % sample は空

        rows_sample = [rows_sample; tmp]; %#ok<AGROW>
    end

    % ========== 2) window 行 ==========
    W = behav.winStats; % 変数: meanAngle_deg, R, RayleighZ, RayleighP, MetaWinIndex, Fish
    tor = behav.tortuosityWin(:);
    sti = behav.straightnessWin(:);

    % WinStartSec を付与（MetaWinIndex 1..N → time_windows）
    winStartCol = time_windows(W.MetaWinIndex)';
 

    PairKey_Fish    = strcat("C", string(cohort_id), "_F", string((cohort_id-1)*num_tracks + W.Fish));
    PairKey_FishWin = strcat(PairKey_Fish, "_W", string(winStartCol));
       % winStartCol = repmat(winStartCol, size(PairKey_Fish));
    % PairKey_FishWin = PairKey_Fish + "_W" + string(winStartCol);

    rows_window = table( ...
        repmat("window", height(W),1), ...
        repmat(string(filename), height(W),1), ...
        repmat(string(condition), height(W),1), ...
        repmat(cohort_id, height(W),1), ...
        repmat(age_group, height(W),1), ...
        (cohort_id-1)*num_tracks + W.Fish, ...
        W.Fish, ...
        repmat(order_in_session, height(W),1), ...
        repmat(nan(height(W),1),1), ...   % Time
        repmat(nan(height(W),1),1), ...   % Speed
        repmat(nan(height(W),1),1), ...   % Heading
        repmat(nan(height(W),1),1), ...   % HeadingMag
        repmat(nan(height(W),1),1), ...   % cosMag
        repmat(nan(height(W),1),1), ...   % sinMag
        repmat(inten_uT, height(W),1), ...
        repmat(nan(height(W),1),1), ...   % TurnRate
        repmat(nan(height(W),1),1), ...   % Curvature
        winStartCol, ...
        tor, sti, ...
        W.meanAngle_deg, W.R, W.RayleighZ, W.RayleighP, ...
        PairKey_Fish, PairKey_FishWin, ...
        'VariableNames', {'Level','Session','Condition','Cohort','AgeGroup','Fish','Track','OrderInSession','Time', ...
                          'Speed','Heading','HeadingMag','cosMag','sinMag','Intensity_uT', ...
                          'TurnRate','Curvature', ...
                          'WinStartSec','Tortuosity','Straightness', ...
                          'meanAngle_deg','R','RayleighZ','RayleighP', ...
                          'PairKey_Fish','PairKey_FishWin'});

    % ========== 3) intensity 行 ==========
    rows_intensity = table();
    for w = 1:numel(time_windows)
        % --- ウィンドウ選択 ---
        idxw = (rows_window.WinStartSec == time_windows(w));

        % --- 安全処理（NaNやサイズ不一致に対応）---
        if ~islogical(idxw)
            % 数値インデックスだった場合
            tmp = false(height(rows_window),1);
            validIdx = idxw(idxw > 0 & idxw <= height(rows_window)); % 範囲内のみ
            tmp(validIdx) = true;
            idxw = tmp;
        elseif numel(idxw) ~= height(rows_window)
            % 論理配列だけどサイズが一致しない場合
            tmp = false(height(rows_window),1);
            ncopy = min(height(rows_window), numel(idxw));
            tmp(1:ncopy) = idxw(1:ncopy);
            idxw = tmp;
        end

        % --- 実際の計算（idxwが空ならNaNを返す）---
        if any(idxw)
            muStraight = mean(rows_window.Straightness(idxw), 'omitnan');
            muTor      = mean(rows_window.Tortuosity(idxw),  'omitnan');
            muR        = mean(rows_window.R(idxw),           'omitnan');
        else
            muStraight = NaN;
            muTor = NaN;
            muR = NaN;
        end

        tmp = table( ...
            "intensity", string(filename), string(condition), cohort_id, age_group, NaN, NaN, ...
            order_in_session, NaN, NaN, NaN, NaN, NaN, NaN, ...
            inten_uT, NaN, NaN, ...
            double(time_windows(w)), muTor, muStraight, ...
            NaN, muR, NaN, NaN, ...
            strcat("C", string(cohort_id), "_FALL"), ...
            strcat("C", string(cohort_id), "_W", string(time_windows(w))), ...
            'VariableNames', {'Level','Session','Condition','Cohort','AgeGroup','Fish','Track','OrderInSession','Time', ...
                              'Speed','Heading','HeadingMag','cosMag','sinMag','Intensity_uT', ...
                              'TurnRate','Curvature','WinStartSec','Tortuosity','Straightness', ...
                              'meanAngle_deg','R','RayleighZ','RayleighP', ...
                              'PairKey_Fish','PairKey_FishWin'});
        rows_intensity = [rows_intensity; tmp]; %#ok<AGROW>
    end

    % ========== 縦積み ==========
    ALL = [ALL; rows_sample; rows_window; rows_intensity]; %#ok<AGROW>

    % 結果構造体（必要最小限）
    if ~exist('results','var') || ~isfield(results, condition)
        results.(condition) = struct(); results.(condition).sessions = [];
    end
    session = struct();
    session.file = filename;
    session.cohort_id = cohort_id;
    session.age_group = age_group;
    session.order_in_session = order_in_session;
    session.fish_global_ids = fish_global_ids;
    session.meta = behav.meta;
    results.(condition).sessions = [results.(condition).sessions; session]; %#ok<AGROW>
end

% ------- 型の整理 -------
ALL.Level           = string(ALL.Level);
ALL.Session         = string(ALL.Session);
ALL.Condition       = string(ALL.Condition);
ALL.AgeGroup        = string(ALL.AgeGroup);
ALL.PairKey_Fish    = string(ALL.PairKey_Fish);
ALL.PairKey_FishWin = string(ALL.PairKey_FishWin);

% ------- 出力 -------
output_csv = sprintf('speed_heading_data2025_%s-all.csv', char(species));
writetable(ALL, output_csv);

if exist('results','var')
    results_stats = results;
else
    results_stats = struct();
end
end
