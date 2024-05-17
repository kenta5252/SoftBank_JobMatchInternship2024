%% 通過するSセルの全滞在時間のデータをとるコード
%% ユーザの移動する線分は深さが一様分布を採用

%% Sセル600

clc;
clear;
tic;


Simulation_kaisu = 1;
S_number = 600;


%% パラメータ値設定

In=0; % 円の中の基地局数(乱数)
Out=0; % 円の外の基地局数(乱数)
R=1; % Mセル半径[km]
r=0.04; %Sセル半径[km]
tani_ugoki=0.0001; % Userの動きの細かさ
User_v=50; % UEの速度[km/h]

Sojourn_Distance_S = zeros; % Sセル滞在時間を確保しておく配列

maxDistance_M = 2*R; % Sセルの最長通過距離
maxTime_M = maxDistance_M / tani_ugoki; % 動きの細かさに対するMセルを通過する最大の時間(回数)

% Sセル内を動く度に計算する変数
% UEを動かすfor文の都合で+2をしている
Cos_a=zeros(In,maxTime_M+2); 
Sin_a=zeros(In,maxTime_M+2);
Cos_b=zeros(In,maxTime_M+2); % b=sita-a
Sin_b=zeros(In,maxTime_M+2);
expect_stay_d=zeros(In,maxTime_M+2); % 予想されるセルの滞在距離

expect_stay_t=0; % UEがSセルに入った時点での予測滞在時間

SojournNum = 1; % 滞在時間を確保しておく配列の順序番号

sum_points = S_number*2;
lam = 1:1:sum_points;

old_x = poissrnd(lam); % 0から200の間でポアソンに従った乱数の300個のデータ
old_y = poissrnd(lam);

x_y = zeros(sum_points, 2);

for number = 1:sum_points
    x_y(number,1) = old_x(1, number) - S_number;
    x_y(number,2) = old_y(1, number) - S_number;
    x_y(number,1) = x_y(number,1) / S_number; % x軸:ー１から+1の間でポアソンに従った乱数が発生
    x_y(number,2) = x_y(number,2) / S_number; %
end


[m, n] = size(x_y);

% 列ごとにデータをシャッフル
shuffled_x_y = zeros(m, n);
for i = 1:n
    idx = randperm(m);  % インデックスをランダムに並び替える
    shuffled_x_y(:, i) = x_y(idx, i);  % データを並び替えたインデックスで更新
end


% 円の中の乱数を探す
for N_in=1:sum_points
    if ((shuffled_x_y(N_in,1)/R)^2+(shuffled_x_y(N_in,2)/R)^2)^(0.5)<=1
        ransu_naka(In+1,1)=shuffled_x_y(N_in,1);
        ransu_naka(In+1,2)=shuffled_x_y(N_in,2);
        In=In+1;
    end
    if In == S_number
        break
    end
end



%% シミュレーション開始(Mセルの角度が一様分布)

for x=1:Simulation_kaisu % シミュレーション回数は何回するか？

    % 半径Rの円から無作為に半径を選ぶ
    selected_radius = R * rand();

    %% 原点 O
    %% 無作為に選んだ点 A
    %% 始点 U

    % 選ばれた半径上の点 A を無作為に選ぶ（極座標形式）
    theta = 2 * pi * rand();  % 0から2πまでのランダムな角度
    OA = [selected_radius * cos(theta), selected_radius * sin(theta)];
    
    % 点 A（任意の点）から原点 O へのベクトル: AO = [-xOA, -yOA]
    AO = -OA;

    % ベクトルAOと直交するベクトル AU を求めるために、ベクトル AO を90度時計回り（ー９０）に回転させる
    AU_vec = [AO(2), -AO(1)];

    % AU の長さを三平方の定理より導出
    AU_length = sqrt(R^2 - norm(OA)^2);
    %fprintf('AUの長さ: (%f)\n', AU_length);

    % ベクトルAUの大きさを本来のAUの長さにするにする
    AU = AU_length * AU_vec / norm(AU_vec);
    %fprintf('AUのノルム: (%f)\n', norm(AU));
    %fprintf('AU: (%f, %f)\n', AU);

    % OU(Uの座標)が知りたい．
    % AU = OA - OU
    % OU = OA - AU
    OU = OA - AU;
    %fprintf('点Uの座標: (%f, %f)\n', OU);
    %fprintf('OUのノルム: (%f)\n', norm(OU));

    % 始点U終点Vの進行角度(OAの直線とx軸に対する角度)
    % 2点を結ぶ直線がx軸となす角度を計算(2点：OUとOA)
    Angle_of_Travel = atan2(OA(2) - OU(2), OA(1) - OU(1));

    % 終点Vの位置をベクトルで導出してもいいけど，計算省略のためにUAの距離の2倍がMセル滞在距離とする．
    Sojourn_Distance_M = 2 * norm(AU); % U-V間距離
    

    kaunto = zeros(In,1);
    x_now(1)= OU(1) ; % ループのために初期状態を配列にいれる
    y_now(1)= OU(2) ; % ループのために初期状態を配列にいれる


    %% UEが移動する
         
    for tt = 2 : 1: maxTime_M+2 % tt回単位当たりで動く

        Sojourn_Distance_S; % 予想される滞在距離[km]

        x_now(tt) = x_now(tt-1) + tani_ugoki * cos(Angle_of_Travel);
        y_now(tt) = y_now(tt-1) + tani_ugoki * sin(Angle_of_Travel);


        if tt > (Sojourn_Distance_M / tani_ugoki) % このifでnnのforループを終わらせるように
            break;
        end

        for n=1:In % スモールセルの数In個をnで全てみる 最短距離のセルと接続
            
            Length_nextScell_user = ((ransu_naka(n,1)-x_now(tt))^2 + (ransu_naka(n,2)-y_now(tt))^2 )^(1/2);

            if kaunto(n,1) == 0
                    
                if Length_nextScell_user <= r
                    Cos_a(n,tt) = (ransu_naka(n,1)-x_now(tt)) / Length_nextScell_user; 
                    Sin_a(n,tt) = (ransu_naka(n,2)-y_now(tt)) / Length_nextScell_user;
                    Cos_b(n,tt) = Cos_a(n,tt) * cos(Angle_of_Travel) + Sin_a(n,tt) * sin(Angle_of_Travel); % 加法定理で接続先候補のセルの滞在距離を導出
                    Sin_b(n,tt) = Sin_a(n,tt) * cos(Angle_of_Travel) - Cos_a(n,tt) * sin(Angle_of_Travel);
                      
                    Sojourn_Distance_S(1, SojournNum) = Length_nextScell_user * Cos_b(n,tt) + (r^2 - (Length_nextScell_user * Sin_b(n,tt))^2) ^ (1/2); % 予想される滞在距離[km]
                    
                    SojournNum = SojournNum + 1;

                    kaunto(n,1) = 1;
                end
                    
            end

        end

    end
    
end % シミュレーション回数分回す

toc;




%% ファイル作成

filename = 'SmallCell600_SojournTime_StatisticalData.txt';

fileID = fopen(filename,'w');

fprintf(fileID, 'Mセル移動経路\n');
fprintf(fileID, '-- Mセルの中心からの深さを一様分布させてMセルを横断．\n\n');

fprintf(fileID, '解析データ　条件詳細\n');
fprintf(fileID, '-- Sセル数：%d\n', S_number);
fprintf(fileID, '-- シミュレーション回数：%d\n', Simulation_kaisu);

[rows, cols] = size(Sojourn_Distance_S(1, :));
fprintf(fileID, '-- 合計Sセル通過数（データ数）：%d\n', cols);

fprintf(fileID, '\n');

fprintf(fileID, '解析データ　共通条件\n');
fprintf(fileID, '-- Mセル数：1\n');
fprintf(fileID, '-- Mセル半径：%d [m]\n', R*1000);
fprintf(fileID, '-- Sセル半径：%d [m]\n', r*1000);

fprintf(fileID, '\n');

fprintf(fileID, '＊＊解析データは「, 」のように，カンマとスペースで区切っており，最後のデータにはカンマはつけていない\n');
fprintf(fileID, '＊＊滞在距離データの単位は[m]\n');

fprintf(fileID, '\n');

fprintf(fileID, '以下，Sセル滞在距離解析データ\n\n');



for j = 1 : 1 : cols

    fprintf(fileID, '%8.6f', Sojourn_Distance_S(1, j)*1000);
    if j == cols
        break;
    end

    fprintf(fileID, ', ');
    
end

fprintf(fileID,'\n\n');

    

fclose(fileID);
