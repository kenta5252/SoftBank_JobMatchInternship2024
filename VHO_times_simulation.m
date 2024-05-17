%% 垂直HO数の解析(滞在時間＝接続時間閾値設けたときのシミュレーション)
%% ユーザの移動する線分は角度が [-pi/2, pi/2) の範囲で一様分布を採用


%% アルゴリズム
% ループ① Sセル数を変える：結果の横軸
%     ループ② 滞在時間閾値を変える：結果の凡例
%         ループ③ シミュレーション開始(Mセルの出発地点を決定)：HO数の平均をとる
%             ループ④ UE(移動端末)がMセルを横断
%                     横断中に垂直HOを加算


%% Sセル個数　100-1000
%% ユーザ速度　50km/h

clc;
clear;
tic;

%% 変えていく変数 %%%%%%%%%%%%%%%%%
Simulation_kaisu = 10;
N_S = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000];
T_threshold = [1.0, 2.0, 3.0, 4.0, 5.0];
User_v=50; % UEの速度[km/h]
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% パラメータ値設定 %%%%%%%%%%%%%%%%%
In=0; % 円の中の基地局数(乱数)
Out=0; % 円の外の基地局数(乱数)
R=1; % Mセル半径[km]
r=0.04; % Sセル半径[km]
tani_ugoki=0.0001; % Userの動きの細かさ
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 変数設定 %%%%%%%%%%%%%%%%%
Sojourn_Distance_S = 0; % Sセル滞在時間を確保しておく配列
Sojourn_Time_S = 0;
maxSojourn_Distance_M = 2*R; % Sセルの最長通過距離
maxSojourn_Time_M = maxSojourn_Distance_M / tani_ugoki; % 動きの細かさに対するMセルを通過する最大の時間(回数)
% Sセル内を動く度に計算する変数
% UEを動かすfor文の都合で+2をしている
Cos_a=zeros(In,maxSojourn_Time_M+2);
Sin_a=zeros(In,maxSojourn_Time_M+2);
Cos_b=zeros(In,maxSojourn_Time_M+2); % b=sita-a
Sin_b=zeros(In,maxSojourn_Time_M+2);

SojournNum = 1; % 滞在したSセルを順番に格納していくためのカウントする変数

VHO = zeros(length(T_threshold), length(N_S)); % 垂直HO数の合計
Average_VHO = zeros(length(T_threshold), length(N_S)); % 垂直HＯ数の平均（これを理論の期待値と比較する）

jyoutai = 0; % 今現在，Ｍセル（1）かＳセル（2）に接続しているかを判断するスカラー値変数

Length_prevScell_user = 0; % 現在接続中のSセルとの距離（jyoutai==2の時のみ使う変数）
Length_nextScell_user = 0; % 次に接続するSセル候補との距離（jyoutai==1,2の時に使う変数）

% SmallCellNum; % 現在のSセル数の順番
Scell = 0; % 現在のSセル数

yesHO(1) = 0; % 一回ループ内においてHOしたかしていないかを判断 していたら「１」
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% ループ①　Sセル数を変える：結果の横軸 %%%%%%%%%%%%%%%%%
for SmallCellNum = 1 : 1 : length(N_S)

    Scell = SmallCellNum * 100; % Sセル数


    %% ポアソン点過程に従ってSセル配置 %%%%%%%%%%%%%%%%%
    sum_points = Scell  * 2;
    lam = 1 : 1 : sum_points;
    
    old_x = poissrnd(lam); % 0から200の間でポアソンに従った乱数の300個のデータ
    old_y = poissrnd(lam);
    
    x_y = zeros(sum_points, 2);
    
    for number = 1 : sum_points
        x_y(number,1) = old_x(1, number) - Scell;
        x_y(number,2) = old_y(1, number) - Scell;
        x_y(number,1) = x_y(number,1) / Scell; % x軸:ー１から+1の間でポアソンに従った乱数が発生
        x_y(number,2) = x_y(number,2) / Scell; %
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
            ransu_naka(In+1,1) = shuffled_x_y(N_in,1);
            ransu_naka(In+1,2) = shuffled_x_y(N_in,2);
            In = In + 1;
        end
        if In == Scell 
            break
        end
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    %% ループ②　滞在時間閾値を変える：結果の凡例 %%%%%%%%%%%%%%%%%
    for TimeThresholdNum = 1 : 1 : length(T_threshold)

        T_th = T_threshold(1, TimeThresholdNum);



        %% ループ③ シミュレーション開始(Mセルの中心とSセルの距離の深さが一様分布)%%%%%%%%%%%%%%%%%
        for x = 1 : Simulation_kaisu % シミュレーション回数は何回するか？
        
            % 半径Rの円から無作為に半径を選ぶ
            selected_radius = R * rand();
        
            % 原点 O
            % 無作為に選んだ点 A
            % 始点 U
        
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

            SelectScell = zeros(1, maxSojourn_Time_M);
        
            %%%%%%%%%　UEの出発始点と角度が決定　%%%%%%%%%
    
    
            kaunto = zeros(In,1);
            x_now(1)= OU(1) ; % ループのために初期状態を配列にいれる
            y_now(1)= OU(2) ; % ループのために初期状態を配列にいれる

            jyoutai = 0; % 動きはじめは，状態を0にしておく
        
            %%%%%%%%%　ループ④ UEが移動するループ　%%%%%%%%%         
            for tt = 2 : 1 : maxSojourn_Time_M + 2 % tt回単位当たりで動く
        
                x_now(tt) = x_now(tt-1) + tani_ugoki * cos(Angle_of_Travel);
                y_now(tt) = y_now(tt-1) + tani_ugoki * sin(Angle_of_Travel);
        
        
                if tt > (Sojourn_Distance_M / tani_ugoki) % このifでnnのforループを終わらせるように
                    break;
                end

                if jyoutai == 2 % Sセルに接続しているとき，そのSセルを出てしまったら，Mセルに接続する
                    Length_prevScell_user = ((ransu_naka(SelectScell(1, tt-1),1)-x_now(tt))^2+(ransu_naka(SelectScell(1, tt-1),2)-y_now(tt))^2)^(1/2);
                    if Length_prevScell_user > r
                        jyoutai = 1; % SセルからMセルにVHO
                        VHO(TimeThresholdNum, SmallCellNum) = VHO(TimeThresholdNum, SmallCellNum) + 1; % M->SのVHO
                    end
                end

                yesHO(1) = 0; % 下のfor文でHOしていないかを判断して，上のprevの計算で利用
        
                for n = 1 : In % Sセルの数In個をnで全てみる 最短距離のセルと接続(簡略化のため，同時に２つとは接続しないという過程で)
                    
                    Length_nextScell_user = ((ransu_naka(n,1)-x_now(tt))^2 + (ransu_naka(n,2)-y_now(tt))^2 )^(1/2);
        
                    % t==3の時点で初期値(jyoutai=0)が変わっていないてことは，出発地点で近くにSセルがいなかったことを意味するからここで状態を変える
                    if tt == 3 % 最初から２つ目
                        if jyoutai == 0
                            jyoutai = 1; % Mセルに接続している
                            yesHO = yesHO + 1; % 知りたいのはSセルに接続し続けているかだからこれもカウントしておく
                        end
                    end
        
                    if kaunto(n,1) == 0
                            
                        if Length_nextScell_user <= r
                            Cos_a(n,tt) = (ransu_naka(n,1)-x_now(tt)) / Length_nextScell_user; 
                            Sin_a(n,tt) = (ransu_naka(n,2)-y_now(tt)) / Length_nextScell_user;
                            Cos_b(n,tt) = Cos_a(n,tt) * cos(Angle_of_Travel) + Sin_a(n,tt) * sin(Angle_of_Travel); % 加法定理で接続先候補のセルの滞在距離を導出
                            Sin_b(n,tt) = Sin_a(n,tt) * cos(Angle_of_Travel) - Cos_a(n,tt) * sin(Angle_of_Travel);

                            Sojourn_Distance_S = Length_nextScell_user * Cos_b(n,tt) + (r^2 - (Length_nextScell_user * Sin_b(n,tt))^2) ^ (1/2); % 予想される滞在距離[km]
                            Sojourn_Time_S = 3600 * Sojourn_Distance_S / User_v;
        
                            %%%%%%%%%　垂直HO数を加算するフェーズ　%%%%%%%%% 
                            if jyoutai == 0 % 最初の場合分け　このif文が入るのは，tt=2の時のみ
                                %UEが移動するときに最初にSセルの位置からいたらそれは垂直HOすることにならないから異動の初めに接続状態を確認して垂直HO数を加算する
                                jyoutai = 2; %Sセルに接続している
                                kaunto(n,1) = 1;
                                SelectScell(1, tt) = n;
                                SojournNum = SojournNum + 1;
                                yesHO = yesHO + 1;
                                break
                            end

                            if Sojourn_Time_S >= T_th % 滞在時間閾値を超えたら

                                if jyoutai == 1 % Mセルに接続していたら
                                    jyoutai = 2; % Sセルに接続セル
                                    VHO(TimeThresholdNum, SmallCellNum) = VHO(TimeThresholdNum, SmallCellNum) + 1; % M->SのVHO
                                    kaunto(n,1) = 1;
                                    SelectScell(1, tt) = n;
                                    SojournNum = SojournNum + 1;
                                    yesHO = yesHO + 1;
                                    break
                                end

                                if jyoutai == 2 % Sセルに接続していたら，違うSセルに接続する
                                    SelectScell(1, tt) = n;
                                    SojournNum = SojournNum + 1;
                                    yesHO = yesHO + 1;
                                    break
                                end

                            end
        
                        end
                            
                    end
        
                end % Sセルの数In個をnで全てみる

                if yesHO == 0
                    SelectScell(1, tt) = SelectScell(1, tt-1);
                end
        
            end % １経路を横断終了
            
        end % シミュレーション回数分回す

        Average_VHO(TimeThresholdNum, SmallCellNum) = VHO(TimeThresholdNum, SmallCellNum) / Simulation_kaisu;

    end % 滞在時間閾値ごとのループ

    toc;


end % Sセル数を変えていく




%% ファイル作成

filename = 'VHO_times_Simulation_with_threshold.txt';

fileID = fopen(filename,'w');


fprintf(fileID, 'シミュレーション回数：%d\n\n', Simulation_kaisu);




for SmallCellNum = 1 : 1 : length(N_S)

    fprintf(fileID, 'Sセル数：%d\n', N_S(1, SmallCellNum));

    for TimeThresholdNum = 1 : 1 : length(T_threshold)

        fprintf(fileID, '滞在時間閾値：%.2f [s]   平均垂直HO数：%.2f [回] \n', T_threshold(1, TimeThresholdNum), Average_VHO(TimeThresholdNum, SmallCellNum));


    end

    fprintf(fileID, '\n');

end

fprintf(fileID,'\n\n');


fclose(fileID);
