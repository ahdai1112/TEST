%% pinhole chessboard projection demo
close all 
clear
clc

%% 設定相機參數

focalLengthX = 600; % 假設焦距為600像素
focalLengthY = 600;
principalPoint = [640, 480]'; % 假設影像中心為(320, 240)像素

% 相機內參數
K = [focalLengthX, 0, principalPoint(1);
     0, focalLengthY, principalPoint(2);
     0, 0, 1];

% 相機初始外參數(位置和姿態)

% 初始姿態(世界到相機)
Rc_w = roty(180);

% 初始位置(世界到相機)
tc_cw = [3;
         3;
         16];

%transformation matrix
pose =rigid3d(Rc_w' ,tc_cw');

T = [Rc_w, tc_cw; 0, 0, 0,1];

%% 設定figure
f = figure("Position",[61,81,1639,874]);
ax1 = axes("Parent",f,"Position", ...
    [0.054037055888621,0.084932230241155,0.44091191191838,0.806029994971516],"NextPlot","add");
axis equal;
box on;
title(ax1,'Virtual Chessboard Projection',"FontSize",15);
xlabel(ax1,'x_{cm}','FontSize',15);
ylabel(ax1,'y_{cm}','FontSize',15);
zlabel(ax1,'z_{cm}','FontSize',15);
axis([-15 20 -15 20 0 25]);
set(gca,'FontSize',15);
grid on;
view(1.560541489891217e+02,32.879971831323658); 

ax2 = axes("Parent",f,"Position",[0.550335570469799,0.206176925391137,0.427089688834655,0.595882571176369],"NextPlot","add");
xlabel(ax2,'u axes');
ylabel(ax2,'v axes');
set(gca,"YDir",'reverse'); %像素座標的v方向朝下
axis([0 1280 0 960]);
set(gca,'FontSize',15);
box on;
grid on;
grid minor;
title(ax2,'Image','FontSize',15);

%% 棋盤格設定

%棋盤格格數和大小
numRows = 15; 
numCols = 11; 
blockSize = 0.5; 

% 初始化頂點和面
vertices = zeros((numRows+1)*(numCols+1), 3);
faces = zeros(numRows*numCols, 4);
colors = zeros(numRows*numCols, 3);

% 計算每一格的頂點和面
for row = 0:numRows-1
    for col = 0:numCols-1
        
        x = (row - 1) * blockSize;
        y = (col - 1) * blockSize;
        z = 0;

        % 計算頂點
        v1 = row * (numCols+1) + col + 1;
        v2 = (row+1) * (numCols+1) + col + 1;
        v3 = (row+1) * (numCols+1) + col + 2;
        v4 = row * (numCols+1) + col + 2;

        % 添加頂點
        vertices([v1 v2 v3 v4], :) = [x y z; x+blockSize y z; x+blockSize y+blockSize z; x y+blockSize z];

        % 添加面
        faces(row*numCols + col + 1, :) = [v1 v2 v3 v4];

        % 用列+行的奇偶來添加顏色
        if mod(row + col, 2) == 0
            colors(row*numCols + col + 1, :) = [0 0 0]; % 黑色
        else
            colors(row*numCols + col + 1, :) = [1 1 1]; % 白色
        end
    end
end

% 使用patch函數繪製世界座標下的棋盤格
chessboard_world = patch(ax1,'Vertices', vertices, 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'flat');

%% plot of chessboard's corner plots
%以紅色點顯示棋盤格頂點
X = vertices(:,1);
Y = vertices(:,2);
Z = vertices(:,3);
plot3(ax1,X,Y,Z, 'r.','MarkerSize',3);

%% 繪製世界座標軸
% 繪製世界座標
% x軸
wordXCoordinate = line(ax1,[0 3], [0 0], [0 0], 'Color', 'r', 'LineWidth', 2);

% y軸
wordYCoordinate = line(ax1,[0 0], [0 3], [0 0], 'Color', 'g', 'LineWidth', 2);

% z軸
wordZCoordinate = line(ax1,[0 0], [0 0], [0 3], 'Color', 'b', 'LineWidth', 2);

%% 繪製camera
cam = plotCamera("AbsolutePose", pose,'Size',1,"AxesVisible",1,"Opacity",0.06,"Parent",ax1);

% movement line
xm = tc_cw(1)+3*cos(0);
ym = tc_cw(2)+3*sin(0);
zm = tc_cw(3);
movement = line(ax1,xm,ym,zm,'Color','r');

%% 初始化相機原點與世界座標中Z = 0的交點
line_leftTop = line(ax1,[0 0], [0 0], [0 0], 'Color', 'r');
line_leftBottom = line(ax1,[0 0], [0 0], [0 0], 'Color', 'r');
line_rightBottom = line(ax1,[0 0], [0 0], [0 0], 'Color', 'r');
line_rightTop = line(ax1,[0 0], [0 0], [0 0], 'Color', 'r');

%% 初始化位姿顯示
txt1 = text(ax1,0,-12,22,['Roll(deg) : ',num2str(0),' Pitch(deg) : ',num2str(0),' Yaw(deg) : ',num2str(0)],'FontWeight','bold');
txt2 = text(ax1,0,-12,20.5,['t_{x} : ',num2str(tc_cw(1)),'  t_{y} : ',num2str(tc_cw(2)),'  t_{z} : ',num2str(tc_cw(3))],'FontWeight','bold');

%% 初始化投影平面
intersection_plane = patch(ax1,'Vertices',zeros(4,3),'Faces',[1,2,3,4],'FaceColor','yellow','FaceAlpha', 0.2);
chessboard_pixel = patch(ax2,'Vertices', zeros(4,2), 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'flat');
chessboard_pixel_point = plot(ax2,0,0,'r.','MarkerSize',8);

%% 動畫錄製
% open the videowriter
% video = VideoWriter("pinhole chessboard projection demo","MPEG-4");
% 
% video.FrameRate = 60;
% open(video)

%% main loop


for i = 0:0.1:60

    r = 3; % 旋轉軌跡半徑
    zr = tc_cw(3)+0.1*i; % z軸移動範圍

    % 更新相機位姿
    updateTranslation = [tc_cw(1)+r*cos(i/3), tc_cw(2)+r*sin(i/3), zr]';
    updateRotation = roty(180)*rotz(20*sin(i/2))*roty(20*cos(i/2))*rotx(20*sin(i/2));
    updatePose = rigid3d(updateRotation, updateTranslation');
    cam.AbsolutePose = updatePose;
    
    % 繪製相機移動軌跡
    xm = [xm, updateTranslation(1)];
    ym = [ym, updateTranslation(2)];
    zm = [zm, updateTranslation(3)];
    
    movement.XData = xm;
    movement.YData = ym;
    movement.ZData = zm;


    T = [updateRotation', updateTranslation; 0, 0, 0,1];
    T = inv(T);  % world to camera
    updateTranslation_c_frame = T(1:3,4);
    % pinhole model

    % 在世界座標上表示的點
    point3D_world = [X, Y, Z]';

    % 世界座標轉換成齊次座標
    point3D_world_homo = [point3D_world;ones(1,length(vertices))];

    % 世界座標的點表示在相機座標上
    point3D_camera = T*point3D_world_homo;

    % 單位從公分轉換到像素
    point_pixel_homo = K*point3D_camera(1:3,:);

    % 去除深度的值，取u和v的值
    point_pixel = round(point_pixel_homo(1:2,:),0)./round(point3D_camera(3,:),0);
    
    %像素平面的四點邊界
    leftTop = [0, 0];
    leftBottom = [0, 960];
    rightTop = [1280, 0];
    rightBottom = [1280, 960];
    
    %解二元一次方程式得到對應像素座標邊界點在世界座標的位置

    %左上
    intersectionMatrix_leftTop = [focalLengthX*updateRotation(1,1)+principalPoint(1)*updateRotation(3,1)-updateRotation(3,1)*leftTop(1),...
                                  focalLengthX*updateRotation(1,2)+principalPoint(1)*updateRotation(3,2)-updateRotation(3,2)*leftTop(1);
                                  focalLengthY*updateRotation(2,1)+principalPoint(2)*updateRotation(3,2)-updateRotation(3,1)*leftTop(2),...
                                  focalLengthY*updateRotation(2,2)+principalPoint(2)*updateRotation(3,2)-updateRotation(3,2)*leftTop(2)];
    constant_leftTop = [-focalLengthX*updateTranslation_c_frame(1)-principalPoint(1)*updateTranslation_c_frame(3)+updateTranslation_c_frame(3)*leftTop(1);
                        -focalLengthY*updateTranslation_c_frame(2)-principalPoint(2)*updateTranslation_c_frame(3)+updateTranslation_c_frame(3)*leftTop(2)];
    sol_leftTop = intersectionMatrix_leftTop\constant_leftTop;
    
    %左下
    intersectionMatrix_leftBottom = [focalLengthX*updateRotation(1,1)+principalPoint(1)*updateRotation(3,1)-updateRotation(3,1)*leftBottom(1),...
                                     focalLengthX*updateRotation(1,2)+principalPoint(1)*updateRotation(3,2)-updateRotation(3,2)*leftBottom(1);
                                     focalLengthY*updateRotation(2,1)+principalPoint(2)*updateRotation(3,2)-updateRotation(3,1)*leftBottom(2),...
                                     focalLengthY*updateRotation(2,2)+principalPoint(2)*updateRotation(3,2)-updateRotation(3,2)*leftBottom(2)];
    constant_leftBottom = [-focalLengthX*updateTranslation_c_frame(1)-principalPoint(1)*updateTranslation_c_frame(3)+updateTranslation_c_frame(3)*leftBottom(1);
                           -focalLengthY*updateTranslation_c_frame(2)-principalPoint(2)*updateTranslation_c_frame(3)+updateTranslation_c_frame(3)*leftBottom(2)];
    sol_leftBottom = intersectionMatrix_leftBottom\constant_leftBottom;

    %右下
    intersectionMatrix_rightTop = [focalLengthX*updateRotation(1,1)+principalPoint(1)*updateRotation(3,1)-updateRotation(3,1)*rightTop(1),...
                                      focalLengthX*updateRotation(1,2)+principalPoint(1)*updateRotation(3,2)-updateRotation(3,2)*rightTop(1);
                                      focalLengthY*updateRotation(2,1)+principalPoint(2)*updateRotation(3,2)-updateRotation(3,1)*rightTop(2),...
                                      focalLengthY*updateRotation(2,2)+principalPoint(2)*updateRotation(3,2)-updateRotation(3,2)*rightTop(2)];
    constant_rightTop = [-focalLengthX*updateTranslation_c_frame(1)-principalPoint(1)*updateTranslation_c_frame(3)+updateTranslation_c_frame(3)*rightTop(1);
                            -focalLengthY*updateTranslation_c_frame(2)-principalPoint(2)*updateTranslation_c_frame(3)+updateTranslation_c_frame(3)*rightTop(2)];
    sol_rightTop = intersectionMatrix_rightTop\constant_rightTop;

    %右上
    intersetrionMatrix_rightBottom = [focalLengthX*updateRotation(1,1)+principalPoint(1)*updateRotation(3,1)-updateRotation(3,1)*rightBottom(1),...
                                   focalLengthX*updateRotation(1,2)+principalPoint(1)*updateRotation(3,2)-updateRotation(3,2)*rightBottom(1);
                                   focalLengthY*updateRotation(2,1)+principalPoint(2)*updateRotation(3,2)-updateRotation(3,1)*rightBottom(2),...
                                   focalLengthY*updateRotation(2,2)+principalPoint(2)*updateRotation(3,2)-updateRotation(3,2)*rightBottom(2)];
    constant_rightBottom = [-focalLengthX*updateTranslation_c_frame(1)-principalPoint(1)*updateTranslation_c_frame(3)+updateTranslation_c_frame(3)*rightBottom(1);
                         -focalLengthY*updateTranslation_c_frame(2)-principalPoint(2)*updateTranslation_c_frame(3)+updateTranslation_c_frame(3)*rightBottom(2)];
    sol_rightBottom = intersetrionMatrix_rightBottom\constant_rightBottom;
    
    % 繪製相機原點至世界座標Z=0的連線
    line_leftTop.XData = [updateTranslation(1,:) sol_leftTop(1)];
    line_leftTop.YData = [updateTranslation(2,:) sol_leftTop(2)];
    line_leftTop.ZData = [updateTranslation(3,:) 0];

    line_leftBottom.XData = [updateTranslation(1,:) sol_leftBottom(1)];
    line_leftBottom.YData = [updateTranslation(2,:) sol_leftBottom(2)];
    line_leftBottom.ZData = [updateTranslation(3,:) 0];

    line_rightBottom.XData = [updateTranslation(1,:) sol_rightBottom(1)];
    line_rightBottom.YData = [updateTranslation(2,:) sol_rightBottom(2)];
    line_rightBottom.ZData = [updateTranslation(3,:) 0];

    line_rightTop.XData = [updateTranslation(1,:) sol_rightTop(1)];
    line_rightTop.YData = [updateTranslation(2,:) sol_rightTop(2)];
    line_rightTop.ZData = [updateTranslation(3,:) 0];

    
    vertices_plane = [sol_leftTop(1) sol_leftTop(2) 0;
                      sol_leftBottom(1) sol_leftBottom(2) 0;
                      sol_rightBottom(1) sol_rightBottom(2) 0;
                      sol_rightTop(1) sol_rightTop(2) 0];

    faces_plane = [1 2 3 4];

    intersection_plane.Vertices = vertices_plane;
    
    % 更新位姿顯示
    roll = atan2(updateRotation(3,2),updateRotation(3,3));
    pitch = asin(-updateRotation(3,1));
    yaw = atan2(updateRotation(2,1),updateRotation(1,1));
    txt1.String = ['Roll(deg) : ',num2str(roll),' Pitch(deg) : ',num2str(pitch),' Yaw(deg) : ',num2str(yaw)];
    txt2.String = ['t_{x} : ',num2str(updateTranslation(1)),'  t_{y} : ',num2str(updateTranslation(2)),'  t_{z} : ',num2str(updateTranslation(3))];

    % 繪製世界點於像素座標上表示的點
    chessboard_pixel.Vertices = point_pixel';
    chessboard_pixel_point.XData = point_pixel(1,:);
    chessboard_pixel_point.YData = point_pixel(2,:);
   
    drawnow

    % frame = getframe(gcf);
    % writeVideo(video,frame);
 
end

%% 結束影片錄製
% close(video)























