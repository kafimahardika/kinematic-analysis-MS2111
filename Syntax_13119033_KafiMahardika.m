%% Identitas
% Nama  : Kafi Mahardika
% NIM   : 13119033
% Kelas : K-03
%% Input
 
% Kita input data-data yang diberikan di awal

L2 = 83; % [50 + 033] mm
L3 = 166; % [2*L2] mm
L4 = 124.5; % [1.5*L2] mm 
c = 141.1; % [1.7*L2] mm
Omega2 = 15; % [rad/s]
IniTetha2 = pi/4 ; % [rad]
IniTetha3 = 0.2 ; % [rad]
IniTetha4 = 4 ; % [rad]
% Analisis pada waktu 0-2 [s] dengan jarak 0.001 [s]
t = 0:0.001:2 ; % [s]
NumData = size(t,2);

%% Matrix coordinate
% Kita tentukan matrix variabel awal q, serta matrix nol posisi, ang.
% velocity, dan ang. acceleration

q = [0;0;0;0;0;IniTetha2;0;0;IniTetha3;0;0;IniTetha4];
q_all = zeros(12,NumData);
v_all = zeros(12,NumData);
a_all = zeros(12,NumData);

%% Calculation

for j = 1:NumData
    for i = 1:3 %Jumlah pengulangan kalkulasi untuk mendapatkan hasil
        % Kita menentukan constraint matrix dari persamaan driving
        % constraint serta joint constraint
        
        % Constraint Matrix
        C = [q(1) ; q(2) ; q(3) ;...
            q(4)-(L2*cos(q(6)))/2;q(5)-(L2*sin(q(6)))/2 ;...
            q(4)+(L2*cos(q(6)))/2-q(7)+(L3*cos(q(9)))/2 ;...
            q(5)+(L2*sin(q(6)))/2-q(8)+(L3*sin(q(9)))/2 ;...
            q(7)+(L3*cos(q(9)))/2-q(10)+(L4*cos(q(12)))/2 ;...
            q(8)+(L3*sin(q(9)))/2-q(11)+(L4*sin(q(12)))/2 ;...
            q(10)+(L4*cos(q(12)))/2-(c*cos(q(3)));...
            q(11)+(L4*sin(q(12)))/2-(c*sin(q(3)));...
            q(6)-IniTetha2-Omega2*t(j)];
            
        % Untuk analisis posisi, kita gunakan deret Taylor yaitu
        % C(q + delta_q,t) = C + Cq*delta_q
        
        % Kita cari Constraint Jacobian Matrix (Cq), turunan C terhadap q
        
        % Jacobian Matrix 
        Cq = [1 0 0 0 0 0 0 0 0 0 0 0 ;
            0 1 0 0 0 0 0 0 0 0 0 0 ;
            0 0 1 0 0 0 0 0 0 0 0 0 ;
            0 0 0 1 0 (L2/2)*sin(q(6)) 0 0 0 0 0 0 ;
            0 0 0 0 1 (-L2/2)*cos(q(6)) 0 0 0 0 0 0 ;
            0 0 0 1 0 (-L2/2)*sin(q(6)) -1 0 (-L3/2)*sin(q(9)) 0 0 0 ;
            0 0 0 0 1 (L2/2)*cos(q(6)) 0 -1 (L3/2)*cos(q(9)) 0 0 0 ;
            0 0 0 0 0 0 1 0 (-L3/2)*sin(q(9)) -1 0 (-L4/2)*sin(q(12)) ;            
            0 0 0 0 0 0 0 1 (L3/2)*cos(q(9)) 0 -1 (L4/2)*cos(q(12)) ;
            1 0 c*sin(q(3)) 0 0 0 0 0 0 1 0 (-L4/2)*sin(q(12)) ;
            0 1 -c*cos(q(3)) 0 0 0 0 0 0 0 1 (L4/2)*cos(q(12)) ;
            0 0 0 0 0 1 0 0 0 0 0 0] ;
        
        % Setelah itu kita tentukan newton differences (delta_q) dan
        % tambahkan ke variabel q pada setiap iterasi untuk mendapatkan 
        % matrix jacobian yang baru dan terus dilakukan hingga iterasi
        % selesai dan mendapatkan posisi masing-masing
        
        delta_q = inv(Cq)*(-C) ;
        q = q + delta_q ;
        
        % Untuk analisis kecepatan, Constraint matrix kita turunkan
        % terhadap waktu (Ct). Kita tau bahwa Cq*qdot + Ct = 0, maka
        
        % Velocity Analysis
        Ct = [0;0;0;0;0;0;0;0;0;0;0;-Omega2];
        qdot = inv(Cq)*(-Ct);
        
        % Untuk Analisis Percepatan, kita gunakan persamaan Cq*qdot2=Qd,
        % dimana Qd=-Cq_qdot*Cq_qdot_q-2*Cqt*qdot-Ctt
        
        %Acceleration Analysis
        Cq_qdot = [qdot(1);qdot(2);qdot(3) ;...
            qdot(4)+(qdot(6)*(L2/2)*sin(q(6))) ;...
            qdot(5)-(qdot(6)*(L2/2)*cos(q(6))) ;...
            qdot(4)-(qdot(6)*(L2/2)*sin(q(6)))-qdot(7)-(qdot(9)*(L3/2)*sin(q(9))) ;...
            qdot(5)+(qdot(6)*(L2/2)*cos(q(6)))-qdot(8)+(qdot(9)*(L3/2)*cos(q(9))) ;...
            qdot(7)-(qdot(9)*(L3/2)*sin(q(9)))-qdot(10)-(qdot(12)*(L4/2)*sin(q(12))) ;...
            qdot(8)+(qdot(9)*(L3/2)*cos(q(9)))-qdot(11)+(qdot(12)*(L4/2)*cos(q(12))) ;...
            qdot(10)-(qdot(12)*(L4/2)*sin(q(12)))+(qdot(3)*c*sin(q(3))) ;...
            qdot(11)+(qdot(12)*(L4/2)*cos(q(12)))-(qdot(3)*c*cos(q(3))) ;...
            qdot(6)];
        
        Cq_qdot_q = [zeros(1,12);
            zeros(1,12);
            zeros(1,12);
            zeros(1,5) -(qdot(6)*(L2/2)*cos(q(6))) zeros(1,6);
            zeros(1,5) -(qdot(6)*(L2/2)*sin(q(6))) zeros(1,6);
            zeros(1,5) -(qdot(6)*(L2/2)*cos(q(6))) zeros(1,2) -(qdot(9)*(L3/2)*cos(q(9))) zeros(1,3);
            zeros(1,5) -(qdot(6)*(L2/2)*sin(q(6))) zeros(1,2) -(qdot(9)*(L3/2)*sin(q(9))) zeros(1,3);
            zeros(1,8) -(qdot(9)*(L3/2)*cos(q(9))) zeros(1,2) -(qdot(12)*(L4/2)*cos(q(12)));
            zeros(1,8) -(qdot(9)*(L3/2)*sin(q(9))) zeros(1,2) -(qdot(12)*(L4/2)*sin(q(12)));
            zeros(1,2) (qdot(3)*c*cos(q(3))) zeros(1,8) -(qdot(12)*(L4/2)*cos(q(12)));
            zeros(1,2) (qdot(3)*c*sin(q(3))) zeros(1,8) -(qdot(12)*(L4/2)*sin(q(12)));
            zeros(1,12)];
        
        % Turunan C terhadap q dan t
        Cqt = zeros(12,12);
        % Turunan Ct terhadap t
        Ctt = zeros(12,1);
        
        % Dari matrix diatas, kita hitung menggunakan persamaan di bawah 
        Qd = -Cq_qdot_q*qdot-2*Cqt*qdot-Ctt;
        % Percepatan didapatkan dengan
        qdot2 = inv(Cq)*Qd;
        
    end
    % Kita masukkan posisi, kecepatan, dan percepatan yang telah dihitung
    % ke suatu matriks, dimana tiap kolom melambangkan tiap detiknya
    q_all(:,j)=q;
    v_all(:,j)=qdot;
    a_all(:,j)=qdot2;
end

% Trajectory of Joint A
figure(1);
plot3(2*q_all(4,:), 2*q_all(5,:), t)
grid on
xlabel('Trajectory (X) [mm]');
ylabel('Trajectory (Y) [mm]');
zlabel('Time [s]');
title('Trajectory of Joint A');

% Trajectory of Joint B
figure(2);
plot3(2*q_all(10,:), 2*q_all(11,:), t)
grid on
xlabel('Trajectory (X) [mm]');
ylabel('Trajectory (Y) [mm]');
zlabel('Time [s]');
title('Trajectory of Joint B');

% Angular Velocity of Link 3 vs Time
figure(3);
plot(t, v_all(9,:));
grid on
xlabel('Time [s]');
ylabel('Angular Velocity [rad/s]');
title('Angular Velocity of Link 3 vs Time')

% Angular Velocity of Link 4 vs Time
figure(4);
plot(t, v_all(12,:));
grid on
xlabel('Time [s]');
ylabel('Angular Velocity [rad/s]');
title('Angular Velocity of Link 4 vs Time')

% Angular Acceleration of Link 3 vs Time
figure(5);
plot(t, a_all(9,:));
grid on
xlabel('Time [s]');
ylabel('Angular Acceleration [rad/s^2]');
title('Angular Acceleration of Link 3 vs Time');

% Angular Acceleration of Link 4 vs Time
figure(6);
plot(t, a_all(12,:));
grid on
xlabel('Time [s]');
ylabel('Angular Acceleration [rad/s^2]');
title('Angular Acceleration of Link 4 vs Time');            