%% CFD Solver of subsonic-supersonic isentropic nozzle flow
%% Method
% Maccormack's Technique

%% Basic parameters
n = 31;             %node numbers
L = 3.0;            %total length
% rou = zeros(1, n);  %density
% v = zeros(1, n);    %velocity
% T = zeros(1, n);    %temperture
% A = zeros(1, n);    %section area
x = linspace(0, L, n);  %node coordinate
dx = L / (n - 1);   %space steplength
dt_l = zeros(1, n); %time steplength list
drdt1 = zeros(1, n);    %dr/dt in predict step
dvdt1 = zeros(1, n);    %dv/dt in predict step
dTdt1 = zeros(1, n);    %dT/dt in predict step
drdt2 = zeros(1, n);    %dr/dt in correct step
dvdt2 = zeros(1, n);    %dv/dt in correct step
dTdt2 = zeros(1, n);    %dT/dt in correct step
C = 0.5;            %Courant number
STEP_N = 1000;      %step number
t = 0;              %total time
gama = 1.4;         %gama

%% Initialize
A = 1 + 2.2 .*(x - 1.5) .^2;
rou = 1 - 0.314 .* x;
T = 1 - 0.2314 .* x;
v = (0.1 + 1.09 .* x) .* sqrt(T);

%% History data list
%Sample point : mid point
rou_hist = zeros(1, STEP_N);
v_hist = zeros(1, STEP_N);
T_hist = zeros(1, STEP_N);
M_hist = zeros(1, STEP_N);

%% Calculation
for iter = 1:STEP_N
    %Record history data
    rou_hist(iter) = rou(ceil(n/2));
    v_hist(iter) = v(ceil(n/2));
    T_hist(iter) = T(ceil(n/2));
    M_hist(iter) = v(ceil(n/2)) / sqrt(T(ceil(n/2)));
    
    %Calculate minimum timespace
    dt_l = C .* dx ./ (abs(v) + sqrt(T));
    dt = min(dt_l);
    t = t + dt;
    
    %Predict step
    for i = 1:n - 1
        drdt1(i) = -v(i) * (rou(i+1) - rou(i)) / dx - rou(i) * (v(i+1) - v(i)) / dx...
            - rou(i) * v(i) * (log(A(i+1)) - log(A(i))) / dx;
        dvdt1(i) = -v(i) * (v(i+1) - v(i)) / dx - ((T(i+1) - T(i)) / dx...
            + T(i) / rou(i) * (rou(i+1) - rou(i)) / dx) / gama;
        dTdt1(i) = -v(i) * (T(i+1) - T(i)) / dx - (gama - 1) * T(i) * ...
            ((v(i+1) - v(i)) / dx + v(i) * (log(A(i+1)) - log(A(i))) / dx);
    end
    
    rou_pre = rou + drdt1 .* dt;
    v_pre = v + dvdt1 .* dt;
    T_pre = T + dTdt1 .* dt;
    
    %Correct step
    for i = 2:n - 1
        drdt2(i) = -v_pre(i) * (rou_pre(i) - rou_pre(i-1)) / dx - ...
            rou_pre(i) * (v_pre(i) - v_pre(i-1)) / dx - ...
            rou_pre(i) * v_pre(i) * (log(A(i)) - log(A(i-1))) / dx;
        dvdt2(i) = -v_pre(i) * (v_pre(i) - v_pre(i-1)) / dx - ...
            ((T_pre(i) - T_pre(i-1)) / dx + T_pre(i) / rou_pre(i) * ...
            (rou_pre(i) - rou_pre(i-1)) / dx) / gama;
        dTdt2(i) = -v_pre(i) * (T_pre(i) - T_pre(i-1)) / dx - (gama - 1) * ...
            T_pre(i) * ((v_pre(i) - v_pre(i-1)) / dx + v_pre(i) * ...
            (log(A(i)) - log(A(i-1))) / dx);
    end
    
    drdt = (drdt1 + drdt2) ./ 2;
    dvdt = (dvdt1 + dvdt2) ./ 2;
    dTdt = (dTdt1 + dTdt2) ./ 2;
    
    rou = rou + drdt .* dt;
    v = v + dvdt .* dt;
    T = T + dTdt .* dt;
    
    v(1) = 2 * v(2) - v(3);
    v(n) = 2 * v(n-1) - v(n-2);
    rou(1) = 1;
    rou(n) = 2 * rou(n-1) - rou(n-2);
    T(1) = 1.0;
    T(n) = 2 * T(n-1) - T(n-2);
    
    %Print debug infomation
%     if mod(i, 10) == 1
%         fprintf('\n    t \t   dt \t  drdt\t  dvdt\t  dTdt\t   rou\t    v \t    T \n');
%     end
%     fprintf('%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t\n',...
%     t, dt, average(drdt), average(dvdt), average(dTdt), rou(n/2), v(n/2), T(n/2));
end

%% Plot final solution
subplot(2,2,1)
plot(1:STEP_N, rou_hist);
subplot(2,2,2)
plot(1:STEP_N, v_hist);
subplot(2,2,3)
plot(1:STEP_N, T_hist);
subplot(2,2,4)
plot(1:STEP_N, M_hist);