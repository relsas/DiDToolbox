clear;

fprintf('=== Test SDID / Synthetic Control ===\n');

% 1. Create Data (SC Style)
% 1 Treated Unit (ID=1), 20 Controls.
% T=20 (10 pre, 10 post).
% Treated Unit = 0.5 * Co1 + 0.5 * Co2 + Noise
% Effect = +5 in post.

T_pre = 30; T_post = 10; T = T_pre + T_post;
% Note: We need T_pre > N_co to ensure the SC weights are uniquely identified.
% With T_pre=10 and N_co=20, the system is underdetermined, and weights are dispersed.
N_co = 20;

time = (1:T)';
% Controls 1-2 are Random Walks (Non-stationary)
% Controls 3-20 are White Noise (Stationary) - SC should ignore them
Y_co = randn(T, N_co);
% Make first 2 cumulative (trends) - STRONG signal
Y_co(:,1:2) = cumsum(randn(T, 2)*1.0) + randn(1, 2);
% Make 3:end just noise (or very weak trend) - WEAK noise
Y_co(:,3:end) = randn(T, N_co-2)*0.1;

Y_trt = 0.5 * Y_co(:,1) + 0.5 * Y_co(:,2) + randn(T,1)*0.01;
Y_trt(T_pre+1:end) = Y_trt(T_pre+1:end) + 5; % Treatment Effect

% Pivot to Long
id_trt = ones(T,1);
time_trt = time;
y_trt = Y_trt;
d_trt = zeros(T,1); d_trt(T_pre+1:end) = 1;

id_co = []; time_co = []; y_co = []; d_co = [];
for i=1:N_co
    id_co = [id_co; repmat(i+1, T, 1)];
    time_co = [time_co; time];
    y_co = [y_co; Y_co(:,i)];
    d_co = [d_co; zeros(T,1)];
end

data = table([id_trt; id_co], [time_trt; time_co], [y_trt; y_co], [d_trt; d_co], ...
    'VariableNames', {'id','time','y','D'});

ds = did.Dataset.fromTable(data,"idVar","id","timeVar","time","yVar","y","dVar","D");


% 2. Run SC
fprintf('\n--- Running Synthetic Control (SC) ---\n');
resSC = did.fit("sc", ds, "Display", true);

fprintf('SC Weight on Unit 2 (Control 1): %.4f\n', resSC.omega(1));
fprintf('SC Weight on Unit 3 (Control 2): %.4f\n', resSC.omega(2));
fprintf('SC Est: %.4f (True: 5.0)\n', resSC.tau(1));

% Test Plotting
try
    fprintf('Testing SDID.plot(resSC)...\n');
    did.estimators.SDID.plot(resSC);
    close(gcf); % Close figure to avoid popup
    fprintf('Plotting successful.\n');
catch ME
    warning('Plotting failed: %s', ME.message);
end

assert(abs(resSC.tau - 5) < 0.5, 'SC Estimate bias large');
assert(abs(sum(resSC.omega(1:2)) - 1) < 0.1, 'SC Weights detection failed');

% 3. Run SDID
fprintf('\n--- Running SDID ---\n');
resSDID = did.fit("sdid", ds, "Display", true);

fprintf('SDID Est: %.4f (True: 5.0)\n', resSDID.tau);
assert(abs(resSDID.tau - 5) < 0.5, 'SDID Estimate bias large');


% 4. Solver Check (Method=Frank-Wolfe)
fprintf('\n--- Running SDID (Frank-Wolfe explicit) ---\n');
resFW = did.fit("sdid", ds, "Display", true, "Solver", "frank-wolfe");
fprintf('SDID(FW) Est: %.4f\n', resFW.tau);

% 5. Test Summary Table
fprintf('\n--- Summary Table Check ---\n');
disp(resSDID.summaryTable);

fprintf('\nSuccess.\n');
