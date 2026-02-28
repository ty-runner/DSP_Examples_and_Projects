fsamp = 20000;
fcuts = [4000 4500];
mags = [1 0]; 
devs = [0.1 0.05];

%% Kaiser Design
[n1,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
n1 = n1 + rem(n1,2);
hh = fir1(n1,Wn,ftype,kaiser(n1+1,beta),"noscale");
fprintf('Kaiser Order %d\n', n1);

[H1,f] = freqz(hh,1,1024,fsamp);

%% Parks-McClellan Design
[n2,fo,ao,w] = firpmord(fcuts,mags,devs,fsamp);
n2 = n2 + rem(n2,2);
b = firpm(n2,fo,ao,w);
fprintf('Parks-McClellan Order %d\n', n2);

[H2,f] = freqz(b,1,1024,fsamp);
disp(H2);
Hmag = abs(H2);

% Passband indices
pb = f <= 4000;

% Stopband indices
sb = f >= 4500;

% Actual ripples
delta_p_actual = max(abs(Hmag(pb) - 1));
delta_s_actual = max(Hmag(sb));

% --- Step 3: Add delta_2 to the center tap in time domain ---
M = n2;                       % FIR order (even)
center = M/2 + 1;             % MATLAB index for n = M/2

delta2 = delta_s_actual;      % your choice (try + or - if needed)

b1 = b;
b1(center) = b1(center) + delta2;

% --- Steps 4 & 5: Modify zeros of b1 ---
tol = 1e-6;                   % tolerance for "on unit circle"
pull = 0.999;                 % if UC zeros not repeated, pull slightly inside

z = roots(b1);                % zeros of polynomial B(z) = b1(1) z^M + ... + b1(end)
r = abs(z);

z_out = z(r > 1 + tol);
z_in  = z(r < 1 - tol);
z_uc  = z(abs(r - 1) <= tol);

% Step 4: reflect outside zeros inside
z_out_ref = 1 ./ conj(z_out);

% Step 5: handle unit-circle zeros
% We'll cluster UC zeros by angle, and within each cluster:
% - if multiplicity is even: keep half (preserving conjugate pairing)
% - if multiplicity is odd: pull them slightly inside instead of deleting
keep_uc = [];
if ~isempty(z_uc)
    ang = angle(z_uc);
    [ang_s, idx] = sort(ang);
    z_uc_s = z_uc(idx);

    clusterStart = 1;
    for k = 2:(length(z_uc_s)+1)
        newCluster = (k > length(z_uc_s)) || (abs(ang_s(k) - ang_s(clusterStart)) > 1e-3);

        if newCluster
            cluster = z_uc_s(clusterStart:k-1);
            m = length(cluster);

            if mod(m,2) == 0 && m > 0
                % keep half of the zeros in this cluster
                keep_uc = [keep_uc; cluster(1:(m/2))];
            else
                % odd multiplicity (or single): pull inside instead of "removing half"
                keep_uc = [keep_uc; pull * cluster];
            end

            clusterStart = k;
        end
    end
end

% Rebuild new zero set and reconstruct FIR
z_new = [z_in; z_out_ref; keep_uc];

b_new = real(poly(z_new));    % poly gives coefficients of ∏(z - z_k)

% Normalize DC gain to match b1 (for lowpass, sum(b) = H(e^{j0}))
b_new = b_new * (sum(b1) / sum(b_new));

% Frequency response of modified filter
[H2_new,f] = freqz(b_new,1,1024,fsamp);
Hmag_new = abs(H2_new);

% Recompute ripple metrics on modified response
pb = f <= 4000;
sb = f >= 4500;

delta_p_after = max(abs(Hmag_new(pb) - 1));
delta_s_after = max(Hmag_new(sb));

fprintf('Actual passband ripple (before) = %f\n', delta_p_actual);
fprintf('Actual stopband ripple (before) = %f\n', delta_s_actual);
fprintf('Actual passband ripple (after)  = %f\n', delta_p_after);
fprintf('Actual stopband ripple (after)  = %f\n', delta_s_after);

% --- Plot both (Kaiser vs modified Parks) ---
figure
plot(f,abs(H1),'LineWidth',1.5)
hold on
plot(f,abs(H2_new),'LineWidth',1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Kaiser','Parks (zero-modified)')
title('Frequency Response Comparison')