% TESTmeetingAug2026.m
run('../src/sandbox/meetingAug2026.m');

assert(exist('kx_est', 'var'), 'kx_est not computed');
assert(exist('ky_est', 'var'), 'ky_est not computed');
assert(abs(kx_est - kx_true) < 0.5, sprintf('kx_est=%g, expected %g', kx_est, kx_true));
assert(abs(ky_est - ky_true) < 0.5, sprintf('ky_est=%g, expected %g', ky_est, ky_true));

disp('TESTmeetingAug2026: ALL PASSED');

