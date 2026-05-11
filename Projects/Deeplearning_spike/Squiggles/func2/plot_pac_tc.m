function plot_pac_tc(tc, tsec, labelStr)
plot(tsec, tc, 'LineWidth', 1.5); grid on
xlabel('time (s) wrt onset'); ylabel('zPAC'); title(labelStr);
end