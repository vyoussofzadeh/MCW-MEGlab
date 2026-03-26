function plot_two_rater_timeline(Tm, r1, r2, base)
    if ~ismember('time_sec', Tm.Properties.VariableNames)
        x = Tm.sample; xlab = 'Sample index';
    else
        x = Tm.time_sec; xlab = 'Time (s)';
    end

    figure; hold on;
    plot(x(Tm.bothAgree),  ones(sum(Tm.bothAgree),1),  'g.', 'MarkerSize', 14);
    plot(x(Tm.bothReject), ones(sum(Tm.bothReject),1), 'r.', 'MarkerSize', 14);
    plot(x(Tm.disagree),   ones(sum(Tm.disagree),1),   'k.', 'MarkerSize', 18);

    ylim([0.5 1.5]);
    yticks(1); yticklabels("Candidates");
    xlabel(xlab);
    title(sprintf('%s: %s vs %s (green both agree, red both reject, black disagree)', base, r1, r2));
    grid on;

    txt = sprintf('bothAgree=%d | bothReject=%d | disagree=%d | overlap=%d', ...
        sum(Tm.bothAgree), sum(Tm.bothReject), sum(Tm.disagree), sum(Tm.bothReviewed));
    text(0.01, 0.92, txt, 'Units','normalized');
end

