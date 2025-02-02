function plotMEGandRsnrDiscordant(MEG_LI, rSNR_left, rSNR_right, discordSubs, subIDs)
% PLOTMEGANDRSNRDISCORDANT  Plots MEG_LI vs. subject index and rSNR (left/right)
%                           only for the discordant samples.
%
%   plotMEGandRsnrDiscordant(MEG_LI, rSNR_left, rSNR_right, discordSubs, subIDs)
%
%   INPUTS:
%       MEG_LI     - (#Subjects x 1) array of MEG laterality index values
%       rSNR_left  - (#Subjects x #TimePoints) numeric array (left hemisphere rSNR)
%       rSNR_right - (#Subjects x #TimePoints) numeric array (right hemisphere rSNR)
%       discordSubs- Vector of subject indices that are "discordant"
%       subIDs     - (#Subjects x 1) numeric or cell array of subject labels (or 1:N)
%
%   This function:
%       1) Creates a figure with two subplots:
%           (a) Subplot 1: a bar or scatter of MEG_LI for discordant subjects
%           (b) Subplot 2: a line plot of rSNR_left vs. rSNR_right
%               averaged across time (or pick a time point) for discordant subjects
%       2) Highlights the subject IDs if provided in subIDs.
%
%   EXAMPLE:
%       plotMEGandRsnrDiscordant(MEG_LI, rSNR_L, rSNR_R, [4,5,10], 1:72);
%
%   Author: (Your Name / Organization)
%   Date: (Date)


if nargin < 5
    error('All 5 inputs are required. Check your usage.');
end

% Validate sizes
nSubjects = length(MEG_LI);
if size(rSNR_left,1) ~= nSubjects || size(rSNR_right,1) ~= nSubjects
    error('rSNR_left/right must have the same #Subjects as MEG_LI');
end
if length(subIDs) ~= nSubjects
    error('subIDs length must match the number of subjects.');
end

% Extract discordant subset for easy reference
discordMEG_LI = MEG_LI(discordSubs);

% For demonstration, let's average rSNR across time for each hemisphere
% or you can select a specific time index if desired
avg_rSNR_left  = mean(rSNR_left,2,'omitnan');
avg_rSNR_right = mean(rSNR_right,2,'omitnan');

discord_rSNR_left  = avg_rSNR_left(discordSubs);
discord_rSNR_right = avg_rSNR_right(discordSubs);

% Create figure
figure('Name','MEG LI & rSNR','Color','w','Position',[200,200,450,450]);

% ---------------- Subplot 1: Plot MEG_LI for discordant subs ----------------
subplot(1,2,1);
hold on; box off; grid on;

% We can do a bar plot or scatter
% Example: bar chart
bar(discordMEG_LI, 'FaceColor',[0.2 0.6 0.8]);
ylabel('MEG LI');
% title('MEG-LI Discordant Sub');
set(gca, 'XTick', 1:length(discordSubs), 'XTickLabel', subIDs(discordSubs));

xtickangle(45);
% optional: add a horizontal y=0 line
yline(0,'k--');

% ---------------- Subplot 2: rSNR Left vs Right (averaged) ----------------
subplot(1,2,2);
hold on; box off; grid on;

% We'll do a scatter comparing left vs. right rSNR
scatter(discord_rSNR_left, discord_rSNR_right, 50, 'filled','MarkerFaceColor','r');
xlabel('Avg rSNR Left');
ylabel('Avg rSNR Right');
% title('rSNR Discordant Subs');
% Optionally label each point with subject ID
for i = 1:length(discordSubs)
    text(discord_rSNR_left(i), discord_rSNR_right(i), num2str(subIDs(discordSubs(i))), ...
        'VerticalAlignment','bottom','HorizontalAlignment','left',...
        'FontSize',8,'Color','k');
end

% reference line
xL = xlim;
plot(xL,xL,'k--','LineWidth',1); % y = x line

hold off;

end
