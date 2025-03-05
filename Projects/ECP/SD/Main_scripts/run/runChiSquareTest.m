function [chi2stat, pVal, df, expected] = runChiSquareTest(observed)
    % runChiSquareTest - Performs a chi-square test of independence.
    %
    % Input:
    %   observed - A contingency table (matrix) of observed frequencies.
    %
    % Output:
    %   chi2stat - The chi-square test statistic.
    %   pVal     - The p-value of the test.
    %   df       - The degrees of freedom.
    %   expected - The expected frequencies under the null hypothesis.
    
    % Calculate the expected frequencies
    rowTotals = sum(observed, 2);
    colTotals = sum(observed, 1);
    total = sum(rowTotals);
    expected = (rowTotals * colTotals) / total;
    
    % Calculate the chi-square statistic
    chi2stat = sum((observed - expected).^2 ./ expected, 'all');
    
    % Calculate degrees of freedom
    [nRows, nCols] = size(observed);
    df = (nRows - 1) * (nCols - 1);
    
    % Calculate the p-value
    pVal = 1 - chi2cdf(chi2stat, df);
end