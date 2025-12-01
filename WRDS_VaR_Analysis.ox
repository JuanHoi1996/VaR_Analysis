#include <oxstd.h>

/**
 * WRDS CRSP Data VaR Analysis
 * Author: Generated Code
 * Description: Calculate VaR for 7 stocks using return data from WRDSCRSPdata.csv
 * Data: Microsoft, Chevron, Apple, Disney, Walmart, Nike, American Express (2015-2020)
 */

// =================== GLOBAL CONSTANTS ===================
decl g_iNumAssets;         // 7 stocks
decl g_iHistoricalPeriod;  // 800 trading days for historical window
decl g_iHoldingPeriod;     // 10 days holding period
decl g_iInvestmentAmount;  // Investment per asset
decl g_iNumSimulations;    // Monte Carlo simulations for historical method
decl g_dZScore99;          // -2.3263 - 99% confidence z-score
decl g_sDataFile;          // Data file name

// Asset names for reference
decl g_asAssetNames;

// Initialize global constants
initializeConstants()
{
    g_iNumAssets = 7;
    g_iHistoricalPeriod = 800;
    g_iHoldingPeriod = 10;
    g_iInvestmentAmount = 5000000;  // Equal investment per asset
    g_iNumSimulations = 1000;
    g_dZScore99 = -2.3263;
    g_sDataFile = "WRDSCRSPdata.csv";
    
    // Asset names for output
    g_asAssetNames = {"Microsoft", "Chevron", "Apple", "Disney", "Walmart", "Nike", "American_Express"};
}

/**
 * Load return data from WRDSCRSPdata.csv
 */
loadReturnData()
{
    println("Loading return data from: ", g_sDataFile);
    
    decl mData = loadmat(g_sDataFile);
    if (mData == <>) {
        println("Error: Could not load ", g_sDataFile);
        return <>;
    }
    
    println("Data loaded: ", rows(mData), " observations x ", columns(mData), " assets");
    
    // Check if we have the expected number of assets
    if (columns(mData) != g_iNumAssets) {
        println("Warning: Expected ", g_iNumAssets, " assets, found ", columns(mData));
    }
    
    // Display first few rows for verification
    println("First 5 observations preview:");
    decl i, j;
    for (i = 0; i < min(5, rows(mData)); i++) {
        print("Day ", i + 1, ": ");
        for (j = 0; j < columns(mData); j++) {
            print(mData[i][j], " ");
        }
        println("");
    }
    
    return mData;
}

/**
 * Calculate portfolio weights (equal investment approach like original codes)
 */
calculatePortfolioWeights()
{
    decl vWeights = ones(g_iNumAssets, 1) / g_iNumAssets;  // Equal weights
    
    println("Portfolio weights (equal investment):");
    decl i;
    for (i = 0; i < g_iNumAssets; i++) {
        print("Asset ", i + 1, " (", g_asAssetNames[i], "): ", 
              vWeights[i][0] * 100, " percent\n");
    }
    
    return vWeights;
}

/**
 * Historical Simulation VaR using return data (adapted from pghs.ox logic)
 * Since we have returns instead of prices, we simulate portfolio returns directly
 */
calculateHistoricalVaR_Returns(const mHistReturns, const vWeights, const i)
{
    decl b, j;
    decl vPortfolioReturns = zeros(g_iNumSimulations, 1);
    
    decl iNumHistoricalDays = rows(mHistReturns);
    
    // Monte Carlo simulation using historical returns
    for (b = 0; b < g_iNumSimulations; b++) {
        decl dCumulativeReturn = 0;
        
        // Simulate g_iHoldingPeriod days using random historical returns
        for (j = 0; j < g_iHoldingPeriod; j++) {
            // Randomly select a historical day
            decl iRandomDay = floor(ranu(1, 1)[0] * iNumHistoricalDays);
            if (iRandomDay >= iNumHistoricalDays) iRandomDay = iNumHistoricalDays - 1;
            
            // Calculate portfolio return for this day
            decl vDayReturns = mHistReturns[iRandomDay][];
            decl dPortfolioReturn = (vWeights' * vDayReturns')[0];
            
            // Add to cumulative return (log returns are additive)
            dCumulativeReturn += dPortfolioReturn;
        }
        
        vPortfolioReturns[b] = dCumulativeReturn;
    }
    
    // Calculate VaR as 1% quantile of losses (negative returns)
    decl vSortedReturns = sortc(vPortfolioReturns);  // Sort returns (lowest first)
    decl iQuantileIndex = floor(0.01 * g_iNumSimulations);
    if (iQuantileIndex >= g_iNumSimulations) iQuantileIndex = g_iNumSimulations - 1;
    
    // VaR is the negative of the 1% quantile (worst loss), expressed as positive number
    decl dVaR = -vSortedReturns[iQuantileIndex] * g_iInvestmentAmount * g_iNumAssets;
    
    return dVaR;
}

/**
 * Parametric VaR using return data (adapted from pgrel1.ox logic)
 */
calculateParametricVaR_Returns(const mHistReturns, const vWeights, const i)
{
    // Calculate portfolio statistics from historical returns
    decl vMeanReturns = meanc(mHistReturns);
    decl mCovMatrix = variance(mHistReturns);
    
    // Portfolio statistics
    decl dPortfolioMean = (vWeights' * vMeanReturns')[0];
    decl dPortfolioVar = (vWeights' * mCovMatrix * vWeights)[0];
    
    // Scale for holding period (returns are already in log form, so we can scale directly)
    decl dScaledMean = dPortfolioMean * g_iHoldingPeriod;
    decl dScaledStd = sqrt(dPortfolioVar * g_iHoldingPeriod);
    
    // Calculate VaR using exact formula from pgrel1.ox:
    // VaR = -(mean + z_score * std) * portfolio_value
    decl dPortfolioValue = g_iInvestmentAmount * g_iNumAssets;
    decl dVaR = -(dScaledMean + g_dZScore99 * dScaledStd) * dPortfolioValue;
    
    return dVaR;
}

/**
 * Calculate actual portfolio loss over holding period using return data
 */
calculateActualLoss_Returns(const mAllReturns, const vWeights, const iStart)
{
    decl dCumulativeReturn = 0;
    decl j;
    
    // Sum returns over holding period
    for (j = 0; j < g_iHoldingPeriod; j++) {
        if (iStart + j >= rows(mAllReturns)) break;
        
        decl vDayReturns = mAllReturns[iStart + j][];
        decl dPortfolioReturn = (vWeights' * vDayReturns')[0];
        dCumulativeReturn += dPortfolioReturn;
    }
    
    // Convert to actual loss in currency units
    decl dPortfolioValue = g_iInvestmentAmount * g_iNumAssets;
    decl dActualLoss = -dCumulativeReturn * dPortfolioValue;
    
    return dActualLoss;
}

/**
 * Save results to CSV files
 */
saveResultsToCSV_WRDS(const mHistoricalResults, const mParametricResults, const mActualLosses)
{
    println("Saving Historical Simulation results to WRDS_pghs_results.csv...");
    decl fp1 = fopen("WRDS_pghs_results.csv", "w");
    fprint(fp1, "Period,Day_Index,Historical_VaR,Actual_Loss,Violation\n");
    
    decl i;
    for (i = 0; i < rows(mHistoricalResults); i++) {
        decl violation = (mActualLosses[i] > mHistoricalResults[i]) ? 1 : 0;
        fprint(fp1, i + 1, ",", g_iHistoricalPeriod + i + 1, ",", 
               mHistoricalResults[i], ",", mActualLosses[i], ",", violation, "\n");
    }
    fclose(fp1);
    
    println("Saving Parametric VaR results to WRDS_pgrel1_results.csv...");
    decl fp2 = fopen("WRDS_pgrel1_results.csv", "w");
    fprint(fp2, "Period,Day_Index,Parametric_VaR,Actual_Loss,Violation\n");
    
    for (i = 0; i < rows(mParametricResults); i++) {
        decl violation = (mActualLosses[i] > mParametricResults[i]) ? 1 : 0;
        fprint(fp2, i + 1, ",", g_iHistoricalPeriod + i + 1, ",", 
               mParametricResults[i], ",", mActualLosses[i], ",", violation, "\n");
    }
    fclose(fp2);
    
    println("Saving comparison results to WRDS_var_comparison.csv...");
    decl fp3 = fopen("WRDS_var_comparison.csv", "w");
    fprint(fp3, "Period,Day_Index,Historical_VaR,Parametric_VaR,Actual_Loss,Hist_Violation,Param_Violation,VaR_Difference\n");
    
    for (i = 0; i < rows(mHistoricalResults); i++) {
        decl hist_violation = (mActualLosses[i] > mHistoricalResults[i]) ? 1 : 0;
        decl param_violation = (mActualLosses[i] > mParametricResults[i]) ? 1 : 0;
        decl var_diff = mHistoricalResults[i] - mParametricResults[i];
        
        fprint(fp3, i + 1, ",", g_iHistoricalPeriod + i + 1, ",", 
               mHistoricalResults[i], ",", mParametricResults[i], ",", mActualLosses[i], ",",
               hist_violation, ",", param_violation, ",", var_diff, "\n");
    }
    fclose(fp3);
    
    println("All WRDS results saved successfully!");
}

/**
 * Display summary statistics for WRDS data
 */
displaySummaryStatistics_WRDS(const mHistoricalResults, const mParametricResults, const mActualLosses, const iTotalPeriods)
{
    println("\n============================================================");
    println("WRDS CRSP DATA - VaR ANALYSIS SUMMARY");
    println("Assets: Microsoft, Chevron, Apple, Disney, Walmart, Nike, American Express");
    println("Historical Period: ", g_iHistoricalPeriod, " days, Holding Period: ", g_iHoldingPeriod, " days");
    println("============================================================");
    
    // Calculate violations
    decl histViolations = 0, paramViolations = 0;
    decl i;
    for (i = 0; i < iTotalPeriods; i++) {
        if (mActualLosses[i] > mHistoricalResults[i]) histViolations++;
        if (mActualLosses[i] > mParametricResults[i]) paramViolations++;
    }
    
    println("Total test periods: ", iTotalPeriods);
    println("Expected violations (1%): ", iTotalPeriods * 0.01);
    println("");
    
    println("Historical Simulation Method:");
    print("- Mean VaR: ", meanc(mHistoricalResults)[0], "\n");
    print("- Min VaR:  ", minc(mHistoricalResults)[0], "\n");
    print("- Max VaR:  ", maxc(mHistoricalResults)[0], "\n");
    print("- Violations: ", histViolations, " (", histViolations / iTotalPeriods * 100, " percent)\n");
    println("");
    
    println("Parametric Method:");
    print("- Mean VaR: ", meanc(mParametricResults)[0], "\n");
    print("- Min VaR:  ", minc(mParametricResults)[0], "\n");
    print("- Max VaR:  ", maxc(mParametricResults)[0], "\n");
    print("- Violations: ", paramViolations, " (", paramViolations / iTotalPeriods * 100, " percent)\n");
    println("");
    
    println("Actual Losses:");
    print("- Mean loss: ", meanc(mActualLosses)[0], "\n");
    print("- Min loss:  ", minc(mActualLosses)[0], "\n");
    print("- Max loss:  ", maxc(mActualLosses)[0], "\n");
    println("");
    
    // Method comparison
    decl vDifferences = mHistoricalResults - mParametricResults;
    println("Method Comparison (Historical - Parametric):");
    print("- Mean difference: ", meanc(vDifferences)[0], "\n");
    
    // Calculate correlation
    decl meanHist = meanc(mHistoricalResults)[0];
    decl meanParam = meanc(mParametricResults)[0];
    decl varHist = variance(mHistoricalResults)[0][0];
    decl varParam = variance(mParametricResults)[0][0];
    decl stdHist = sqrt(varHist);
    decl stdParam = sqrt(varParam);
    decl covar = meanc((mHistoricalResults - meanHist) .* (mParametricResults - meanParam))[0];
    decl corr = covar / (stdHist * stdParam);
    print("- Correlation: ", corr, "\n");
    
    // Model performance assessment
    println("\nModel Performance Assessment:");
    decl dExpectedRate = 0.01;
    decl dHistRate = histViolations / iTotalPeriods;
    decl dParamRate = paramViolations / iTotalPeriods;
    
    print("Historical Simulation: ");
    if (fabs(dHistRate - dExpectedRate) < 0.01) {
        println("GOOD performance");
    } else if (dHistRate > dExpectedRate + 0.02) {
        println("POOR performance (too many violations)");
    } else if (dHistRate < dExpectedRate - 0.02) {
        println("CONSERVATIVE (too few violations)");
    } else {
        println("ACCEPTABLE performance");
    }
    
    print("Parametric Method: ");
    if (fabs(dParamRate - dExpectedRate) < 0.01) {
        println("GOOD performance");
    } else if (dParamRate > dExpectedRate + 0.02) {
        println("POOR performance (too many violations)");
    } else if (dParamRate < dExpectedRate - 0.02) {
        println("CONSERVATIVE (too few violations)");
    } else {
        println("ACCEPTABLE performance");
    }
}

// =================== MAIN FUNCTION ===================
main()
{
    // Initialize constants
    initializeConstants();
    
    // Set same random seed for reproducibility
    ranseed(967537412);
    
    println("============================================================");
    println("WRDS CRSP DATA - VaR CALCULATION");
    println("7 Stocks: Microsoft, Chevron, Apple, Disney, Walmart, Nike, American Express");
    println("Period: 2015-2020, Historical Window: 800 days, Holding Period: 10 days");
    println("============================================================");
    
    // Step 1: Load return data
    decl mReturns = loadReturnData();
    if (mReturns == <>) return;
    
    decl iTotalDays = rows(mReturns);
    println("Total trading days available: ", iTotalDays);
    
    // Check if we have enough data
    if (iTotalDays < g_iHistoricalPeriod + g_iHoldingPeriod) {
        println("Error: Insufficient data! Need at least ", g_iHistoricalPeriod + g_iHoldingPeriod, " days.");
        return;
    }
    
    // Step 2: Calculate portfolio weights
    decl vWeights = calculatePortfolioWeights();
    decl dTotalInvestment = g_iInvestmentAmount * g_iNumAssets;
    print("Total portfolio value: ", dTotalInvestment, " currency units\n");
    
    // Step 3: Calculate number of backtesting periods
    decl iBacktestPeriods = iTotalDays - g_iHistoricalPeriod - g_iHoldingPeriod + 1;
    println("Number of backtesting periods: ", iBacktestPeriods);
    
    if (iBacktestPeriods <= 0) {
        println("Error: Insufficient data for backtesting");
        return;
    }
    
    // Step 4: VaR calculations and backtesting
    println("\n============================================================");
    println("CALCULATING VaR FOR ", iBacktestPeriods, " PERIODS");
    println("Using ", g_iHistoricalPeriod, " days history for each prediction");
    println("============================================================");
    
    decl mHistoricalResults = zeros(iBacktestPeriods, 1);
    decl mParametricResults = zeros(iBacktestPeriods, 1);
    decl mActualLosses = zeros(iBacktestPeriods, 1);
    
    decl i;
    for (i = 0; i < iBacktestPeriods; i++) {
        // Progress indicator
        if ((i - (i / 50) * 50) == 0 || i == iBacktestPeriods - 1) {
            print("Progress: ", (i + 1), "/", iBacktestPeriods, " (", 
                  (i + 1) * 100.0 / iBacktestPeriods, " percent)\n");
        }
        
        // Get historical data window: i to i+g_iHistoricalPeriod-1
        decl mHistReturns = mReturns[i : i + g_iHistoricalPeriod - 1][];
        
        // Calculate Historical VaR
        mHistoricalResults[i] = calculateHistoricalVaR_Returns(mHistReturns, vWeights, i);
        
        // Calculate Parametric VaR
        mParametricResults[i] = calculateParametricVaR_Returns(mHistReturns, vWeights, i);
        
        // Calculate actual loss: from i+g_iHistoricalPeriod for g_iHoldingPeriod days
        decl iLossStart = i + g_iHistoricalPeriod;
        mActualLosses[i] = calculateActualLoss_Returns(mReturns, vWeights, iLossStart);
    }
    
    println("\nCalculations completed!");
    
    // Step 5: Save results and display statistics
    saveResultsToCSV_WRDS(mHistoricalResults, mParametricResults, mActualLosses);
    displaySummaryStatistics_WRDS(mHistoricalResults, mParametricResults, mActualLosses, iBacktestPeriods);
    
    println("\n============================================================");
    println("WRDS ANALYSIS COMPLETED SUCCESSFULLY");
    println("Results saved to:");
    println("- WRDS_pghs_results.csv (Historical Simulation)");
    println("- WRDS_pgrel1_results.csv (Parametric Method)");
    println("- WRDS_var_comparison.csv (Method Comparison)");
    println("============================================================");
} 